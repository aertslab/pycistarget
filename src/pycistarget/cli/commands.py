from __future__ import annotations

import contextlib
import logging
import multiprocessing
import os
import sys
import time
from collections.abc import Callable, Iterator
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Literal

import pandas as pd
import pyranges as pr

from pycistarget.motif_enrichment_cistarget import cisTarget, cisTargetDatabase
from pycistarget.motif_enrichment_dem import (
    DEM,
    DEMDatabase,
    DEMScores,
    get_foreground_and_background_regions,
)
from pycistarget.motif_enrichment_result import MotifEnrichmentResult, Species

log = logging.getLogger(__name__)


def _bed_name(path: str) -> str:
    """File name of a bed file without .bed / .bed.gz extension."""
    return Path(path).name.replace(".gz", "").replace(".bed", "")


def _write_result(
    result: MotifEnrichmentResult,
    output_folder: str,
    prefix: str,
    output_mode: Literal["tsv", "hdf5", "hdf5+"],
    write_html: bool,
):
    if output_mode not in ("tsv", "hdf5", "hdf5+"):
        raise ValueError(f"Output mode: {output_mode} is not supported!")
    out = Path(output_folder)
    out.mkdir(parents=True, exist_ok=True)
    if write_html:
        (out / f"{prefix}_{result.name}.html").write_text(result.to_html())
    if output_mode == "tsv":
        result.motif_enrichment.to_csv(
            out / f"{prefix}_{result.name}.tsv", sep="\t", header=True, index=True
        )
    else:
        result.write_hdf5(
            out / f"{prefix}_{result.name}.hdf5", mode="w" if output_mode == "hdf5" else "a"
        )


# Loggers with per region set messages, which are replaced by the progress bar.
_REGION_SET_LOGGERS = [
    "pycistarget.motif_enrichment_cistarget",
    "pycistarget.motif_enrichment_dem",
    "pycistarget.motif_enrichment_result",
]


@contextlib.contextmanager
def _progress(total: int, description: str) -> Iterator[Callable[[str, int], None]]:
    """
    Progress of analysing ``total`` region sets: yields a function to call when one is finished.

    Shows a progress bar on a terminal and logs a line per region set otherwise (e.g. in a log
    file). The INFO messages for each region set are hidden while the progress bar is shown.
    """
    from tqdm.auto import tqdm

    try:
        # A terminal without a width (e.g. some pseudo terminals) does not show a progress bar.
        show_bar = sys.stderr.isatty() and os.get_terminal_size(sys.stderr.fileno()).columns > 0
    except OSError:
        show_bar = False
    bar = tqdm(
        total=total,
        desc=description,
        disable=not show_bar,
        bar_format="{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} region sets "
        "[{elapsed}<{remaining}{postfix}]",
    )
    start = time.perf_counter()
    loggers = [logging.getLogger(name) for name in _REGION_SET_LOGGERS] if show_bar else []
    levels = [logger.level for logger in loggers]

    n_done = 0

    def done(name: str, n_motifs: int):
        nonlocal n_done
        n_done += 1
        bar.set_postfix_str(f"{name}: {n_motifs} enriched motifs", refresh=False)
        bar.update()
        if not show_bar:
            log.info(
                f"{description}: finished {name} ({n_motifs} enriched motifs) "
                f"[{n_done}/{total}, {time.perf_counter() - start:.0f} s]"
            )

    for logger in loggers:
        logger.setLevel(logging.WARNING)
    try:
        yield done
    finally:
        for logger, level in zip(loggers, levels, strict=True):
            logger.setLevel(level)
        bar.close()


# State shared with forked worker processes (the database is not copied).
_worker_state: dict[str, Any] = {}


def _init_worker(n_threads: int):
    import numba

    numba.set_num_threads(n_threads)


def _run_cistarget_region_set(i: int) -> tuple[str, int]:
    state = _worker_state
    result = cisTarget(
        region_set=state["region_sets"][i], name=state["names"][i], **state["kwargs"]
    )
    result.run_ctx(state["ctx_db"])
    _write_result(result, **state["output"])
    return result.name, len(result.motif_enrichment)


def _run_in_parallel(
    run: Callable[[int], tuple[str, int]], n: int, n_jobs: int, done: Callable[[str, int], None]
):
    """
    Run ``run(i)`` for i in 0..n-1 in ``n_jobs`` forked processes (sharing ``_worker_state``).

    The numba threads are divided over the processes.
    """
    import numba
    import numba.np.ufunc.parallel

    n_jobs = max(1, min(n_jobs, n))
    if n_jobs > 1 and numba.np.ufunc.parallel._is_initialized:
        # The OpenMP threading layer can deadlock in forked processes once it is used.
        log.warning(
            "Parallel numba code already ran in this process: analysing region sets one by one."
        )
        n_jobs = 1
    if n_jobs == 1:
        for i in range(n):
            done(*run(i))
        return
    # numba.get_num_threads() would start the threading layer in this process.
    n_threads = max(1, numba.config.NUMBA_NUM_THREADS // n_jobs)  # type: ignore[attr-defined]
    with ProcessPoolExecutor(
        n_jobs,
        mp_context=multiprocessing.get_context("fork"),
        initializer=_init_worker,
        initargs=(n_threads,),
    ) as executor:
        futures = [executor.submit(run, i) for i in range(n)]
        try:
            for future in as_completed(futures):
                done(*future.result())
        except BaseException:
            executor.shutdown(wait=True, cancel_futures=True)
            raise


def run_cistarget_command(
    cistarget_db_fname: str,
    bed_fname: str | list[str],
    output_folder: str,
    fraction_overlap_w_cistarget_database: float,
    auc_threshold: float,
    nes_threshold: float,
    rank_threshold: float,
    path_to_motif_annotations: str | None,
    annotation_version: str,
    motif_similarity_fdr: float,
    orthologous_identity_threshold: float,
    species: Species,
    annotations_to_use: list[str],
    name: str | None = None,
    output_mode: Literal["tsv", "hdf5", "hdf5+"] = "hdf5",
    write_html: bool = True,
    force: bool = False,
    n_jobs: int = 1,
):
    """
    Run pycisTarget on a region set.

    Parameters
    ----------
    cistarget_db_fname
        Path to the cisTarget rankings database.
    bed_fname
        Path to a bed file, or a list of paths to bed files, on which to run motif enrichment
        analysis. The analysis is run for each bed file separately (one output file per bed file),
        the database is read once for the regions of all bed files.
    output_folder
        Path to the folder where to store the results.
    fraction_overlap_w_cistarget_database
        Fraction of nucleotides, of regions in the bed file, that should overlap with regions
        in the cistarget database in order for them to be included in the analysis.
    auc_threshold
        Threshold on the AUC value for calling significant motifs.
    nes_threshold
        Threshold on the NES value for calling significant motifs.
        NES - Normalised Enrichment Score - is defined as (AUC - Avg(AUC)) / sd(AUC).
    rank_threshold
        The fraction of ranked regions to take into account when creating recovery curves.
    path_to_motif_annotations
        Path to the motif-to-TF annotations. Set to None to download these from the
        cisTarget webservers.
    annotation_version
        Version of the motif-to-TF annotation to use. This parameter is used
        to download the correct motif-to-TF data from the cisTarget webservers.
    motif_similarity_fdr
        Threshold on motif similarity scores for calling similar motifs.
    orthologous_identity_threshold
        Threshold on the protein-protein orthology score for calling orthologous motifs.
    species
        Species used for the analysis. This parameter is used to download the correct
        motif-to-TF annotations from the cisTarget webservers.
    annotations_to_use
        Which annotations to use for annotation motifs to TFs.
    name
        Name of this analysis. This name is appended to the output file name.
        Defaults to the file name of the bed file. Can only be set for a single bed file.
    output_mode
        tsv: tab-separated text file containing which motifs are enriched.
        hdf5: a new hdf5 file will be created containing the full results.
        hdf5+: an existing hdf5 file will be appended with the full motifs enrichment results.
    write_html
        Whether or not to save the results as an html file.
    force
        If the motif-to-TF annotations can not be loaded, continue without annotations and
        cistromes (with a warning) instead of failing.
    n_jobs
        Number of region sets to analyse in parallel (in forked processes which share the
        database in memory). The numba threads are divided over the processes.
    """
    if n_jobs < 1:
        raise ValueError("n_jobs should be at least 1.")
    bed_fnames = [bed_fname] if isinstance(bed_fname, str) else list(bed_fname)
    if name is not None and len(bed_fnames) > 1:
        raise ValueError("A name can only be given when analysing a single bed file.")
    names = [name] if name is not None else [_bed_name(f) for f in bed_fnames]
    if len(set(names)) != len(names):
        raise ValueError("Bed files should have different file names.")
    region_sets = [pr.read_bed(bed) for bed in bed_fnames]
    # Read the database once, for the regions of all region sets (each region once, so the mapping
    # of each region set to the database is the same as when it is analysed on its own).
    all_regions = pr.PyRanges(
        pd.concat(
            [
                rs.df[["Chromosome", "Start", "End"]].astype({"Chromosome": str})
                for rs in region_sets
            ]
        ).drop_duplicates(ignore_index=True)
    )
    ctx_db = cisTargetDatabase(
        fname=cistarget_db_fname,
        region_sets=all_regions,
        fraction_overlap=fraction_overlap_w_cistarget_database,
    )
    kwargs: dict[str, Any] = dict(
        species=species,
        auc_threshold=auc_threshold,
        nes_threshold=nes_threshold,
        rank_threshold=rank_threshold,
        path_to_motif_annotations=path_to_motif_annotations,
        annotation_version=annotation_version,
        annotation_to_use=annotations_to_use,
        motif_similarity_fdr=motif_similarity_fdr,
        orthologous_identity_threshold=orthologous_identity_threshold,
        force=force,
    )
    # Load the motif-to-TF annotations once (also for the worker processes), before the analysis.
    cisTarget(region_set=region_sets[0], name=names[0], **kwargs).load_motif_annotation()
    _worker_state.update(
        region_sets=region_sets,
        names=names,
        kwargs=kwargs,
        ctx_db=ctx_db,
        output=dict(
            output_folder=output_folder,
            prefix="motif_enrichment_cistarget",
            output_mode=output_mode,
            write_html=write_html,
        ),
    )
    try:
        with _progress(len(names), "cisTarget") as done:
            _run_in_parallel(_run_cistarget_region_set, len(names), n_jobs, done)
    finally:
        _worker_state.clear()


def run_dem_command(
    dem_db_fname: str,
    paths_to_foreground_bed: list[str],
    paths_to_background_bed: list[str],
    output_folder: str,
    species: Species,
    fraction_overlap_w_dem_database: float = 0.4,
    max_bg_regions: int | None = None,
    path_to_genome_annotation: str | None = None,
    balance_number_of_promoters: bool = False,
    promoter_space: int = 1_000,
    adjpval_thr: float = 0.05,
    log2fc_thr: float = 1.0,
    mean_fg_thr: float = 0.0,
    motif_hit_thr: float | None = None,
    path_to_motif_annotations: str | None = None,
    annotation_version: str = "v10nr_clust",
    annotations_to_use: list[str] | None = None,
    motif_similarity_fdr: float = 0.001,
    orthologous_identity_threshold: float = 0.0,
    seed: int = 555,
    name: str | None = None,
    output_mode: Literal["tsv", "hdf5", "hdf5+"] = "hdf5",
    write_html: bool = True,
    force: bool = False,
):
    """
    Run DEM analysis between foreground and background regions.

    See ``pycistarget dem --help`` for a description of the parameters.
    """
    if name is None:
        name = (
            "__".join(_bed_name(f) for f in paths_to_foreground_bed)
            + "__VS__"
            + "__".join(_bed_name(f) for f in paths_to_background_bed)
        )

    genome_annotation = (
        pd.read_table(path_to_genome_annotation) if path_to_genome_annotation is not None else None
    )
    dem_db = DEMDatabase(dem_db_fname, fraction_overlap=fraction_overlap_w_dem_database)
    _run_dem_contrast(
        dem_db=dem_db,
        scores=None,
        foreground_region_sets=[pr.read_bed(f) for f in paths_to_foreground_bed],
        background_region_sets=[pr.read_bed(f) for f in paths_to_background_bed],
        name=name,
        output_folder=output_folder,
        species=species,
        max_bg_regions=max_bg_regions,
        genome_annotation=genome_annotation,
        balance_number_of_promoters=balance_number_of_promoters,
        promoter_space=promoter_space,
        adjpval_thr=adjpval_thr,
        log2fc_thr=log2fc_thr,
        mean_fg_thr=mean_fg_thr,
        motif_hit_thr=motif_hit_thr,
        path_to_motif_annotations=path_to_motif_annotations,
        annotation_version=annotation_version,
        annotations_to_use=annotations_to_use,
        motif_similarity_fdr=motif_similarity_fdr,
        orthologous_identity_threshold=orthologous_identity_threshold,
        seed=seed,
        output_mode=output_mode,
        write_html=write_html,
        force=force,
    )


def run_dem_one_vs_all_command(
    dem_db_fname: str,
    paths_to_bed: list[str],
    output_folder: str,
    species: Species,
    fraction_overlap_w_dem_database: float = 0.4,
    max_bg_regions: int | None = None,
    path_to_genome_annotation: str | None = None,
    balance_number_of_promoters: bool = False,
    promoter_space: int = 1_000,
    adjpval_thr: float = 0.05,
    log2fc_thr: float = 1.0,
    mean_fg_thr: float = 0.0,
    motif_hit_thr: float | None = None,
    path_to_motif_annotations: str | None = None,
    annotation_version: str = "v10nr_clust",
    annotations_to_use: list[str] | None = None,
    motif_similarity_fdr: float = 0.001,
    orthologous_identity_threshold: float = 0.0,
    seed: int = 555,
    output_mode: Literal["tsv", "hdf5", "hdf5+"] = "hdf5",
    write_html: bool = True,
    force: bool = False,
):
    """
    Run DEM analysis for each bed file against all other bed files.

    The database is read once. The analysis of bed file ``x.bed`` is named ``x_vs_all``, giving
    the same result as ``run_dem_command`` with ``x.bed`` as foreground, all other bed files as
    background and ``name="x_vs_all"``.

    See ``pycistarget dem --help`` for a description of the parameters.
    """
    names = [_bed_name(f) for f in paths_to_bed]
    if len(paths_to_bed) < 2:
        raise ValueError("At least two bed files are needed for a one vs all analysis.")
    if len(set(names)) != len(names):
        raise ValueError("Bed files should have different file names.")
    genome_annotation = (
        pd.read_table(path_to_genome_annotation) if path_to_genome_annotation is not None else None
    )
    region_sets = [pr.read_bed(f) for f in paths_to_bed]
    dem_db = DEMDatabase(dem_db_fname, fraction_overlap=fraction_overlap_w_dem_database)
    # All database regions overlapping with any of the region sets.
    db_regions: set[str] = set()
    for region_set in region_sets:
        db_regions.update(dem_db.map_regions(region_set)["Query"])
    log.info(f"Reading the scores of {len(db_regions)} regions from the DEM database")
    scores = dem_db.load_scores(db_regions)
    for i, name in enumerate(names):
        _run_dem_contrast(
            dem_db=dem_db,
            scores=scores,
            foreground_region_sets=[region_sets[i]],
            background_region_sets=region_sets[:i] + region_sets[i + 1 :],
            name=f"{name}_vs_all",
            output_folder=output_folder,
            species=species,
            max_bg_regions=max_bg_regions,
            genome_annotation=genome_annotation,
            balance_number_of_promoters=balance_number_of_promoters,
            promoter_space=promoter_space,
            adjpval_thr=adjpval_thr,
            log2fc_thr=log2fc_thr,
            mean_fg_thr=mean_fg_thr,
            motif_hit_thr=motif_hit_thr,
            path_to_motif_annotations=path_to_motif_annotations,
            annotation_version=annotation_version,
            annotations_to_use=annotations_to_use,
            motif_similarity_fdr=motif_similarity_fdr,
            orthologous_identity_threshold=orthologous_identity_threshold,
            seed=seed,
            output_mode=output_mode,
            write_html=write_html,
            force=force,
        )


def _run_dem_contrast(
    dem_db: DEMDatabase,
    scores: DEMScores | None,
    foreground_region_sets: list[pr.PyRanges],
    background_region_sets: list[pr.PyRanges],
    name: str,
    output_folder: str,
    species: Species,
    max_bg_regions: int | None,
    genome_annotation: pd.DataFrame | None,
    balance_number_of_promoters: bool,
    promoter_space: int,
    adjpval_thr: float,
    log2fc_thr: float,
    mean_fg_thr: float,
    motif_hit_thr: float | None,
    path_to_motif_annotations: str | None,
    annotation_version: str,
    annotations_to_use: list[str] | None,
    motif_similarity_fdr: float,
    orthologous_identity_threshold: float,
    seed: int,
    output_mode: Literal["tsv", "hdf5", "hdf5+"],
    write_html: bool,
    force: bool,
):
    foreground_regions, background_regions = get_foreground_and_background_regions(
        foreground_region_sets=foreground_region_sets,
        background_region_sets=background_region_sets,
        max_bg_regions=max_bg_regions,
        genome_annotation=genome_annotation,
        balance_number_of_promoters=balance_number_of_promoters,
        promoter_space=promoter_space,
        seed=seed,
    )
    dem_result = DEM(
        foreground_regions=foreground_regions,
        background_regions=background_regions,
        name=name,
        species=species,
        adjpval_thr=adjpval_thr,
        log2fc_thr=log2fc_thr,
        mean_fg_thr=mean_fg_thr,
        motif_hit_thr=motif_hit_thr,
        path_to_motif_annotations=path_to_motif_annotations,
        annotation_version=annotation_version,
        annotation_to_use=annotations_to_use,
        motif_similarity_fdr=motif_similarity_fdr,
        orthologous_identity_threshold=orthologous_identity_threshold,
        force=force,
    )
    dem_result.run(dem_db, scores)
    _write_result(dem_result, output_folder, "motif_enrichment_dem", output_mode, write_html)
