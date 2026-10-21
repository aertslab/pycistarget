"""Utility functions: region name conversion, region mapping, motif-to-TF annotation and cistromes."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence

import pandas as pd
import pyranges as pr

ANNOTATION_COLUMNS = [
    "Direct_annot",
    "Motif_similarity_annot",
    "Orthology_annot",
    "Motif_similarity_and_Orthology_annot",
]

_MOTIF2TF_URL = (
    "https://resources.aertslab.org/cistarget/motif2tf/motifs-{version}-nr.{name}-m0.001-o0.0.tbl"
)
_SPECIES_TO_ANNOTATION_NAME = {
    "homo_sapiens": "hgnc",
    "mus_musculus": "mgi",
    "drosophila_melanogaster": "flybase",
}


def coord_to_region_names(coord: pr.PyRanges | pd.DataFrame) -> list[str]:
    """Convert coordinates to region names (``chrom:start-end``)."""
    if isinstance(coord, pr.PyRanges):
        coord = coord.df
    if len(coord) == 0:
        return []
    return [
        f"{chrom}:{start}-{end}"
        for chrom, start, end in zip(
            coord["Chromosome"].astype(str).tolist(),
            coord["Start"].tolist(),
            coord["End"].tolist(),
            strict=True,
        )
    ]


def region_set_motif_hits(
    regions_to_db: pd.DataFrame, motif_hits: dict[str, list[str]]
) -> dict[str, list[str]]:
    """
    Convert motif hits (database regions) to the regions of the region set mapping to them.

    ``regions_to_db`` maps region set regions ("Target") to database regions ("Query").
    """
    targets_per_query: dict[str, set[str]] = {}
    for query, target in zip(
        regions_to_db["Query"].tolist(), regions_to_db["Target"].tolist(), strict=True
    ):
        targets_per_query.setdefault(query, set()).add(target)
    return {
        motif: sorted(set().union(*(targets_per_query[hit] for hit in hits)))
        for motif, hits in motif_hits.items()
    }


def region_names_to_coordinates(region_names: Iterable[str]) -> pd.DataFrame:
    """Convert region names (``chrom:start-end``) to a dataframe with coordinates, indexed by region name."""
    region_names = [name for name in region_names if ":" in name]
    chromosomes, starts, ends = [], [], []
    for name in region_names:
        chrom, coord = name.split(":", 1)
        start, end = coord.split("-", 1)
        chromosomes.append(chrom)
        starts.append(int(start))
        ends.append(int(end))
    return pd.DataFrame(
        {"Chromosome": chromosomes, "Start": starts, "End": ends},
        index=pd.Index(region_names, dtype=object),
    )


def _to_pyranges(regions: pr.PyRanges | Sequence[str] | str) -> pr.PyRanges:
    if isinstance(regions, pr.PyRanges):
        return regions
    if isinstance(regions, str):
        return pr.read_bed(regions)
    return pr.PyRanges(region_names_to_coordinates(regions))


def target_to_query(
    target: pr.PyRanges | Sequence[str] | str,
    query: pr.PyRanges | Sequence[str] | str,
    fraction_overlap: float = 0.4,
) -> pd.DataFrame:
    """
    Map target regions to query regions (e.g. regions of a cisTarget database).

    A target and query region are matched when the overlap is more than ``fraction_overlap``
    of the length of either region.

    Returns a dataframe with columns "Target" and "Query" (region names).
    """
    join_pr = _to_pyranges(target).join(_to_pyranges(query), report_overlap=True)
    if len(join_pr) == 0:
        raise ValueError(
            "None of the regions in the bed file overlap with the cistarget database. "
            "Check whether the chromosome names of the input regions are correctly "
            "formatted and that you are using the correct database."
        )
    joined = join_pr.df
    overlap_query = joined["Overlap"] / (joined["End_b"] - joined["Start_b"])
    overlap_target = joined["Overlap"] / (joined["End"] - joined["Start"])
    joined = joined[(overlap_query > fraction_overlap) | (overlap_target > fraction_overlap)]
    chrom = joined["Chromosome"].astype(str)
    return pd.DataFrame(
        {
            "Target": (
                chrom + ":" + joined["Start"].astype(str) + "-" + joined["End"].astype(str)
            ).to_numpy(dtype=object),
            "Query": (
                chrom + ":" + joined["Start_b"].astype(str) + "-" + joined["End_b"].astype(str)
            ).to_numpy(dtype=object),
        }
    )


def _join_tfs(tfs: pd.Series) -> str:
    return ", ".join(sorted(set(tfs)))


def load_motif_annotations(
    specie: str,
    version: str = "v10nr_clust",
    fname: str | None = None,
    column_names: Sequence[str] = (
        "#motif_id",
        "gene_name",
        "motif_similarity_qvalue",
        "orthologous_identity",
        "description",
    ),
    motif_similarity_fdr: float = 0.001,
    orthologous_identity_threshold: float = 0.0,
) -> pd.DataFrame:
    """
    Load motif-to-TF annotations from a motif2TF snapshot.

    Parameters
    ----------
    specie
        Species to download annotations for (homo_sapiens, mus_musculus or drosophila_melanogaster).
        Only used when ``fname`` is None.
    version
        Motif collection version. Only used when ``fname`` is None.
    fname
        Path (or URL) to the motif2TF snapshot.
    column_names
        The names of the columns in the snapshot to load.
    motif_similarity_fdr
        The maximum FDR on motif similarity for annotations by motif similarity.
    orthologous_identity_threshold
        The minimum orthologous identity for annotations by orthology.

    Returns
    -------
    A dataframe indexed by motif ID with one column per annotation type, containing comma separated TF names.
    """
    if fname is None:
        if specie not in _SPECIES_TO_ANNOTATION_NAME:
            raise ValueError(
                f"Can not download motif annotations for species {specie!r}, choose one of "
                f"{list(_SPECIES_TO_ANNOTATION_NAME)} or provide a motif annotation file."
            )
        fname = _MOTIF2TF_URL.format(version=version, name=_SPECIES_TO_ANNOTATION_NAME[specie])
    df = pd.read_csv(fname, sep="\t", usecols=list(column_names))
    df = df.rename(
        columns={
            "#motif_id": "MotifID",
            "gene_name": "TF",
            "motif_similarity_qvalue": "MotifSimilarityQvalue",
            "orthologous_identity": "OrthologousIdentity",
            "description": "Annotation",
        }
    )
    df = df[
        (df["MotifSimilarityQvalue"] <= motif_similarity_fdr)
        & (df["OrthologousIdentity"] >= orthologous_identity_threshold)
    ]
    is_similar = df["Annotation"].str.contains("similar")
    is_orthologous = df["Annotation"].str.contains("orthologous")
    selections = {
        "Direct_annot": df["Annotation"] == "gene is directly annotated",
        "Motif_similarity_annot": is_similar & ~is_orthologous,
        "Orthology_annot": ~is_similar & is_orthologous,
        "Motif_similarity_and_Orthology_annot": is_similar & is_orthologous,
    }
    return pd.concat(
        [
            df[selection].groupby("MotifID")["TF"].agg(_join_tfs).rename(column)
            for column, selection in selections.items()
        ],
        axis=1,
        sort=False,
    )


def _split_tfs(annotation: pd.Series) -> set[str]:
    return {tf for tfs in annotation.dropna() for tf in tfs.split(", ")}


def get_TF_list(
    motif_enrichment_table: pd.DataFrame, annotation: Sequence[str] = ANNOTATION_COLUMNS
) -> list[str]:
    """Get all TFs annotated to the motifs in a motif enrichment table."""
    tfs: set[str] = set()
    for column in annotation:
        if column in motif_enrichment_table:
            tfs |= _split_tfs(motif_enrichment_table[column])
    return sorted(tfs)


def get_motifs_per_TF(
    motif_enrichment_table: pd.DataFrame,
    tf: str,
    motif_column: str,
    annotation: Sequence[str] = ANNOTATION_COLUMNS,
) -> list[str]:
    """
    Get the motifs annotated to a TF in a motif enrichment table.

    ``motif_column`` is the column containing the motif names, or "Index" when motifs are the index.
    """
    motifs: set[str] = set()
    for column in annotation:
        if column in motif_enrichment_table:
            hits = motif_enrichment_table[column].map(
                lambda tfs: isinstance(tfs, str) and tf in tfs.split(", ")
            )
            selected = motif_enrichment_table[hits.astype(bool)]
            motifs |= set(selected.index if motif_column == "Index" else selected[motif_column])
    return sorted(motifs)


def _get_motifs_per_tf(
    motif_enrichment_table: pd.DataFrame, annotation: Sequence[str]
) -> dict[str, set[str]]:
    """Get the motifs (index of the table) annotated to each TF by any of the annotation columns."""
    motifs_per_tf: dict[str, set[str]] = defaultdict(set)
    for column in annotation:
        if column in motif_enrichment_table:
            for motif, tfs in motif_enrichment_table[column].items():
                if isinstance(tfs, str):
                    for tf in tfs.split(", "):
                        motifs_per_tf[tf].add(str(motif))
    return motifs_per_tf


def get_cistrome_per_TF(
    motif_hits_dict: Mapping[str, Iterable[str]], motifs: Iterable[str]
) -> list[str]:
    """Get the union of the motif hits of a set of motifs."""
    regions: set[str] = set()
    for motif in motifs:
        if motif in motif_hits_dict:
            regions.update(motif_hits_dict[motif])
    return sorted(regions)


def get_cistromes_per_region_set(
    motif_enrichment_region_set: pd.DataFrame,
    motif_hits_regions_set: Mapping[str, Iterable[str]],
    annotation: Sequence[str] = ANNOTATION_COLUMNS,
) -> dict[str, list[str]]:
    """
    Get (direct and extended) cistromes for all TFs in a motif enrichment table.

    Direct cistromes contain the motif hits of motifs directly annotated to a TF.
    Extended cistromes (suffix "_extended") contain the motif hits of motifs annotated
    to a TF by any of the annotations in ``annotation``.
    Cistrome names get the number of regions appended, e.g. "TF_(123r)".
    """
    cistromes: dict[str, list[str]] = {}
    if "Direct_annot" in annotation:
        motifs_per_tf = _get_motifs_per_tf(motif_enrichment_region_set, ["Direct_annot"])
        for tf in sorted(motifs_per_tf):
            cistromes[tf] = get_cistrome_per_TF(motif_hits_regions_set, motifs_per_tf[tf])
    if "Direct_annot" not in annotation or len(annotation) > 1:
        motifs_per_tf = _get_motifs_per_tf(motif_enrichment_region_set, annotation)
        for tf in get_TF_list(motif_enrichment_region_set):
            cistromes[f"{tf}_extended"] = get_cistrome_per_TF(
                motif_hits_regions_set, motifs_per_tf.get(tf, set())
            )
    return {f"{name}_({len(regions)}r)": regions for name, regions in cistromes.items()}
