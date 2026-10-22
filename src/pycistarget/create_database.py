"""
Create cisTarget motif databases from a FASTA file and Cluster-Buster motif files.

This is a faster replacement for ``create_cistarget_motif_databases.py`` from
https://github.com/aertslab/create_cisTarget_databases. Motifs are scored with a reimplementation
of Cluster-Buster (``cbust -f 4 -c 0.0 -r 10000 -b <bg_padding> -t 1 [-l]``) which gives the same
scores, and the rankings are created with the same seeded tie-breaking, so the databases are the
same as those of ``create_cistarget_motif_databases.py``.
"""

from __future__ import annotations

import logging
import os
import random
import time
from collections import Counter, deque
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from os import PathLike
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.feather as pf

from pycistarget._cbust import Sequences, score_motifs

log = logging.getLogger(__name__)

#: Number of motif files scored per call to the scoring code (for progress messages).
_MOTIFS_PER_CHUNK = 32


def read_motif_list(motifs_dir: str | PathLike, motifs_list: str | PathLike) -> dict[str, Path]:
    """
    Motif IDs and their Cluster-Buster motif file (``<motifs_dir>/<motif ID>.cb``).

    The motif list has one motif ID per line (optionally with a ``.cb`` extension); empty lines and
    lines starting with ``#`` are ignored.
    """
    motifs: dict[str, Path] = {}
    with open(motifs_list) as fh:
        for line in fh:
            motif_id = line.rstrip()
            if not motif_id or motif_id.startswith("#"):
                continue
            motif_id = motif_id.removesuffix(".cb")
            path = Path(motifs_dir) / f"{motif_id}.cb"
            if not path.is_file():
                raise FileNotFoundError(
                    f'Cluster-Buster motif file "{path}" does not exist for motif {motif_id}.'
                )
            motifs[motif_id] = path
    if not motifs:
        raise ValueError(f'No motifs in "{motifs_list}".')
    return motifs


def _chunks(n: int, size: int) -> Iterator[slice]:
    for start in range(0, n, size):
        yield slice(start, min(start + size, n))


def score_regions(
    fasta: str | PathLike,
    motif_files: list[Path],
    bg_padding: int = 0,
    mask: bool = False,
    n_threads: int = 0,
) -> tuple[np.ndarray, list[str]]:
    """
    Score each sequence of a FASTA file with Cluster-Buster motif files.

    Parameters
    ----------
    fasta
        FASTA file with the regions to score. Sequence names (the first word of the title line)
        are the region IDs.
    motif_files
        Cluster-Buster motif files.
    bg_padding
        Number of bases at the start and end of each sequence which are only used for calculating
        the background nucleotide frequencies, not for scoring motifs.
    mask
        Score lowercase (masked) bases as N.
    n_threads
        Number of threads (0: all cores).

    Returns
    -------
    The top cluster score of each sequence (columns) for each motif file (rows), as a float32 array,
    and the region IDs (in FASTA order).
    """
    start_time = time.monotonic()
    sequences = Sequences(str(fasta), mask=mask)
    region_ids = sequences.names
    if len(set(region_ids)) != len(region_ids):
        duplicated = sorted(r for r, n in Counter(region_ids).items() if n > 1)
        raise ValueError(f'Region IDs are not unique in "{fasta}": {", ".join(duplicated[:5])}')
    log.info(
        f"Read {len(region_ids)} sequences from {fasta} in {time.monotonic() - start_time:.1f} s."
    )

    scores = np.empty((len(motif_files), len(region_ids)), dtype=np.float32)
    start_time = time.monotonic()
    for chunk in _chunks(len(motif_files), _MOTIFS_PER_CHUNK):
        scores[chunk] = score_motifs(
            sequences,
            [str(path) for path in motif_files[chunk]],
            bg_padding=bg_padding,
            n_threads=n_threads,
        )
        elapsed = time.monotonic() - start_time
        log.info(f"Scored {chunk.stop} of {len(motif_files)} motifs ({elapsed:.0f} s).")
    return scores, region_ids


def rankings_from_scores(
    scores: np.ndarray, seed: int, n_threads: int = 0, transposed: bool = False
) -> np.ndarray:
    """
    Rank regions by score (highest score: rank 0) for each motif (row of ``scores``).

    Ties are broken randomly, exactly like create_cisTarget_databases does: for each motif in turn,
    a permutation of the regions is drawn from ``numpy.random.default_rng(seed)``.

    With ``transposed``, the rankings are returned as a regions x motifs array.
    """
    n_motifs, n_regions = scores.shape
    dtype = np.int16 if n_regions <= 2**15 else np.int32
    rankings = np.empty((n_regions, n_motifs) if transposed else (n_motifs, n_regions), dtype=dtype)
    order = np.arange(n_regions, dtype=dtype)
    rng = np.random.default_rng(seed=seed)

    def rank(i: int, permutation: np.ndarray) -> None:
        ranking = np.empty(n_regions, dtype=dtype)
        ranking[permutation[(-scores[i])[permutation].argsort()]] = order
        if transposed:
            rankings[:, i] = ranking
        else:
            rankings[i] = ranking

    # The permutations are drawn in motif order (as the random state is shared), the sorting
    # (which releases the GIL) is done in parallel.
    n_workers = n_threads or os.cpu_count() or 1
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        pending: deque = deque()
        for i in range(n_motifs):
            if len(pending) == 2 * n_workers:
                pending.popleft().result()
            pending.append(executor.submit(rank, i, rng.permutation(n_regions)))
        for future in pending:
            future.result()
    return rankings


def database_filename(db_prefix: str | PathLike, columns: str, rows: str, kind: str) -> Path:
    """File name used by create_cisTarget_databases, e.g. ``<prefix>.regions_vs_motifs.rankings.feather``."""
    return Path(f"{db_prefix}.{columns}_vs_{rows}.{kind}.feather")


def _write_feather(
    path: Path, values: np.ndarray, column_names: list[str], index_name: str, index: list[str]
) -> None:
    """Write the rows of ``values`` as columns, followed by an index column (like create_cisTarget_databases)."""
    table = pa.Table.from_arrays(
        [pa.array(row) for row in values] + [pa.array(index, type=pa.string())],
        names=[*column_names, index_name],
    )
    pf.write_feather(table, path, compression="zstd", compression_level=6)


def create_cistarget_motif_databases(
    fasta: str | PathLike,
    motifs_dir: str | PathLike,
    motifs_list: str | PathLike,
    db_prefix: str | PathLike,
    bg_padding: int = 0,
    mask: bool = False,
    seed: int | None = None,
    n_threads: int = 0,
) -> list[Path]:
    """
    Create cisTarget motif databases (scores and rankings), like create_cistarget_motif_databases.py.

    Writes ``<db_prefix>.motifs_vs_regions.scores.feather``,
    ``<db_prefix>.regions_vs_motifs.scores.feather`` and
    ``<db_prefix>.regions_vs_motifs.rankings.feather``. Regions and motifs are sorted by name.

    Parameters
    ----------
    fasta
        FASTA file with the regions (and their background padding).
    motifs_dir
        Directory with Cluster-Buster motif files (``<motif ID>.cb``).
    motifs_list
        File with the motif IDs to score, one per line.
    db_prefix
        Prefix of the database files.
    bg_padding
        Number of bases at the start and end of each sequence which are only used for calculating
        the background nucleotide frequencies, not for scoring motifs.
    mask
        Score lowercase (masked) bases as N.
    seed
        Random seed for breaking ties when ranking. By default a random seed is used (and logged).
    n_threads
        Number of threads (0: all cores).

    Returns
    -------
    The paths of the written database files.
    """
    if seed is None:
        seed = random.randint(0, 2**64)
    motifs = read_motif_list(motifs_dir, motifs_list)
    motif_ids = sorted(motifs)
    log.info(f"Scoring {len(motif_ids)} motifs (rankings seed: {seed}).")

    scores, region_ids = score_regions(
        fasta, [motifs[m] for m in motif_ids], bg_padding, mask, n_threads
    )
    region_order = sorted(range(len(region_ids)), key=region_ids.__getitem__)
    if region_order != list(range(len(region_ids))):
        scores = scores[:, region_order]
        region_ids = [region_ids[i] for i in region_order]

    paths = []
    start_time = time.monotonic()
    path = database_filename(db_prefix, "motifs", "regions", "scores")
    _write_feather(path, scores, motif_ids, "regions", region_ids)
    paths.append(path)

    path = database_filename(db_prefix, "regions", "motifs", "scores")
    _write_feather(path, np.ascontiguousarray(scores.T), region_ids, "motifs", motif_ids)
    paths.append(path)
    log.info(f"Wrote the scores databases in {time.monotonic() - start_time:.0f} s.")

    start_time = time.monotonic()
    rankings = rankings_from_scores(scores, seed, n_threads, transposed=True)
    del scores
    path = database_filename(db_prefix, "regions", "motifs", "rankings")
    _write_feather(path, rankings, region_ids, "motifs", motif_ids)
    paths.append(path)
    log.info(f"Created and wrote the rankings database in {time.monotonic() - start_time:.0f} s.")
    return paths
