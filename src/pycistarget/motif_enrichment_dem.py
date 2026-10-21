"""Differentially Enriched Motif (DEM) analysis: compare motif scores of foreground and background regions."""

from __future__ import annotations

import logging
import math
import random
from collections.abc import Iterable
from os import PathLike

import numba
import numpy as np
import pandas as pd
import pyranges as pr
from sklearn.metrics import roc_curve

from pycistarget.database import FeatherCisTargetDatabase
from pycistarget.motif_enrichment_result import MotifEnrichmentResult, Species
from pycistarget.utils import (
    coord_to_region_names,
    region_names_to_coordinates,
    region_set_motif_hits,
    target_to_query,
)

log = logging.getLogger(__name__)


@numba.njit(parallel=True, cache=True)
def mean_axis1(arr):
    """
    Calculate the mean of each row of a 2D-numpy matrix with numba.

    Values are summed in float64. This is spelled out instead of using `np.mean`, as the
    accumulator precision of numba's `np.mean` differs between numba versions and between
    parallel and non-parallel code (results would change slightly).
    """
    mean_axis1_array = np.empty(arr.shape[0], dtype=np.float64)
    for i in numba.prange(arr.shape[0]):
        row_sum = 0.0
        for j in range(arr.shape[1]):
            row_sum += np.float64(arr[i, j])
        mean_axis1_array[i] = row_sum / arr.shape[1]
    return mean_axis1_array


@numba.njit(cache=True)
def get_log2_fc(fg_mat, bg_mat):
    """
    Calculate log2 fold change between the row means of a foreground and background matrix.
    """
    if fg_mat.shape[0] != bg_mat.shape[0]:
        raise ValueError(
            "Foreground matrix and background matrix have a different first dimensions!"
        )
    return np.log2((mean_axis1(fg_mat) + 10**-12) / (mean_axis1(bg_mat) + 10**-12))


@numba.njit(cache=True)
def rankdata_average_numba(arr: np.ndarray):
    """
    Rank data, assigning average ranks to ties (based on scipy.stats.rankdata).
    """
    sorter = np.argsort(arr, kind="quicksort")
    inv = np.empty(sorter.size, dtype=np.intp)
    inv[sorter] = np.arange(sorter.size, dtype=np.intp)
    arr = arr[sorter]
    obs = np.empty(arr.shape, dtype=np.intp)
    obs[0] = True
    obs[1:] = arr[1:] != arr[:-1]
    dense = obs.cumsum()[inv]
    non_zero = np.nonzero(obs)[0]
    count = np.empty(non_zero.shape[0] + 1)
    count[0:-1] = non_zero
    count[-1] = len(obs)
    result = 0.5 * (count[dense] + count[dense - 1] + 1)
    return result


@numba.njit(cache=True)
def norm_sf(z):
    return (1.0 + math.erf(-z / math.sqrt(2.0))) / 2.0


@numba.njit(cache=True)
def ranksums_numba(x: np.ndarray, y: np.ndarray):
    """
    Wilcoxon rank-sum test (based on scipy.stats.ranksums).
    """
    n1 = len(x)
    n2 = len(y)
    alldata = np.concatenate((x, y))
    ranked = rankdata_average_numba(alldata)
    x = ranked[:n1]
    s = np.sum(x, axis=0)
    expected = n1 * (n1 + n2 + 1) / 2.0
    z = (s - expected) / np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    prob = 2 * norm_sf(abs(z))
    return z, prob


@numba.njit(parallel=True, cache=True)
def ranksums_numba_multiple(X: np.ndarray, Y: np.ndarray):
    """
    Wilcoxon rank-sum test between each row of X and the corresponding row of Y.
    """
    n = X.shape[0]
    if Y.shape[0] != n:
        raise ValueError("X and Y should have the same shape on dimension 0")
    ranksums_z = np.zeros((n), dtype=np.float64)
    ranksums_p = np.zeros((n), dtype=np.float64)
    for i in numba.prange(n):
        z, p = ranksums_numba(X[i], Y[i])
        ranksums_z[i] = z
        ranksums_p[i] = p
    return ranksums_z, ranksums_p


# Size of the hash table used to count the distinct scores of a motif. Cluster-Buster scores
# in cisTarget databases have a few thousand distinct values per motif; motifs with more
# distinct values than half the table size use the sorting based rank-sum test.
_HASH_BITS = 14
_HASH_TABLE_SIZE = 1 << _HASH_BITS
_HASH_MAX_DISTINCT = _HASH_TABLE_SIZE // 2


@numba.njit(cache=True)
def _count_distinct_scores(row_bits, fg_idx, bg_idx, keys, used, fg_counts, bg_counts):
    """
    Count how often each distinct score occurs in the foreground and background columns of a row.

    Scores are hashed on their float32 bit pattern. Returns the number of distinct scores, or -1
    when there are too many distinct scores or NaNs.
    """
    shift = np.uint64(32 - _HASH_BITS)
    mask = keys.shape[0] - 1
    n_distinct = 0
    for which in range(2):
        idx = fg_idx if which == 0 else bg_idx
        counts = fg_counts if which == 0 else bg_counts
        for j in range(idx.shape[0]):
            bits = row_bits[idx[j]]
            if (bits & np.uint32(0x7FFFFFFF)) > np.uint32(0x7F800000):
                return -1  # NaN
            if bits == np.uint32(0x80000000):
                bits = np.uint32(0)  # -0.0 is tied with 0.0
            # Fibonacci hashing of the 32 bits.
            slot = np.int64(
                ((np.uint64(bits) * np.uint64(2654435761)) & np.uint64(0xFFFFFFFF)) >> shift
            )
            while used[slot] and keys[slot] != bits:
                slot = (slot + 1) & mask
            if not used[slot]:
                used[slot] = True
                keys[slot] = bits
                n_distinct += 1
                if n_distinct > _HASH_MAX_DISTINCT:
                    return -1
            counts[slot] += 1
    return n_distinct


@numba.njit(cache=True)
def _ranksums_gather(row, fg_idx, bg_idx):
    """
    Wilcoxon rank-sum test between ``row[fg_idx]`` and ``row[bg_idx]`` (float32 scores).

    Gives exactly the same result as ``ranksums_numba``, but counts distinct scores instead of
    sorting all scores: twice the rank sum of the foreground is an integer, so it is computed
    exactly and the statistic is identical to the one from summing (half integer) average ranks.
    """
    keys = np.empty(_HASH_TABLE_SIZE, dtype=np.uint32)
    used = np.zeros(_HASH_TABLE_SIZE, dtype=np.bool_)
    fg_counts = np.zeros(_HASH_TABLE_SIZE, dtype=np.int64)
    bg_counts = np.zeros(_HASH_TABLE_SIZE, dtype=np.int64)
    if (
        _count_distinct_scores(
            row.view(np.uint32), fg_idx, bg_idx, keys, used, fg_counts, bg_counts
        )
        < 0
    ):
        return ranksums_numba(row[fg_idx], row[bg_idx])
    slots = np.nonzero(used)[0]
    slots = slots[np.argsort(keys[slots].view(np.float32))]
    n_before = 0
    twice_rank_sum = 0
    for slot in slots:
        n_tied = fg_counts[slot] + bg_counts[slot]
        # Average rank (1-based) of the tied scores is n_before + (n_tied + 1) / 2.
        twice_rank_sum += fg_counts[slot] * (2 * n_before + n_tied + 1)
        n_before += n_tied
    n1 = len(fg_idx)
    n2 = len(bg_idx)
    s = twice_rank_sum / 2.0
    expected = n1 * (n1 + n2 + 1) / 2.0
    z = (s - expected) / np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    prob = 2 * norm_sf(abs(z))
    return z, prob


@numba.njit(cache=True)
def _mean_gather(row, idx):
    """Mean of ``row[idx]``, summed in float64 in the order of ``idx`` (as ``mean_axis1``)."""
    total = 0.0
    for j in range(idx.shape[0]):
        total += np.float64(row[idx[j]])
    return total / idx.shape[0]


@numba.njit(parallel=True, cache=True)
def _means_gather(scores, fg_idx, bg_idx, motif_idx):
    mean_fg = np.empty(motif_idx.shape[0], dtype=np.float64)
    mean_bg = np.empty(motif_idx.shape[0], dtype=np.float64)
    for i in numba.prange(motif_idx.shape[0]):
        mean_fg[i] = _mean_gather(scores[motif_idx[i]], fg_idx)
        mean_bg[i] = _mean_gather(scores[motif_idx[i]], bg_idx)
    return mean_fg, mean_bg


@numba.njit(parallel=True, cache=True)
def _ranksums_counting_pvalues(scores, fg_idx, bg_idx, motif_idx):
    pvalues = np.empty(motif_idx.shape[0], dtype=np.float64)
    for i in numba.prange(motif_idx.shape[0]):
        _, pvalues[i] = _ranksums_gather(scores[motif_idx[i]], fg_idx, bg_idx)
    return pvalues


@numba.njit(parallel=True, cache=True)
def _ranksums_sorting_pvalues(scores, fg_idx, bg_idx, motif_idx):
    pvalues = np.empty(motif_idx.shape[0], dtype=np.float64)
    for i in numba.prange(motif_idx.shape[0]):
        row = scores[motif_idx[i]]
        _, pvalues[i] = ranksums_numba(row[fg_idx], row[bg_idx])
    return pvalues


def _dem_statistics(
    scores: np.ndarray, fg_idx: np.ndarray, bg_idx: np.ndarray, motif_idx: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Mean foreground score, mean background score and rank-sum p-value of each motif.

    ``scores`` is a C-contiguous motifs x regions matrix, the foreground and background
    are columns ``fg_idx`` and ``bg_idx`` and the motifs are rows ``motif_idx``.
    """
    mean_fg, mean_bg = _means_gather(scores, fg_idx, bg_idx, motif_idx)
    # Counting distinct scores relies on the float32 bit pattern.
    ranksums_pvalues = (
        _ranksums_counting_pvalues if scores.dtype == np.float32 else _ranksums_sorting_pvalues
    )
    return mean_fg, mean_bg, ranksums_pvalues(scores, fg_idx, bg_idx, motif_idx)


@numba.njit(cache=True)
def _log2_fc_from_means(mean_fg, mean_bg):
    return np.log2((mean_fg + 10**-12) / (mean_bg + 10**-12))


@numba.njit(cache=True)
def _optimal_threshold_roc_gather(row, fg_idx, bg_idx):
    """
    Same as ``get_optimal_threshold_roc(row[fg_idx], row[bg_idx])``: the score with the highest
    tpr - fpr on the ROC curve of ``sklearn.metrics.roc_curve`` (``drop_intermediate=True``),
    only considering scores above zero.

    Returns NaN when the foreground or background has no scores above zero.
    """
    n_fg = 0
    for j in range(fg_idx.shape[0]):
        if row[fg_idx[j]] > 0:
            n_fg += 1
    n_bg = 0
    for j in range(bg_idx.shape[0]):
        if row[bg_idx[j]] > 0:
            n_bg += 1
    if n_fg == 0 or n_bg == 0:
        return np.nan
    n = n_fg + n_bg
    scores = np.empty(n, dtype=np.float64)
    is_fg = np.zeros(n, dtype=np.int64)
    k = 0
    for j in range(fg_idx.shape[0]):
        if row[fg_idx[j]] > 0:
            scores[k] = row[fg_idx[j]]
            is_fg[k] = 1
            k += 1
    for j in range(bg_idx.shape[0]):
        if row[bg_idx[j]] > 0:
            scores[k] = row[bg_idx[j]]
            k += 1
    order = np.argsort(scores)[::-1]
    # True and false positives at each distinct score, from high to low scores.
    thresholds = np.empty(n, dtype=np.float64)
    tps = np.empty(n, dtype=np.int64)
    fps = np.empty(n, dtype=np.int64)
    m = 0
    tp = 0
    for i in range(n):
        tp += is_fg[order[i]]
        if i == n - 1 or scores[order[i + 1]] != scores[order[i]]:
            thresholds[m] = scores[order[i]]
            tps[m] = tp
            fps[m] = i + 1 - tp
            m += 1
    # Points that are dropped by roc_curve(drop_intermediate=True).
    keep = np.ones(m, dtype=np.bool_)
    if m > 2:
        for i in range(1, m - 1):
            keep[i] = (fps[i + 1] - 2 * fps[i] + fps[i - 1] != 0) or (
                tps[i + 1] - 2 * tps[i] + tps[i - 1] != 0
            )
    # roc_curve starts at (0, 0) with threshold inf; argmax returns the first maximum.
    best_j = 0.0
    best_threshold = np.inf
    for i in range(m):
        if keep[i]:
            j = tps[i] / np.float64(tps[m - 1]) - fps[i] / np.float64(fps[m - 1])
            if j > best_j:
                best_j = j
                best_threshold = thresholds[i]
    return best_threshold


@numba.njit(parallel=True, cache=True)
def _optimal_thresholds_roc(scores, fg_idx, bg_idx, motif_idx):
    thresholds = np.empty(motif_idx.shape[0], dtype=np.float64)
    for i in numba.prange(motif_idx.shape[0]):
        thresholds[i] = _optimal_threshold_roc_gather(scores[motif_idx[i]], fg_idx, bg_idx)
    return thresholds


def p_adjust_bh(p: np.ndarray) -> np.ndarray:
    """
    Benjamini-Hochberg p-value correction for multiple hypothesis testing.
    """
    p = np.asarray(p, dtype=np.float64)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def get_optimal_threshold_roc(
    foreground_scores_motif: np.ndarray, background_scores_motif: np.ndarray
) -> float:
    """
    Get the motif score threshold that best separates foreground from background regions (Youden's J).

    Only regions with a score above zero are considered.
    """
    foreground_scores_motif = foreground_scores_motif[foreground_scores_motif > 0]
    background_scores_motif = background_scores_motif[background_scores_motif > 0]
    all_scores = np.concatenate([foreground_scores_motif, background_scores_motif]).astype(float)
    labels = np.zeros_like(all_scores, dtype=int)
    labels[0 : foreground_scores_motif.shape[0]] = 1
    fpr, tpr, thresholds = roc_curve(labels, all_scores)
    return thresholds[np.argmax(tpr - fpr)]


class DEMScores:
    """
    Motif scores of a set of database regions, loaded in memory.

    Attributes
    ----------
    scores
        C-contiguous matrix with motifs as rows and regions as columns.
    motifs
        Motif names (rows).
    regions
        Database region names (columns).
    """

    def __init__(self, scores: np.ndarray, motifs: pd.Index, regions: list[str]):
        self.scores = np.ascontiguousarray(scores)
        self.motifs = motifs
        self.regions = regions
        self._column = {region: i for i, region in enumerate(regions)}

    def columns(self, regions: Iterable[str]) -> np.ndarray:
        """Column indices of ``regions`` (deduplicated and in increasing order)."""
        try:
            return np.array(sorted({self._column[r] for r in regions}), dtype=np.int64)
        except KeyError as e:
            raise ValueError(f"Region {e} was not loaded from the DEM database.") from None


class DEMDatabase(FeatherCisTargetDatabase):
    """
    cisTarget scores database (regions_vs_motifs) used for DEM analysis.

    Parameters
    ----------
    fname
        Path to a ``*.regions_vs_motifs.scores.feather`` file.
    fraction_overlap
        Minimal overlap between query regions and regions in the database for the mapping.
    """

    def __init__(self, fname: str | PathLike, fraction_overlap: float = 0.4):
        super().__init__(fname)
        self.db_regions = pr.PyRanges(region_names_to_coordinates(self.region_ids))
        self.fraction_overlap = fraction_overlap

    def map_regions(self, regions: pr.PyRanges) -> pd.DataFrame:
        """Mapping of ``regions`` ("Target") to database regions ("Query")."""
        return target_to_query(
            target=regions, query=self.db_regions, fraction_overlap=self.fraction_overlap
        )

    def load_scores(self, region_ids: Iterable[str]) -> DEMScores:
        """Load the scores of the database regions ``region_ids`` in memory."""
        scores, motifs, regions = self.load_array(region_ids)
        return DEMScores(scores, motifs, regions)

    def get_scores(self, regions: pr.PyRanges) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Get the scores for the database regions overlapping with ``regions``.

        Returns the mapping of regions to database regions and the scores (motifs x database regions).
        """
        target_to_db = self.map_regions(regions)
        return target_to_db, self.load(set(target_to_db["Query"]))


class DEM(MotifEnrichmentResult):
    """
    Differentially Enriched Motif analysis between foreground and background regions.

    Motif scores of foreground and background regions are compared with a Wilcoxon rank-sum test,
    motifs are enriched when the Benjamini-Hochberg adjusted p-value, log2 fold change and
    mean foreground score pass the thresholds.
    """

    _analysis_type = "DEM"

    def __init__(
        self,
        foreground_regions: pr.PyRanges,
        background_regions: pr.PyRanges,
        name: str,
        species: Species,
        adjpval_thr: float = 0.05,
        log2fc_thr: float = 1.0,
        mean_fg_thr: float = 0.0,
        motif_hit_thr: float | None = None,
        path_to_motif_annotations: str | None = None,
        annotation_version: str = "v10nr_clust",
        annotation_to_use: list[str] | None = None,
        motif_similarity_fdr: float = 0.001,
        orthologous_identity_threshold: float = 0.0,
        motifs_to_use: list[str] | None = None,
    ):
        self.foreground_regions = foreground_regions
        self.background_regions = background_regions
        self.adjpval_thr = adjpval_thr
        self.log2fc_thr = log2fc_thr
        self.mean_fg_thr = mean_fg_thr
        self.motif_hit_thr = motif_hit_thr
        super().__init__(
            name=name,
            species=species,
            path_to_motif_annotations=path_to_motif_annotations,
            annotation_version=annotation_version,
            annotation_to_use=(
                annotation_to_use
                if annotation_to_use is not None
                else [
                    "Direct_annot",
                    "Motif_similarity_annot",
                    "Orthology_annot",
                    "Motif_similarity_and_Orthology_annot",
                ]
            ),
            motif_similarity_fdr=motif_similarity_fdr,
            orthologous_identity_threshold=orthologous_identity_threshold,
            motifs_to_use=motifs_to_use,
        )

    def run(self, dem_db: DEMDatabase, scores: DEMScores | None = None):
        """
        Run DEM.

        Parameters
        ----------
        dem_db
            The DEM (scores) database.
        scores
            Scores loaded with ``dem_db.load_scores``, for (at least) the database regions
            overlapping with the foreground and background regions. Pass these to run several
            analyses without reading the database again. By default only the database regions
            needed for this analysis are read.
        """
        log.info(f"Running DEM for {self.name}")
        foreground_region_to_db = dem_db.map_regions(self.foreground_regions)
        background_region_to_db = dem_db.map_regions(self.background_regions)
        self.regions_to_db = (
            pd.concat([foreground_region_to_db, background_region_to_db])
            .drop_duplicates()
            .reset_index(drop=True)
        )
        if scores is None:
            scores = dem_db.load_scores(
                set(foreground_region_to_db["Query"]) | set(background_region_to_db["Query"])
            )
        fg_idx = scores.columns(foreground_region_to_db["Query"])
        bg_idx = scores.columns(background_region_to_db["Query"])
        if self.motifs_to_use is not None:
            motif_idx = np.flatnonzero(scores.motifs.isin(self.motifs_to_use))
        else:
            motif_idx = np.arange(len(scores.motifs))
        motif_names = scores.motifs[motif_idx]

        # Wilcoxon rank-sum test, log2 fold change and Benjamini-Hochberg correction.
        mean_fg, mean_bg, pvalues = _dem_statistics(scores.scores, fg_idx, bg_idx, motif_idx)
        result = pd.DataFrame(
            data={
                "Log2FC": _log2_fc_from_means(mean_fg, mean_bg),
                "Adjusted_pval": p_adjust_bh(pvalues),
                "Mean_fg": mean_fg,
                "Mean_bg": mean_bg,
            },
            index=pd.Index(motif_names, dtype=object, name=motif_names.name),
        )
        result = result.loc[
            (result["Adjusted_pval"] <= self.adjpval_thr)
            & (result["Log2FC"] >= self.log2fc_thr)
            & (result["Mean_fg"] >= self.mean_fg_thr)
        ]
        self.motif_enrichment = result.sort_values(
            ["Log2FC", "Adjusted_pval"], ascending=[False, True]
        )
        log.info("Adding motif-to-TF annotation")
        self.add_motif_annotation()

        # Motif hits: foreground regions with a score above the motif hit threshold.
        enriched_idx = scores.motifs.get_indexer(self.motif_enrichment.index)
        if self.motif_hit_thr is None:
            motif_hit_thr = _optimal_thresholds_roc(scores.scores, fg_idx, bg_idx, enriched_idx)
            for i in np.flatnonzero(np.isnan(motif_hit_thr)):
                # No scores above zero in the foreground or background.
                row = scores.scores[enriched_idx[i]]
                motif_hit_thr[i] = get_optimal_threshold_roc(row[fg_idx], row[bg_idx])
        else:
            motif_hit_thr = np.full(len(enriched_idx), self.motif_hit_thr, dtype=np.float64)
        fg_regions = np.array(scores.regions, dtype=object)[fg_idx]
        motif_hits: dict[str, list[str]] = {
            motif: fg_regions[scores.scores[row, fg_idx] > thr].tolist()
            for motif, row, thr in zip(
                self.motif_enrichment.index, enriched_idx, motif_hit_thr, strict=True
            )
        }
        self.motif_enrichment["Motif_hit_thr"] = motif_hit_thr
        self.motif_enrichment["Motif_hits"] = [
            len(motif_hits[m]) for m in self.motif_enrichment.index
        ]

        rs_motif_hits = region_set_motif_hits(self.regions_to_db, motif_hits)
        self.motif_hits = {"Database": motif_hits, "Region_set": rs_motif_hits}
        self.get_cistromes()


def get_foreground_and_background_regions(
    foreground_region_sets: list[pr.PyRanges],
    background_region_sets: list[pr.PyRanges],
    max_bg_regions: int | None = None,
    genome_annotation: pd.DataFrame | None = None,
    balance_number_of_promoters: bool = True,
    promoter_space: int = 1_000,
    seed: int = 555,
) -> tuple[pr.PyRanges, pr.PyRanges]:
    """
    Get foreground and background regions for DEM analysis.

    Foreground regions are removed from the background. Optionally the background is
    subsampled to ``max_bg_regions`` regions and/or the fraction of promoter regions in
    the background is matched to the fraction of promoters in the foreground.
    """
    foreground_region_names = set(
        coord_to_region_names(pd.concat([ranges.df for ranges in foreground_region_sets]))
    )
    background_region_names_set = (
        set(coord_to_region_names(pd.concat([ranges.df for ranges in background_region_sets])))
        - foreground_region_names
    )
    if max_bg_regions is not None:
        max_bg_regions = min(max_bg_regions, len(background_region_names_set))

    # Sort regions, so shuffling with the same seed results in the same regions.
    foreground_regions_sorted = sorted(foreground_region_names)
    background_region_names = sorted(background_region_names_set)

    if balance_number_of_promoters:
        if genome_annotation is None:
            raise ValueError(
                "Please provide genome_annotation if you wish to balance the number of promoters"
            )
        tmp_pr_foreground = pr.PyRanges(region_names_to_coordinates(foreground_regions_sorted))
        tmp_pr_background = pr.PyRanges(region_names_to_coordinates(background_region_names))

        # Define promoter locations
        promoter_annotation = genome_annotation.copy()
        promoter_annotation["End"] = promoter_annotation["Start"] + promoter_space
        promoter_annotation["Start"] = promoter_annotation["Start"] - promoter_space
        pr_promoter_annotation = pr.PyRanges(
            promoter_annotation[["Chromosome", "Start", "End"]].drop_duplicates()
        )

        # Calculate fraction of promoters in foreground
        nr_promoter_in_fg = len(
            tmp_pr_foreground.count_overlaps(pr_promoter_annotation).df.query("NumberOverlaps > 0")
        )
        fraction_promoter_fg = nr_promoter_in_fg / len(foreground_regions_sorted)

        # Calculate number of allowed promoter regions in background
        n_background = (
            max_bg_regions if max_bg_regions is not None else len(background_region_names)
        )
        nr_allowed_promoters_background = int(n_background * fraction_promoter_fg)

        # Separate promoter and non promoter regions in background
        tmp_background_promoter_overlap = tmp_pr_background.count_overlaps(
            pr_promoter_annotation
        ).df
        background_region_names_no_promoter = sorted(
            coord_to_region_names(tmp_background_promoter_overlap.query("NumberOverlaps == 0"))
        )
        background_region_names_promoter = sorted(
            coord_to_region_names(tmp_background_promoter_overlap.query("NumberOverlaps > 0"))
        )

        # Sample promoter regions from background
        nr_allowed_promoters_background = min(
            len(background_region_names_promoter), nr_allowed_promoters_background
        )
        random.Random(seed).shuffle(background_region_names_promoter)
        background_region_names_promoter = background_region_names_promoter[
            0:nr_allowed_promoters_background
        ]

        # Sample non promoter regions from background
        nr_non_promoter_regions_background = n_background - nr_allowed_promoters_background
        random.Random(seed).shuffle(background_region_names_no_promoter)
        background_region_names_no_promoter = background_region_names_no_promoter[
            0:nr_non_promoter_regions_background
        ]

        background_region_names = [
            *background_region_names_promoter,
            *background_region_names_no_promoter,
        ]
    elif max_bg_regions is not None:
        random.Random(seed).shuffle(background_region_names)
        background_region_names = background_region_names[0:max_bg_regions]

    foreground_regions = pr.PyRanges(region_names_to_coordinates(foreground_regions_sorted))
    background_regions = pr.PyRanges(region_names_to_coordinates(background_region_names))
    return foreground_regions, background_regions
