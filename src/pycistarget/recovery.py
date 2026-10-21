"""
Recovery curve, AUC and leading edge calculations for cisTarget.

All region sets are unweighted, so recovery curves are cumulative hit counts and
the AUC numerator is an exact integer. This keeps results identical to the
(R)cisTarget / ctxcore definitions:

- rank cutoff for the AUC: ``round(auc_threshold * total_regions) - 1`` (0-based ranks,
  ranks equal to the cutoff are excluded);
- AUC: area under the recovery curve up to the cutoff, divided by
  ``(rank_cutoff + 1) * n_regions``;
- leading edge: regions ranked at or before the rank at which the recovery curve
  exceeds the mean + 2 sd of all recovery curves the most.
"""

from __future__ import annotations

import numba
import numpy as np


def derive_rank_cutoff(auc_threshold: float, total_regions: int) -> int:
    if not 0.0 < auc_threshold <= 1.0:
        raise ValueError("AUC threshold must be a fraction between 0.0 and 1.0.")
    rank_cutoff = int(round(auc_threshold * total_regions))
    if not 0 < rank_cutoff < total_regions:
        raise ValueError(
            f"An AUC threshold of {auc_threshold:f} corresponds to {rank_cutoff:d} top ranked "
            "regions in the database. Please increase the AUC threshold."
        )
    # The region at the cutoff itself is not included in the AUC (same as R cisTarget).
    return rank_cutoff - 1


def calculate_aucs(rankings: np.ndarray, total_regions: int, auc_threshold: float) -> np.ndarray:
    """
    Calculate the normalised AUC of the recovery curve of each row of ``rankings``.

    Parameters
    ----------
    rankings
        Array of shape (n_motifs, n_regions) with the 0-based rank of each region for each motif.
    total_regions
        Total number of regions in the database.
    auc_threshold
        Fraction of the ranked database to take into account.
    """
    rank_cutoff = derive_rank_cutoff(auc_threshold, total_regions)
    max_auc = float((rank_cutoff + 1) * rankings.shape[1])
    return _auc_areas(rankings, rank_cutoff) / max_auc


@numba.njit(parallel=True, cache=True)
def _auc_areas(rankings: np.ndarray, rank_cutoff: int) -> np.ndarray:
    # Each region ranked at r < rank_cutoff contributes (rank_cutoff - r) to the area under
    # the (unnormalised) recovery curve. Summed as integers, so the result is exact.
    areas = np.empty(rankings.shape[0], dtype=np.float64)
    for i in numba.prange(rankings.shape[0]):
        area = 0
        for j in range(rankings.shape[1]):
            if rankings[i, j] < rank_cutoff:
                area += rank_cutoff - np.int64(rankings[i, j])
        areas[i] = area
    return areas


def recovery_curves(rankings: np.ndarray, rank_threshold: int) -> np.ndarray:
    """
    Calculate recovery curves: the number of regions recovered at each rank below ``rank_threshold``.

    Returns a float64 array of shape (n_motifs, rank_threshold).
    """
    rccs = np.zeros((rankings.shape[0], rank_threshold), dtype=np.float64)
    _fill_recovery_curves(rankings, rccs)
    return rccs


@numba.njit(parallel=True, cache=True)
def _fill_recovery_curves(rankings: np.ndarray, rccs: np.ndarray):
    rank_threshold = rccs.shape[1]
    for i in numba.prange(rankings.shape[0]):
        for j in range(rankings.shape[1]):
            if rankings[i, j] < rank_threshold:
                rccs[i, rankings[i, j]] += 1.0
        for r in range(1, rank_threshold):
            rccs[i, r] += rccs[i, r - 1]


def leading_edge(
    rcc: np.ndarray, avg2stdrcc: np.ndarray, ranking: np.ndarray, regions: np.ndarray
) -> tuple[list[str], int]:
    """
    Get the regions in the leading edge of a recovery curve.

    Returns the regions ranked at or before the rank at maximum difference between
    ``rcc`` and ``avg2stdrcc`` (sorted by rank), and that rank.
    """
    rank_at_max = int(np.argmax(rcc - avg2stdrcc))
    order = np.argsort(ranking)
    in_leading_edge = ranking[order] <= rank_at_max
    return regions[order][in_leading_edge].tolist(), rank_at_max
