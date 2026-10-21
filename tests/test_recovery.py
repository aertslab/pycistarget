import numpy as np
import pytest

from pycistarget.recovery import calculate_aucs, derive_rank_cutoff, leading_edge, recovery_curves


def reference_auc(ranking: np.ndarray, rank_cutoff: int) -> float:
    """AUC as the area under the step-wise recovery curve (R cisTarget definition)."""
    x = np.concatenate([np.sort(ranking[ranking < rank_cutoff]), [rank_cutoff]])
    y = np.arange(1, x.size, dtype=np.float64)
    return float(np.sum(np.diff(x) * y))


@pytest.fixture
def rankings():
    rng = np.random.default_rng(1)
    total = 10_000
    return np.stack([rng.permutation(total)[:300] for _ in range(50)]).astype(np.int32), total


def test_derive_rank_cutoff():
    assert derive_rank_cutoff(0.005, 488_108) == round(0.005 * 488_108) - 1
    with pytest.raises(ValueError):
        derive_rank_cutoff(0.0, 100)
    with pytest.raises(ValueError):
        derive_rank_cutoff(0.001, 100)


def test_aucs_match_reference(rankings):
    ranks, total = rankings
    auc_threshold = 0.05
    rank_cutoff = derive_rank_cutoff(auc_threshold, total)
    max_auc = (rank_cutoff + 1) * ranks.shape[1]
    expected = np.array([reference_auc(r, rank_cutoff) for r in ranks]) / max_auc
    np.testing.assert_array_equal(calculate_aucs(ranks, total, auc_threshold), expected)


def test_recovery_curves_match_reference(rankings):
    ranks, _ = rankings
    rank_threshold = 500
    rccs = recovery_curves(ranks, rank_threshold)
    expected = np.stack([(r[:, None] <= np.arange(rank_threshold)).sum(axis=0) for r in ranks])
    np.testing.assert_array_equal(rccs, expected.astype(np.float64))
    # Works for non C-contiguous input as well.
    np.testing.assert_array_equal(recovery_curves(np.asfortranarray(ranks), rank_threshold), rccs)


def test_leading_edge():
    regions = np.array(["a", "b", "c", "d"], dtype=object)
    ranking = np.array([5, 0, 9, 2])
    rcc = np.array([1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0])
    avg2stdrcc = np.array([0.5, 0.5, 0.5, 2.0, 2.0, 2.0, 2.0])
    hits, rank_at_max = leading_edge(rcc, avg2stdrcc, ranking, regions)
    assert rank_at_max == 2
    assert hits == ["b", "d"]
