import numpy as np
import pandas as pd
import pyranges as pr
import pytest

from pycistarget.database import FeatherCisTargetDatabase
from pycistarget.motif_enrichment_dem import (
    _dem_statistics,
    _log2_fc_from_means,
    _optimal_thresholds_roc,
    get_log2_fc,
    get_optimal_threshold_roc,
    mean_axis1,
    p_adjust_bh,
    ranksums_numba_multiple,
)
from pycistarget.utils import (
    coord_to_region_names,
    get_cistromes_per_region_set,
    get_motifs_per_TF,
    get_TF_list,
    load_motif_annotations,
    region_names_to_coordinates,
    region_set_motif_hits,
    target_to_query,
)


def test_region_name_round_trip():
    names = ["chr1:100-600", "chrX:5-10"]
    coords = region_names_to_coordinates(names)
    assert coords.to_dict("list") == {
        "Chromosome": ["chr1", "chrX"],
        "Start": [100, 5],
        "End": [600, 10],
    }
    assert list(coords.index) == names
    assert coord_to_region_names(coords) == names
    assert coord_to_region_names(pr.PyRanges(coords)) == names


def test_target_to_query_fraction_overlap():
    target = pr.PyRanges(region_names_to_coordinates(["chr1:100-200", "chr1:1000-1100"]))
    query = ["chr1:150-250", "chr1:1090-1190", "chr2:100-200"]
    mapping = target_to_query(target, query, fraction_overlap=0.4)
    # 50% overlap passes, 10% overlap does not.
    assert mapping.to_dict("list") == {"Target": ["chr1:100-200"], "Query": ["chr1:150-250"]}
    with pytest.raises(ValueError, match="None of the regions"):
        target_to_query(target, ["chr3:1-10"])


def test_database_load_in_database_order(synthetic):
    db = FeatherCisTargetDatabase(synthetic.rankings_db)
    assert db.total_regions == len(synthetic.regions)
    # regions[1] (chr2) comes before regions[3] (chr1) in the database.
    requested = [synthetic.regions[3], synthetic.regions[1], "chr9:1-2", synthetic.regions[3]]
    expected_columns = [synthetic.regions[1], synthetic.regions[3]]
    df = db.load(requested)
    assert list(df.columns) == expected_columns
    assert list(df.index) == synthetic.motifs
    np.testing.assert_array_equal(df[synthetic.regions[1]].to_numpy(), synthetic.rankings[:, 1])
    values, motifs, regions = db.load_array(requested, chunk_size=1)
    assert values.flags.c_contiguous
    assert regions == expected_columns
    assert list(motifs) == synthetic.motifs
    np.testing.assert_array_equal(values, df.to_numpy())


def test_load_motif_annotations(synthetic):
    annot = load_motif_annotations("homo_sapiens", fname=str(synthetic.annotation))
    # TF lists are sorted, so results do not depend on hash randomisation.
    assert annot.loc["motif_0", "Direct_annot"] == "TFA, TFB"
    assert annot.loc["motif_2", "Motif_similarity_annot"] == "TFC"
    assert annot.loc["motif_3", "Orthology_annot"] == "TFD"
    strict = load_motif_annotations(
        "homo_sapiens", fname=str(synthetic.annotation), orthologous_identity_threshold=0.9
    )
    assert "motif_3" not in strict.index


def test_cistromes():
    motif_enrichment = pd.DataFrame(
        {
            "Direct_annot": ["TFA, TFB", "TFA", np.nan],
            "Motif_similarity_annot": [np.nan, "TFC", "TFA"],
        },
        index=["m1", "m2", "m3"],
    )
    motif_hits = {"m1": ["r1", "r2"], "m2": ["r2", "r3"], "m3": ["r4"]}
    assert get_TF_list(motif_enrichment) == ["TFA", "TFB", "TFC"]
    assert get_motifs_per_TF(motif_enrichment, "TFA", "Index", ["Direct_annot"]) == ["m1", "m2"]
    cistromes = get_cistromes_per_region_set(
        motif_enrichment, motif_hits, ["Direct_annot", "Motif_similarity_annot"]
    )
    assert cistromes == {
        "TFA_(3r)": ["r1", "r2", "r3"],
        "TFB_(2r)": ["r1", "r2"],
        "TFA_extended_(4r)": ["r1", "r2", "r3", "r4"],
        "TFB_extended_(2r)": ["r1", "r2"],
        "TFC_extended_(2r)": ["r2", "r3"],
    }


def test_mean_axis1_accumulates_in_float64():
    rng = np.random.default_rng(0)
    x = np.round(rng.gamma(1.0, 2.0, size=(20, 50_001)), 1).astype(np.float32)
    expected = x.astype(np.float64).sum(axis=1) / x.shape[1]
    np.testing.assert_array_equal(mean_axis1(x), expected)
    np.testing.assert_array_equal(mean_axis1(np.asfortranarray(x)), expected)


def test_p_adjust_bh():
    from scipy.stats import false_discovery_control

    p = np.random.default_rng(0).random(1_000) ** 3
    np.testing.assert_allclose(p_adjust_bh(p), false_discovery_control(p, method="bh"), rtol=1e-12)


def test_motifs_per_tf_matches_whole_tf_names():
    # "YOD1" is a suffix of "MYOD1" and "SP1" of "MESP1": these should not match.
    motif_enrichment = pd.DataFrame(
        {"Direct_annot": ["MYOD1, TCF3", "YOD1", "MESP1", "SP1, SP2", "E(spl)m3-HLH"]},
        index=["m1", "m2", "m3", "m4", "m5"],
    )
    assert get_motifs_per_TF(motif_enrichment, "YOD1", "Index", ["Direct_annot"]) == ["m2"]
    assert get_motifs_per_TF(motif_enrichment, "SP1", "Index", ["Direct_annot"]) == ["m4"]
    assert get_motifs_per_TF(motif_enrichment, "E(spl)m3-HLH", "Index", ["Direct_annot"]) == ["m5"]


N_DEM_TEST = 20_000  # more than _HASH_MAX_DISTINCT distinct values for the sorting path


def _dem_test_matrix():
    rng = np.random.default_rng(1)
    # Few distinct values (many ties, counting path), as cisTarget scores.
    few = np.round(rng.gamma(1.0, 2.0, size=(6, N_DEM_TEST)), 1)
    few[rng.random(few.shape) < 0.4] = 0.0
    few[0, :10] = -0.0
    few[1, :1_000] += 3.0
    # More distinct values than the hash table holds (sorting path).
    many = rng.gamma(1.0, 2.0, size=(2, N_DEM_TEST))
    nan = few[:1].copy()
    nan[0, 5] = np.nan
    return np.vstack([few, many, nan]).astype(np.float32)


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_dem_statistics_identical_to_full_matrix_functions(dtype):
    scores = _dem_test_matrix().astype(dtype)
    fg_idx = np.arange(0, N_DEM_TEST, 3, dtype=np.int64)
    bg_idx = np.setdiff1d(np.arange(N_DEM_TEST), fg_idx)
    motif_idx = np.arange(scores.shape[0])
    mean_fg, mean_bg, pvalues = _dem_statistics(scores, fg_idx, bg_idx, motif_idx)
    fg, bg = scores[:, fg_idx], scores[:, bg_idx]
    np.testing.assert_array_equal(mean_fg, mean_axis1(fg))
    np.testing.assert_array_equal(mean_bg, mean_axis1(bg))
    np.testing.assert_array_equal(pvalues, ranksums_numba_multiple(fg, bg)[1])
    np.testing.assert_array_equal(_log2_fc_from_means(mean_fg, mean_bg), get_log2_fc(fg, bg))


def test_optimal_thresholds_identical_to_sklearn():
    scores = _dem_test_matrix()[:-1]  # roc_curve does not accept NaN
    fg_idx = np.arange(0, N_DEM_TEST, 3, dtype=np.int64)
    bg_idx = np.setdiff1d(np.arange(N_DEM_TEST), fg_idx)
    thresholds = _optimal_thresholds_roc(scores, fg_idx, bg_idx, np.arange(scores.shape[0]))
    expected = [get_optimal_threshold_roc(row[fg_idx], row[bg_idx]) for row in scores]
    np.testing.assert_array_equal(thresholds, expected)
    # No scores above zero in the background.
    row = np.zeros((1, N_DEM_TEST), dtype=np.float32)
    row[0, fg_idx] = 1.0
    assert np.isnan(_optimal_thresholds_roc(row, fg_idx, bg_idx, np.array([0]))[0])


def test_region_set_motif_hits():
    regions_to_db = pd.DataFrame(
        {"Target": ["t1", "t2", "t3", "t3"], "Query": ["q1", "q1", "q2", "q3"]}
    )
    hits = {"m1": ["q1", "q3"], "m2": ["q2"], "m3": []}
    assert region_set_motif_hits(regions_to_db, hits) == {
        "m1": ["t1", "t2", "t3"],
        "m2": ["t3"],
        "m3": [],
    }
