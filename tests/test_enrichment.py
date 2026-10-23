import h5py
import pandas as pd
import pyranges as pr
import pytest

from pycistarget.cli.pycistarget import main
from pycistarget.input_output import read_hdf5
from pycistarget.motif_enrichment_cistarget import cisTarget, cisTargetDatabase
from pycistarget.motif_enrichment_dem import DEM, DEMDatabase, get_foreground_and_background_regions
from pycistarget.utils import region_names_to_coordinates

from .conftest import ENRICHED_MOTIFS

ANNOTATIONS = ["Direct_annot", "Motif_similarity_annot", "Orthology_annot"]


def _pyranges(regions):
    return pr.PyRanges(region_names_to_coordinates(regions))


@pytest.fixture(scope="module")
def cistarget_result(synthetic):
    regions = pr.read_bed(str(synthetic.foreground_bed))
    ctx_db = cisTargetDatabase(synthetic.rankings_db, region_sets=regions)
    result = cisTarget(
        region_set=regions,
        name="fg",
        species="homo_sapiens",
        auc_threshold=0.05,
        rank_threshold=0.2,
        path_to_motif_annotations=str(synthetic.annotation),
        annotation_to_use=ANNOTATIONS,
    )
    result.run_ctx(ctx_db)
    return result


@pytest.fixture(scope="module")
def dem_result(synthetic):
    foreground, background = get_foreground_and_background_regions(
        [pr.read_bed(str(synthetic.foreground_bed))],
        [pr.read_bed(str(synthetic.background_bed))],
        balance_number_of_promoters=False,
    )
    result = DEM(
        foreground_regions=foreground,
        background_regions=background,
        name="fg_vs_bg",
        species="homo_sapiens",
        path_to_motif_annotations=str(synthetic.annotation),
        annotation_to_use=ANNOTATIONS,
    )
    result.run(DEMDatabase(synthetic.scores_db))
    return result


def test_cistarget_finds_planted_motifs(cistarget_result, synthetic):
    me = cistarget_result.motif_enrichment
    assert set(me.index) == set(ENRICHED_MOTIFS)
    assert list(me.columns) == [
        "Logo", "Region_set", "NES", "AUC", "Rank_at_max", *ANNOTATIONS, "Motif_hits",
    ]  # fmt: skip
    assert me["NES"].is_monotonic_decreasing
    assert me.loc["motif_1", "Direct_annot"] == "TFA"
    for motif, hits in cistarget_result.motif_hits["Database"].items():
        assert me.loc[motif, "Motif_hits"] == len(hits)
        assert set(hits) <= set(synthetic.foreground)
    assert "TFA_(200r)" in cistarget_result.cistromes["Region_set"]


def test_cistarget_without_enriched_motifs(synthetic):
    regions = _pyranges(synthetic.background[:200])
    result = cisTarget(
        region_set=regions,
        name="bg",
        species="homo_sapiens",
        nes_threshold=100.0,
        path_to_motif_annotations=str(synthetic.annotation),
    )
    result.run_ctx(cisTargetDatabase(synthetic.rankings_db, region_sets=regions))
    assert len(result.motif_enrichment) == 0
    assert result.cistromes == {"Database": {}, "Region_set": {}}


def test_dem_finds_planted_motifs(dem_result):
    me = dem_result.motif_enrichment
    assert set(me.index) == set(ENRICHED_MOTIFS)
    assert (me["Log2FC"] >= 1.0).all() and (me["Adjusted_pval"] <= 0.05).all()
    # Motif hit thresholds and counts are part of the result table.
    for motif in me.index:
        assert me.loc[motif, "Motif_hits"] == len(dem_result.motif_hits["Database"][motif])
        assert me.loc[motif, "Motif_hit_thr"] > 0


def test_max_bg_regions_is_respected(synthetic):
    _, background = get_foreground_and_background_regions(
        [_pyranges(synthetic.foreground)],
        [_pyranges(synthetic.background)],
        max_bg_regions=100,
        balance_number_of_promoters=False,
    )
    assert len(background) == 100


def test_balancing_promoters_requires_genome_annotation(synthetic):
    with pytest.raises(ValueError, match="genome_annotation"):
        get_foreground_and_background_regions(
            [_pyranges(synthetic.foreground)], [_pyranges(synthetic.background)]
        )


@pytest.mark.parametrize("result_fixture", ["cistarget_result", "dem_result"])
def test_hdf5_round_trip(result_fixture, request, tmp_path):
    result = request.getfixturevalue(result_fixture)
    path = tmp_path / "result.hdf5"
    result.write_hdf5(path)
    loaded = read_hdf5(path)[result.name]
    assert type(loaded) is type(result)
    assert loaded.species == "homo_sapiens"
    assert loaded.annotation_to_use == ANNOTATIONS
    pd.testing.assert_frame_equal(
        loaded.motif_enrichment, result.motif_enrichment, check_dtype=False, check_index_type=False
    )
    # Keys in the file are lowercase (relied upon by SCENIC+).
    assert loaded.motif_hits["region_set"] == result.motif_hits["Region_set"]
    assert loaded.cistromes["database"] == result.cistromes["Database"]
    # Several results can be stored in one file.
    renamed = request.getfixturevalue(result_fixture)
    renamed.name, original_name = "other", renamed.name
    try:
        renamed.write_hdf5(path, mode="a")
    finally:
        renamed.name = original_name
    assert set(read_hdf5(path)) == {result.name, "other"}


def test_cli(synthetic, tmp_path):
    common = ["--species", "homo_sapiens", "--path_to_motif_annotations", str(synthetic.annotation)]
    main(
        [
            "cistarget",
            "--cistarget_db_fname", str(synthetic.rankings_db),
            "--bed_fname", str(synthetic.foreground_bed),
            "--output_folder", str(tmp_path / "ctx"),
            "--auc_threshold", "0.05",
            "--rank_threshold", "0.2",
            "--write_html",
            *common,
        ]
    )  # fmt: skip
    assert (tmp_path / "ctx" / "motif_enrichment_cistarget_foreground.hdf5").exists()
    assert (tmp_path / "ctx" / "motif_enrichment_cistarget_foreground.html").exists()
    # dem without --genome_annotation, into a folder that does not exist yet.
    main(
        [
            "dem",
            "--dem_db_fname", str(synthetic.scores_db),
            "--foreground_beds", str(synthetic.foreground_bed),
            "--background_beds", str(synthetic.background_bed),
            "--output_folder", str(tmp_path / "dem"),
            "--output_mode", "tsv",
            *common,
        ]
    )  # fmt: skip
    tsv = pd.read_table(
        tmp_path / "dem" / "motif_enrichment_dem_foreground__VS__background.tsv", index_col=0
    )
    assert set(tsv.index) == set(ENRICHED_MOTIFS)
    with h5py.File(tmp_path / "ctx" / "motif_enrichment_cistarget_foreground.hdf5") as h5:
        assert h5["foreground/metadata/type"][()] == b"cisTarget"


def _assert_same_result(a, b):
    pd.testing.assert_frame_equal(a.motif_enrichment, b.motif_enrichment)
    assert a.motif_hits == b.motif_hits
    assert a.cistromes == b.cistromes


def test_cli_batch_modes_equal_single_runs(synthetic, tmp_path):
    common = ["--species", "homo_sapiens", "--path_to_motif_annotations", str(synthetic.annotation)]
    beds = [tmp_path / "set_a.bed", tmp_path / "set_b.bed", tmp_path / "set_c.bed"]
    regions = [
        synthetic.foreground[:150] + synthetic.background[:100],
        synthetic.background[100:900],
    ]
    regions.append(synthetic.foreground[100:] + synthetic.background[900:])
    for bed, bed_regions in zip(beds, regions, strict=True):
        bed.write_text("".join(r.replace(":", "\t").replace("-", "\t") + "\n" for r in bed_regions))
    ctx = ["cistarget", "--cistarget_db_fname", str(synthetic.rankings_db), "--auc_threshold", "0.05",
           "--rank_threshold", "0.2", *common]  # fmt: skip
    dem = ["dem", "--dem_db_fname", str(synthetic.scores_db), *common]
    main([*ctx, "--bed_fname", *map(str, beds), "--output_folder", str(tmp_path / "batch")])
    main([*dem, "--one_vs_all_beds", *map(str, beds), "--output_folder", str(tmp_path / "batch")])
    for i, bed in enumerate(beds):
        name = bed.stem
        main([*ctx, "--bed_fname", str(bed), "--output_folder", str(tmp_path / "single")])
        background = [str(b) for b in beds if b != bed]
        main([*dem, "--foreground_beds", str(bed), "--background_beds", *background,
              "--name", f"{name}_vs_all", "--output_folder", str(tmp_path / "single")])  # fmt: skip
        for analysis, result_name in [("cistarget", name), ("dem", f"{name}_vs_all")]:
            fname = f"motif_enrichment_{analysis}_{result_name}.hdf5"
            batch = read_hdf5(tmp_path / "batch" / fname)[result_name]
            single = read_hdf5(tmp_path / "single" / fname)[result_name]
            _assert_same_result(batch, single)
            if analysis == "dem" and i == 0:
                assert len(batch.motif_enrichment) > 0


def test_cli_rejects_invalid_batch_arguments(synthetic, tmp_path):
    dem = ["dem", "--dem_db_fname", str(synthetic.scores_db), "--species", "homo_sapiens",
           "--output_folder", str(tmp_path)]  # fmt: skip
    for extra in [
        ["--one_vs_all_beds", str(synthetic.foreground_bed), str(synthetic.background_bed),
         "--foreground_beds", str(synthetic.foreground_bed)],
        ["--one_vs_all_beds", str(synthetic.foreground_bed)],
        ["--foreground_beds", str(synthetic.foreground_bed)],
    ]:  # fmt: skip
        with pytest.raises(SystemExit):
            main([*dem, *extra])
    with pytest.raises(ValueError, match="single bed file"):
        main(["cistarget", "--cistarget_db_fname", str(synthetic.rankings_db), "--species",
              "homo_sapiens", "--output_folder", str(tmp_path), "--name", "x",
              "--bed_fname", str(synthetic.foreground_bed), str(synthetic.background_bed)])  # fmt: skip


def test_force_runs_without_motif_annotations(synthetic, cistarget_result, tmp_path, caplog):
    """With --force, a missing motif-to-TF annotation gives results without annotations/cistromes."""
    missing = tmp_path / "missing_annotations.tbl"
    common = ["--species", "homo_sapiens", "--path_to_motif_annotations", str(missing)]
    ctx = [
        "cistarget", "--cistarget_db_fname", str(synthetic.rankings_db),
        "--auc_threshold", "0.05", "--rank_threshold", "0.2", *common,
    ]  # fmt: skip
    beds = [synthetic.foreground_bed, synthetic.background_bed]
    with pytest.raises(RuntimeError, match="Unable to load motif-to-TF annotations"):
        main([*ctx, "--bed_fname", str(beds[0]), "--output_folder", str(tmp_path / "fail")])

    with caplog.at_level("WARNING"):
        main(
            [*ctx, "--bed_fname", *map(str, beds), "--output_folder", str(tmp_path / "ctx"),
             "--write_html", "--force"]
        )  # fmt: skip
        main(
            ["dem", "--dem_db_fname", str(synthetic.scores_db),
             "--foreground_beds", str(beds[0]), "--background_beds", str(beds[1]),
             "--output_folder", str(tmp_path / "dem"), "--name", "fg_vs_bg", "--force", *common]
        )  # fmt: skip
    # One warning, the annotations are not loaded again for each region set.
    assert caplog.text.count("Continuing without motif-to-TF annotations") == 1

    result = read_hdf5(tmp_path / "ctx" / "motif_enrichment_cistarget_foreground.hdf5")[
        "foreground"
    ]
    expected = cistarget_result.motif_enrichment
    pd.testing.assert_series_equal(
        result.motif_enrichment["NES"], expected["NES"], check_names=False, check_index_type=False
    )
    assert result.motif_enrichment[ANNOTATIONS[0]].isna().all()
    # read_hdf5 returns lower case keys.
    assert result.cistromes == {"database": {}, "region_set": {}}
    assert result.motif_hits["region_set"] == cistarget_result.motif_hits["Region_set"]
    dem = read_hdf5(tmp_path / "dem" / "motif_enrichment_dem_fg_vs_bg.hdf5")["fg_vs_bg"]
    assert set(dem.motif_enrichment.index) == set(ENRICHED_MOTIFS)
    assert dem.cistromes == {"database": {}, "region_set": {}}


def test_force_without_annotations_for_species(synthetic):
    regions = pr.read_bed(str(synthetic.foreground_bed))
    result = cisTarget(
        region_set=regions,
        name="fg",
        species="danio_rerio",  # type: ignore[arg-type]
        auc_threshold=0.05,
        rank_threshold=0.2,
        force=True,
    )
    result.run_ctx(cisTargetDatabase(synthetic.rankings_db, region_sets=regions))
    assert len(result.motif_enrichment) > 0
    assert result.motif_enrichment["Direct_annot"].isna().all()
