"""Logo plots (pycistarget plot barplot / heatmap)."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pyranges as pr
import pytest

pytest.importorskip("logomaker")

from pycistarget.cli.pycistarget import main
from pycistarget.input_output import read_hdf5
from pycistarget.motif_enrichment_cistarget import cisTarget
from pycistarget.motif_enrichment_dem import DEM
from pycistarget.plotting import (
    LogoPlot,
    draw_logo,
    enrichment_scores,
    information_matrix,
    load_logo,
    plot_logo_barplot,
    plot_logo_heatmap,
    read_cluster_buster_matrices,
    read_enrichment_scores,
)

MOTIFS_DIR = Path(__file__).parent / "data" / "cbust"


def _cistarget_result(name: str, nes: dict[str, float]) -> cisTarget:
    result = cisTarget(region_set=pr.PyRanges(), name=name, species="homo_sapiens")
    result.motif_enrichment = pd.DataFrame(
        {"Region_set": name, "NES": list(nes.values()), "AUC": 0.1, "Rank_at_max": 10.0},
        index=pd.Index(list(nes), name="motifs"),
    )
    return result


def _dem_result(name: str, log2fc: dict[str, float]) -> DEM:
    result = DEM(
        foreground_regions=pr.PyRanges(),
        background_regions=pr.PyRanges(),
        name=name,
        species="homo_sapiens",
    )
    result.motif_enrichment = pd.DataFrame(
        {"Log2FC": list(log2fc.values()), "Adjusted_pval": 0.0},
        index=pd.Index(list(log2fc), name="motifs"),
    )
    return result


def test_read_cluster_buster_matrices():
    matrices = read_cluster_buster_matrices(MOTIFS_DIR / "metacluster_52.4.cb")
    assert [m.shape for m in matrices] == [(25, 4), (26, 4)]
    # WEIGHT / GAP lines and comments are skipped.
    assert len(read_cluster_buster_matrices(MOTIFS_DIR / "weighted.cb")) == 2


def test_information_matrix():
    info = information_matrix(np.array([[1.0, 1.0, 1.0, 1.0], [0.0, 10.0, 0.0, 0.0], [1, 1, 0, 0]]))
    assert list(info.columns) == list("ACGT")
    np.testing.assert_allclose(info.to_numpy().sum(axis=1), [0.0, 2.0, 1.0])
    np.testing.assert_allclose(info.iloc[1], [0.0, 2.0, 0.0, 0.0])


def test_load_logo(tmp_path, caplog):
    logo = load_logo(MOTIFS_DIR, "metacluster_52.4")
    assert logo is not None
    # The matrix with the highest information content of the metacluster.
    infos = [
        information_matrix(m)
        for m in read_cluster_buster_matrices(MOTIFS_DIR / "metacluster_52.4.cb")
    ]
    best = max(infos, key=lambda i: i.to_numpy().sum())
    pd.testing.assert_frame_equal(logo, best)

    (tmp_path / "bad.cb").write_text("1 2 3 4\n")
    with caplog.at_level(logging.WARNING):
        assert load_logo(MOTIFS_DIR, "missing_motif") is None
        assert load_logo(tmp_path, "bad") is None
    assert "No logo for motif missing_motif" in caplog.text
    assert "No logo for motif bad" in caplog.text


def test_read_enrichment_scores(tmp_path):
    path = tmp_path / "results.hdf5"
    _cistarget_result("a", {"m1": 5.0, "m2": 3.0}).write_hdf5(path)
    _cistarget_result("b", {"m2": 4.0}).write_hdf5(path, mode="a")
    scores = read_enrichment_scores(path)
    assert list(scores) == ["a", "b"]
    assert scores["a"][0] == "cisTarget"
    assert scores["a"][1].to_dict() == {"m1": 5.0, "m2": 3.0}

    tsv = tmp_path / "motif_enrichment_dem_x_vs_all.tsv"
    _dem_result("x_vs_all", {"m1": 1.5}).motif_enrichment.to_csv(tsv, sep="\t")
    assert {k: (t, s.to_dict()) for k, (t, s) in read_enrichment_scores(tsv).items()} == {
        "x_vs_all": ("DEM", {"m1": 1.5})
    }


def test_plot_barplot_command(tmp_path, caplog):
    path = tmp_path / "motif_enrichment_cistarget_a.hdf5"
    _cistarget_result(
        "a", {"jaspar__MA0130.1": 6.0, "metacluster_52.4": 4.0, "not_in_collection": 3.5}
    ).write_hdf5(path)
    output = tmp_path / "plots" / "a.png"
    with caplog.at_level(logging.WARNING):
        main(
            [
                "plot",
                "barplot",
                "--result",
                str(path),
                "--motifs_dir",
                str(MOTIFS_DIR),
                "--output",
                str(output),
            ]
        )
    assert output.stat().st_size > 0
    assert "No logo for motif not_in_collection" in caplog.text

    # A threshold that excludes all motifs.
    with pytest.raises(ValueError, match="No motifs to plot"):
        main(
            ["plot", "barplot", "--result", str(path), "--motifs_dir", str(MOTIFS_DIR),
             "--output", str(output), "--threshold", "10"]
        )  # fmt: skip


def test_plot_barplot_needs_name_for_several_results(tmp_path):
    path = tmp_path / "results.hdf5"
    _cistarget_result("a", {"weighted": 5.0}).write_hdf5(path)
    _cistarget_result("b", {"weighted": 4.0}).write_hdf5(path, mode="a")
    args = ["plot", "barplot", "--result", str(path), "--motifs_dir", str(MOTIFS_DIR)]
    with pytest.raises(ValueError, match="Several results to plot"):
        main([*args, "--output", str(tmp_path / "x.png")])
    main([*args, "--output", str(tmp_path / "b.pdf"), "--name", "b"])
    assert (tmp_path / "b.pdf").stat().st_size > 0


def test_plot_heatmap_command(tmp_path):
    paths = []
    for name, log2fc in {
        "a_vs_all": {"jaspar__MA0130.1": 2.0, "weighted": 1.2},
        "b_vs_all": {"metacluster_52.4": 1.5, "weighted": 1.1},
        "c_vs_all": {},  # no enriched motifs
    }.items():
        paths.append(tmp_path / f"motif_enrichment_dem_{name}.tsv")
        _dem_result(name, log2fc).motif_enrichment.to_csv(paths[-1], sep="\t")
    output = tmp_path / "heatmap.svg"
    main(
        ["plot", "heatmap", "--results", *map(str, paths), "--motifs_dir", str(MOTIFS_DIR),
         "--output", str(output), "--threshold", "1.15"]
    )  # fmt: skip
    assert output.stat().st_size > 0

    cistarget_path = tmp_path / "ctx.hdf5"
    _cistarget_result("ctx", {"weighted": 4.0}).write_hdf5(cistarget_path)
    with pytest.raises(ValueError, match="can not be plotted together"):
        main(
            ["plot", "heatmap", "--results", str(paths[0]), str(cistarget_path),
             "--motifs_dir", str(MOTIFS_DIR), "--output", str(output)]
        )  # fmt: skip


@pytest.fixture
def close_figures():
    import matplotlib.pyplot as plt

    yield
    plt.close("all")


def test_heatmap_motif_selection_and_order(close_figures):
    scores = pd.DataFrame(
        {
            "x": [5.0, np.nan, 3.0, 9.0, 3.1],
            "y": [6.0, 8.0, np.nan, np.nan, 3.2],
        },
        index=["m1", "m2", "m3", "m4", "m5"],
    )
    # Grouped by the region set with the highest score, sorted by that score.
    assert plot_logo_heatmap(scores, MOTIFS_DIR).motifs == ["m4", "m3", "m2", "m1", "m5"]
    assert plot_logo_heatmap(scores, MOTIFS_DIR, threshold=4.0).motifs == ["m4", "m2", "m1"]
    assert plot_logo_heatmap(scores, MOTIFS_DIR, top=1).motifs == ["m4", "m2"]


def test_enrichment_scores_inputs(tmp_path):
    a = _cistarget_result("a", {"m1": 5.0, "m2": 3.0})
    b = _cistarget_result("b", {"m2": 4.0})
    path = tmp_path / "results.hdf5"
    a.write_hdf5(path)
    b.write_hdf5(path, mode="a")
    expected = pd.DataFrame({"a": [5.0, 3.0], "b": [np.nan, 4.0]}, index=["m1", "m2"])

    for results in (path, str(path), read_hdf5(path), [a, b], {"a": a, "b": b}, [path]):
        table, score_name = enrichment_scores(results)
        assert score_name == "NES"
        pd.testing.assert_frame_equal(table, expected, check_names=False)

    # A column of a motif enrichment table: no region set name, the score name from the column.
    table, score_name = enrichment_scores(a.motif_enrichment["NES"])
    assert (list(table.columns), score_name) == ([""], "NES")
    # Named scores.
    table, score_name = enrichment_scores({"x": pd.Series({"m1": 1.0})})
    assert (list(table.columns), score_name) == (["x"], "Score")

    with pytest.raises(ValueError, match="Several results are named: a"):
        enrichment_scores([a, a])
    with pytest.raises(TypeError, match="Can not plot results of type int"):
        enrichment_scores(1)  # type: ignore[arg-type]


def test_plot_api(tmp_path, close_figures):
    result = _cistarget_result(
        "a", {"jaspar__MA0130.1": 6.0, "metacluster_52.4": 4.0, "weighted": 3.5}
    )
    plot = plot_logo_barplot(result, MOTIFS_DIR, top=2, color="black")
    assert isinstance(plot, LogoPlot)
    assert plot.motifs == ["jaspar__MA0130.1", "metacluster_52.4"]
    assert len(plot.logo_axes) == 2 and plot.colorbar is None
    assert plot.figure._suptitle.get_text() == "a"
    assert [bar.get_width() for bar in plot.ax.patches] == [6.0, 4.0]
    assert plot.ax.get_xlabel() == "NES"
    # No labels next to the logos.
    assert not any(label.get_text() for label in plot.ax.get_yticklabels())
    plot.save(tmp_path / "out" / "a.png")
    assert (tmp_path / "out" / "a.png").stat().st_size > 0

    assert plot_logo_barplot(result, MOTIFS_DIR, title="").figure._suptitle is None
    with pytest.raises(ValueError, match="No result named b"):
        plot_logo_barplot(result, MOTIFS_DIR, name="b")

    other = _cistarget_result("b", {"weighted": 5.0})
    plot = plot_logo_heatmap(
        [result, other], MOTIFS_DIR, cmap="viridis", vmin=0, vmax=10, title="topics"
    )
    assert plot.motifs == ["jaspar__MA0130.1", "metacluster_52.4", "weighted"]
    assert plot.colorbar is not None and plot.colorbar.ax.get_ylabel() == "NES"
    assert [t.get_text() for t in plot.ax.get_xticklabels()] == ["a", "b"]
    assert plot.ax.images[0].get_clim() == (0, 10)
    assert plot.ax.images[0].get_cmap().name.startswith("viridis")


def test_draw_logo(close_figures):
    import matplotlib.pyplot as plt

    _, ax = plt.subplots()
    logo = load_logo(MOTIFS_DIR, "jaspar__MA0130.1")
    draw_logo(ax, logo, width=30)
    assert ax.get_xlim() == (-0.5, 29.5)
    assert not ax.axison
    draw_logo(ax, None)  # missing logo: nothing drawn
