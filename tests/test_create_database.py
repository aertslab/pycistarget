"""Database creation: scores must be identical to Cluster-Buster's (expected outputs in data/cbust)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.feather as pf
import pytest

from pycistarget._cbust import Sequences, score_motifs
from pycistarget.cli.pycistarget import main
from pycistarget.create_database import rankings_from_scores, read_motif_list
from pycistarget.database import FeatherCisTargetDatabase

DATA = Path(__file__).parent / "data" / "cbust"
MOTIFS = ["jaspar__MA0130.1", "metacluster_52.4", "weighted"]
BG_PADDING = 1000


def _expected(motif: str, bg_range: int, mask: bool, names: list[str]) -> np.ndarray:
    """Cluster-Buster's output, read like create_cisTarget_databases does (missing regions: 0)."""
    path = DATA / f"expected.{motif}.r{bg_range}.{'mask' if mask else 'nomask'}.txt"
    df = pd.read_csv(
        path,
        sep="\t",
        header=0,
        names=["seq_name", "crm_score"],
        index_col="seq_name",
        usecols=["seq_name", "crm_score"],
        dtype={"seq_name": str, "crm_score": np.float32},
        engine="c",
    )
    return df["crm_score"].reindex(names, fill_value=0).to_numpy(np.float32)


@pytest.mark.parametrize("bg_range, mask", [(10_000, True), (1_500, True), (10_000, False)])
@pytest.mark.parametrize("n_threads", [1, 3])
def test_scores_identical_to_cluster_buster(bg_range, mask, n_threads):
    sequences = Sequences(str(DATA / "regions.fa"), mask=mask)
    scores = score_motifs(
        sequences,
        [str(DATA / f"{m}.cb") for m in MOTIFS],
        bg_padding=BG_PADDING,
        bg_range=bg_range,
        n_threads=n_threads,
    )
    assert scores.dtype == np.float32
    assert scores.shape == (len(MOTIFS), len(sequences))
    for motif, row in zip(MOTIFS, scores, strict=True):
        np.testing.assert_array_equal(row, _expected(motif, bg_range, mask, sequences.names))


def test_sequences():
    sequences = Sequences(str(DATA / "regions.fa"))
    assert len(sequences) == 43
    assert sequences.names[-4:] == ["long", "odd_a", "odd_b", "odd_c"]
    assert sequences.lengths[-4:] == [12_000, 2_300, 2_300, 2_700]


def test_score_errors(tmp_path):
    sequences = Sequences(str(DATA / "regions.fa"))
    with pytest.raises(ValueError, match="should be at least 2401 bp long"):
        score_motifs(sequences, [str(DATA / "weighted.cb")], bg_padding=1200)
    with pytest.raises(OSError, match="missing.cb"):
        score_motifs(sequences, [str(tmp_path / "missing.cb")])
    (tmp_path / "bad.cb").write_text(">bad\n1 2 3\nx\n")
    with pytest.raises(ValueError, match="Could not read motif file"):
        score_motifs(sequences, [str(tmp_path / "bad.cb")])


def _rankings_create_cistarget_databases(scores: np.ndarray, seed: int) -> np.ndarray:
    """The ranking code of create_cisTarget_databases (CisTargetDatabase.convert_scores_db_to_rankings_db)."""
    rng = np.random.default_rng(seed=seed)
    dtype = np.int16 if scores.shape[1] <= 2**15 else np.int32
    rankings = np.empty(scores.shape, dtype=dtype)
    for i, motif_scores in enumerate(scores):
        permutation = rng.permutation(motif_scores.shape[0])
        ranking = np.empty(motif_scores.shape[0], dtype=dtype)
        ranking[permutation[(-motif_scores)[permutation].argsort()]] = np.arange(
            motif_scores.shape[0], dtype=dtype
        )
        rankings[i] = ranking
    return rankings


@pytest.mark.parametrize("n_regions", [1_000, 40_000])
def test_rankings_like_create_cistarget_databases(n_regions):
    rng = np.random.default_rng(1)
    # Many ties, as in real databases (scores have 3 significant digits and many are 0).
    scores = (
        np.round(rng.exponential(1.0, size=(30, n_regions)) - 1.0, 1).clip(0).astype(np.float32)
    )
    expected = _rankings_create_cistarget_databases(scores, seed=666)
    rankings = rankings_from_scores(scores, seed=666, n_threads=4)
    assert rankings.dtype == expected.dtype == (np.int16 if n_regions <= 2**15 else np.int32)
    np.testing.assert_array_equal(rankings, expected)
    np.testing.assert_array_equal(
        rankings_from_scores(scores, seed=666, n_threads=2, transposed=True), expected.T
    )
    assert not np.array_equal(rankings_from_scores(scores, seed=667), expected)


def test_read_motif_list(tmp_path):
    (tmp_path / "motifs.txt").write_text("# comment\nweighted.cb\n\njaspar__MA0130.1\n")
    assert read_motif_list(DATA, tmp_path / "motifs.txt") == {
        "weighted": DATA / "weighted.cb",
        "jaspar__MA0130.1": DATA / "jaspar__MA0130.1.cb",
    }
    (tmp_path / "missing.txt").write_text("weighted\nmissing\n")
    with pytest.raises(FileNotFoundError, match="missing"):
        read_motif_list(DATA, tmp_path / "missing.txt")


def test_create_db_command(tmp_path):
    (tmp_path / "motifs.txt").write_text("\n".join(reversed(MOTIFS)) + "\n")
    prefix = tmp_path / "db"
    main(
        [
            "create_db",
            "--fasta",
            str(DATA / "regions.fa"),
            "--motifs_dir",
            str(DATA),
            "--motifs",
            str(tmp_path / "motifs.txt"),
            "--output_prefix",
            str(prefix),
            "--bg_padding",
            str(BG_PADDING),
            "--mask",
            "--seed",
            "123",
            "--threads",
            "2",
        ]
    )
    regions = sorted(Sequences(str(DATA / "regions.fa")).names)
    expected = np.stack([_expected(m, 10_000, True, regions) for m in MOTIFS])

    table = pf.read_table(f"{prefix}.motifs_vs_regions.scores.feather")
    assert table.column_names == [*MOTIFS, "regions"]
    assert table.column("regions").to_pylist() == regions
    np.testing.assert_array_equal(np.stack([table.column(m).to_numpy() for m in MOTIFS]), expected)

    scores_db = FeatherCisTargetDatabase(f"{prefix}.regions_vs_motifs.scores.feather")
    assert list(scores_db.region_ids) == regions
    scores, motifs, _ = scores_db.load_array(regions)
    assert list(motifs) == MOTIFS
    np.testing.assert_array_equal(scores, expected)

    rankings_db = FeatherCisTargetDatabase(f"{prefix}.regions_vs_motifs.rankings.feather")
    rankings, motifs, _ = rankings_db.load_array(regions)
    assert list(motifs) == MOTIFS
    np.testing.assert_array_equal(rankings, _rankings_create_cistarget_databases(expected, 123))


def test_create_db_rejects_duplicate_regions(tmp_path):
    fasta = tmp_path / "dup.fa"
    fasta.write_text(">a\n" + "A" * 50 + "\n>a\n" + "C" * 50 + "\n")
    (tmp_path / "motifs.txt").write_text("weighted\n")
    with pytest.raises(ValueError, match="not unique"):
        main(
            [
                "create_db",
                "--fasta",
                str(fasta),
                "--motifs_dir",
                str(DATA),
                "--motifs",
                str(tmp_path / "motifs.txt"),
                "--output_prefix",
                str(tmp_path / "db"),
            ]
        )
