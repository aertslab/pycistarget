"""Small synthetic cisTarget databases, region sets and motif annotations for tests."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.feather as pf
import pytest

N_REGIONS = 2_000
N_MOTIFS = 40
REGION_LENGTH = 500
# Motifs (rows) for which the foreground regions get high scores.
ENRICHED_MOTIFS = ["motif_0", "motif_1", "motif_2"]


@dataclass
class SyntheticData:
    rankings_db: Path
    scores_db: Path
    foreground_bed: Path
    background_bed: Path
    annotation: Path
    regions: list[str]
    motifs: list[str]
    scores: np.ndarray  # motifs x regions
    rankings: np.ndarray  # motifs x regions
    foreground: list[str]
    background: list[str]


def _write_db(path: Path, values: np.ndarray, regions: list[str], motifs: list[str]):
    table = pa.table({**{r: values[:, i] for i, r in enumerate(regions)}, "motifs": motifs})
    pf.write_feather(table, path, compression="lz4")


def _write_bed(path: Path, regions: list[str]):
    rows = [r.replace(":", "\t").replace("-", "\t") for r in regions]
    path.write_text("\n".join(rows) + "\n")


@pytest.fixture(scope="session")
def synthetic(tmp_path_factory) -> SyntheticData:
    tmp = tmp_path_factory.mktemp("synthetic")
    rng = np.random.default_rng(0)
    regions = [
        f"chr{1 + i % 3}:{1_000 * (i // 3) + 10_000}-{1_000 * (i // 3) + 10_000 + REGION_LENGTH}"
        for i in range(N_REGIONS)
    ]
    motifs = [f"motif_{i}" for i in range(N_MOTIFS)]
    foreground_idx = np.arange(0, N_REGIONS, 10)  # 200 regions
    background_idx = np.setdiff1d(np.arange(N_REGIONS), foreground_idx)

    # Scores rounded to 1 decimal (as Cluster-Buster), with many ties and zeros.
    scores = np.round(rng.gamma(1.0, 2.0, size=(N_MOTIFS, N_REGIONS)), 1).astype(np.float32)
    scores[rng.random(scores.shape) < 0.3] = 0.0
    for i, motif in enumerate(motifs):
        if motif in ENRICHED_MOTIFS:
            scores[i, foreground_idx] += np.float32(6.0)
    # Rankings: 0-based, highest score first, ties broken at random.
    rankings = np.empty(scores.shape, dtype=np.int32)
    for i in range(N_MOTIFS):
        order = np.lexsort((rng.random(N_REGIONS), -scores[i]))
        rankings[i, order] = np.arange(N_REGIONS, dtype=np.int32)

    data = SyntheticData(
        rankings_db=tmp / "synthetic.regions_vs_motifs.rankings.feather",
        scores_db=tmp / "synthetic.regions_vs_motifs.scores.feather",
        foreground_bed=tmp / "foreground.bed",
        background_bed=tmp / "background.bed",
        annotation=tmp / "motif_annotation.tbl",
        regions=regions,
        motifs=motifs,
        scores=scores,
        rankings=rankings,
        foreground=[regions[i] for i in foreground_idx],
        background=[regions[i] for i in background_idx],
    )
    _write_db(data.rankings_db, rankings, regions, motifs)
    _write_db(data.scores_db, scores, regions, motifs)
    _write_bed(data.foreground_bed, data.foreground)
    _write_bed(data.background_bed, data.background)
    annotation = pd.DataFrame(
        {
            "#motif_id": ["motif_0", "motif_0", "motif_1", "motif_2", "motif_2", "motif_3"],
            "gene_name": ["TFB", "TFA", "TFA", "MYOD1", "TFC", "TFD"],
            "motif_similarity_qvalue": [0.0, 0.0, 0.0, 0.0, 0.0001, 0.0],
            "orthologous_identity": [1.0, 1.0, 1.0, 1.0, 1.0, 0.5],
            "description": [
                "gene is directly annotated",
                "gene is directly annotated",
                "gene is directly annotated",
                "gene is directly annotated",
                "motif is annotated for similar motif X",
                "gene is orthologous to Y which is directly annotated",
            ],
        }
    )
    annotation.to_csv(data.annotation, sep="\t", index=False)
    return data
