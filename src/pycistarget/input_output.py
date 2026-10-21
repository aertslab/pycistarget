"""Read cistarget and DEM results from hdf5 files."""

from __future__ import annotations

from os import PathLike
from typing import Any

import h5py
import numpy as np
import pandas as pd
import pyranges as pr

from pycistarget.motif_enrichment_cistarget import cisTarget
from pycistarget.motif_enrichment_dem import DEM

_REQUIRED_METADATA = {
    "cisTarget": [
        "annotation_to_use",
        "annotation_version",
        "auc_threshold",
        "motif_similarity_fdr",
        "nes_threshold",
        "orthologous_identity_threshold",
        "rank_threshold",
        "species",
    ],
    "DEM": [
        "annotation_to_use",
        "annotation_version",
        "log2fc_thr",
        "adjpval_thr",
        "motif_similarity_fdr",
        "mean_fg_thr",
        "orthologous_identity_threshold",
        "species",
    ],
}
_OPTIONAL_METADATA = ["motifs_to_use", "path_to_motif_annotations", "motif_hit_thr"]


def _decode(value: Any) -> Any:
    """Convert values read from hdf5 to python types (bytes -> str, arrays -> lists)."""
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.ndarray):
        return [_decode(v) for v in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _read_region_lists(group: h5py.Group) -> dict[str, dict[str, list[str]]]:
    return {
        database_or_regionset: {
            name: [region.decode() for region in group[database_or_regionset][name][:]]
            for name in group[database_or_regionset].keys()
        }
        for database_or_regionset in group.keys()
    }


def _decode_strings(df: Any) -> pd.DataFrame:
    assert isinstance(df, pd.DataFrame)
    for column in df.columns:
        if df[column].dtype == object:
            df[column] = df[column].map(lambda x: x.decode() if isinstance(x, bytes) else x)
    return df


def read_hdf5(path: str | PathLike) -> dict[str, cisTarget | DEM]:
    """
    Read motif enrichment results written with ``write_hdf5``.

    Returns a dictionary with the analysis name as key and a :class:`cisTarget` or :class:`DEM`
    object as value.
    """
    return_dict: dict[str, cisTarget | DEM] = {}
    with h5py.File(path, mode="r") as h5:
        for name in h5.keys():
            metadata_grp = h5[name]["metadata"]
            analysis_type = _decode(metadata_grp["type"][()])
            if analysis_type not in _REQUIRED_METADATA:
                raise ValueError(f"Unrecognized type {analysis_type}")
            metadata_kwargs = {
                field: _decode(metadata_grp[field][()])
                for field in _REQUIRED_METADATA[analysis_type] + _OPTIONAL_METADATA
                if field in metadata_grp and (analysis_type == "DEM" or field != "motif_hit_thr")
            }
            result: cisTarget | DEM
            if analysis_type == "cisTarget":
                result = cisTarget(region_set=pr.PyRanges(), name=name, **metadata_kwargs)
            else:
                result = DEM(
                    foreground_regions=pr.PyRanges(),
                    background_regions=pr.PyRanges(),
                    name=name,
                    **metadata_kwargs,
                )
            result.cistromes = _read_region_lists(h5[name]["cistromes"])
            result.motif_hits = _read_region_lists(h5[name]["motif_hits"])
            return_dict[name] = result

    for name, result in return_dict.items():
        result.motif_enrichment = _read_table(path, f"{name}/motif_enrichment")
        result.regions_to_db = _read_table(path, f"{name}/regions_to_db")
    return return_dict


def _read_table(path: str | PathLike, key: str) -> pd.DataFrame:
    try:
        return _decode_strings(pd.read_hdf(path, key=key))
    except KeyError:
        # Empty tables (e.g. no enriched motifs) are not written by pandas.
        return pd.DataFrame()
