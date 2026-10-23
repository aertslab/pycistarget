"""Base class for motif enrichment results (cisTarget and DEM) and (de)serialisation to hdf5."""

from __future__ import annotations

import functools
import logging
import os
from os import PathLike
from typing import Any, Literal

import h5py
import numpy as np
import pandas as pd

import pycistarget
from pycistarget.utils import get_cistromes_per_region_set, load_motif_annotations

log = logging.getLogger(__name__)

Species = Literal["homo_sapiens", "mus_musculus", "drosophila_melanogaster"]


def _mtime(path: str | None) -> float | None:
    try:
        return os.stat(path).st_mtime if path is not None else None
    except OSError:  # e.g. a URL
        return None


@functools.lru_cache(maxsize=4)
def _load_motif_annotations_cached(
    species: str,
    version: str,
    fname: str | None,
    motif_similarity_fdr: float,
    orthologous_identity_threshold: float,
    mtime: float | None,
) -> pd.DataFrame:
    """
    Motif annotations are loaded once when analysing several region sets (do not modify the result).

    ``mtime`` is the modification time of a local annotation file, so changes to it are picked up.
    """
    return load_motif_annotations(
        species,
        version=version,
        fname=fname,
        motif_similarity_fdr=motif_similarity_fdr,
        orthologous_identity_threshold=orthologous_identity_threshold,
    )


# Motif-to-TF annotations that failed to load (key: arguments of _load_motif_annotations_cached).
_failed_annotations: dict[tuple, Exception] = {}


_METADATA_FIELDS = {
    "cisTarget": [
        "species",
        "auc_threshold",
        "nes_threshold",
        "rank_threshold",
        "annotation_version",
        "annotation_to_use",
        "path_to_motif_annotations",
        "motif_similarity_fdr",
        "orthologous_identity_threshold",
        "motifs_to_use",
    ],
    "DEM": [
        "species",
        "adjpval_thr",
        "log2fc_thr",
        "mean_fg_thr",
        "motif_hit_thr",
        "path_to_motif_annotations",
        "annotation_version",
        "annotation_to_use",
        "motif_similarity_fdr",
        "orthologous_identity_threshold",
        "motifs_to_use",
    ],
}

# Strings are kept as object columns: PyTables stores them as (utf-8 encoded) fixed width strings.
_DTYPE_MAPPING: dict[str, Any] = {
    "Logo": object,
    "Region_set": object,
    "Direct_annot": object,
    "Motif_similarity_annot": object,
    "Orthology_annot": object,
    "Motif_similarity_and_Orthology_annot": object,
    "NES": np.float64,
    "AUC": np.float64,
    "Rank_at_max": np.float64,
    "Motif_hits": np.int64,
    "Motif_hit_thr": np.float64,
    "Target": object,
    "Query": object,
    "Log2FC": np.float64,
    "Adjusted_pval": np.float64,
    "Mean_fg": np.float64,
    "Mean_bg": np.float64,
}


def _with_hdf5_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    """Cast columns to the dtypes used in the hdf5 file."""
    df = df.copy()
    for column, dtype in _DTYPE_MAPPING.items():
        if column in df.columns:
            df[column] = df[column].astype(dtype)
    if df.index.dtype != object:
        df.index = df.index.astype(object)
    return df


class MotifEnrichmentResult:
    def __init__(
        self,
        name: str,
        species: Species,
        path_to_motif_annotations: str | None = None,
        annotation_version: str = "v10nr_clust",
        annotation_to_use: list[str] | None = None,
        motif_similarity_fdr: float = 0.001,
        orthologous_identity_threshold: float = 0.0,
        motifs_to_use: list[str] | None = None,
        force: bool = False,
    ):
        self.name = name
        self.species = species
        self.path_to_motif_annotations = path_to_motif_annotations
        self.annotation_version = annotation_version
        self.annotation_to_use = (
            list(annotation_to_use)
            if annotation_to_use is not None
            else ["Direct_annot", "Orthology_annot"]
        )
        self.motif_similarity_fdr = motif_similarity_fdr
        self.orthologous_identity_threshold = orthologous_identity_threshold
        self.motifs_to_use = list(motifs_to_use) if motifs_to_use is not None else None
        # Run without motif-to-TF annotations (and cistromes) when they can not be loaded.
        self.force = force
        self.logo_url = f"https://motifcollections.aertslab.org/{self.annotation_version}/logos/"
        self.motif_enrichment = pd.DataFrame()
        self.motif_hits: dict[str, dict[str, list[str]]] = {"Database": {}, "Region_set": {}}
        self.cistromes: dict[str, dict[str, list[str]]] = {"Database": {}, "Region_set": {}}
        self.regions_to_db = pd.DataFrame()

    def add_motif_annotation(self, add_logo: bool = True):
        """
        Add motif-to-TF annotation columns (and a motif logo column) to the motif enrichment table.

        If the annotations can not be loaded, an error is raised, unless ``force`` is set: then the
        annotation columns are left empty (and no cistromes are created), with a warning.
        """
        key = (
            self.species,
            self.annotation_version,
            self.path_to_motif_annotations,
            self.motif_similarity_fdr,
            self.orthologous_identity_threshold,
            _mtime(self.path_to_motif_annotations),
        )
        try:
            if key in _failed_annotations:
                # Do not try to load (download) them again for each region set.
                raise _failed_annotations[key]
            annot_df = _load_motif_annotations_cached(*key)
        except Exception as e:
            message = (
                "Unable to load motif-to-TF annotations"
                + (
                    f" from {self.path_to_motif_annotations}"
                    if self.path_to_motif_annotations
                    else " from the cisTarget webserver; provide them with path_to_motif_annotations"
                )
                + f": {e}"
            )
            if not self.force:
                raise RuntimeError(message) from e
            if key not in _failed_annotations:
                _failed_annotations[key] = e
                log.warning(f"{message}. Continuing without motif-to-TF annotations and cistromes.")
            annot_df = pd.DataFrame(columns=self.annotation_to_use, dtype=object)
        motif_enrichment_w_annot = self.motif_enrichment.merge(
            annot_df[self.annotation_to_use], how="left", left_index=True, right_index=True
        )
        if add_logo:
            logos = [
                f'<img src="{self.logo_url}{motif_name}.png" width="200" >'
                for motif_name in motif_enrichment_w_annot.index
            ]
            motif_enrichment_w_annot.insert(0, "Logo", logos)
        self.motif_enrichment = motif_enrichment_w_annot

    def get_cistromes(self):
        self.cistromes = {
            "Database": get_cistromes_per_region_set(
                self.motif_enrichment, self.motif_hits["Database"], self.annotation_to_use
            ),
            "Region_set": get_cistromes_per_region_set(
                self.motif_enrichment, self.motif_hits["Region_set"], self.annotation_to_use
            ),
        }

    def to_html(self) -> str:
        with pd.option_context("display.max_colwidth", None):
            return self.motif_enrichment.to_html(escape=False, col_space=80)

    def show_result(self):
        """Show the motif enrichment table (with logos) in a notebook."""
        try:
            from IPython.display import HTML
        except ImportError as e:
            raise ImportError(
                "IPython is required to show results, install it with: pip install pycistarget[notebook]"
            ) from e
        return HTML(self.to_html())

    @property
    def _analysis_type(self) -> str:
        raise NotImplementedError

    def write_hdf5(self, path: str | PathLike, mode: Literal["w", "a"] = "w"):
        """
        Write the motif enrichment result to a hdf5 file, under a group with the name of this analysis.

        Use ``mode="a"`` to add the result to an existing file.
        """
        analysis_type = self._analysis_type
        with h5py.File(path, mode) as h5:
            h5_root = h5.create_group(self.name)

            h5_metadata_grp = h5_root.create_group("metadata")
            for metadata_field in _METADATA_FIELDS[analysis_type]:
                data = getattr(self, metadata_field)
                if data is not None:
                    h5_metadata_grp.create_dataset(name=metadata_field, data=data)
            h5_metadata_grp.create_dataset(name="version", data=pycistarget.__version__)
            h5_metadata_grp.create_dataset(name="type", data=analysis_type)

            for group_name, data in (
                ("cistromes", self.cistromes),
                ("motif_hits", self.motif_hits),
            ):
                h5_grp = h5_root.create_group(group_name)
                for key, h5_sub_grp_name in (
                    ("Database", "database"),
                    ("Region_set", "region_set"),
                ):
                    h5_sub_grp = h5_grp.create_group(h5_sub_grp_name)
                    for name, regions in data[key].items():
                        h5_sub_grp.create_dataset(name=name, data=np.array(regions, dtype="S"))

        _with_hdf5_dtypes(self.motif_enrichment).to_hdf(
            path, key=f"{self.name}/motif_enrichment", append=True, mode="r+"
        )
        _with_hdf5_dtypes(self.regions_to_db).to_hdf(
            path, key=f"{self.name}/regions_to_db", append=True, mode="r+"
        )
