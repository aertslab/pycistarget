import pycistarget
from pycistarget.utils import load_motif_annotations, get_cistromes_per_region_set
from typing import Optional, List, Dict, Literal
import h5py
import numpy as np
import pandas as pd
import ssl
from IPython.display import HTML
ssl._create_default_https_context = ssl._create_unverified_context
pd.set_option('display.max_colwidth', None)

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
        "motifs_to_use"
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
        "motifs_to_use"
    ]
}

_DTYPE_MAPPING = {
    "Logo": "S",
    "Region_set": "S",
    "Direct_annot": "S",
    "Motif_similarity_annot": "S",
    "Orthology_annot": "S",
    "Motif_similarity_and_Orthology_annot": "S",
    "NES": np.float64,
    "AUC": np.float64,
    "Rank_at_max": np.float64,
    "Motif_hits": np.int64,
    "Target": "S",
    "Query": "S",
    "Log2FC": np.float64,
    "Adjusted_pval": np.float64,
    "Mean_fg": np.float64,
    "Mean_bg": np.float64
}

class MotifEnrichmentResult:
    def __init__(
            self,
            name: str,
            species: Literal[
                "homo_sapiens", "mus_musculus", "drosophila_melanogaster"],
            path_to_motif_annotations: Optional[str] = None,
            annotation_version: str = "v10nr_clust",
            annotation_to_use: List[str] = [
                'Direct_annot',  'Orthology_annot'],
            motif_similarity_fdr: float = 0.001,
            orthologous_identity_threshold: float = 0.0,
            motifs_to_use: Optional[List[str]] = None):
        self.name = name
        self.species = species
        self.path_to_motif_annotations = path_to_motif_annotations
        self.annotation_version = annotation_version
        self.annotation_to_use = annotation_to_use
        self.motif_similarity_fdr = motif_similarity_fdr
        self.orthologous_identity_threshold = orthologous_identity_threshold
        self.motifs_to_use = motifs_to_use
        self.logo_url = f"https://motifcollections.aertslab.org/{self.annotation_version}/logos/"
        self.motif_enrichment = pd.DataFrame()
        self.motif_hits: Dict[str, Dict[str, List[str]]] = {
            "Database": dict(), "Region_set": dict()}
        self.cistromes: Dict[str, Dict[str, List[str]]] = {
            "Database": dict(), "Region_set": dict()}
        self.regions_to_db = pd.DataFrame()

    def add_motif_annotation(
            self,
            add_logo: bool = True):
        try:
            annot_df = load_motif_annotations(
                self.species,
                version = self.annotation_version,
                fname=self.path_to_motif_annotations,
                motif_similarity_fdr = self.motif_similarity_fdr,
                orthologous_identity_threshold = self.orthologous_identity_threshold)
            annot_df = annot_df[self.annotation_to_use]
            motif_enrichment_w_annot = self.motif_enrichment.merge(
                annot_df, how = "left", left_index = True, right_index = True)
        except:
            Warning("Unable to load motif-to-TF annotations!")
            motif_enrichment_w_annot = self.motif_enrichment
        if add_logo:
            motif_enrichment_w_annot['Logo']=[
                '<img src="' + self.logo_url + motif_name + '.png' + '" width="200" >'
                for motif_name in motif_enrichment_w_annot.index]
            # Put the logo at the first position
            motif_enrichment_w_annot = motif_enrichment_w_annot[
                [
                    *motif_enrichment_w_annot.columns[-1:],
                    *motif_enrichment_w_annot.columns[:-1]
                ]
            ]
        self.motif_enrichment = motif_enrichment_w_annot
    
    def get_cistromes(self):
        cistromes_db = get_cistromes_per_region_set(
            self.motif_enrichment, 
            self.motif_hits['Database'], 
            self.annotation_to_use)
        cistromes_rs = get_cistromes_per_region_set(
            self.motif_enrichment,
            self.motif_hits['Region_set'],
            self.annotation_to_use)
        self.cistromes = {"Database": cistromes_db, "Region_set": cistromes_rs}
    
    def show_result(self):
        return HTML(self.motif_enrichment.to_html(escape = False, col_space = 80))
    
    def write_hdf5(
            self,
            path: str,
            mode: Literal["w", "a"] = "w"):
        from pycistarget.motif_enrichment_cistarget import cisTarget
        from pycistarget.motif_enrichment_dem import DEM
        # Get analysis type
        if isinstance(self, DEM):
            analysis_type = "DEM"
        elif isinstance(self, cisTarget):
            analysis_type = "cisTarget"
        else:
            raise ValueError(f"Unrecognized analysis type: {type(self)}")
        h5 = h5py.File(path, mode)
        # Set root to name of cistarget run
        h5_root = h5.create_group(self.name)

        # Save metadata
        h5_metadata_grp = h5_root.create_group("metadata")
        for metadata_field in _METADATA_FIELDS[analysis_type]:
            data = getattr(self, metadata_field)
            if data is not None:
                h5_metadata_grp.create_dataset(name=metadata_field, data=data)
        h5_metadata_grp.create_dataset(name="version", data=pycistarget.__version__)
        h5_metadata_grp.create_dataset(name="type", data=analysis_type)

        # Save cistromes
        h5_cistromes_grp = h5_root.create_group("cistromes")

        # Save database coordinates of cistromes
        h5_cistromes_database_grp = h5_cistromes_grp.create_group("database")
        for cistrome in self.cistromes["Database"].keys():
            h5_cistromes_database_grp.create_dataset(
                name=cistrome,
                data=np.array(self.cistromes["Database"][cistrome], dtype="S"),
            )

        # Save region_set coordinates of cistromes
        h5_cistromes_regionset_grp = h5_cistromes_grp.create_group("region_set")
        for cistrome in self.cistromes["Region_set"].keys():
            h5_cistromes_regionset_grp.create_dataset(
                name=cistrome,
                data=np.array(self.cistromes["Region_set"][cistrome], dtype="S"),
            )

        # Save motif_hits
        h5_motifhits_grp = h5_root.create_group("motif_hits")

        # Save database coordinates of motif hits
        h5_motifhits_database_grp = h5_motifhits_grp.create_group("database")
        for motif in self.motif_hits["Database"].keys():
            h5_motifhits_database_grp.create_dataset(
                name=motif,
                data=np.array(self.motif_hits["Database"][motif], dtype="S"),
            )

        # Save region_set coordinates of motif hits
        h5_motifhits_regionset_grp = h5_motifhits_grp.create_group("region_set")
        for motif in self.motif_hits["Region_set"].keys():
            h5_motifhits_regionset_grp.create_dataset(
                name=motif,
                data=np.array(self.motif_hits["Region_set"][motif], dtype="S"),
            )

        # Make copy of motif enrichment dataframe and set data types
        motif_enrichment = self.motif_enrichment.copy()
        motif_enrichment = motif_enrichment.astype(
            {
                key: _DTYPE_MAPPING[key]
                for key in _DTYPE_MAPPING.keys()
                if key in motif_enrichment.columns
            }
        )

        # Close file handle.
        h5.close()

        # Write motif enrichment to disk
        motif_enrichment.to_hdf(
            path, key=f"{self.name}/motif_enrichment", append=True, mode="r+"
        )

        # Make a copy of regions_to_db dataframe
        # remove index, which is a copy of the target columns
        # and set data types.
        regions_to_db = self.regions_to_db.copy()
        regions_to_db.reset_index(drop=True)
        regions_to_db = regions_to_db.astype(
            {
                key: _DTYPE_MAPPING[key]
                for key in _DTYPE_MAPPING.keys()
                if key in regions_to_db.columns
            }
        )
        # Write regions_to_db to disk
        regions_to_db.to_hdf(
            path, key=f"{self.name}/regions_to_db", append=True, mode="r+"
        )