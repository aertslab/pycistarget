"""Read/write cistarget and DEM results to/from hdf5 files."""

import h5py
from typing import  Dict, List, Union

import pyranges as pr
from pycistarget.motif_enrichment_cistarget import cisTarget
from pycistarget.motif_enrichment_dem import DEM
import pandas as pd

def read_hdf5(path: str) -> Dict[str, Union[cisTarget, DEM]]:
    # Open hdf5 file
    h5 = h5py.File(path, mode="r")

    return_dict: Dict[str, Union[cisTarget, DEM]] = {}

    # Loop over root keys of hdf5 file.
    # In case multiple cisTarget results are saved in a single file.
    for name in h5.keys():
        analysis_type = h5[name]["metadata"]["type"][()].decode()
        # Read cistromes
        cistromes: Dict[str, Dict[str, List[str]]] = {}
        for database_or_regionset in h5[name]["cistromes"].keys():
            cistromes[database_or_regionset] = {}
            for cistrome_name in h5[name]["cistromes"][database_or_regionset].keys():
                cistromes[database_or_regionset][cistrome_name] = list(map(
                    bytes.decode,
                    h5[name]["cistromes"][database_or_regionset][cistrome_name][:]))
    
        # Read motif_hits
        motif_hits: Dict[str, Dict[str, List[str]]] = {}
        for database_or_regionset in h5[name]["motif_hits"].keys():
            motif_hits[database_or_regionset] = {}
            for motif_name in h5[name]["motif_hits"][database_or_regionset].keys():
                motif_hits[database_or_regionset][motif_name] = list(map(
                            bytes.decode,
                            h5[name]["motif_hits"][database_or_regionset][motif_name][:]))
        # Read metadata
        if analysis_type == "cisTarget":
            metadata_kwargs = dict(
                annotation_to_use = map(bytes.decode, h5[name]["metadata"]["annotation_to_use"][()]),
                annotation_version = h5[name]["metadata"]["annotation_version"][()].decode(),
                auc_threshold = h5[name]["metadata"]["auc_threshold"][()],
                motif_similarity_fdr = h5[name]["metadata"]["motif_similarity_fdr"][()],
                nes_threshold = h5[name]["metadata"]["nes_threshold"][()],
                orthologous_identity_threshold =  h5[name]["metadata"]["orthologous_identity_threshold"][()],
                rank_threshold = h5[name]["metadata"]["rank_threshold"][()],
                species = h5[name]["metadata"]["species"][()],
            )
        elif analysis_type == "DEM":
            metadata_kwargs = dict(
                annotation_to_use = map(bytes.decode, h5[name]["metadata"]["annotation_to_use"][()]),
                annotation_version = h5[name]["metadata"]["annotation_version"][()].decode(),
                log2fc_thr = h5[name]["metadata"]["log2fc_thr"][()],
                adjpval_thr = h5[name]["metadata"]["adjpval_thr"][()],
                motif_similarity_fdr = h5[name]["metadata"]["motif_similarity_fdr"][()],
                mean_fg_thr = h5[name]["metadata"]["mean_fg_thr"][()],
                orthologous_identity_threshold =  h5[name]["metadata"]["orthologous_identity_threshold"][()],
                species = h5[name]["metadata"]["species"][()],
            )
        else:
            raise ValueError(f"Unrecognized type {analysis_type}")

        #optional arguments
        if "motifs_to_use" in h5[name]["metadata"].keys():
            metadata_kwargs["motifs_to_use"] = map(bytes.decode, h5[name]["metadata"]["motifs_to_use"][()])
        if "path_to_motif_annotations" in h5[name]["metadata"].keys():
            metadata_kwargs["path_to_motif_annotations"] = h5[name]["metadata"]["path_to_motif_annotations"][()].decode()
        if "motif_hit_thr" in h5[name]["metadata"].keys():
            metadata_kwargs["motif_hit_thr"] = h5[name]["metadata"]["motif_hit_thr"][()]

        # Read motif enrichment and regions_to_db
        motif_enrichment = pd.read_hdf(
            path, key = f"{name}/motif_enrichment"
        )

        regions_to_db = pd.read_hdf(
            path, key = f"{name}/regions_to_db"
        )

        if analysis_type == "cisTarget":
            return_dict[name] = cisTarget(
                region_set=pr.PyRanges(), # Empty pyranges
                name = name,
                **metadata_kwargs)
        elif analysis_type == "DEM":
            return_dict[name] = DEM(
                foreground_regions=pr.PyRanges(), # Empty pyranges
                background_regions=pr.PyRanges(), # Empty pyranges TODO: save foreground and background regions?
                name = name,
                **metadata_kwargs)
        return_dict[name].cistromes = cistromes
        return_dict[name].motif_hits = motif_hits
        return_dict[name].motif_enrichment = motif_enrichment
        return_dict[name].regions_to_db = regions_to_db

    h5.close()
    return return_dict
