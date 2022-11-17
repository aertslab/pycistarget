from mudata import AnnData, MuData
from scipy.sparse import csr_matrix
from typing import Tuple, Iterator, Mapping, Union
import pandas as pd
import numpy as np
from pycistarget.motif_enrichment_cistarget import cisTarget
from pycistarget.motif_enrichment_dem import DEM

def _motif_hits_to_matrix(motif_hits: dict, features: list) -> Tuple[csr_matrix, list, list]:
    motif_hit_matrix = pd.DataFrame(index = features, columns = motif_hits.keys()).fillna(0)
    for motif in motif_hits.keys():
        motif_hit_matrix.loc[motif_hits[motif], motif] = 1
    matrix = csr_matrix(motif_hit_matrix.to_numpy())
    index = motif_hit_matrix.index
    columns = motif_hit_matrix.columns
    del(motif_hit_matrix)
    return (matrix, index, columns)

def cisTarget_to_AnnData(cistarget_obj: cisTarget) -> AnnData:
    """
    Convert cisTarget class object to AnnData object

    Parameters
    ----------
    cistarget_obj: :class:`cisTarget`
        A cisTarget class object
    
    Returns
    -------
    :class:`AnnData`
    """
    if not isinstance(cistarget_obj, cisTarget):
        raise ValueError("cistarget_obj should be an instance of the class cisTarget")
    #convert motif hits from dict of lists to sparse matrix, will become X in AnnData
    region_names = cistarget_obj.regions_to_db['Target'].tolist()
    motif_hit_matrix, region_names, motif_names = _motif_hits_to_matrix(cistarget_obj.motif_hits['Region_set'], region_names)
    #extract AnnData obs field, which is the mapping between region_set regions and cistarget database regions
    obs = cistarget_obj.regions_to_db.copy()
    obs.index = obs['Target'].to_numpy()
    obs = obs.loc[region_names]
    obs = obs.infer_objects()
    #extract AnnData var field, which is the motif enrichment dataframe
    var = cistarget_obj.motif_enrichment.copy()
    var = var.loc[motif_names]
    var = var.infer_objects()
    #create AnnData
    adata = AnnData(
        X = motif_hit_matrix, dtype = np.uint8,
        obs = obs,
        var = var)
    #add unstructured data
    adata.uns['name'] = cistarget_obj.name
    adata.uns['species'] = cistarget_obj.specie
    adata.uns['auc_threshold'] = cistarget_obj.auc_threshold
    adata.uns['nes_threshold'] = cistarget_obj.nes_threshold
    adata.uns['rank_threshold'] = cistarget_obj.rank_threshold
    adata.uns['annotation_version'] = cistarget_obj.annotation_version
    adata.uns['annotation'] = cistarget_obj.annotation
    adata.uns['path_to_motif_annotations'] = cistarget_obj.path_to_motif_annotations
    adata.uns['motif_similarity_fdr'] = cistarget_obj.motif_similarity_fdr
    adata.uns['orthologous_identity_threshold'] = cistarget_obj.orthologous_identity_threshold
    adata.uns['motifs_to_use'] = cistarget_obj.motifs_to_use
    return adata

def _dict_of_cisTarget_to_AnnData(cistarget_obj_dict: Mapping[str, cisTarget]) -> Tuple[str, AnnData]:
    for key in cistarget_obj_dict.keys():
        yield (key, cisTarget_to_AnnData(cistarget_obj_dict[key]))

def DEM_to_AnnData(dem_obj: DEM) -> Iterator[Tuple[str, AnnData]]:
    """
    Convert DEM class object to an iterator of AnnData objects, one for each contrast.

    Parameters
    ----------
    dem_obj: :class:`DEM`
        A DEM class object
    
    Returns
    -------
    An iterator of tuples for each contrast containing the name of the contrast and an AnnData object for that contrast.
    """
    if not isinstance(dem_obj, DEM):
        raise ValueError("dem_obj should be an instance of the class DEM")
    for contrast in dem_obj.region_sets.keys():
        #convert motif hits from dict of lists to sparse matrix, will become X in AnnData
        region_names = dem_obj.regions_to_db[contrast]['Target'].tolist()
        motif_hit_matrix, region_names, motif_names = _motif_hits_to_matrix(dem_obj.motif_hits['Region_set'][contrast], region_names)
        #extract AnnData obs field, which is the mapping between region_set regions and cistarget database regions
        obs = dem_obj.regions_to_db[contrast].copy()
        obs.index = obs['Target'].to_numpy()
        obs = obs.loc[region_names]
        obs = obs.infer_objects()
        #extract AnnData var field, which is the motif enrichment dataframe
        var = dem_obj.motif_enrichment[contrast].copy()
        var = var.loc[motif_names]
        var = var.infer_objects()
        #create AnnData
        adata = AnnData(
            X = motif_hit_matrix, dtype = np.uint8,
            obs = obs,
            var = var)
        adata.uns['specie'] = dem_obj.specie
        adata.uns['subset_motifs'] = dem_obj.subset_motifs
        adata.uns['contrasts'] = dem_obj.contrasts
        adata.uns['name'] = dem_obj.name
        adata.uns['max_bg_regions'] = dem_obj.max_bg_regions
        adata.uns['adjpval_thr'] = dem_obj.adjpval_thr
        adata.uns['log2fc_thr'] = dem_obj.log2fc_thr
        adata.uns['mean_fg_thr'] = dem_obj.mean_fg_thr
        adata.uns['motif_hit_thr'] = dem_obj.motif_hit_thr
        adata.uns['n_cpu'] = dem_obj.n_cpu
        adata.uns['cluster_buster_path'] = dem_obj.cluster_buster_path
        adata.uns['path_to_genome_fasta'] = dem_obj.path_to_genome_fasta
        adata.uns['path_to_motifs'] = dem_obj.path_to_motifs
        #adata.uns['genome_annotation'] = dem_obj.genome_annotation can not be stored
        adata.uns['promoter_space'] = dem_obj.promoter_space
        adata.uns['annotation_version'] = dem_obj.annotation_version
        adata.uns['motif_annotation'] = dem_obj.motif_annotation
        adata.uns['path_to_motif_annotations'] = dem_obj.path_to_motif_annotations
        adata.uns['motif_similarity_fdr'] = dem_obj.motif_similarity_fdr
        adata.uns['orthologous_identity_threshold'] = dem_obj.orthologous_identity_threshold
        adata.uns['tmp_dir'] = dem_obj.tmp_dir
        yield (contrast, adata)

def dict_motif_enrichment_results_to_mudata(menr: Mapping[str, Union[cisTarget, DEM]]) -> MuData:
    """
    Convert a dictionary of motif enrichment results to a single MuData

    Parameters
    ----------
    menr: Mapping[str, Union[cisTarget, DEM]]
        A dictionary of motif enrichment results
    
    Returns
    -------
    A single MuData containing all motif enrichment results
    """
    if not type(menr) == dict:
        raise ValueError('menr should be a dict!')
    mudata_constructor = {}
    for key in menr.keys():
        if type(menr[key]) == dict and all([isinstance(x, cisTarget) for x in menr[key].values()]):
            mudata_constructor.update(
                {f"{key}_{region_set_name}": adata for region_set_name, adata in _dict_of_cisTarget_to_AnnData(menr[key])})
        elif isinstance(menr[key], DEM):
            mudata_constructor.update(
                {f"{key}_{region_set_name}": adata for region_set_name, adata in DEM_to_AnnData(menr[key])})
        else:
            raise ValueError(f'unknown datatype for {key}: {type(menr[key])}')
    return MuData(mudata_constructor)
