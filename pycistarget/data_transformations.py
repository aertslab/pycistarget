from mudata import MuData, AnnData
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from typing import List, Mapping
from pycistarget.motif_enrichment_cistarget import cisTargetDatabase
from pycistarget.utils import region_names_to_coordinates
import pyranges as pr

_flatten_list = lambda l: [item for sublist in l for item in sublist]

def _group_colnames_by_sample(columns: pd.Index, split_pattern: str) -> dict:
    #turn dataframe columns in the form of {sample}split_pattern{column_name} into a dictionary mapping column_names to samples
    column_names_across_samples = pd.DataFrame(columns.str.split(split_pattern).tolist()).groupby(1).agg(list).to_dict()[0]
    #add column_name to sample
    column_names_across_samples = {key: [f'{x}:{key}' for x in column_names_across_samples[key]] for key in column_names_across_samples.keys()}
    return column_names_across_samples

def _merge_columns_across_samples_if_possible(df, split_pattern = ':'):
    column_names_across_samples = _group_colnames_by_sample(df.columns, split_pattern)
    uniq = lambda x: list(set(x[~pd.isnull(x)]))
    for column_name in column_names_across_samples.keys():
        #check wether across the sample the same value are in each column
        uniq_values_columns = df[column_names_across_samples[column_name]].apply(uniq, axis = 1).tolist()
        if all([len(x) <= 1 for x in uniq_values_columns]):
            print(f"Merging columns for {column_name}")
            #create a new column containing the unique values
            df[column_name] = [x[0] if len(x) == 1 else np.nan for x in uniq_values_columns]
            #drop the old duplicate columns
            df.drop(column_names_across_samples[column_name], axis = 1, inplace = True)

def merge(mdata: MuData, split_pattern:str = ':') -> AnnData:
    """
    Merge `MuData` containing multiple motif enrichment results into a single `AnnData`

    Parameters
    ----------
    mdata: MuData
        A `MuData` containing multiple motif enrichment results
    split_patter: str
        Pattern used to seperate motif enrichment result name and column names in `.var` if merging the columns is not possible.
    
    Returns
    -------
    A single AnnData containing all motif enrichment results.

    Examples
    --------
    >>> mdata_motifs = dict_motif_enrichment_results_to_mudata(scplus_obj.menr)
    >>> adata_motifs = merge(mdata_motifs)
    >>> adata_motifs
        AnnData object with n_obs x n_vars = 139020 x 1337
            obs: 'Query', 'Target'
            var: 'DARs_cellLine_All_CTX_MM001:Region_set', 'DARs_cellLine_All_CTX_MM001:NES', 'DARs_cellLine_All_CTX_MM001:AUC', ...
    """
    var = pd.concat(
        objs = [mdata[k].var.add_prefix(f'{k}{split_pattern}') for k in mdata.mod.keys()],
        axis = 1)
    _merge_columns_across_samples_if_possible(var)
    obs = mdata.obs.copy()
    _merge_columns_across_samples_if_possible(obs)
    x = (pd.concat([mdata[key].to_df() for key in mdata.mod.keys()], axis = 1).T.groupby(level = 0).sum().T > 0) * 1
    return AnnData(
        X = csr_matrix(x.to_numpy()), dtype = np.uint8,
        var = var.loc[x.columns],
        obs = obs.loc[x.index])

def _get_TF_to_motif_mapping(adata: AnnData, annotations_to_use = ['Direct_annot']) -> Mapping[str, List[str]]:
    annotations_not_found = list(set(annotations_to_use) - set(adata.var.columns))
    if len(annotations_not_found) > 0:
        raise ValueError(f"The following annotations were not found in adata.var: {', '.join(annotations_not_found)}")
    annotations_to_list = lambda x: sorted(x.split(', ')) if not pd.isnull(x) else []
    return adata.var[annotations_to_use].applymap(annotations_to_list).sum(1).map(set).map(list) \
         .explode().reset_index().rename({'index': 'motif', 0: 'TF'}, axis = 1) \
         .dropna().groupby('TF').agg(list).to_dict()['motif']

def get_max_rank_for_TF_to_region(
    adata_motifs: AnnData,
    ctx_db_fname: str,
    annotations_to_use: List[str] = ['Orthology_annot', 'Direct_annot']) -> AnnData:
    """
    Create an `AnnData` object containing the max rank for each region across motifs of each TF

    Parameters
    ----------
    adata_motifs: AnnData,
        An AnnData containing motif enrichment results.
    ctx_db_fname: str,
        File path to a cistarget database (has to be the same database as used to generate the motif enrichment results).
    annotations_to_use: List[str] = ['Orthology_annot', 'Direct_annot']
        Which evidence of motif-to-TF annotations to use

    Returns
    -------
        An `AnnData` object containing the max rank for each region across motifs of each TF.

    Examples
    --------
    >>> mdata_motifs = dict_motif_enrichment_results_to_mudata(scplus_obj.menr)
    >>> adata_motifs = merge(mdata_motifs)
    >>> db_fname = 'cluster_SCREEN.regions_vs_motifs.rankings.v2.feather'
    >>> adata_max_rank = get_max_rank_for_TF_to_region(adata_motifs, db_fname)
    >>> adata_max_rank
        AnnData object with n_obs × n_vars = 139020 × 902
    
    """
    ctx_db = cisTargetDatabase(
        fname=ctx_db_fname, region_sets = pr.PyRanges(region_names_to_coordinates(adata_motifs.obs_names)))
    TF_to_motifs = _get_TF_to_motif_mapping(adata_motifs, annotations_to_use = annotations_to_use)
    l_motifs = list(TF_to_motifs.values())
    #convert motif names to numerical indices to allow indexing on the numpy array
    #this is 10x faster compared to indexing on the dataframe itself
    l_motifs_idx = [[ctx_db.db_rankings.index.get_loc(x) for x in m] for m in l_motifs]
    rankings = ctx_db.db_rankings.to_numpy()
    max_rank = np.array([rankings[x].min(0) for x in l_motifs_idx])
    #convert region names from database regions names in adata_motifs
    #this is important in case the database was not generated in the consensus peaks
    #can be skipped in case the databse is generated on consensus peaks
    if not all(ctx_db.regions_to_db['Target'] == ctx_db.regions_to_db['Query']):
        df_max_rank = pd.DataFrame(max_rank.T, index = ctx_db.db_rankings.columns, columns = list(TF_to_motifs.keys()))
        df_max_rank = ctx_db.regions_to_db \
            .merge(right = df_max_rank, how = 'left', left_on = 'Query', right_on = 'regions') \
            .drop('Query', axis = 1).set_index('Target')
        X = csr_matrix(df_max_rank.to_numpy())
        obs = pd.DataFrame(index = df_max_rank.index)
        var = pd.DataFrame(index = df_max_rank.columns)
    else:
        X = csr_matrix(max_rank.T)
        obs = pd.DataFrame(index = ctx_db.db_rankings.columns)
        var = pd.DataFrame(index = pd.Index(list(TF_to_motifs.keys())))
    #create and return AnnData
    return AnnData(X = X, dtype = np.uint32, obs = obs, var = var)

def motif_adata_to_TF_mudata(
    adata: AnnData,
    annotations_to_use_as_direct: List[str] = ['Direct_annot'],
    annotations_to_use_as_exteded: List[str] = ['Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot'],
    split_pattern = ':') -> MuData:
    """
    !THIS FUNCTION IS STILL EXPERIMENTAL!
    """
    column_names_across_samples = _group_colnames_by_sample(adata.var.columns, split_pattern)
    annotations_not_found = [x for x in [*annotations_to_use_as_direct, *annotations_to_use_as_exteded] if x not in column_names_across_samples.keys()]
    if len(annotations_not_found) > 0:
        raise ValueError(f"The following annotations were not found in adata.var: {', '.join(annotations_not_found)}")
    
    annotations_to_list = lambda x: sorted(x.split(', ')) if not pd.isnull(x) else []
    #1 get annotations, 2 convert string of fromat "TF1, TF2" to list and merge across columns, 3 make long table (explode), 4 generate a dict mapping TFs to a list of motifs
    if len(annotations_to_use_as_direct) > 0:
        TFs_to_motifs_direct = adata.var[
            _flatten_list([column_names_across_samples[annot] for annot in annotations_to_use_as_direct])] \
                .applymap(annotations_to_list).sum(1).map(set).map(list) \
                .explode().reset_index().rename({'index': 'motif', 0: 'TF'}, axis = 1) \
                .dropna().groupby('TF').agg(list).to_dict()['motif']
    if len(annotations_to_use_as_exteded) > 0:
        TFs_to_motifs_extended = adata.var[
            _flatten_list([column_names_across_samples[annot] for annot in annotations_to_use_as_exteded])] \
                .applymap(annotations_to_list).sum(1).map(set).map(list) \
                .explode().reset_index().rename({'index': 'motif', 0: 'TF'}, axis = 1) \
                .dropna().groupby('TF').agg(list).to_dict()['motif']
    
    mudata_constructor = {}

    X = adata.to_df()
    if len(annotations_to_use_as_direct) > 0:
        X_direct = pd.DataFrame(index = adata.obs_names, columns = TFs_to_motifs_direct.keys()).fillna(0)
        for TF in TFs_to_motifs_direct.keys():
            #merge the hits for all motifs linked to the TF
            X_direct[TF] = ((X[TFs_to_motifs_direct[TF]].sum(1) > 0) * 1).to_numpy()
        obs = adata.obs.loc[X_direct.index]
        var = pd.DataFrame(index = X_direct.columns, data = {'n_target_regions': [sum(X_direct[TF] != 0) for TF in X_direct.columns], 'evidence': 'direct'})
        adata_direct = AnnData(
            X = csr_matrix(X_direct.to_numpy()), dtype = np.uint8,
            var = var, obs = obs)
        mudata_constructor['direct'] = adata_direct

    if len(annotations_to_use_as_exteded) > 0:
        X_extended = pd.DataFrame(index = adata.obs_names, columns = TFs_to_motifs_extended.keys()).fillna(0)
        for TF in TFs_to_motifs_extended.keys():
            #merge the hits for all motifs linked to the TF
            X_extended[TF] = ((X[TFs_to_motifs_extended[TF]].sum(1) > 0) * 1).to_numpy()
        obs = adata.obs.loc[X_extended.index]
        var = pd.DataFrame(index = X_extended.columns, data = {'n_target_regions': [sum(X_extended[TF] != 0) for TF in X_extended.columns], 'evidence': 'extended'})
        adata_extended = AnnData(
            X = csr_matrix(X_extended.to_numpy()), dtype = np.uint8,
            var = var, obs = obs)
        mudata_constructor['extended'] = adata_extended
    return MuData(mudata_constructor)