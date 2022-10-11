from mudata import MuData, AnnData
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from typing import List

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

def merge(mdata: MuData, split_pattern = ':') -> AnnData:
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

def motif_adata_to_TF_mudata(
    adata: AnnData,
    annotations_to_use_as_direct: List[str] = ['Direct_annot'],
    annotations_to_use_as_exteded: List[str] = ['Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot'],
    split_pattern = ':') -> MuData:
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