import pyranges as pr
import pandas as pd
import numpy as np
import math
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
from PIL import Image
import pyranges as pr
import re
from scipy import sparse
from pyscenic.genesig import Regulon
from typing import Dict, Sequence
import ssl

def coord_to_region_names(coord):
    if isinstance(coord, pr.PyRanges):
        coord = coord.as_df()
        return list(coord['Chromosome'].astype(str) + ':' + coord['Start'].astype(str) + '-' + coord['End'].astype(str))

def region_names_to_coordinates(region_names):
    chrom=pd.DataFrame([i.split(':', 1)[0] for i in region_names])
    coor=[i.split(':', 1)[1] for i in region_names]
    start=pd.DataFrame([int(i.split('-', 1)[0]) for i in coor])
    end=pd.DataFrame([int(i.split('-', 1)[1]) for i in coor])
    regiondf=pd.concat([chrom, start, end], axis=1, sort=False)
    regiondf.index=region_names
    regiondf.columns=['Chromosome', 'Start', 'End']
    return(regiondf)

def regions_overlap(target, query):
    # Read input
    if isinstance(target, str):
        target_pr=pr.read_bed(target)
    if isinstance(target, list):
        target_pr=pr.PyRanges(region_names_to_coordinates(target))
    if isinstance(target, pr.PyRanges):
        target_pr=target
    # Read input
    if isinstance(query, str):
        query_pr=pr.read_bed(query)
    if isinstance(query, list):
        query_pr=pr.PyRanges(region_names_to_coordinates(query))
    if isinstance(query, pr.PyRanges):
        query_pr=query
    
    target_pr = target_pr.overlap(query_pr)
    selected_regions = [str(chrom) + ":" + str(start) + '-' + str(end) for chrom, start, end in zip(list(target_pr.Chromosome), list(target_pr.Start), list(target_pr.End))]
    return selected_regions

def region_sets_to_signature(pr_region_set: pr.PyRanges,region_set_name:str, weights_col: str = None) -> Regulon:
    """
    generates a gene signature object from a dict of PyRanges objects
    :param pr_region_set: PyRanges object to be converted in genesignature object
    :param region_set_name: name of the regions set
    :param weights_col: if set uses this column of the pyranges object as gene2weight
    :return gene signature object of input region dict
    """
    
    if weights_col in pr_region_set.columns and weights_col != None:
        weights = pr_region_set.as_df()[weights_col]
    else:
        weights = np.ones(len(pr_region_set))
    regions_name = coord_to_region_names(pr_region_set)
    signature = Regulon(
                    name                 = region_set_name,
                    gene2weight          = dict(zip(regions_name, weights)),
                    transcription_factor = region_set_name,
                    gene2occurrence      = [])
    
    return signature
    
ssl._create_default_https_context = ssl._create_unverified_context

def load_motif_annotations(specie: str,
                           version: str = 'v9',
                           fname: str = None,
                           column_names=('#motif_id', 'gene_name',
                                         'motif_similarity_qvalue', 'orthologous_identity', 'description'),
                           motif_similarity_fdr: float = 0.001,
                           orthologous_identity_threshold: float = 0.0) -> pd.DataFrame:
    """
    Load motif annotations from a motif2TF snapshot.
    :param fname: the snapshot taken from motif2TF.
    :param column_names: the names of the columns in the snapshot to load.
    :param motif_similarity_fdr: The maximum False Discovery Rate to find factor annotations for enriched motifs.
    :param orthologuous_identity_threshold: The minimum orthologuous identity to find factor annotations
        for enriched motifs.
    :return: A dataframe.
    """
    # Create a MultiIndex for the index combining unique gene name and motif ID. This should facilitate
    # later merging.
    if fname is None:
        if specie == 'mus_musculus':
            name='mgi'
        elif specie == 'homo_sapiens':
            name='hgnc'
        elif specie == 'drosophila_melanogaster':
            name='flybase'
        fname = 'https://resources.aertslab.org/cistarget/motif2tf/motifs-'+version+'-nr.'+name+'-m0.001-o0.0.tbl'
    df = pd.read_csv(fname, sep='\t', usecols=column_names)
    df.rename(columns={'#motif_id':"MotifID",
                       'gene_name':"TF",
                       'motif_similarity_qvalue': "MotifSimilarityQvalue",
                       'orthologous_identity': "OrthologousIdentity",
                       'description': "Annotation" }, inplace=True)
    df = df[(df["MotifSimilarityQvalue"] <= motif_similarity_fdr) &
            (df["OrthologousIdentity"] >= orthologous_identity_threshold)]
    
    # Direct annotation
    df_direct_annot = df[(df["MotifSimilarityQvalue"]<= 0) & (df["OrthologousIdentity"] >= 1)]
    df_direct_annot = df_direct_annot.groupby(['MotifID'])['TF'].apply(lambda x: ', '.join(x)).reset_index()
    df_direct_annot.index = df_direct_annot['MotifID']
    df_direct_annot = pd.DataFrame(df_direct_annot['TF'])
    df_direct_annot.columns = ['Direct_annot']
    # Indirect annotation - by motif similarity
    motif_similarity_annot = df[(df["MotifSimilarityQvalue"]> 0) & (df["OrthologousIdentity"] >= 1)]
    motif_similarity_annot = motif_similarity_annot.groupby(['MotifID'])['TF'].apply(lambda x: ', '.join(x)).reset_index()
    motif_similarity_annot.index =  motif_similarity_annot['MotifID']
    motif_similarity_annot = pd.DataFrame(motif_similarity_annot['TF'])
    motif_similarity_annot.columns = ['Motif_similarity_annot']
    # Indirect annotation - by orthology
    orthology_annot = df[(df["MotifSimilarityQvalue"]<= 0) & (df["OrthologousIdentity"] < 1)]
    orthology_annot = orthology_annot.groupby(['MotifID'])['TF'].apply(lambda x: ', '.join(x)).reset_index()
    orthology_annot.index = orthology_annot['MotifID']
    orthology_annot = pd.DataFrame(orthology_annot['TF'])
    orthology_annot.columns = ['Orthology_annot']
    # Indirect annotation - by orthology
    motif_similarity_and_orthology_annot = df[(df["MotifSimilarityQvalue"]> 0) & (df["OrthologousIdentity"] < 1)]
    motif_similarity_and_orthology_annot = motif_similarity_and_orthology_annot.groupby(['MotifID'])['TF'].apply(lambda x: ', '.join(x)).reset_index()
    motif_similarity_and_orthology_annot.index = motif_similarity_and_orthology_annot['MotifID']
    motif_similarity_and_orthology_annot = pd.DataFrame(motif_similarity_and_orthology_annot['TF'])
    motif_similarity_and_orthology_annot.columns = ['Motif_similarity_and_Orthology_annot']
    # Combine
    df = pd.concat([df_direct_annot, motif_similarity_annot, orthology_annot, motif_similarity_and_orthology_annot], axis=1, sort=False)
    return df
    
def add_motif_url(df: pd.DataFrame, motif_column_name: str = None, motif_names: list = None, base_url: str = "http://motifcollections.aertslab.org/v9/logos/", key_to_add: str = 'Motif_url', verbose: bool = True) :
    """
    Add motif urls to dataframe with motif names.
    :param df: Dataframe containing motif names.
    :param base_url: url to motif collections containing images for each motif in the dataframe
    :param motif_column_name: column name containing the motif names.
    :param motif_names: list of motif names, must be in the same order as the dataframe rows.
    :param key_to_add: column name where the urls should be stored.
    :param verbose: if set to True prints warning messages.
    :return: DataFrame with motif urls.
    """
    if motif_column_name != None and motif_names == None:
        df[key_to_add] = [urljoin(base = base_url, url = motif) for motif in df[motif_column_name]]
    elif isinstance(motif_names, Sequence):
        if verbose:
            Warning('Using a list of motif names, this function assumes that this list has the same order as rows in the dataframe!')
        df[key_to_add] = [urljoin(base = base_url, url = motif) for motif in motif_names]
    else:
        raise Exception('Either provide a column name or a list of motif names.')
    return df

def flatten_list(l):
    return [item for sublist in l for item in sublist]

def invert_dict_one_to_many(d):
    keys = list(d.keys())
    values = list(d.values())
    return {v: [keys[idx] for idx, l in enumerate(values) if v in l ] for v in set(flatten_list(values))}