from ctxcore.genesig import Regulon
import numpy as np
import os
import pandas as pd
import pyranges as pr
import re
import subprocess
import ssl
from typing import Dict, List, Sequence, Union


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

def region_sets_to_signature(region_set: list,
                             region_set_name:str) -> Regulon:
    """
    generates a gene signature object from a dict of PyRanges objects
    :param pr_region_set: PyRanges object to be converted in genesignature object
    :param region_set_name: name of the regions set
    :param weights_col: if set uses this column of the pyranges object as gene2weight
    :return gene signature object of input region dict
    """
    
    weights = np.ones(len(region_set))
    regions_name = region_set
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
    
       
# Only implemented for Homer motif at the moment, but we can easily adapt it for other type of motifs as long as we
# write the corresponding conversion to meme 
def tomtom(homer_motif_path: str,
          meme_path: str,
          meme_collection_path: str):
    homer2meme(homer_motif_path)
    meme_motif_path = homer_motif_path.replace('.motif', '.meme')
    motif_name = os.path.splitext(os.path.basename(meme_motif_path))[0]
    os.makedirs(os.path.join(os.path.dirname(homer_motif_path), 'tomtom', motif_name), exist_ok=True)
    cmd = os.path.join(meme_path, 'tomtom') + ' -thresh 0.3 -oc %s %s %s'
    cmd = cmd % (os.path.join(os.path.dirname(homer_motif_path), 'tomtom', motif_name), meme_motif_path, meme_collection_path)
    try:
        subprocess.check_output(args=cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))
    tomtom_results = pd.read_csv(os.path.join(os.path.dirname(homer_motif_path), 'tomtom', motif_name, 'tomtom.tsv'), sep='\t')
    if not tomtom_results.dropna().shape[0] == 0:
        tomtom_results = tomtom_results.sort_values(by=['E-value'])
        homer_motif_name = tomtom_results.iloc[0,0]
        homer_motif_name = homer_motif_name.split('BestGuess:')[1]
        best_match_motif_name = tomtom_results.iloc[0,1]
        evalue = tomtom_results.iloc[0,4]
        return pd.DataFrame([homer_motif_name, best_match_motif_name, evalue], index=['Best Match/Details', 'Best Match/Tomtom', 'E-value/Tomtom']).transpose()
    
def homer2meme(homer_motif_path: str):
    out_file = open(homer_motif_path.replace('.motif', '.meme'), 'w')
    with open(homer_motif_path) as f:
        data = f.readlines()
    motif_name = data[0].split()[1]
    motif_id = data[0].split()[0][1:]
    out_file.write('MEME version 4.4\n\nALPHABET= ACGT\n\nstrands: + -\n\n' +
                   'Background letter frequencies (from uniform background):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n' +
                   'MOTIF '+ motif_name + ' ' + motif_id  + '\n')
    out_file.write('letter-probability matrix: nsites= 20 alength= 4 w= '+str(len(data)-1)+' E= 0 \n')
    for line in data[1:]:
        out_file.write('  ' + line)
    out_file.write('\n')
    out_file.close()
    
def get_TF_list(motif_enrichment_table: pd.DataFrame,
               annotation: List[str] = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']):
    tf_list = []
    for name in annotation:
        if name in motif_enrichment_table:
            tfs = motif_enrichment_table.loc[:,name].tolist()
            tfs = [x for x in tfs if str(x) != 'nan']
            tfs = list(set([element for item in tfs for element in item.split(', ') if str(re.search(',', str(item)))]))
            tf_list = tf_list + tfs
            
    return list(set(tf_list))

def get_motifs_per_TF(motif_enrichment_table: pd.DataFrame,
                    tf: str, 
                    motif_column: str,
                    annotation: List[str] = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']):
        motifs= []
        for name in annotation:
            if name in motif_enrichment_table:
                    if motif_column != 'Index':
                        motifs = motifs + motif_enrichment_table[motif_enrichment_table[name].str.contains(tf, na=False)][motif_column].tolist()
                    else:
                        motifs = motifs + motif_enrichment_table[motif_enrichment_table[name].str.contains(tf, na=False)].index.tolist()
        return list(set(motifs))
        
def get_cistrome_per_TF(motif_hits_dict,
                       motifs):
    return list(set(sum([motif_hits_dict[x] for x in motifs if x in motif_hits_dict.keys()],[])))
    
def inplace_change(filename, old_string, new_string):
    # Safely read the input filename using 'with'
    with open(filename) as f:
        s = f.read()
        if old_string not in s:
            return
    # Safely write the changed content, if found in the file
    with open(filename, 'w') as f:
        s = s.replace(old_string, new_string)
        f.write(s)
        
def get_position_index(query_list, target_list):
    d = {k: v for v, k in enumerate(target_list)}
    index = (d[k] for k in query_list)
    return list(index)

def target_to_query(target: Union[pr.PyRanges, List[str]],
         query: Union[pr.PyRanges, List[str]],
         fraction_overlap: float = 0.4):
    #Read input
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
    
    join_pr = target_pr.join(query_pr, report_overlap = True)
    join_pr.Overlap_query =  join_pr.Overlap/(join_pr.End_b - join_pr.Start_b)
    join_pr.Overlap_target =  join_pr.Overlap/(join_pr.End - join_pr.Start)
    join_pr = join_pr[(join_pr.Overlap_query > fraction_overlap) | (join_pr.Overlap_target > fraction_overlap)]
    target_regions = [str(chrom) + ":" + str(start) + '-' + str(end) for chrom, start, end in zip(list(join_pr.Chromosome), list(join_pr.Start), list(join_pr.End))]
    query_regions = [str(chrom) + ":" + str(start) + '-' + str(end) for chrom, start, end in zip(list(join_pr.Chromosome), list(join_pr.Start_b), list(join_pr.End_b))]
    target_to_query = pd.DataFrame([target_regions, query_regions], index=['Target', 'Query']).T
    return target_to_query
    
def get_cistromes_per_region_set(motif_enrichment_region_set,
                  motif_hits_regions_set,
                  annotation: List[str] = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']):
    tfs = get_TF_list(motif_enrichment_region_set)
    cistromes_per_region_set = {tf: get_cistrome_per_TF(motif_hits_regions_set,  get_motifs_per_TF(motif_enrichment_region_set, tf, motif_column = 'Index', annotation=annotation)) for tf in tfs}
    return cistromes_per_region_set
