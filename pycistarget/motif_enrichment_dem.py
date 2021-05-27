from ctxcore.rnkdb import FeatherRankingDatabase
from ctxcore.genesig import GeneSignature
import io
import logging
import numpy as np
import os
import pandas as pd
import pyranges as pr
import random
import ray
import re
from scipy.stats import ranksums
import scipy.sparse as sparse
from sklearn.metrics import roc_curve, auc
import ssl
import subprocess
import sys
from typing import Union, Dict, Sequence, Optional

from IPython.display import HTML
ssl._create_default_https_context = ssl._create_unverified_context
pd.set_option('display.max_colwidth', None)

from .cluster_buster import *
from .utils import *

# DEM database
class DEMDatabase: 
    def __init__(self, 
                fname: str,
                region_sets: Dict[str, pr.PyRanges] = None,
                name: str = None,
                fraction_overlap: float = 0.4):
        self.regions_to_db, self.db_scores = self.load_db(fname,
                                                          region_sets,
                                                          name,
                                                          fraction_overlap)
    
    def load_db(self,
                fname: str,
                region_sets: Dict[str, pr.PyRanges] = None,
                name: str = None,
                fraction_overlap: float = 0.4):
        #Create logger
        level    = logging.INFO
        format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
        handlers = [logging.StreamHandler(stream=sys.stdout)]
        logging.basicConfig(level = level, format = format, handlers = handlers)
        log = logging.getLogger('DEM')
        
        log.info('Reading DEM database')
        
        if name is None:
            name = os.path.basename(fname)
        db = FeatherRankingDatabase(fname, name=name)
        db_regions = db.genes
        if region_sets is not None:
            target_to_db_dict = {x: target_to_query(region_sets[x], list(db_regions)) for x in region_sets.keys()}
            target_regions_in_db = list(set(sum([target_to_db_dict[x]['Query'].tolist() for x in target_to_db_dict.keys()],[])))
            target_regions_in_db = GeneSignature(name=name, gene2weight=target_regions_in_db)
            db_scores = db.load(target_regions_in_db)
        else:
            target_to_db_dict = None
            db_scores = db.load_full()
        return target_to_db_dict, db_scores

class DEM():
    def __init__(self,
                 dem_db: Union[str, 'DEMDatabase'],
                 region_sets: Dict[str, pr.PyRanges],
                 specie: str,
                 subset_motifs: Optional[List[str]] = None,
                 contrasts: Optional[Union[str,List]] = 'Other',
                 name: Optional[str] = 'DEM',
                 max_bg_regions: int = None,
                 adjpval_thr: Optional[float] = 0.05,
                 log2fc_thr: Optional[float] = 1,
                 mean_fg_thr: Optional[float] = 0,
                 motif_hit_thr : Optional[float] = None,
                 n_cpu: Optional[int] = 1,
                 fraction_overlap: float = 0.4,
                 cluster_buster_path: str = None,
                 path_to_genome_fasta: str = None,
                 path_to_motifs: str = None,
                 path_to_motif_annotations: str = None,
                 annotation_version: str = 'v9',
                 annotation: list = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot'],
                 tmp_dir: int = None,
                 **kwargs):
        
        # Load database
        if isinstance(dem_db, str):
            dem_db = DEMDatabase(dem_db,
                                 region_sets,
                                 name = name,
                                 fraction_overlap = fraction_overlap)

        self.regions_to_db = dem_db.regions_to_db
        # Other params
        self.region_sets = region_sets
        self.specie = specie
        self.subset_motifs= subset_motifs
        if subset_motifs is not None:
            dem_db_scores = dem_db_scores.loc[subset_motifs,:]
        self.contrasts = contrasts
        self.name = name
        self.max_bg_regions = max_bg_regions
        self.adjpval_thr = adjpval_thr
        self.log2fc_thr = log2fc_thr
        self.mean_fg_thr = mean_fg_thr
        self.motif_hit_thr = motif_hit_thr
        self.n_cpu = n_cpu
        # For Shuffle background options
        self.cluster_buster_path = cluster_buster_path
        self.path_to_genome_fasta =  path_to_genome_fasta
        self.path_to_motifs = path_to_motifs
        # For annotation
        self.annotation_version = annotation_version
        self.annotation = annotation
        self.path_to_motif_annotations = path_to_motif_annotations
        # Tmp
        self.tmp_dir = tmp_dir
        # Info
        self.motif_enrichment = None
        self.motif_hits = None
        self.cistromes = None
        self.run(dem_db.db_scores)
        
    def run(self, dem_db_scores, **kwargs):
        # Create logger
        level    = logging.INFO
        format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
        handlers = [logging.StreamHandler(stream=sys.stdout)]
        logging.basicConfig(level = level, format = format, handlers = handlers)
        log = logging.getLogger('DEM')
        
        contrast_keys=[x for x in self.region_sets.keys()]
        
        region_sets_names = {key: self.regions_to_db[key]['Query'].tolist() for key in self.regions_to_db.keys()}

        if self.contrasts == 'Other':
            if len(self.region_sets) > 1:
                levels=list(self.region_sets.keys())
                contrasts=[[[x], levels[:levels.index(x)] + levels[levels.index(x)+1:]] for x in levels]
                contrasts_names=levels
            else:
                contrasts = [[self.region_sets.keys()[0]],['Shuffle']]
                contrasts_names = [self.region_sets.keys()[0] + '_VS_Shuffle']  
                
        elif isinstance(self.contrasts, list):
            contrasts = self.contrasts
            contrasts_names=['_'.join(contrasts[i][0]) + '_VS_' +'_'.join(contrasts[i][1]) for i in range(len(contrasts))]
            for i in range(len(contrasts)):
                self.regions_to_db[contrasts_names[i]] = pd.concat([self.regions_to_db[x] for x in contrasts[i][0]])
                    
        elif self.contrasts == 'Shuffle':
            levels=list(self.region_sets.keys())
            contrasts=[[[x], 'Shuffle'] for x in levels]
            contrasts_names=levels
        
        # Get region groups
        log.info('Creating contrast groups')
        region_groups = [create_groups(contrast = contrasts[x],
                                   region_sets_names = region_sets_names,
                                   max_bg_regions = self.max_bg_regions,
                                   path_to_genome_fasta = self.path_to_genome_fasta,
                                   path_to_regions_fasta = os.path.join(self.tmp_dir, contrasts_names[x] +'.fa'),  
                                   cbust_path = self.cluster_buster_path,
                                   path_to_motifs = self.path_to_motifs,
                                   motifs = dem_db_scores.index.tolist(),
                                   n_cpu = self.n_cpu) for x in range(len(contrasts))]

        # Compute p-val and log2FC
        if self.n_cpu > len(region_groups):
            self.n_cpu = len(region_groups)
        
        if self.n_cpu > 1:
            ray.init(num_cpus=self.n_cpu, **kwargs)
            DEM_list = ray.get([DEM_internal_ray.remote(dem_db_scores,
                                             region_groups[i],
                                             contrasts_names[i],
                                             adjpval_thr = self.adjpval_thr,
                                             log2fc_thr = self.log2fc_thr,
                                             mean_fg_thr = self.mean_fg_thr,
                                             motif_hit_thr = self.motif_hit_thr) for i in range(len(contrasts))])
            ray.shutdown()
        else:
            DEM_list = [DEM_internal(dem_db_scores,
                                region_groups[i],
                                contrasts_names[i],
                                adjpval_thr = self.adjpval_thr,
                                log2fc_thr = self.log2fc_thr,
                                mean_fg_thr = self.mean_fg_thr,
                                motif_hit_thr = self.motif_hit_thr) for i in range(len(contrasts))]
        self.motif_enrichment = {contrasts_names[i]: DEM_list[i][0] for i in range(len(DEM_list))} 
        db_motif_hits =  {contrasts_names[i]: DEM_list[i][1] for i in range(len(DEM_list))} 
        # Add annotation and logo
        self.add_motif_annotation_dem()
        # Format motif hits
        rs_motif_hits = {name_1: {name_2: list(set(self.regions_to_db[name_1].loc[self.regions_to_db[name_1]['Query'].isin(db_motif_hits[name_1][name_2]), 'Target'].tolist())) for name_2 in db_motif_hits[name_1].keys()} for name_1 in contrasts_names}
        self.motif_hits = {'Database': db_motif_hits, 'Region_set': rs_motif_hits}
        # TF cistromes
        log.info('Forming cistromes')
        cistromes_db = {key: get_cistromes_per_region_set(self.motif_enrichment[key], self.motif_hits['Database'][key], self.annotation) for key in self.motif_enrichment.keys()}
        cistromed_rs = {key: get_cistromes_per_region_set(self.motif_enrichment[key], self.motif_hits['Region_set'][key], self.annotation) for key in self.motif_enrichment.keys()}
        self.cistromes = {'Database': cistromes_db, 'Region_set': cistromed_rs}
        log.info('Done!')
        
    def add_motif_annotation_dem(self,
                       motif_similarity_fdr: Optional[float] = 0.001,
                       orthologous_identity_threshold: Optional[float] = 0.0,
                       add_logo: Optional[bool] = True):
        # Create DEM logger
        level = logging.INFO
        format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
        handlers = [logging.StreamHandler(stream=sys.stdout)]
        logging.basicConfig(level=level, format=format, handlers=handlers)
        log = logging.getLogger('DEM')

        # Read motif annotation. 
        try:
            annot_df = load_motif_annotations(self.specie,
                                          self.annotation_version,
                                          fname=self.path_to_motif_annotations,
                                          motif_similarity_fdr = motif_similarity_fdr,
                                          orthologous_identity_threshold = orthologous_identity_threshold)
            motif_enrichment_dict_w_annot = {key: pd.concat([self.motif_enrichment[key], annot_df], axis=1, sort=False).loc[self.motif_enrichment[key].index.tolist(),:] for key in self.motif_enrichment.keys()}
        except:
            log.info('Unable to load annotation for ' + self.specie)
            annot_df = None
            motif_enrichment_dict_w_annot = self.motif_enrichment
        # Add info to elements in dict
        if add_logo == True:
            for key in self.motif_enrichment.keys():
                motif_enrichment_dict_w_annot[key]['Logo']=['<img src="' +'https://motifcollections.aertslab.org/' + self.annotation_version + '/logos/'+ motif_enrichment_dict_w_annot[key].index.tolist()[i] + '.png' + '" width="200" >' for i in range(motif_enrichment_dict_w_annot[key].shape[0])]
            if annot_df is not None:
                motif_enrichment_dict_w_annot = {key: motif_enrichment_dict_w_annot[key].loc[:,['Logo', 'Contrast', 'Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot', 'Log2FC', 'Adjusted_pval', 'Mean_fg', 'Mean_bg', 'Motif_hit_thr', 'Motif_hits']] for key in motif_enrichment_dict_w_annot.keys()}
            else:
                motif_enrichment_dict_w_annot = {key: motif_enrichment_dict_w_annot[key].loc[:,['Logo', 'Contrast', 'Log2FC', 'Adjusted_pval', 'Mean_fg', 'Mean_bg', 'Motif_hit_thr', 'Motif_hits']] for key in motif_enrichment_dict_w_annot.keys()}
        else:
            if annot_df is not None:
                motif_enrichment_dict_w_annot = {key: motif_enrichment_dict_w_annot[key].loc[:,['Contrast', 'Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot', 'Log2FC', 'Adjusted_pval', 'Mean_fg', 'Mean_bg', 'Motif_hit_thr', 'Motif_hits']] for key in motif_enrichment_dict_w_annot.keys()}
            else:
                motif_enrichment_dict_w_annot = {key: motif_enrichment_dict_w_annot[key].loc[:,['Contrast', 'Log2FC', 'Adjusted_pval', 'Mean_fg', 'Mean_bg', 'Motif_hit_thr', 'Motif_hits']] for key in motif_enrichment_dict_w_annot.keys()}
        
        self.motif_enrichment = motif_enrichment_dict_w_annot 
    
        
    def DEM_results(self,
                    name: Optional[str] = None):
        motif_enrichment_dict = self.motif_enrichment
        if name is None:
            motif_enrichment_table=pd.concat([motif_enrichment_dict[key] for key in motif_enrichment_dict.keys()], axis=0, sort=False)
        else:
            motif_enrichment_table=motif_enrichment_dict[name]
        return HTML(motif_enrichment_table.to_html(escape=False, col_space=80))
        

# Utils
## Shuffle sequences for shuffle background
def shuffle_sequence(sequence: str):
    shuffled_sequence = np.frombuffer(sequence.encode('utf-8'), dtype='uint8')
    np.random.shuffle(shuffled_sequence)
    return shuffled_sequence.tobytes().decode('utf-8')
## Create groups to compare
def create_groups(contrast: list,
                  region_sets_names: list,
                  max_bg_regions: int,
                  path_to_genome_fasta: str,
                  path_to_regions_fasta: str,
                  cbust_path: str,
                  path_to_motifs: str,
                  motifs: list = None,
                  n_cpu: int = 1):
    # Create DEM logger
    level = logging.INFO
    format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level=level, format=format, handlers=handlers)
    log = logging.getLogger('DEM')
    foreground = list(set(sum([region_sets_names[key] for key in contrast[0]],[])))
    if contrast[1] != 'Shuffle':
        background = list(set(sum([region_sets_names[key] for key in contrast[1]],[])))
        if max_bg_regions is not None:
            random.Random(555).shuffle(background)
            background = background[0:max_bg_regions]
    else:
        log.info('Generating and scoring shuffled background')
        if max_bg_regions is None:
            background_pr = pr.PyRanges(region_names_to_coordinates(foreground))
            background_sequences = pd.DataFrame([foreground, pr.get_fasta(background_pr, path_to_genome_fasta).tolist()], index=['Name', 'Sequence']).T
        else:
            foreground_subset = foreground.copy()
            random.Random(555).shuffle(foreground_subset)
            foreground_subset = foreground_subset[0:max_bg_regions]
            background_pr = pr.PyRanges(region_names_to_coordinates(foreground_subset))
            background_sequences = pd.DataFrame([foreground_subset, pr.get_fasta(background_pr, path_to_genome_fasta).tolist()], index=['Name', 'Sequence']).T
        background_sequences['Sequence'] = [shuffle_sequence(x) for x in background_sequences['Sequence']] 
        background_sequences['Name'] = '>' + background_sequences['Name'] 
        background_sequences.to_csv(path_to_regions_fasta, header=False, index=False, sep='\n')
        # Motifs should include .cb
        motifs = [motif + '.cb' for motif in motifs]
        background  = cluster_buster(cbust_path, path_to_motifs, path_to_regions_fasta=path_to_regions_fasta, n_cpu=n_cpu, motifs=motifs)
    return [foreground, background]
 
## Ray function for getting LogFC, pAdj and motif hits between groups   
@ray.remote
def DEM_internal_ray(dem_db_scores: pd.DataFrame,
            region_group: List[List[str]],
            contrast_name: str,
            adjpval_thr: Optional[float] = 0.05,
            log2fc_thr: Optional[float] = 1,
            mean_fg_thr: Optional[float] = 0,
            motif_hit_thr: Optional[float] = None):
            
    return DEM_internal(dem_db_scores, region_group, contrast_name, adjpval_thr, log2fc_thr, mean_fg_thr, motif_hit_thr)



def DEM_internal(dem_db_scores: pd.DataFrame,
            region_group: List[List[str]],
            contrast_name: str,
            adjpval_thr: Optional[float] = 0.05,
            log2fc_thr: Optional[float] = 1,
            mean_fg_thr: Optional[float] = 0,
            motif_hit_thr: Optional[float] = None):
    """
    Find differentiallly enriched motifs.

    Parameters
    ---------
    input_mat: :class:`pd.DataFrame`
        A data frame including CRM scores.
    barcode_group: List
        List of length 2, including foreground cells on the first slot and background on the second (or shuffled scores).
    contrast_name: str
        Name of the contrast
    adjpval_thr: float, optional
        Adjusted p-values threshold. Default: 0.05
    log2fc_thr: float, optional
        Log2FC threshold. Default: np.log2(1.5)

    Return
    ------
    List
        `class::pd.DataFrame` with the selected motifs and logFC and adjusted p-values.
    """
    # Create DEM logger
    level = logging.INFO
    format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level=level, format=format, handlers=handlers)
    log = logging.getLogger('DEM')

    motifs = dem_db_scores.index.tolist()
    fg_mat = sparse.csr_matrix(dem_db_scores.loc[:,region_group[0]].values).toarray()
    
    if isinstance(region_group[1], list):
        bg_mat = sparse.csr_matrix(dem_db_scores.loc[:,region_group[1]].values).toarray()
    else:
        bg_mat = sparse.csr_matrix(region_group[1].values).toarray()
    
    # Delete db
    del dem_db_scores

    # Wilcox
    log.info('Computing DEM for ' + contrast_name)
    wilcox_test = [ranksums(fg_mat[x], y=bg_mat[x]) for x in range(fg_mat.shape[0])]
    # Log2FC
    mean_fg = fg_mat.mean(axis=1)
    mean_bg = bg_mat.mean(axis=1)
    logFC = [np.log2((mean_fg[x] + 10**-12) / 
            (mean_bg[x] + 10**-12)) for x in range(fg_mat.shape[0])]
    # P-value correction
    pvalue = [wilcox_test[x].pvalue for x in range(len(wilcox_test))]
    adj_pvalue = p_adjust_bh(pvalue)
    name = [contrast_name] * len(adj_pvalue)
    # Motif df
    motif_df = pd.DataFrame([logFC, adj_pvalue, mean_fg, mean_bg, name], index=[
                                     'Log2FC', 'Adjusted_pval', 'Mean_fg', 'Mean_bg', 'Contrast'], columns=motifs).transpose()
    motif_df = motif_df.loc[motif_df['Adjusted_pval']
                                              <= adjpval_thr, :]
    motif_df = motif_df.loc[motif_df['Log2FC']
                                              >= log2fc_thr, :]
    motif_df = motif_df.loc[motif_df['Mean_fg']
                                              >= mean_fg_thr, :]
    motif_df = motif_df.sort_values(
        ['Log2FC', 'Adjusted_pval'], ascending=[False, True])
    
    # Motif hits versus background
    keep_motifs = motif_df.index.tolist()
    keep_motifs_index = get_position_index(keep_motifs, motifs)
    scores_mat = sparse.vstack([fg_mat[keep_motifs_index,].T, bg_mat[keep_motifs_index,].T], format='csr').T
    regions = region_group[0] + ['Bg']*bg_mat.shape[1]
    labels = [1]*fg_mat.shape[1] + [0]*bg_mat.shape[1]
    motif_hits_list = [get_motif_hits(scores_mat[i], regions, labels, motif_hit_thr) for i in range(len(keep_motifs))]
    motif_hits = {keep_motifs[i]: motif_hits_list[i][0] for i in range(len(keep_motifs))}
    motif_df['Motif_hit_thr'] = [motif_hits_list[i][1] for i in range(len(keep_motifs))]
    motif_df['Motif_hits'] = [len(motif_hits_list[i][0]) for i in range(len(keep_motifs))]
    motif_df['Motif_hits'] = motif_df['Motif_hits'].astype(int)
    return motif_df, motif_hits

# Helper function to adjust p-value
def p_adjust_bh(p: float):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

# Helper function to determine optimal threshold for motif hits
def get_motif_hits(scores, regions, labels, optimal_threshold=None):
    df = pd.DataFrame([labels, scores.toarray()[0]], columns=regions, index=['Label', 'Score']).sample(frac=1).T
    if optimal_threshold is None:
        df = df[df['Score'] > 0]
        fpr, tpr, thresholds = roc_curve(df['Label'], df['Score'])
        optimal_idx = np.argmax(tpr - fpr)
        optimal_threshold = thresholds[optimal_idx]
    motif_hits = df[(df['Score'] > optimal_threshold) & (df['Label'] == 1)].index.to_list()
    return motif_hits, optimal_threshold
    

