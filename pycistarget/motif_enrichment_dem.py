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

# Set stderr to null when using ray.init to avoid ray printing Broken pipe million times
_stderr = sys.stderr                                                         
null = open(os.devnull,'wb') 

from .cluster_buster import *
from .utils import *

# DEM database
class DEMDatabase: 
    """
    DEM Database class.
    :class:`DEMDatabase` contains a dataframe with motifs as rows, regions as columns and CRM scores as
    values. In addition, is contains a slot to map query regions to regions in the database. For more
    information on how to generate databases, please visit: https://github.com/aertslab/create_cisTarget_databases
    
    Attributes
    ---------
    regions_to_db: pd.DataFrame
        A dataframe containing the mapping between query regions and regions in the database.
    db_scores: pd.DataFrame
        A dataframe with motifs as rows, regions as columns and CRM scores as values.
    total_regions: int
        Total number of regions in the database
    """
    def __init__(self, 
                fname: str,
                region_sets: Dict[str, pr.PyRanges] = None,
                name: str = None,
                fraction_overlap: float = 0.4):
        """
        Initialize DEMDatabase
        
        Parameters
        ---------
        fname: str
            Path to feather file containing the DEM database (regions_vs_motifs)
        region_sets: Dict or pr.PyRanges, optional
            Dictionary or pr.PyRanges that are going to be analyzed with DEM. Default: None.
        name: str, optional
            Name for the DEM database. Default: None
        fraction_overlap: float, optional
            Minimal overlap between query and regions in the database for the mapping.     
        """
        self.regions_to_db, self.db_scores = self.load_db(fname,
                                                          region_sets,
                                                          name,
                                                          fraction_overlap)
    
    def load_db(self,
                fname: str,
                region_sets: Dict[str, pr.PyRanges] = None,
                name: str = None,
                fraction_overlap: float = 0.4):
        """
        Load DEMDatabase
        
        Parameters
        ---------
        fname: str
            Path to feather file containing the DEM database (regions_vs_motifs)
        region_sets: Dict or pr.PyRanges, optional
            Dictionary or pr.PyRanges that are going to be analyzed with DEM. Default: None.
        name: str, optional
            Name for the DEM database. Default: None
        fraction_overlap: float, optional
            Minimal overlap between query and regions in the database for the mapping.     
            
        Return
        ---------
        target_to_db_dict: pd.DataFrame
            A dataframe containing the mapping between query regions and regions in the database.
        db_rankings: pd.DataFrame
            A dataframe with motifs as rows, regions as columns and CRM scores as values.
        total_regions: int
            Total number of regions in the database
        """
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
    """
    DEM class.
    :class:`DEM` contains DEM method for motif enrichment analysis on sets of regions. 
    
    Attributes
    ---------
    regions_to_db: pd.DataFrame
        A dataframe containing the mapping between query regions and regions in the database.
    region_sets: Dict
        A dictionary of PyRanges containing region coordinates for the regions to be analyzed.
    specie: str
        Specie from which genomic coordinates come from.
    subset_motifs: List, optional
        List of motifs to disregard in the analysis. Default: None
    contrasts: List, optional
        List of contrasts to perform. Default: None (Each group versus all the rest)
    name: str
        Analysis name
    max_bg_regions: int, optional
        Maximum number of regions to use as background. Default: None (All)
    adjpval_thr: float, optional
        Adjusted p-value threshold to consider a motif enriched. Default: 0.05
    log2fc_thr: float, optional
        Log2 Fold-change threshold to consider a motif enriched. Default: 1
    mean_fg_thr: float, optional
        Minimul mean signal in the foreground to consider a motif enriched. Default: 0
    motif_hit_thr: float, optional
        Minimal CRM score to consider a region enriched for a motif. Default: None (It will be automatically
        calculated based on precision-recall).
    n_cpu: int, optional
        Number of cores to use. Default: 1
    fraction_overlap: float, optional
        Minimal overlap between query and regions in the database for the mapping.
    cluster_buster_path: str, optional
        Path to cluster buster bin. Only required if using a shuffled background. Default: None
    path_to_genome_fasta: str, optional.
         Path to genome fasta file. Only required if using a shuffled background. Default: None
    path_to_motifs: str, optional.
         Path to motif collection folder (in .cb format). Only required if using a shuffled background. 
         Default: None
    genome_annotation: pr.PyRanges, optional.
         Pyranges containing genome annotation (e.g. biomart). Only required if using promoter balance. 
         Default: None
    genome_annotation: pr.PyRanges, optional.
        A pyRanges containing transcription start sites for each gene, with 'Chromosome', 'Start' and 
        'Strand' as columns (additional columns will be ignored). This data frame can be easily obtained
        via pybiomart:
            # Get TSS annotations
            import pybiomart as pbm
            # For mouse
            dataset = pbm.Dataset(name='mmusculus_gene_ensembl',  host='http://www.ensembl.org')
            # For human
            dataset = pbm.Dataset(name='hsapiens_gene_ensembl',  host='http://www.ensembl.org')
            # For fly
            dataset = pbm.Dataset(name='dmelanogaster_gene_ensembl',  host='http://www.ensembl.org')
            # Query TSS list and format
            annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
            filter = annot['Chromosome/scaffold name'].str.contains('CHR|GL|JH|MT')
            annot = annot[~filter]
            annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].str.replace(r'(\b\\S)', r'chr\1')
            annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
            # Select TSSs of protein coding genes
            annot = annot[annot.Transcript_type == 'protein_coding']
         Only required if using promoter balance. Default: None
    promoter_space: int, optional
        Space around TSS to consider a region promoter. Only used if using promoter balance.
        Default: 1000
    path_to_motif_annotations: str, optional
        Path to motif annotations. If not provided, they will be downloaded from 
        https://resources.aertslab.org based on the specie name provided (only possible for mus_musculus,
        homo_sapiens and drosophila_melanogaster). Default: None
    motif_similarity_fdr: float, optional
        Minimal motif similarity value to consider two motifs similar. Default: 0.001
    orthologous_identity_threshold: float, optional
        Minimal orthology value for considering two TFs orthologous. Default: 0.0
    motifs_to_use: List, optional
        A subset of motifs to use for the analysis. Default: None (All)
    tmp_dir: str, optional
        Temp directory to use if running cluster_buster. Default: None (\tmp)
    motif_enrichment: pd.DataFrame
        A dataframe containing motif enrichment results
    motif_hits: Dict
        A dictionary containing regions that are considered enriched for each motif.
    cistromes: Dict
        A dictionary containing TF cistromes. Cistromes with no extension contain regions linked to directly
        annotated motifs, while '_extended' cistromes can contain regions linked to motifs annotated by 
        similarity or orthology.
    """
    def __init__(self,
                 dem_db_fname,
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
                 genome_annotation: pr.PyRanges = None,
                 promoter_space: int = 1000,
                 path_to_motif_annotations: str = None,
                 annotation_version: str = 'v9',
                 motif_annotation: list = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot'],
                 motif_similarity_fdr: float = 0.001,
                 orthologous_identity_threshold: float = 0.0,
                 tmp_dir: int = None,
                 **kwargs):
        """
        Initialize DEM 
    
        Parameters
        ---------
        dem_db: :class:`DEMDatabase`
            A DEM database object.
        region_sets: Dict
            A dictionary of PyRanges containing region coordinates for the regions to be analyzed.
        specie: str
            Specie from which genomic coordinates come from.
        subset_motifs: List, optional
            List of motifs to disregard in the analysis. Default: None
        contrasts: List, optional
            List of contrasts to perform. Default: None (Each group versus all the rest)
        name: str
            Analysis name
        max_bg_regions: int, optional
            Maximum number of regions to use as background. Default: None (All)
        adjpval_thr: float, optional
            Adjusted p-value threshold to consider a motif enriched. Default: 0.05
        log2fc_thr: float, optional
            Log2 Fold-change threshold to consider a motif enriched. Default: 1
        mean_fg_thr: float, optional
            Minimul mean signal in the foreground to consider a motif enriched. Default: 0
        motif_hit_thr: float, optional
            Minimal CRM score to consider a region enriched for a motif. Default: None (It will be automatically
            calculated based on precision-recall).
        n_cpu: int, optional
            Number of cores to use. Default: 1
        fraction_overlap: float, optional
            Minimal overlap between query and regions in the database for the mapping.
        cluster_buster_path: str, optional
            Path to cluster buster bin. Only required if using a shuffled background. Default: None
        path_to_genome_fasta: str, optional.
            Path to genome fasta file. Only required if using a shuffled background. Default: None
        path_to_motifs: str, optional.
            Path to motif collection folder (in .cb format). Only required if using a shuffled background. 
            Default: None
        genome_annotation: pr.PyRanges, optional.
            Pyranges containing genome annotation (e.g. biomart). Only required if using promoter balance. 
            Default: None
        genome_annotation: pr.PyRanges, optional.
            A pyRanges containing transcription start sites for each gene, with 'Chromosome', 'Start' and 
            'Strand' as columns (additional columns will be ignored). This data frame can be easily obtained
            via pybiomart:
                # Get TSS annotations
                import pybiomart as pbm
                # For mouse
                dataset = pbm.Dataset(name='mmusculus_gene_ensembl',  host='http://www.ensembl.org')
                # For human
                dataset = pbm.Dataset(name='hsapiens_gene_ensembl',  host='http://www.ensembl.org')
                # For fly
                dataset = pbm.Dataset(name='dmelanogaster_gene_ensembl',  host='http://www.ensembl.org')
                # Query TSS list and format
                annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
                filter = annot['Chromosome/scaffold name'].str.contains('CHR|GL|JH|MT')
                annot = annot[~filter]
                annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].str.replace(r'(\b\\S)', r'chr\1')
                annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
                # Select TSSs of protein coding genes
                annot = annot[annot.Transcript_type == 'protein_coding']
            Only required if using promoter balance. Default: None
        promoter_space: int, optional
            Space around TSS to consider a region promoter. Only used if using promoter balance.
            Default: 1000
        path_to_motif_annotations: str, optional
            Path to motif annotations. If not provided, they will be downloaded from 
            https://resources.aertslab.org based on the specie name provided (only possible for mus_musculus,
            homo_sapiens and drosophila_melanogaster). Default: None
        motif_similarity_fdr: float, optional
            Minimal motif similarity value to consider two motifs similar. Default: 0.001
        orthologous_identity_threshold: float, optional
            Minimal orthology value for considering two TFs orthologous. Default: 0.0
        motifs_to_use: List, optional
            A subset of motifs to use for the analysis. Default: None (All)
        tmp_dir: str, optional
            Temp directory to use if running cluster_buster. Default: None (\tmp)
        **kwargs:
            Additional parameters to pass to `ray.init()`
        """
        self.dem_db_fname = dem_db_fname
        # Other params
        self.fraction_overlap = fraction_overlap
        self.region_sets = region_sets
        self.specie = specie
        self.subset_motifs= subset_motifs
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
        # For promoter awareness
        self.genome_annotation = genome_annotation
        self.promoter_space = promoter_space
        # For annotation
        self.annotation_version = annotation_version
        self.motif_annotation = motif_annotation
        self.path_to_motif_annotations = path_to_motif_annotations
        self.motif_similarity_fdr = motif_similarity_fdr
        self.orthologous_identity_threshold = orthologous_identity_threshold
        # Tmp
        self.tmp_dir = tmp_dir
        # Info
        self.motif_enrichment = None
        self.motif_hits = None
        self.cistromes = None
        
    def run(self, **kwargs):
        """
        Run DEM
    
        Parameters
        ---------
        dem_db_scores: pd.DataFrame
            A dataframe containing maximum CRM score for each motif in each regions.
        **kwargs
            Additional parameters to pass to `ray.init()`
        """
        # Create logger
        level    = logging.INFO
        format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
        handlers = [logging.StreamHandler(stream=sys.stdout)]
        logging.basicConfig(level = level, format = format, handlers = handlers)
        log = logging.getLogger('DEM')
        
        dem_db = DEMDatabase(self.dem_db_fname,
                            self.region_sets,
                            name = self.name,
                            fraction_overlap = self.fraction_overlap)
        self.regions_to_db = dem_db.regions_to_db
        if self.subset_motifs is not None:
            dem_db.db_scores = dem_db.db_scores.loc[self.subset_motifs,:]

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
                                   annotation = self.genome_annotation,
                                   promoter_space = self.promoter_space,
                                   motifs = dem_db.db_scores.index.tolist(),
                                   n_cpu = self.n_cpu) for x in range(len(contrasts))]

        # Compute p-val and log2FC
        if self.n_cpu > len(region_groups):
            self.n_cpu = len(region_groups)
        
        if self.n_cpu > 1:
            ray.init(num_cpus=self.n_cpu, **kwargs)
            sys.stderr = null
            DEM_list = ray.get([DEM_internal_ray.remote(dem_db.db_scores,
                                             region_groups[i],
                                             contrasts_names[i],
                                             adjpval_thr = self.adjpval_thr,
                                             log2fc_thr = self.log2fc_thr,
                                             mean_fg_thr = self.mean_fg_thr,
                                             motif_hit_thr = self.motif_hit_thr) for i in range(len(contrasts))])
            ray.shutdown()
            sys.stderr = sys.__stderr__ 
        else:
            DEM_list = [DEM_internal(dem_db.db_scores,
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
        cistromes_db = {key: get_cistromes_per_region_set(self.motif_enrichment[key], self.motif_hits['Database'][key], self.motif_annotation) for key in self.motif_enrichment.keys()}
        cistromed_rs = {key: get_cistromes_per_region_set(self.motif_enrichment[key], self.motif_hits['Region_set'][key], self.motif_annotation) for key in self.motif_enrichment.keys()}
        self.cistromes = {'Database': cistromes_db, 'Region_set': cistromed_rs}
        log.info('Done!')
        
    def add_motif_annotation_dem(self,
                       add_logo: Optional[bool] = True):
        """
        Add motif annotation

        Parameters
        ---------
        add_logo: boolean, optional
            Whether to add the motif logo to the motif enrichment dataframe
        """
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
                                          motif_similarity_fdr = self.motif_similarity_fdr,
                                          orthologous_identity_threshold = self.orthologous_identity_threshold)
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
                motif_enrichment_dict_w_annot = {key: motif_enrichment_dict_w_annot[key][sum([['Logo', 'Contrast'], self.motif_annotation, ['Log2FC', 'Adjusted_pval', 'Mean_fg', 'Mean_bg', 'Motif_hit_thr', 'Motif_hits']], [])] for key in motif_enrichment_dict_w_annot.keys()}
            else:
                motif_enrichment_dict_w_annot = {key: motif_enrichment_dict_w_annot[key][['Logo', 'Contrast', 'Log2FC', 'Adjusted_pval', 'Mean_fg', 'Mean_bg', 'Motif_hit_thr', 'Motif_hits']] for key in motif_enrichment_dict_w_annot.keys()}
        else:
            if annot_df is not None:
                motif_enrichment_dict_w_annot = {key: motif_enrichment_dict_w_annot[key][sum([['Contrast'], self.motif_annotation, ['Log2FC', 'Adjusted_pval', 'Mean_fg', 'Mean_bg', 'Motif_hit_thr', 'Motif_hits']],[])] for key in motif_enrichment_dict_w_annot.keys()}
            else:
                motif_enrichment_dict_w_annot = {key: motif_enrichment_dict_w_annot[key][['Contrast', 'Log2FC', 'Adjusted_pval', 'Mean_fg', 'Mean_bg', 'Motif_hit_thr', 'Motif_hits']] for key in motif_enrichment_dict_w_annot.keys()}
        
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
    """
    Shuffle given sequence
    """
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
                  annotation: pr.PyRanges = None,
                  promoter_space: int = 1000,
                  motifs: list = None,
                  n_cpu: int = 1):
    """"
    Format contrast groups
    """
    # Create DEM logger
    level = logging.INFO
    format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level=level, format=format, handlers=handlers)
    log = logging.getLogger('DEM')
    foreground = list(set(sum([region_sets_names[key] for key in contrast[0]],[])))
    if contrast[1] != 'Shuffle':
        background = list(set(sum([region_sets_names[key] for key in contrast[1]],[])))
        background = list(np.setdiff1d(background,foreground))
        if max_bg_regions is not None:
            if annotation is None:
                random.Random(555).shuffle(background)
                background = background[0:max_bg_regions]
            else:
                # Annotation of promoters
                annotation['End'] = annotation['Start']+promoter_space
                annotation['Start'] = annotation['Start']-promoter_space
                annotation = pr.PyRanges(annotation[['Chromosome', 'Start', 'End']])
                # Nr of promoters in the foreground
                fg_pr_overlap = pr.PyRanges(region_names_to_coordinates(foreground)).count_overlaps(annotation)
                fg_pr = coord_to_region_names(fg_pr_overlap[fg_pr_overlap.NumberOverlaps != 0])
                if len(fg_pr) == len(foreground):
                    nr_pr = max_bg_regions
                elif len(fg_pr) == 0:
                    nr_pr = 0
                else:
                    fg_pr = coord_to_region_names(fg_pr_overlap[fg_pr_overlap.NumberOverlaps != 0][['Chromosome', 'Start', 'End']])
                    fg_no_pr =  coord_to_region_names(fg_pr_overlap[fg_pr_overlap.NumberOverlaps == 0][['Chromosome', 'Start', 'End']])
                    nr_pr = int(max_bg_regions*(len(fg_pr)/(len(fg_pr) + len(fg_no_pr))))
                
                nr_no_pr = max_bg_regions-nr_pr
                # Nr of promoters in the background
                bg_pr_overlap = pr.PyRanges(region_names_to_coordinates(background)).count_overlaps(annotation)
                bg_pr = coord_to_region_names(bg_pr_overlap[bg_pr_overlap.NumberOverlaps != 0][['Chromosome', 'Start', 'End']])
                bg_no_pr = coord_to_region_names(bg_pr_overlap[bg_pr_overlap.NumberOverlaps == 0][['Chromosome', 'Start', 'End']])
                if len(bg_pr) < nr_pr:
                    nr_pr = len(bg_pr)
                nr_no_pr = max_bg_regions-nr_pr
                if len(bg_no_pr) < nr_no_pr:
                    nr_no_pr = len(bg_no_pr)
                # For promoters   
                bg_pr_subset = bg_pr.copy()
                random.Random(555).shuffle(bg_pr_subset)
                bg_pr = bg_pr_subset[0:nr_pr]
                # For others        
                bg_no_pr_subset = bg_no_pr.copy()
                random.Random(555).shuffle(bg_no_pr_subset)
                bg_no_pr = bg_no_pr_subset[0:nr_no_pr]         
                background = bg_no_pr+bg_pr 
        
    else:
        log.info('Generating and scoring shuffled background')
        if max_bg_regions is None:
            background_pr = pr.PyRanges(region_names_to_coordinates(foreground))
            background_sequences = pd.DataFrame([foreground, pr.get_fasta(background_pr, path_to_genome_fasta).tolist()], index=['Name', 'Sequence']).T
        else:
            if annotation is not None:
                annotation['End'] = annotation['Start']+promoter_space
                annotation['Start'] = annotation['Start']-promoter_space
                annotation = pr.PyRanges(annotation[['Chromosome', 'Start', 'End']])
                fg_pr_overlap = pr.PyRanges(region_names_to_coordinates(foreground)).count_overlaps(annotation)
                if len(fg_pr_overlap) == len(foreground):
                    nr_pr = len(foreground)
                elif len(fg_pr_overlap) == 0:
                    nr_pr = 0
                else:
                    fg_pr = coord_to_region_names(fg_pr_overlap[fg_pr_overlap.NumberOverlaps != 0][['Chromosome', 'Start', 'End']])
                    fg_no_pr =  coord_to_region_names(fg_pr_overlap[fg_pr_overlap.NumberOverlaps == 0][['Chromosome', 'Start', 'End']])
                    nr_pr = int(max_bg_regions*(len(fg_pr)/(len(fg_pr) + len(fg_no_pr))))
                nr_no_pr = max_bg_regions-nr_pr
                # For promoters           
                fg_pr_subset = fg_pr.copy()
                random.Random(555).shuffle(fg_pr_subset)
                bg_pr = fg_pr_subset[0:nr_pr]
                # For others        
                fg_no_pr_subset = fg_no_pr.copy()
                random.Random(555).shuffle(fg_no_pr_subset)
                bg_no_pr = fg_no_pr_subset[0:nr_no_pr] 
                region_names = bg_no_pr+bg_pr
                background_pr = pr.PyRanges(region_names_to_coordinates(region_names))  
            else:
                foreground_subset = foreground.copy()
                random.Random(555).shuffle(foreground_subset)
                region_names = foreground_subset[0:max_bg_regions]
                background_pr = pr.PyRanges(region_names_to_coordinates(region_names))   
                                          
            background_sequences = pd.DataFrame([region_names, pr.get_fasta(background_pr, path_to_genome_fasta).tolist()], index=['Name', 'Sequence']).T
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
    """"
    Internal DEM function to use with ray.
    """
            
    return DEM_internal(dem_db_scores, region_group, contrast_name, adjpval_thr, log2fc_thr, mean_fg_thr, motif_hit_thr)



def DEM_internal(dem_db_scores: pd.DataFrame,
            region_group: List[List[str]],
            contrast_name: str,
            adjpval_thr: Optional[float] = 0.05,
            log2fc_thr: Optional[float] = 1,
            mean_fg_thr: Optional[float] = 0,
            motif_hit_thr: Optional[float] = None):
    """
    Internal operations for DEM.
    """
    # Create DEM logger
    level = logging.INFO
    format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level=level, format=format, handlers=handlers)
    log = logging.getLogger('DEM')

    motifs = dem_db_scores.index.tolist()
    fg_mat = sparse.csc_matrix(dem_db_scores.loc[:,region_group[0]].values).toarray()
    
    if isinstance(region_group[1], list):
        bg_mat = sparse.csc_matrix(dem_db_scores.loc[:,region_group[1]].values).toarray()
    else:
        bg_mat = sparse.csc_matrix(region_group[1].values).toarray()
    
    # Delete db
    del dem_db_scores

    # Wilcox
    log.info('Computing DEM for ' + contrast_name)
    wilcox_test = [ranksums(fg_mat[x], y=bg_mat[x]) for x in range(fg_mat.shape[0])]
    # Log2FC
    mean_fg = fg_mat.mean(axis=1)
    mean_bg = bg_mat.mean(axis=1)
    logFC = np.log2((mean_fg + 10 ** -12) / (mean_bg + 10 ** -12)).tolist()

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
    """
    Benjamini-Hochberg p-value correction for multiple hypothesis testing.
    """
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

# Helper function to determine optimal threshold for motif hits
def get_motif_hits(scores, regions, labels, optimal_threshold=None):
    """
    Determine optimal score threshold based on precision-recall.
    """
    df = pd.DataFrame([labels, scores.toarray()[0]], columns=regions, index=['Label', 'Score']).sample(frac=1).T
    if optimal_threshold is None:
        df = df[df['Score'] > 0]
        fpr, tpr, thresholds = roc_curve(df['Label'], df['Score'])
        optimal_idx = np.argmax(tpr - fpr)
        optimal_threshold = thresholds[optimal_idx]
    motif_hits = df[(df['Score'] > optimal_threshold) & (df['Label'] == 1)].index.to_list()
    return motif_hits, optimal_threshold
    


    

