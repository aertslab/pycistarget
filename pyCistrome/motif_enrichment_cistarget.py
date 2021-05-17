import os
import pyranges as pr
import pyCistrome.utils as utils
import numpy as np
from pyscenic.recovery import recovery, aucs as calc_aucs
from pyscenic.genesig import Regulon
from pyarrow.feather import FeatherReader
import pandas as pd
from itertools import repeat
from pyscenic.recovery import leading_edge4row
from functools import partial
from typing import Union, Dict, Sequence
from tqdm import tqdm
import logging
import ray
import sys


class cisTargetDatabase:
    def __init__(self, ctx_db_fname: str, index_name: str):
        if not os.path.isfile(ctx_db_fname):
            raise FileNotFoundError('{} does not exist!'.format(ctx_db_fname))
        if not ctx_db_fname.endswith('.feather'):
            raise NotImplementedError('Currently only feather formatted databases are supported!')
        
        self.ctx_db_fname = ctx_db_fname
        self.index_name = index_name

        #set name of the database based on file name
        self.name = os.path.basename(ctx_db_fname).split('.')[0]
        
        #try to open feather database
        try:
            reader = FeatherReader(ctx_db_fname)
        except:
            print('{} could not be read by FeatherReader'.format(ctx_db_fname))
        
        # Do not count column 1 as it contains the index with the name of the features.
        self.total_regions = reader.num_columns - 1

        #load regions from ctx database
        self.ctx_regions = tuple(
            reader.get_column_name(idx) for idx in range(reader.num_columns)
            if reader.get_column_name(idx) != index_name
        )
        #convert ctx_regions to pyranges obect
        self.pr_ctx_regions = pr.PyRanges(utils.region_names_to_coordinates(self.ctx_regions))
        del(reader)

    def load_rankings_for_input_regions(self, input_regions: Union[Sequence[str], pr.PyRanges]) -> None:
        """
        Loads rankings from cistarget database for the input regions
        :params input_regions: either a sequence of region names in the form of "chr:start-end" strings or a PyRanges object
        """
        if isinstance(input_regions, pr.PyRanges):
            #convert pyranges input regions to tuple of region names
            input_regions = tuple(utils.coord_to_region_names(input_regions))
        elif isinstance(input_regions, Sequence):
            if not all([type(v) == str for v in input_regions]):
                raise ValueError('input_regions should be either a PyRanges object or a sequence of region names in the form of "chr:start-end"')
            input_regions = tuple(input_regions)
        else:
            raise ValueError('input_regions should be either a PyRanges object or a sequence of region names in the form of "chr:start-end"')
        reader = FeatherReader(self.ctx_db_fname)
        ctx_regions2idx = {region: idx for idx, region in enumerate(self.ctx_regions)}
        self.rankings_db = reader.read_pandas(
            columns = (self.index_name, ) + tuple(sorted(input_regions, key=lambda region: ctx_regions2idx[region]))
        )
        if self.rankings_db.empty:
            raise Warning('Database is empty for the provided input regions, check input! Returning empty dataframe.')
            return pd.DataFrame()
        # Avoid copying the whole dataframe by replacing the index in place.
        # This makes loading a database twice as fast in case the database file is already in the filesystem cache.
        self.rankings_db.set_index(self.index_name, inplace = True)

def overlap_regions_w_cistarget_db(pr_regions: pr.PyRanges, 
                                   ctx_db: cisTargetDatabase,
                                   min_fraction_overlap: float = 0.4,
                                   n_cpu: int = 1, 
                                   log = None, 
                                   region_set_name = 'input') -> pr.PyRanges:
    """
    Overlap regions from the cistarget database with input regions taking into account a minimal fractional overlap.
    :param pr_regions: pyranges object with input regions
    :param ctx_db: cistarget database
    :param min_fraction_overlap: minimum fraction of overlap needed between a cistarget region and an input region.
    :param n_cpu: number of threads to use.
    :param log: logging object to log to
    :param region_set_name: name of the regions set
    :param return_type: what to return. One of: "forward", "reverse" or "both". Returns cistarget regions which overlap with input, input regions which overlap with cistarget regions or both respectively
    :return: regions from the cistarget database which overlap with the input regions 
    """
    if log == None:
        # Create logger
        level    = logging.INFO
        format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
        handlers = [logging.StreamHandler(stream=sys.stdout)]
        logging.basicConfig(level = level, format = format, handlers = handlers)
        log = logging.getLogger('pyCistrome')

    pr_ctx_regions = ctx_db.pr_ctx_regions
    input_regions = pr_regions
    if isinstance(input_regions, pr.PyRanges):
        #overlap ctx regions with topic regions taking into account fraction of overlap
        joined_regions = input_regions.join(other = pr_ctx_regions, strandedness = False, how = None, report_overlap=True, suffix = '_ctx')
        joined_regions = joined_regions.insert(pd.Series(name = 'Length', data = joined_regions.as_df()['End'] - joined_regions.as_df()['Start']))
        joined_regions = joined_regions.insert(pd.Series(name = 'Length_ctx', data = joined_regions.as_df()['End_ctx'] - joined_regions.as_df()['Start_ctx']))
        joined_regions = joined_regions.insert(pd.Series(name = 'FractionOverlaps', data = joined_regions.as_df()['Overlap'] / joined_regions.as_df()['Length']))
        joined_regions = joined_regions.insert(pd.Series(name = 'FractionOverlaps_ctx', data = joined_regions.as_df()['Overlap'] / joined_regions.as_df()['Length_ctx']))
        
        #subset so only regions with a minimal fractional overlap of min_fraction_overlap are kept
        joined_regions = joined_regions.subset(
            lambda region: np.logical_and(region.FractionOverlaps >= min_fraction_overlap, region.FractionOverlaps_ctx >= min_fraction_overlap),
            strand = False,
            nb_cpu = n_cpu)

        df_ctx_input_overlap = joined_regions.as_df()[['Chromosome', 'Start_ctx', 'End_ctx', 'FractionOverlaps', 'FractionOverlaps_ctx']]
        df_ctx_input_overlap.columns = ['Chromosome', 'Start', 'End', 'FractionOverlaps', 'FractionOverlaps_ctx']
        pr_ctx_input_overlap = pr.PyRanges(df_ctx_input_overlap)

        fraction_input_regions_in_database = len(pr_ctx_input_overlap) / len(input_regions)

    else:
        raise ValueError('input_regions has to be a single PyRanges object!')
    
    
    log.info('{}: {} % of regions overlap with regions in the cistarget database with a minmal fractional overlap of {}.'.
        format(region_set_name,
                str(np.round(fraction_input_regions_in_database, 3) * 100),
                str(min_fraction_overlap)))
    return joined_regions, pr_ctx_input_overlap 

class cisTarget:
    def __init__(self, pr_regions: pr.PyRanges, name: str, species: str, log):
        self.input_regions = pr_regions
        self.name = name
        self.species = species
        self.log = log
    
    def prepare_regions(self, ctx_db, min_fraction_overlap: float = 0.4, n_cpu: int = 1):
        self.pr_joined_regions, self.pr_ctx_input_overlap = overlap_regions_w_cistarget_db(  pr_regions = self.input_regions, 
                                                                                                ctx_db = ctx_db,
                                                                                                min_fraction_overlap = min_fraction_overlap,
                                                                                                n_cpu = n_cpu, 
                                                                                                log = self.log, 
                                                                                                region_set_name = self.name)
        self.region_set_signature = utils.region_sets_to_signature(pr_region_set = self.pr_ctx_input_overlap, region_set_name = self.name)

    def run(self,
            ctx_db: cisTargetDatabase,
            weighted_recovery: bool = False,
            auc_threshold: float = 0.005,
            nes_threshold: float = 3.0,
            rank_threshold: int = 20000) -> pd.DataFrame:
        """
        Finds features of which the rankings are enriched in the input region set
        :param ctx_db: cistarget database object with loaded regions
        :param weighted_recovery: wether or not to use weighted recovery in the analysis
        :param auc_threshold: the cut-off for fraction of ranked genomic regions at which to calculate AUC. This measure is then used for comparing all the features.
        :param nes_threshold: only the enriched regions with normalized enrichment score (NES) higher than the threshold will be returned
        :param rank_threshold: The total number of ranked regions to take into account when creating a recovery curve
        :return: a pandas dataframe with enriched features
        REFERENCES:
        ----------
        Van de Sande B., Flerin C., et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis.
        Nat Protoc. June 2020:1-30. doi:10.1038/s41596-020-0336-2
        """
        #if not isinstance(ctx_db, cisTargetDatabase):
        #    raise ValueError('ctx_db should be a cisTargetDatabase!')
        #if not isinstance(region_set_signature, Regulon):
        #    raise ValueError('region_set_signature should be a gene signature object!')
        log = self.log
        #hardcoded values
        COLUMN_NAME_NES = "NES"
        COLUMN_NAME_AUC = "AUC"
        COLUMN_NAME_GRP = "GROUP"
        COLUMN_NAME_MOTIF_ID = "MotifID"
        COLUMN_NAME_TARGET_GENES = "TargetRegions"
        COLUMN_NAME_RANK_AT_MAX = "RankAtMax"
        COLUMN_NAME_TYPE = "Type"

        region_set_signature = self.region_set_signature

        if not all(np.isin(region_set_signature.genes, ctx_db.ctx_regions)):
            raise ValueError('Not all input regions are in the cistarget database, overlap your input regions first!')
        
        #load rankings for input regions
        #print('Loading cistarget database ... ')
        #ctx_db.load_rankings_for_input_regions(region_set_signature.genes)

        #assertion just for debuging
        assert all(np.isin(region_set_signature.genes, ctx_db.rankings_db.columns.values)), 'all regions in input should be in the cistarget database after loading!'
        regions = np.array(list(region_set_signature.genes))
        features, rankings = ctx_db.rankings_db.index.values, ctx_db.rankings_db[regions].values
        weights = np.asarray([region_set_signature[region] for region in regions]) if weighted_recovery else np.ones(len(regions))

        # Calculate recovery curves, AUC and NES values.
        log.info('{}: Calculating AUC ...'.format(region_set_signature.transcription_factor))
        aucs = calc_aucs(ctx_db.rankings_db[regions], ctx_db.total_regions, weights, auc_threshold)
        ness = (aucs - aucs.mean()) / aucs.std()

        # Keep only features that are enriched, i.e. NES sufficiently high.
        enriched_features_idx = ness >= nes_threshold

        enriched_features = pd.DataFrame(index=pd.Index(features[enriched_features_idx], name = COLUMN_NAME_MOTIF_ID),
                                    data={COLUMN_NAME_NES: ness[enriched_features_idx],
                                        COLUMN_NAME_AUC: aucs[enriched_features_idx],
                                        COLUMN_NAME_GRP: repeat(region_set_signature.transcription_factor, sum(enriched_features_idx))})
        
        log.info('{}: Recovery analysis ...'.format(region_set_signature.transcription_factor))
        rccs, _ = recovery(ctx_db.rankings_db[regions], ctx_db.total_regions, weights, rank_threshold, auc_threshold, no_auc=True)  
        avgrcc = rccs.mean(axis=0)        
        avg2stdrcc = avgrcc + 2.0 * rccs.std(axis=0)

        rccs = rccs[enriched_features_idx, :]
        rankings = rankings[enriched_features_idx, :]

        enriched_features.columns = pd.MultiIndex.from_tuples(list(zip(repeat("Enrichment"),
                                                                        enriched_features.columns)))
        df_rnks = pd.DataFrame(index=enriched_features.index,
                            columns=pd.MultiIndex.from_tuples(list(zip(repeat("Ranking"), regions))),
                            data=rankings)
        df_rccs = pd.DataFrame(index=enriched_features.index,
                            columns=pd.MultiIndex.from_tuples(list(zip(repeat("Recovery"), np.arange(rank_threshold)))),
                            data=rccs)
        enriched_features = pd.concat([enriched_features, df_rccs, df_rnks], axis=1)
        # Calculate the leading edges for each row. Always return importance from gene inference phase.
        weights = np.asarray([region_set_signature[region] for region in regions])
        enriched_features[[("Enrichment", COLUMN_NAME_TARGET_GENES), ("Enrichment", COLUMN_NAME_RANK_AT_MAX)]] = enriched_features.apply(
            partial(leading_edge4row, avg2stdrcc=avg2stdrcc, genes=regions, weights=weights), axis=1)
        log.info('{}: Done!'.format(region_set_signature.transcription_factor))
        self.enriched_features = enriched_features['Enrichment']

    def cistarget_get_cistrome(self, 
                               type_annotation_to_include: list = ['Direct_annot'],
                               motif2tf_fname: str = None,
                               version: str = 'v9',
                               motif_similarity_fdr: float = 0.001,
                               orthologous_identity_threshold: float = 0.0):

        motif_annotation = utils.load_motif_annotations(specie = self.species, 
                                                        fname = motif2tf_fname, 
                                                        version = version, 
                                                        motif_similarity_fdr = motif_similarity_fdr, 
                                                        orthologous_identity_threshold = orthologous_identity_threshold)
        
        annotated = pd.merge(self.enriched_features, motif_annotation, on = 'MotifID', how = 'left')
        motif_to_tf_dict = {motif: utils.flatten_list([v.split(', ') for v in annotated.loc[motif][type_annotation_to_include].dropna().values]) for motif in annotated.index}
        tf_to_motif_dict = utils.invert_dict_one_to_many(motif_to_tf_dict)
        tf_to_target_regions_dict = {TF: np.concatenate(annotated.loc[tf_to_motif_dict[TF]]['TargetRegions'].dropna())[:, 0] for TF in tf_to_motif_dict.keys()}
        df_joined_regions = self.pr_joined_regions.as_df()
        input_region_names = [chrom + ':' + str(start) + '-' + str(end) 
                                for chrom, start, end 
                                in zip(df_joined_regions['Chromosome'].values.to_list(), df_joined_regions['Start'].values, df_joined_regions['End'].values)]
        ctx_region_names = [chrom + ':' + str(start) + '-' + str(end) 
                                for chrom, start, end 
                                in zip(df_joined_regions['Chromosome'].values.to_list(), df_joined_regions['Start_ctx'].values, df_joined_regions['End_ctx'].values)]
        ctx_to_input = {ctx_region: input_region for ctx_region, input_region in zip(ctx_region_names, input_region_names)}
        self.cistromes = {TF: [ctx_to_input[region] for region in tf_to_target_regions_dict[TF]] for TF in tf_to_target_regions_dict.keys()}

@ray.remote
def ctx_ray(ctx_db: cisTargetDatabase,
            pr_regions: pr.PyRanges,
            name: str,
            species: str,
            weighted_recovery: bool = False,
            auc_threshold: float = 0.005,
            nes_threshold: float = 3.0,
            rank_threshold: int = 20000) -> pd.DataFrame:
    """
    Finds features of which the rankings are enriched in the input region set
    :param ctx_db: cistarget database object with loaded regions
    :param pr_regions: pyranges object with input regions
    :param name: name of the region set
    :param species: species from which the regions originate
    :param weighted_recovery: wether or not to use weighted recovery in the analysis
    :param auc_threshold: the cut-off for fraction of ranked genomic regions at which to calculate AUC. This measure is then used for comparing all the features.
    :param nes_threshold: only the enriched regions with normalized enrichment score (NES) higher than the threshold will be returned
    :param rank_threshold: The total number of ranked regions to take into account when creating a recovery curve
    :return: a pandas dataframe with enriched features
    REFERENCES:
    ----------
    Van de Sande B., Flerin C., et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis.
    Nat Protoc. June 2020:1-30. doi:10.1038/s41596-020-0336-2
    """
    #if not isinstance(ctx_db, cisTargetDatabase):
    #    raise ValueError('ctx_db should be a cisTargetDatabase!')
    #if not isinstance(region_set_signature, Regulon):
    #    raise ValueError('region_set_signature should be a gene signature object!')
    
    # Create logger
    level    = logging.INFO
    format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level = level, format = format, handlers = handlers)
    log = logging.getLogger('pyCistrome')
    log.info('Running '+ name)
    ctx_result = cisTarget(pr_regions, name, species, log)
    log.info('Overlapping regions from {} with cistarget database: {}'.format(name, ctx_db.name))
    ctx_result.prepare_regions(ctx_db)
    ctx_result.run(ctx_db)
    log.info('getting cistromes for {}'.format(name))
    ctx_result.cistarget_get_cistrome()
    return ctx_result
    
def cistarget_motif_enrichment(ctx_db: cisTargetDatabase,
                               pr_regions_dict: dict,
                               species:str,
                               weighted_recovery: bool = False,
                               auc_threshold: float = 0.005,
                               nes_threshold: float = 3.0,
                               rank_threshold: int = 20000,
                               n_cpu : int = 5,
                               **kwargs):
    """
    Finds features of which the rankings are enriched in the input region set
    :param ctx_db: cistarget database
    :param pr_regions_dict: dict of pyranges objects with input regions
    :param species: species from which the regions originate
    :param weighted_recovery: wether or not to use weighted recovery in the analysis
    :param auc_threshold: the cut-off for fraction of ranked genomic regions at which to calculate AUC. This measure is then used for comparing all the features.
    :param nes_threshold: only the enriched regions with normalized enrichment score (NES) higher than the threshold will be returned
    :param rank_threshold: The total number of ranked regions to take into account when creating a recovery curve
    :param n_cpu: number of cores to use
    :param **kwargs: additional parameters to pass to ray.init.
    :return: a dictionary of pandas data frames with enriched features
    """
    # Create logger
    level    = logging.INFO
    format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level = level, format = format, handlers = handlers)
    log = logging.getLogger('pyCistrome')
    
    #load rankings for all regions
    pr_all_regions = pr.PyRanges(utils.region_names_to_coordinates(
        list(set(utils.flatten_list([utils.coord_to_region_names(pr_regions_dict[key]) 
        for key in pr_regions_dict.keys()])))))
    _, pr_all_regions_ctx_overlap = overlap_regions_w_cistarget_db(pr_all_regions, ctx_db)
    all_regions_ctx_overlap = utils.coord_to_region_names(pr_all_regions_ctx_overlap)
    log.info('Loading rankings for {} regions from cistarget database: {}. This can take a while!'.format(len(all_regions_ctx_overlap), ctx_db.name))
    ctx_db.load_rankings_for_input_regions(all_regions_ctx_overlap)
    log.info('Done loading rankings for {} regions'.format(len(all_regions_ctx_overlap)))
    #run cistarget analysis in parallel
    ray.init(num_cpus=n_cpu, **kwargs); ctx_dict = ray.get([ctx_ray.remote(ctx_db = ctx_db, 
                                                                 pr_regions = pr_regions_dict[key], 
                                                                 name = key,  
                                                                 species = species,
                                                                 weighted_recovery = weighted_recovery,
                                                                 auc_threshold = auc_threshold, 
                                                                 nes_threshold = nes_threshold, 
                                                                 rank_threshold = rank_threshold) for key in list(pr_regions_dict.keys())]); ray.shutdown()
    
    ctx_dict = {key: ctx_result for key, ctx_result in zip(list(pr_regions_dict.keys()), ctx_dict)}
    return ctx_dict
