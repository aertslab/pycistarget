from ctxcore.genesig import Regulon, GeneSignature
from ctxcore.recovery import recovery, aucs as calc_aucs
from ctxcore.recovery import leading_edge4row
from ctxcore.rnkdb import FeatherRankingDatabase
from itertools import repeat
from functools import partial
import logging
import os
import numpy as np
import pandas as pd
import pyranges as pr
import ray
import ssl
import sys
from typing import Union, Dict, Sequence, Optional

from IPython.display import HTML
ssl._create_default_https_context = ssl._create_unverified_context
pd.set_option('display.max_colwidth', None)

from .utils import *

class cisTargetDatabase: 
    def __init__(self, 
                fname: str,
                region_sets: Union[Dict[str, pr.PyRanges], pr.PyRanges] = None,
                name: str = None,
                fraction_overlap: float = 0.4):
        self.regions_to_db, self.db_rankings, self.total_regions = self.load_db(fname,
                                                          region_sets,
                                                          name,
                                                          fraction_overlap)
    def load_db(self,
                fname: str,
                region_sets: Union[Dict[str, pr.PyRanges], pr.PyRanges] = None,
                name: str = None,
                fraction_overlap: float = 0.4):
        #Create logger
        level    = logging.INFO
        format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
        handlers = [logging.StreamHandler(stream=sys.stdout)]
        logging.basicConfig(level = level, format = format, handlers = handlers)
        log = logging.getLogger('cisTarget')
        
        log.info('Reading cisTarget database')
        
        if name is None:
            name = os.path.basename(fname)
        db = FeatherRankingDatabase(fname, name=name)
        total_regions = db.total_genes
        db_regions = db.genes
        prefix = None
        if '__' in db_regions[0]:
            prefix = db_regions[0].split('__')[0]
            db_regions = [x.split('__')[1] for x in db_regions]
        if region_sets is not None:
            if type(region_sets) == dict:
                target_to_db_dict = {x: target_to_query(region_sets[x], list(db_regions), fraction_overlap = fraction_overlap) for x in region_sets.keys()}
                target_regions_in_db = list(set(sum([target_to_db_dict[x]['Query'].tolist() for x in target_to_db_dict.keys()],[])))
            elif type(region_sets) == pr.PyRanges:
                target_to_db = target_to_query(region_sets, list(db_regions), fraction_overlap = fraction_overlap)
                target_to_db.index = target_to_db['Target']
                target_to_db_dict = target_to_db #for return purposes
                target_regions_in_db = list(set(target_to_db['Query'].tolist()))
            else:
                raise ValueError('region_sets should be either a dict of PyRanges objects or a single PyRanges object, not {}'.format(type(region_sets)))
            name='test'
            if prefix is not None:
                target_regions_in_db = [prefix + '__' + x for x in target_regions_in_db]
            target_regions_in_db = GeneSignature(name=name, gene2weight=target_regions_in_db)
            db_rankings = db.load(target_regions_in_db)
            if prefix is not None:
                db_rankings.columns = [x.split('__')[1] for x in db_rankings.columns]
        else:
            log.warn('Loading complete cistarget database, this can take a long time and consumes a lot of memory!')
            target_to_db_dict = None
            db_rankings = db.load_full()
        return target_to_db_dict, db_rankings, total_regions

# cisTarget class
class cisTarget:
    def __init__(self, 
                 ctx_db, 
                 region_set: pr.PyRanges,
                 name: str,
                 specie: str,
                 auc_threshold: float = 0.005,
                 nes_threshold: float = 3.0,
                 rank_threshold: float = 0.05,
                 path_to_motif_annotations: str = None,
                 annotation_version: str = 'v9',
                 annotation: list = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot'],
                 motif_similarity_fdr: float = 0.001,
                 orthologous_identity_threshold: float = 0.0):
        self.regions_to_db = ctx_db.regions_to_db[name] if type(ctx_db.regions_to_db) == dict else ctx_db.regions_to_db.loc[coord_to_region_names(region_set)]
        self.region_set = region_set
        self.name = name
        self.specie = specie
        self.auc_threshold = auc_threshold
        self.nes_threshold = nes_threshold
        self.rank_threshold = rank_threshold
        self.annotation_version = annotation_version
        self.annotation = annotation
        self.path_to_motif_annotations = path_to_motif_annotations
        self.motif_similarity_fdr = motif_similarity_fdr
        self.orthologous_identity_threshold = orthologous_identity_threshold
        # Run ctx
        self.run_ctx(ctx_db)

    def run_ctx(self,
            ctx_db: cisTargetDatabase) -> pd.DataFrame:
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
        
        # Create logger
        level    = logging.INFO
        format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
        handlers = [logging.StreamHandler(stream=sys.stdout)]
        logging.basicConfig(level = level, format = format, handlers = handlers)
        log = logging.getLogger('cisTarget')

        #Hardcoded values
        COLUMN_NAME_NES = "NES"
        COLUMN_NAME_AUC = "AUC"
        COLUMN_NAME_GRP = "GROUP"
        COLUMN_NAME_MOTIF_ID = "MotifID"
        COLUMN_NAME_TARGET_GENES = "TargetRegions"
        COLUMN_NAME_RANK_AT_MAX = "RankAtMax"
        COLUMN_NAME_TYPE = "Type"

        # Log
        log.info("Running cisTarget for {} which has {} regions".format(self.name, len(self.regions_to_db['Query'].tolist())))
        # Load signature as Regulon
        region_set_signature = region_sets_to_signature(self.regions_to_db['Query'].tolist(), region_set_name = self.name)
        # Get regions
        regions = np.array(list(region_set_signature.genes))
        db_rankings_regions = ctx_db.db_rankings[regions]
        #Get features, rankings and weights
        features, rankings = ctx_db.db_rankings.index.values, db_rankings_regions.values
        weights = np.asarray(np.ones(len(regions)))
        # Calculate recovery curves, AUC and NES values.
        aucs = calc_aucs(db_rankings_regions, ctx_db.total_regions, weights, self.auc_threshold)
        ness = (aucs - aucs.mean()) / aucs.std()
        # Keep only features that are enriched, i.e. NES sufficiently high.
        enriched_features_idx = ness >= self.nes_threshold
        # Make dataframe
        enriched_features = pd.DataFrame(index=pd.Index(features[enriched_features_idx], name = COLUMN_NAME_MOTIF_ID),
                                    data={COLUMN_NAME_NES: ness[enriched_features_idx],
                                        COLUMN_NAME_AUC: aucs[enriched_features_idx],
                                        COLUMN_NAME_GRP: repeat(region_set_signature.transcription_factor, sum(enriched_features_idx))})
        # Recovery analysis
        rccs, _ = recovery(db_rankings_regions, ctx_db.total_regions, weights, int(self.rank_threshold*ctx_db.total_regions), self.auc_threshold, no_auc=True)  
        avgrcc = rccs.mean(axis=0)        
        avg2stdrcc = avgrcc + 2.0 * rccs.std(axis=0)
        # Select features
        rccs = rccs[enriched_features_idx, :]
        rankings = rankings[enriched_features_idx, :]
        # Format df
        enriched_features.columns = pd.MultiIndex.from_tuples(list(zip(repeat("Enrichment"),
                                                                        enriched_features.columns)))
        df_rnks = pd.DataFrame(index=enriched_features.index,
                            columns=pd.MultiIndex.from_tuples(list(zip(repeat("Ranking"), regions))),
                            data=rankings)
        df_rccs = pd.DataFrame(index=enriched_features.index,
                            columns=pd.MultiIndex.from_tuples(list(zip(repeat("Recovery"), np.arange(int(self.rank_threshold*ctx_db.total_regions))))),
                            data=rccs)
        enriched_features = pd.concat([enriched_features, df_rccs, df_rnks], axis=1)
        # Calculate the leading edges for each row. Always return importance from gene inference phase.
        weights = np.asarray([region_set_signature[region] for region in regions])
        enriched_features[[("Enrichment", COLUMN_NAME_TARGET_GENES), ("Enrichment", COLUMN_NAME_RANK_AT_MAX)]] = enriched_features.apply(
            partial(leading_edge4row, avg2stdrcc=avg2stdrcc, genes=regions, weights=weights), axis=1)
        enriched_features = enriched_features['Enrichment'].rename_axis(None)
        # Format enriched features
        enriched_features.columns = ['NES', 'AUC', 'Region_set', 'Motif_hits', 'Rank_at_max']
        enriched_features = enriched_features.sort_values('NES', ascending=False)
        self.motif_enrichment = enriched_features[['Region_set', 'NES', 'AUC', 'Rank_at_max']]
        # Annotation
        log.info("Annotating motifs for " + self.name)
        self.add_motif_annotation_cistarget()
        # Motif hits
        db_motif_hits = {key: [enriched_features.loc[key, 'Motif_hits'][i][0] for i in range(len(enriched_features.loc[key, 'Motif_hits']))] for key in enriched_features.index}
        rs_motif_hits = {key: list(set(self.regions_to_db.loc[self.regions_to_db['Query'].isin(db_motif_hits[key]), 'Target'].tolist())) for key in db_motif_hits.keys()}
        self.motif_hits = {'Database': db_motif_hits, 'Region_set': rs_motif_hits}
        self.motif_enrichment['Motif_hits'] = [len(db_motif_hits[i]) for i in db_motif_hits.keys()]
        # Cistromes
        log.info("Getting cistromes for " + self.name)
        cistromes_db = get_cistromes_per_region_set(self.motif_enrichment, self.motif_hits['Database'], self.annotation)
        cistromes_rs = get_cistromes_per_region_set(self.motif_enrichment, self.motif_hits['Region_set'], self.annotation)
        self.cistromes = {'Database': cistromes_db, 'Region_set': cistromes_rs}
        
    def add_motif_annotation_cistarget(self,
                       add_logo: Optional[bool] = True):
        # Create cisTarget logger
        level = logging.INFO
        format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
        handlers = [logging.StreamHandler(stream=sys.stdout)]
        logging.basicConfig(level=level, format=format, handlers=handlers)
        log = logging.getLogger('cisTarget')

        # Read motif annotation. 
        try:
            annot_df = load_motif_annotations(self.specie,
                                          version = self.annotation_version,
                                          fname=self.path_to_motif_annotations,
                                          motif_similarity_fdr = self.motif_similarity_fdr,
                                          orthologous_identity_threshold = self.orthologous_identity_threshold)
            motif_enrichment_w_annot = pd.concat([self.motif_enrichment, annot_df], axis=1, sort=False).loc[self.motif_enrichment.index.tolist(),:]
        except:
            log.info('Unable to load annotation for ' + self.specie)
            annot_df = None
            motif_enrichment_w_annot = self.motif_enrichment
        # Add info to elements in dict
        if add_logo == True:
            motif_enrichment_w_annot['Logo']=['<img src="' +'https://motifcollections.aertslab.org/' + self.annotation_version + '/logos/'+ motif_enrichment_w_annot.index.tolist()[i] + '.png' + '" width="200" >' for i in range(motif_enrichment_w_annot.shape[0])]
            if annot_df is not None:
                motif_enrichment_w_annot = motif_enrichment_w_annot[sum([['Logo', 'Region_set'], self.annotation, ['NES', 'AUC', 'Rank_at_max']], [])]
            else:
                motif_enrichment_w_annot = motif_enrichment_w_annot[['Logo', 'Region_set', 'NES', 'AUC', 'Rank_at_max']]
        else:
            if annot_df is not None:
                motif_enrichment_w_annot = motif_enrichment_w_annot[sum([['Region_set'], self.annotation, ['NES', 'AUC', 'Rank_at_max']], [])]
            else:
                motif_enrichment_w_annot = motif_enrichment_w_annot[['Region_set', 'NES', 'AUC', 'Rank_at_max']]
        self.motif_enrichment = motif_enrichment_w_annot 

# Run cisTarget            
def run_cistarget(ctx_db: cisTargetDatabase,
                               region_sets: dict,
                               specie: str,
                               name: str = 'cisTarget',
                               fraction_overlap: float = 0.4,
                               auc_threshold: float = 0.005,
                               nes_threshold: float = 3.0,
                               rank_threshold: float = 0.05,
                               path_to_motif_annotations: str = None,
                               annotation_version: str = 'v9',
                               annotation: list = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot'],
                               motif_similarity_fdr: float = 0.001,
                               orthologous_identity_threshold: float = 0.0,
                               n_cpu : int = 1,
                               **kwargs):
    """
    Finds features of which the rankings are enriched in the input region set
    :param ctx_db: cistarget database
    :param pr_regions_dict: dict of pyranges objects with input regions
    :param species: species from which the regions originate
    :param auc_threshold: the cut-off for fraction of ranked genomic regions at which to calculate AUC. This measure is then used for comparing all the features.
    :param nes_threshold: only the enriched regions with normalized enrichment score (NES) higher than the threshold will be returned
    :param rank_threshold: The total number of ranked regions to take into account when creating a recovery curve
    :param n_cpu: number of cores to use
    :param **kwargs: additional parameters to pass to ray.init.
    :return: a dictionary of pandas data frames with enriched features
    """
    # Create cisTarget logger
    level = logging.INFO
    format = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level=level, format=format, handlers=handlers)
    log = logging.getLogger('cisTarget')
    
    # Load database
    if isinstance(ctx_db, str):
        ctx_db = cisTargetDatabase(ctx_db,
                             region_sets,
                             name = name,
                             fraction_overlap = fraction_overlap)
    
    # Run cistarget analysis in parallel
    if n_cpu > 1:
        ray.init(num_cpus=n_cpu, **kwargs)
        ctx_dict = ray.get([ctx_internal_ray.remote(ctx_db = ctx_db, 
                                            region_set = region_sets[key], 
                                            name = key,  
                                            specie = specie,
                                            auc_threshold = auc_threshold, 
                                            nes_threshold = nes_threshold, 
                                            rank_threshold = rank_threshold,
                                            path_to_motif_annotations = path_to_motif_annotations,
                                            annotation_version = annotation_version,
                                            annotation = annotation,
                                            motif_similarity_fdr = motif_similarity_fdr,
                                            orthologous_identity_threshold = orthologous_identity_threshold) for key in list(region_sets.keys())])
        ray.shutdown()
    else:
        ctx_dict = [ctx_internal(ctx_db = ctx_db, 
                                            region_set = region_sets[key], 
                                            name = key,  
                                            specie = specie,
                                            auc_threshold = auc_threshold, 
                                            nes_threshold = nes_threshold, 
                                            rank_threshold = rank_threshold,
                                            path_to_motif_annotations = path_to_motif_annotations,
                                            annotation_version = annotation_version,
                                            annotation = annotation,
                                            motif_similarity_fdr = motif_similarity_fdr,
                                            orthologous_identity_threshold = orthologous_identity_threshold) for key in list(region_sets.keys())]
    ctx_dict = {key: ctx_result for key, ctx_result in zip(list(region_sets.keys()), ctx_dict)}
    log.info('Done!')
    return ctx_dict
        
@ray.remote
def ctx_internal_ray(ctx_db: cisTargetDatabase,
            region_set: pr.PyRanges,
            name: str,
            specie: str,
            auc_threshold: float = 0.005,
            nes_threshold: float = 3.0,
            rank_threshold: float = 0.05,
            path_to_motif_annotations: str = None,
            annotation_version: str = 'v9',
            annotation: list = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot'],
            motif_similarity_fdr: float = 0.001,
            orthologous_identity_threshold: float = 0.0) -> pd.DataFrame:
    
    return ctx_internal(ctx_db = ctx_db,
                        region_set = region_set,
                        name = name,
                        specie = specie,
                        auc_threshold = auc_threshold,
                        nes_threshold = nes_threshold,
                        rank_threshold = rank_threshold,
                        path_to_motif_annotations = path_to_motif_annotations,
                        annotation_version = annotation_version,
                        annotation = annotation,
                        motif_similarity_fdr = motif_similarity_fdr,
                        orthologous_identity_threshold = orthologous_identity_threshold)


def ctx_internal(ctx_db: cisTargetDatabase,
            region_set: pr.PyRanges,
            name: str,
            specie: str,
            auc_threshold: float = 0.005,
            nes_threshold: float = 3.0,
            rank_threshold: float = 0.05,
            path_to_motif_annotations: str = None,
            annotation_version: str = 'v9',
            annotation: list = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot'],
            motif_similarity_fdr: float = 0.001,
            orthologous_identity_threshold: float = 0.0) -> pd.DataFrame:
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
    ctx_result = cisTarget(ctx_db,
                           region_set, 
                           name, 
                           specie,
                           auc_threshold,
                           nes_threshold,
                           rank_threshold,
                           path_to_motif_annotations,
                           annotation_version,
                           annotation,
                           motif_similarity_fdr,
                           orthologous_identity_threshold)
    return ctx_result
    
## Show results 
def cistarget_results(cistarget_dict,
                    name: Optional[str] = None):
        motif_enrichment_dict = {key: cistarget_dict[key].motif_enrichment for key in cistarget_dict.keys()}
        if name is None:
            motif_enrichment_table=pd.concat([motif_enrichment_dict[key] for key in motif_enrichment_dict.keys()], axis=0, sort=False)
        else:
            motif_enrichment_table=motif_enrichment_dict[name]
        return HTML(motif_enrichment_table.to_html(escape=False, col_space=80))