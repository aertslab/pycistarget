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
import h5py
from collections.abc import Mapping
from utils import is_iterable_not_string

from IPython.display import HTML
ssl._create_default_https_context = ssl._create_unverified_context
pd.set_option('display.max_colwidth', None)

# Set stderr to null when using ray.init to avoid ray printing Broken pipe million times
_stderr = sys.stderr                                                         
null = open(os.devnull,'wb') 


from .utils import *

class cisTargetDatabase: 
    """
    cisTarget Database class.
    :class:`cisTargetDatabase` contains a dataframe with motifs as rows, regions as columns and rank as
    values. In addition, is contains a slot to map query regions to regions in the database. For more
    information on how to generate databases, please visit: https://github.com/aertslab/create_cisTarget_databases
    
    Attributes
    ---------
    regions_to_db: pd.DataFrame
        A dataframe containing the mapping between query regions and regions in the database.
    db_rankings: pd.DataFrame
        A dataframe with motifs as rows, regions as columns and rank as values.
    total_regions: int
        Total number of regions in the database
    """
    def __init__(self, 
                fname: str,
                region_sets: Union[Dict[str, pr.PyRanges], pr.PyRanges] = None,
                name: str = None,
                fraction_overlap: float = 0.4):
        """
        Initialize cisTargetDatabase
        
        Parameters
        ---------
        fname: str
            Path to feather file containing the cisTarget database (regions_vs_motifs)
        region_sets: Dict or pr.PyRanges, optional
            Dictionary or pr.PyRanges that are going to be analyzed with cistarget. Default: None.
        name: str, optional
            Name for the cistarget database. Default: None
        fraction_overlap: float, optional
            Minimal overlap between query and regions in the database for the mapping.     
        """
        self.regions_to_db, self.db_rankings, self.total_regions = self.load_db(fname,
                                                          region_sets,
                                                          name,
                                                          fraction_overlap)
    def load_db(self,
                fname: str,
                region_sets: Union[Dict[str, pr.PyRanges], pr.PyRanges] = None,
                name: str = None,
                fraction_overlap: float = 0.4):
        """
        Load cisTargetDatabase
        
        Parameters
        ---------
        fname: str
            Path to feather file containing the cisTarget database (regions_vs_motifs)
        region_sets: Dict or pr.PyRanges, optional
            Dictionary or pr.PyRanges that are going to be analyzed with cistarget. Default: None.
        name: str, optional
            Name for the cistarget database. Default: None
        fraction_overlap: float, optional
            Minimal overlap between query and regions in the database for the mapping.     
            
        Return
        ---------
        target_to_db_dict: pd.DataFrame
            A dataframe containing the mapping between query regions and regions in the database.
        db_rankings: pd.DataFrame
            A dataframe with motifs as rows, regions as columns and rank as values.
        total_regions: int
            Total number of regions in the database
        """
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
    """
    cisTarget class.
    :class:`cisTarget` contains method for motif enrichment analysis on sets of regions. 
    
    Attributes
    ---------
    regions_to_db: pd.DataFrame
        A dataframe containing the mapping between query regions and regions in the database.
    region_set: pr.PyRanges
        A PyRanges containing region coordinates for the regions to be analyzed.
    name: str
        Analysis name
    specie: str
        Specie from which genomic coordinates come from
    auc_threshold: float, optional
        The fraction of the ranked genome to take into account for the calculation of the
        Area Under the recovery Curve. Default: 0.005
    nes_threshold: float, optional
        The Normalized Enrichment Score (NES) threshold to select enriched features. Default: 3.0
    rank_threshold: float, optional
        The total number of ranked genes to take into account when creating a recovery curve.
        Default: 0.05
    path_to_motif_annotations: str, optional
        Path to motif annotations. If not provided, they will be downloaded from 
        https://resources.aertslab.org based on the specie name provided (only possible for mus_musculus,
        homo_sapiens and drosophila_melanogaster). Default: None
    annotation_version: str, optional
        Motif collection version. Default: v9
    annotation: List, optional
        Annotation to use for forming cistromes. It can be 'Direct_annot' (direct evidence that the motif is 
        linked to that TF), 'Motif_similarity_annot' (based on tomtom motif similarity), 'Orthology_annot'
        (based on orthology with a TF that is directly linked to that motif) or 'Motif_similarity_and_Orthology_annot'.
        Default: ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']
    motif_similarity_fdr: float, optional
        Minimal motif similarity value to consider two motifs similar. Default: 0.001
    orthologous_identity_threshold: float, optional
        Minimal orthology value for considering two TFs orthologous. Default: 0.0
    motifs_to_use: List, optional
        A subset of motifs to use for the analysis. Default: None (All)
    motif_enrichment: pd.DataFrame
        A dataframe containing motif enrichment results
    motif_hits: Dict
        A dictionary containing regions that are considered enriched for each motif.
    cistromes: Dict
        A dictionary containing TF cistromes. Cistromes with no extension contain regions linked to directly
        annotated motifs, while '_extended' cistromes can contain regions linked to motifs annotated by 
        similarity or orthology.
        
    References
    ---------
    Van de Sande B., Flerin C., et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis.
    Nat Protoc. June 2020:1-30. doi:10.1038/s41596-020-0336-2
    """
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
                 orthologous_identity_threshold: float = 0.0,
                 motifs_to_use: list = None):
        """
        Initialize cisTarget class.

        Parameters
        ---------
        ctx_db: :class:`cisTargetDatabase`
            A cistarget database object.
        region_set: pr.PyRanges
            A PyRanges containing region coordinates for the regions to be analyzed.
        name: str
            Analysis name
        specie: str
            Specie from which genomic coordinates come from
        auc_threshold: float, optional
            The fraction of the ranked genome to take into account for the calculation of the
            Area Under the recovery Curve. Default: 0.005
        nes_threshold: float, optional
            The Normalized Enrichment Score (NES) threshold to select enriched features. Default: 3.0
        rank_threshold: float, optional
            The total number of ranked genes to take into account when creating a recovery curve.
            Default: 0.05
        path_to_motif_annotations: str, optional
            Path to motif annotations. If not provided, they will be downloaded from 
            https://resources.aertslab.org based on the specie name provided (only possible for mus_musculus,
            homo_sapiens and drosophila_melanogaster). Default: None
        annotation_version: str, optional
            Motif collection version. Default: v9
        annotation: List, optional
            Annotation to use for forming cistromes. It can be 'Direct_annot' (direct evidence that the motif is 
            linked to that TF), 'Motif_similarity_annot' (based on tomtom motif similarity), 'Orthology_annot'
            (based on orthology with a TF that is directly linked to that motif) or 'Motif_similarity_and_Orthology_annot'.
            Default: ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']
        motif_similarity_fdr: float, optional
            Minimal motif similarity value to consider two motifs similar. Default: 0.001
        orthologous_identity_threshold: float, optional
            Minimal orthology value for considering two TFs orthologous. Default: 0.0
        motifs_to_use: List, optional
            A subset of motifs to use for the analysis. Default: None (All)
    
        References
        ---------
        Van de Sande B., Flerin C., et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis.
        Nat Protoc. June 2020:1-30. doi:10.1038/s41596-020-0336-2
        """
        self.regions_to_db = ctx_db.regions_to_db[name] if type(ctx_db.regions_to_db) == dict else ctx_db.regions_to_db.loc[set(coord_to_region_names(region_set)) & set(ctx_db.regions_to_db['Target'])]
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
        self.motifs_to_use = motifs_to_use
        # Run ctx
        self.run_ctx(ctx_db)

    def run_ctx(self,
            ctx_db: cisTargetDatabase) -> pd.DataFrame:
        """
        Run cisTarget

        Parameters
        ---------
        ctx_db: :class:`cisTargetDatabase`
            A cistarget database object.
    
        References
        ---------
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
        
        #subset rankings database on motifs and regions
        if self.motifs_to_use is not None:
            log.info('Using only user provided motifs')
            motifs_not_in_db = set.difference(set(self.motifs_to_use), set(ctx_db.db_rankings.index.values))
            if len(motifs_not_in_db) > 0:
                log.info('Some motifs provided by the parameter <motifs_to_use> are not in the rankings database: {}'.format(motifs_not_in_db))
            motifs_to_use = set(self.motifs_to_use) & set(ctx_db.db_rankings.index.values)
            db_rankings_regions = ctx_db.db_rankings.loc[motifs_to_use, regions]
        else:
            db_rankings_regions = ctx_db.db_rankings[regions]

        #Get features, rankings and weights
        features, rankings = db_rankings_regions.index.values, db_rankings_regions.values
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
        """
        Add motif annotation

        Parameters
        ---------
        add_logo: boolean, optional
            Whether to add the motif logo to the motif enrichment dataframe
    
        References
        ---------
        Van de Sande B., Flerin C., et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis.
        Nat Protoc. June 2020:1-30. doi:10.1038/s41596-020-0336-2
        """
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

    def to_h5ad(self,
                grp_or_fname: Union[h5py.Group, str],
                compression_type = 'gzip',
                compression_opts = 4,
                verbose = False):
        if not isinstance(grp_or_fname, h5py.Group):
            #create a hdf5 file with a new group
            had5_file = h5py.File(grp_or_fname, "w")
            if verbose:
                print(f"Creating group: {self.name}")
            grp = had5_file.create_group(name = self.name)
        else:
            grp = grp_or_fname
        try:
            #save each attribute
            for attr_name in self.__dict__.keys():
                if getattr(self, attr_name) is not None:
                    if not is_iterable_not_string(getattr(self, attr_name)):
                        #numbers and strings can be saved directly as attributes
                        if verbose:
                            print(f"Adding attribute {attr_name} to group {grp.name}")
                        grp.attrs[attr_name] = getattr(self, attr_name)
                    elif not isinstance(getattr(self, attr_name), Mapping):
                        #array like objects can be saved as datasets
                        if verbose:
                            print(f"Creating dataset {attr_name} in group {grp.name}")
                        if isinstance(getattr(self, attr_name), pr.PyRanges):
                            grp.create_dataset(name = attr_name,
                                            data = getattr(self, attr_name).as_df().to_numpy(dtype = 'S'),
                                            compression = compression_type,
                                            compression_opts = compression_opts)
                        elif type(getattr(self, attr_name)) == list:
                            grp.create_dataset(name = attr_name,
                                            data = np.array(getattr(self, attr_name), dtype = 'S'),
                                            compression = compression_type,
                                            compression_opts = compression_opts)
                        else:
                            grp.create_dataset(name = attr_name,
                                            data = getattr(self, attr_name).to_numpy(dtype = 'S'),
                                            compression = compression_type,
                                            compression_opts = compression_opts)
                    else:
                        #attribute is a mappable
                        attr = getattr(self, attr_name)
                        def _walk_down_mappable(x, grp, name = None, verbose = False):
                            if isinstance(x, Mapping):
                                for key in x.keys():
                                    if verbose:
                                        print(f"Creating group {key} in group {grp.name}")
                                    new_grp = grp.create_group(name = key)
                                    _walk_down_mappable(x[key], new_grp, key, verbose)
                            else:
                                if verbose:
                                    print(f"Creating dataset {name} in group {grp.name}")
                                if isinstance(x, pr.PyRanges):
                                    grp.create_dataset(name = name,
                                                    data = x.as_df().to_numpy(dtype = 'S'),
                                                    compression = compression_type,
                                                    compression_opts  = compression_opts)
                                elif type(x) == list:
                                    grp.create_dataset(name = attr_name,
                                            data = np.array(x, dtype = 'S'),
                                            compression = compression_type,
                                            compression_opts = compression_opts)
                                else:
                                    grp.create_dataset(name = name,
                                                    data = x.to_numpy(dtype = 'S'),
                                                    compression = compression_type,
                                                    compression_opts  = compression_opts)
                        if verbose:
                            print(f"Creating group {attr_name} in group {grp.name}")
                        grp = grp.create_group(name = attr_name)
                        _walk_down_mappable(attr, grp, name = attr_name, verbose = verbose) 
        except Exception as e:
            print(e)
        finally:
            #close had5_file if it was opened before
            if not isinstance(grp_or_fname, h5py.Group):
                had5_file.close()

def write_h5ad(ctx_result: Union[cisTarget, dict],
               had5_file_or_file_path: Union[h5py.File, str],
               compression_type = 'gzip',
               compression_opts = 4,
               verbose = False):
    if isinstance(ctx_result, cisTarget):
        ctx_result.to_h5ad( 
            grp_or_fname = had5_file_or_file_path, 
            compression_type = compression_type, 
            compression_opts = compression_opts, 
            verbose = verbose)
    elif type(ctx_result, dict):
        if isinstance(had5_file_or_file_path, h5py.File):
            had5_file = had5_file_or_file_path
        elif type(had5_file_or_file_path) == str:
            had5_file = h5py.File(had5_file_or_file_path, "w")
        try:
            for key in ctx_result.keys():
                if verbose:
                    print(f"Creating group: {key}")
                grp = had5_file.create_group(name = key)
                ctx_result[key].to_h5ad(
                    grp_or_fname = grp,
                    compression_type = compression_type, 
                    compression_opts = compression_opts, 
                    verbose = verbose)
        except Exception as e:
            raise Exception(e)
        finally:
            if type(had5_file_or_file_path) == str:
                had5_file.close()
    else:
        raise ValueError('ctx_result should be an instance of class cisTarget or a dict of instances.')

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
                               motifs_to_use: list = None,
                               **kwargs):
    """
    Run cisTarget.

    Parameters
    ---------
    ctx_db: :class:`cisTargetDatabase`
        A cistarget database object.
    region_sets: Dict
        A dictionary of PyRanges containing region coordinates for the region sets to be analyzed.
    specie: str
        Specie from which genomic coordinates come from
    name: str, optional
        Analysis name. Default: cisTarget
    fraction_overlap: float, optional
        Minimal overlap between query and regions in the database for the mapping. Default: 0.4
    auc_threshold: float, optional
        The fraction of the ranked genome to take into account for the calculation of the
        Area Under the recovery Curve. Default: 0.005
    nes_threshold: float, optional
        The Normalized Enrichment Score (NES) threshold to select enriched features. Default: 3.0
    rank_threshold: float, optional
        The total number of ranked genes to take into account when creating a recovery curve.
        Default: 0.05
    path_to_motif_annotations: str, optional
        Path to motif annotations. If not provided, they will be downloaded from 
        https://resources.aertslab.org based on the specie name provided (only possible for mus_musculus,
        homo_sapiens and drosophila_melanogaster). Default: None
    annotation_version: str, optional
        Motif collection version. Default: v9
    annotation: List, optional
        Annotation to use for forming cistromes. It can be 'Direct_annot' (direct evidence that the motif is 
        linked to that TF), 'Motif_similarity_annot' (based on tomtom motif similarity), 'Orthology_annot'
        (based on orthology with a TF that is directly linked to that motif) or 'Motif_similarity_and_Orthology_annot'.
        Default: ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']
    motif_similarity_fdr: float, optional
        Minimal motif similarity value to consider two motifs similar. Default: 0.001
    orthologous_identity_threshold: float, optional
        Minimal orthology value for considering two TFs orthologous. Default: 0.0
    n_cpu: int, optional
        Number of cores to use. If 1, ray will not be used. Default: 1
    motifs_to_use: List, optional
        A subset of motifs to use for the analysis. Default: None (All)
    **kwargs:
        Extra parameters to pass to `ray.init()`
        
    Return
    ---------
        A dictionary with one :class:`cisTarget` object per analysis.
        
    References
    ---------
    Van de Sande B., Flerin C., et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis.
    Nat Protoc. June 2020:1-30. doi:10.1038/s41596-020-0336-2
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
        sys.stderr = null
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
                                            orthologous_identity_threshold = orthologous_identity_threshold,
                                            motifs_to_use = motifs_to_use) for key in list(region_sets.keys())])
        ray.shutdown()
        sys.stderr = sys.__stderr__
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
                                            orthologous_identity_threshold = orthologous_identity_threshold,
                                            motifs_to_use = motifs_to_use) for key in list(region_sets.keys())]
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
            orthologous_identity_threshold: float = 0.0,
            motifs_to_use: list = None) -> pd.DataFrame:
            
    """
    Internal function to run cistarget in parallel with Ray.
    
    Parameters
    ---------
    ctx_db: :class:`cisTargetDatabase`
        A cistarget database object.
    name: str
        Analysis name
    specie: str
        Specie from which genomic coordinates come from
    auc_threshold: float, optional
        The fraction of the ranked genome to take into account for the calculation of the
        Area Under the recovery Curve. Default: 0.005
    nes_threshold: float, optional
        The Normalized Enrichment Score (NES) threshold to select enriched features. Default: 3.0
    rank_threshold: float, optional
        The total number of ranked genes to take into account when creating a recovery curve.
        Default: 0.05
    path_to_motif_annotations: str, optional
        Path to motif annotations. If not provided, they will be downloaded from 
        https://resources.aertslab.org based on the specie name provided (only possible for mus_musculus,
        homo_sapiens and drosophila_melanogaster). Default: None
    annotation_version: str, optional
        Motif collection version. Default: v9
    annotation: List, optional
        Annotation to use for forming cistromes. It can be 'Direct_annot' (direct evidence that the motif is 
        linked to that TF), 'Motif_similarity_annot' (based on tomtom motif similarity), 'Orthology_annot'
        (based on orthology with a TF that is directly linked to that motif) or 'Motif_similarity_and_Orthology_annot'.
        Default: ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']
    motif_similarity_fdr: float, optional
        Minimal motif similarity value to consider two motifs similar. Default: 0.001
    orthologous_identity_threshold: float, optional
        Minimal orthology value for considering two TFs orthologous. Default: 0.0
    motifs_to_use: List, optional
        A subset of motifs to use for the analysis. Default: None (All)
        
    Return
    ---------
        A :class:`cisTarget` object.
        
    References
    ---------
    Van de Sande B., Flerin C., et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis.
    Nat Protoc. June 2020:1-30. doi:10.1038/s41596-020-0336-2
    """
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
                        orthologous_identity_threshold = orthologous_identity_threshold,
                        motifs_to_use = motifs_to_use)


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
            orthologous_identity_threshold: float = 0.0,
            motifs_to_use: list = None ):
    """
    Internal function to run cistarget.
    
    Parameters
    ---------
    ctx_db: :class:`cisTargetDatabase`
        A cistarget database object.
    name: str
        Analysis name
    specie: str
        Specie from which genomic coordinates come from
    auc_threshold: float, optional
        The fraction of the ranked genome to take into account for the calculation of the
        Area Under the recovery Curve. Default: 0.005
    nes_threshold: float, optional
        The Normalized Enrichment Score (NES) threshold to select enriched features. Default: 3.0
    rank_threshold: float, optional
        The total number of ranked genes to take into account when creating a recovery curve.
        Default: 0.05
    path_to_motif_annotations: str, optional
        Path to motif annotations. If not provided, they will be downloaded from 
        https://resources.aertslab.org based on the specie name provided (only possible for mus_musculus,
        homo_sapiens and drosophila_melanogaster). Default: None
    annotation_version: str, optional
        Motif collection version. Default: v9
    annotation: List, optional
        Annotation to use for forming cistromes. It can be 'Direct_annot' (direct evidence that the motif is 
        linked to that TF), 'Motif_similarity_annot' (based on tomtom motif similarity), 'Orthology_annot'
        (based on orthology with a TF that is directly linked to that motif) or 'Motif_similarity_and_Orthology_annot'.
        Default: ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot']
    motif_similarity_fdr: float, optional
        Minimal motif similarity value to consider two motifs similar. Default: 0.001
    orthologous_identity_threshold: float, optional
        Minimal orthology value for considering two TFs orthologous. Default: 0.0
    motifs_to_use: List, optional
        A subset of motifs to use for the analysis. Default: None (All)
        
    Return
    ---------
        A :class:`cisTarget` object.
        
    References
    ---------
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
                           orthologous_identity_threshold,
                           motifs_to_use)
    return ctx_result
    
## Show results 
def cistarget_results(cistarget_dict,
                    name: Optional[str] = None):
    """
    A function to show cistarget results in jupyter notebooks.
    
    Parameters
    ---------
    cistarget_dict: Dict
        A dictionary with one :class:`cisTarget` object per slot.
    name: str
        Dictionary key of the analysis result to show. Default: None (All)
    """
    motif_enrichment_dict = {key: cistarget_dict[key].motif_enrichment for key in cistarget_dict.keys()}
    if name is None:
        motif_enrichment_table=pd.concat([motif_enrichment_dict[key] for key in motif_enrichment_dict.keys()], axis=0, sort=False)
    else:
        motif_enrichment_table=motif_enrichment_dict[name]
    return HTML(motif_enrichment_table.to_html(escape=False, col_space=80))