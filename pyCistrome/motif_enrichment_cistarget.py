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

def overlap_regions_w_cistarget_db( ctx_db: cisTargetDatabase, 
                                    input_regions: Union[pr.PyRanges, Dict[str, pr.PyRanges]], 
                                    min_fraction_overlap: float = 0.4,
                                    n_cpu: int = 5) -> Union[pr.PyRanges, Dict[str, pr.PyRanges]]:
    """
    Overlap regions from the cistarget database with input regions taking into account a minimal fractional overlap.
    :param input_regions: a single PyRanges object or a dictionary of PyRanges objects.
    :param min_fraction_overlap: minimum fraction of overlap needed between a cistarget region and an input region.
    :param n_cpy: number of threads to use.
    :return: regions from the cistarget database which overlap with the input regions 
    """
    #convert ctx_regions name to PyRanges object
    pr_ctx_regions = pr.PyRanges(utils.region_names_to_coordinates(ctx_db.ctx_regions))

    if isinstance(input_regions, pr.PyRanges):
        #overlap ctx regions with topic regions taking into account fraction of overlap
        pr_ctx_input_overlap = pr_ctx_regions.coverage( other               = input_regions,
                                                        strandedness        = False,
                                                        keep_nonoverlapping = False,              #Eventhough keep_nonoverlapping = False, still non-overlapping are returned ... bug?
                                                        overlap_col         = 'NumberOverlaps',
                                                        fraction_col        = 'FractionOverlaps',
                                                        nb_cpu              = n_cpu)
        
        #subset so only regions with a minimal fractional overlap of min_fraction_overlap are kept
        pr_ctx_input_overlap = pr_ctx_input_overlap.subset(
            lambda region: region.FractionOverlaps >= min_fraction_overlap,
            strand = False,
            nb_cpu = n_cpu)
        
        #calculate the reverse overlap to see how many regions of the input regions overlap with ctx regions
        pr_input_ctx_overlap = input_regions.coverage(  other               = pr_ctx_regions,
                                                        strandedness        = False,
                                                        keep_nonoverlapping = False,              #Eventhough keep_nonoverlapping = False, still non-overlapping are returned ... bug?
                                                        overlap_col         = 'NumberOverlaps',
                                                        fraction_col        = 'FractionOverlaps',
                                                        nb_cpu              = n_cpu)
        
        #subset so only regions with a minimal fractional overlap of min_fraction_overlap are kept
        pr_input_ctx_overlap = pr_input_ctx_overlap.subset(
            lambda region: region.FractionOverlaps >= min_fraction_overlap,
            strand = False,
            nb_cpu = n_cpu)
        
        fraction_input_regions_ctx_overlap = len(pr_input_ctx_overlap) / len(input_regions)

        pr_ctx_input_overlap = {'input': pr_ctx_input_overlap}
        fraction_input_regions_in_database = {'input': fraction_input_regions_ctx_overlap}
        
    elif isinstance(input_regions, dict):
        if not all([isinstance(v, pr.PyRanges) for v in input_regions.values()]):
            raise ValueError('A dictionary with PyRanges objects as values has to be provided for input_regions!')
        pr_ctx_input_overlap = {}
        fraction_input_regions_in_database = {}
        #generate progress bar
        pbar = tqdm(input_regions.keys(), total = len(input_regions.keys()))
        for region_set_name in pbar:
            #set progress bar description
            pbar.set_description(region_set_name, refresh=True)

            #overlap ctx regions with topic regions taking into account fraction of overlap
            pr_ctx_input_overlap[region_set_name] = pr_ctx_regions.coverage( 
                                                            other               = input_regions[region_set_name],
                                                            strandedness        = False,
                                                            keep_nonoverlapping = False,              #Eventhough keep_nonoverlapping = False, still non-overlapping are returned ... bug?
                                                            overlap_col         = 'NumberOverlaps',
                                                            fraction_col        = 'FractionOverlaps',
                                                            nb_cpu              = n_cpu)
            #subset so only regions with a minimal fractional overlap of min_fraction_overlap are kept
            pr_ctx_input_overlap[region_set_name] = pr_ctx_input_overlap[region_set_name].subset(
                lambda region: region.FractionOverlaps >= min_fraction_overlap,
                strand = False,
                nb_cpu = n_cpu)
            #calculate the reverse overlap to see how many regions of the input regions overlap with ctx regions
            pr_input_ctx_overlap = input_regions[region_set_name].coverage(  
                                                            other               = pr_ctx_regions,
                                                            strandedness        = False,
                                                            keep_nonoverlapping = False,              #Eventhough keep_nonoverlapping = False, still non-overlapping are returned ... bug?
                                                            overlap_col         = 'NumberOverlaps',
                                                            fraction_col        = 'FractionOverlaps',
                                                            nb_cpu              = n_cpu)
            #subset so only regions with a minimal fractional overlap of min_fraction_overlap are kept
            pr_input_ctx_overlap = pr_input_ctx_overlap.subset(
                lambda region: region.FractionOverlaps >= min_fraction_overlap,
                strand = False,
                nb_cpu = n_cpu)
            fraction_input_regions_in_database[region_set_name] = len(pr_input_ctx_overlap) / len(input_regions[region_set_name])

    else:
        raise ValueError('input_regions has to be a single PyRanges object or a dictionary of PyRanges objects!')
    
    #report fraction of input regions in cistarget database
    for region_set_name in fraction_input_regions_in_database.keys():
        print('{}: {} % of regions overlap with regions in the cistarget database with a minmal fractional overlap of {}.'.
        format(region_set_name,
                str(np.round(fraction_input_regions_in_database[region_set_name], 3) * 100),
                str(min_fraction_overlap)))
    
    return pr_ctx_input_overlap

def annotate_features(features_df: pd.DataFrame, motif_annotations_fname: str, index_col: str = 'MotifID') -> pd.DataFrame:
    """
    annotates enriched featues
    :param features_df: pandas dataframe with enriched features
    :param motif_annotations_fname: path to motif2tf dataframe
    :param index_col: column name to use for indexing
    return pandas dataframe with annoted features
    """
    feature_annotations = utils.load_motif_annotations(motif_annotations_fname = motif_annotations_fname, index_col = index_col)
    annotated = pd.merge(features_df, feature_annotations, left_on = index_col, right_on = index_col, how = 'left')
    annotated['gene_name'] = annotated['gene_name'].fillna('None')
    annotated['Annotation'] = annotated['Annotation'].fillna('None')
    annotated = annotated.drop(labels = ['MotifSimilarityQvalue', 'OrthologousIdentity'], axis = 1).groupby(index_col).agg(
        {
            'NES'       : lambda x: x.unique(),
            'AUC'       : lambda x: x.unique(),
            'gene_name' : lambda x: ';'.join(x.unique()),
            'Annotation': lambda x: ';'.join(x.unique())
        }
    )
    return annotated

def find_enriched_features( ctx_db: cisTargetDatabase,
                            region_set_signature: Regulon,
                            motif_annotations_fname: str,
                            weighted_recovery: bool = False,
                            auc_threshold: float = 0.005,
                            nes_threshold: float = 3.0,
                            rank_threshold: int = 20000,
                            get_cistrome: bool = True) -> pd.DataFrame:
    """
    Finds features of which the rankings are enriched in the input region set
    :param ctx_db: cistarget database object
    :param region_set_signature: gene signature object containing the input region set, regions should be present in cistarget database
    :param weighted_recovery: wether or not to use weighted recovery in the analysis
    :param auc_threshold: the cut-off for fraction of ranked genomic regions at which to calculate AUC. This measure is then used for comparing all the features.
    :param nes_threshold: only the enriched regions with normalized enrichment score (NES) higher than the threshold will be returned
    :param rank_threshold: The total number of ranked regions to take into account when creating a recovery curve
    :param motif_annotations_fname: path to motif2tf dataframe
    :param get_cistrome: when True returns for each motif sets of regions where that motif is enriched
    :return: a pandas dataframe with enriched features
    REFERENCES:
    ----------
    Van de Sande B., Flerin C., et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis.
    Nat Protoc. June 2020:1-30. doi:10.1038/s41596-020-0336-2
    """
    if not isinstance(ctx_db, cisTargetDatabase):
        raise ValueError('ctx_db should be a cisTargetDatabase!')
    if not isinstance(region_set_signature, Regulon):
        raise ValueError('region_set_signature should be a gene signature object!')
    
    #hardcoded values
    COLUMN_NAME_NES = "NES"
    COLUMN_NAME_AUC = "AUC"
    COLUMN_NAME_TF = "TF"
    COLUMN_NAME_MOTIF_ID = "MotifID"
    COLUMN_NAME_TARGET_GENES = "TargetGenes"
    COLUMN_NAME_RANK_AT_MAX = "RankAtMax"
    COLUMN_NAME_TYPE = "Type"

    if not all(np.isin(region_set_signature.genes, ctx_db.ctx_regions)):
        raise ValueError('Not all input regions are in the cistarget database, overlap your input regions first!')
    
    #load rankings for input regions
    print('Loading cistarget database ... ')
    ctx_db.load_rankings_for_input_regions(region_set_signature.genes)

    #assertion just for debuging
    assert all(np.isin(region_set_signature.genes, ctx_db.rankings_db.columns.values)) and len(region_set_signature.genes) == len(ctx_db.rankings_db.columns.values), 'all regions in input should be in the cistarget database after loading!'
    
    features, regions, rankings = ctx_db.rankings_db.index.values, ctx_db.rankings_db.columns.values, ctx_db.rankings_db.values
    weights = np.asarray([region_set_signature[region] for region in regions]) if weighted_recovery else np.ones(len(regions))

    # Calculate recovery curves, AUC and NES values.
    print('Calculating AUC ...')
    aucs = calc_aucs(ctx_db.rankings_db, ctx_db.total_regions, weights, auc_threshold)
    ness = (aucs - aucs.mean()) / aucs.std()

    # Keep only features that are enriched, i.e. NES sufficiently high.
    enriched_features_idx = ness >= nes_threshold

    enriched_features = pd.DataFrame(index=pd.MultiIndex.from_tuples(list(zip(repeat(region_set_signature.transcription_factor),
                                                                              features[enriched_features_idx])),
                                                                     names=[COLUMN_NAME_TF, COLUMN_NAME_MOTIF_ID]),
                                data={COLUMN_NAME_NES: ness[enriched_features_idx],
                                      COLUMN_NAME_AUC: aucs[enriched_features_idx]})
    
    #annotate features
    print('Annotating features ...')
    enriched_features = annotate_features(enriched_features, motif_annotations_fname)

    if get_cistrome:
        print('Recovery analysis ...')
        rccs, _ = recovery(ctx_db.rankings_db, ctx_db.total_regions, weights, rank_threshold, auc_threshold, no_auc=True)  
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
        #delete unecessary data
        print('Cleaning up ...')
        del enriched_features['Ranking']
    return enriched_features




