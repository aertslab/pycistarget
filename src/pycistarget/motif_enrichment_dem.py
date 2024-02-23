from ctxcore.rnkdb import FeatherRankingDatabase
from ctxcore.genesig import GeneSignature
import logging
import numpy as np
import pandas as pd
import pyranges as pr
import random
from sklearn.metrics import roc_curve
import numba
import sys
from typing import Optional, Tuple, Literal, List
from pycistarget.utils import (
    target_to_query,
    coord_to_region_names,
    region_names_to_coordinates)
from pycistarget.motif_enrichment_result import MotifEnrichmentResult
import math

# TODO: Add these generic functions to another package
@numba.njit(parallel=True)
def mean_axis1(arr):
    """
    Calculate column wise mean of 2D-numpy matrix with numba, mimicking `np.mean(x, axis=1)`.

    Parameters
    ----------
    arr
        2D-numpy array to calculate the mean per column for.
    """

    mean_axis1_array = np.empty(arr.shape[0], dtype=np.float64)
    for i in numba.prange(arr.shape[0]):
        mean_axis1_array[i] = np.mean(arr[i, :])
    return mean_axis1_array

@numba.njit
def get_log2_fc(fg_mat, bg_mat):
    """
    Calculate log2 fold change between foreground and background matrix.

    Parameters
    ----------
    fg_mat
        2D-numpy foreground matrix.
    bg_mat
        2D-numpy background matrix.
    """

    if fg_mat.shape[0] != bg_mat.shape[0]:
        _ERROR_MESSAGE =  "Foreground matrix and background matrix have a different first dimensions!"
        raise ValueError(_ERROR_MESSAGE)

    # Calculate log2 fold change between foreground and background matrix with numba in
    # a similar way as the following numpy code:
    #    np.log2(
    #        (np.mean(fg_mat, axis=1) + 10**-12) / (np.mean(bg_mat, axis=1) + 10**-12)
    #    )
    return np.log2(
        (mean_axis1(fg_mat) + 10**-12) / (mean_axis1(bg_mat) + 10**-12)
    )

@numba.jit(nopython=True)
def rankdata_average_numba(arr: np.ndarray):
    """
    based on code from scipy
    """
    sorter = np.argsort(arr, kind = "quicksort")
    inv = np.empty(sorter.size, dtype = np.intp)
    inv[sorter] = np.arange(sorter.size, dtype=np.intp)
    arr = arr[sorter]
    obs = np.empty(arr.shape, dtype=np.intp)
    obs[0] = True
    obs[1:] = arr[1:] != arr[:-1]
    dense = obs.cumsum()[inv]
    non_zero = np.nonzero(obs)[0]
    count = np.empty(non_zero.shape[0] + 1)
    count[0:-1] = non_zero
    count[-1] = len(obs)
    result = .5 * (count[dense] + count[dense - 1] + 1)
    return result

@numba.jit(nopython=True)
def norm_sf(z):
    return ((1.0 + math.erf(-z / math.sqrt(2.0))) / 2.0)

@numba.jit(nopython=True)
def ranksums_numba(x: np.ndarray, y:np.ndarray):
    """
    based on code from scipy
    """
    n1 = len(x)
    n2 = len(y)
    alldata = np.concatenate((x, y))
    ranked = rankdata_average_numba(alldata)
    x = ranked[:n1]
    s = np.sum(x, axis=0)
    expected = n1 * (n1+n2+1) / 2.0
    z = (s - expected) / np.sqrt(n1*n2*(n1+n2+1)/12.0)
    prob = 2 * norm_sf(abs(z))
    return z, prob

@numba.jit(nopython=True, parallel = True)
def ranksums_numba_multiple(X: np.ndarray, Y: np.ndarray):
    n = X.shape[0]
    if Y.shape[0] != n:
        raise ValueError("X and Y should have the same shape on dimension 0")
    ranksums_z = np.zeros((n), dtype=np.float64)
    ranksums_p = np.zeros((n), dtype=np.float64)
    for i in numba.prange(n):
        z, p = ranksums_numba(X[i], Y[i])
        ranksums_z[i] = z
        ranksums_p[i] = p
    return ranksums_z, ranksums_p

def p_adjust_bh(p: np.ndarray):
    """
    Benjamini-Hochberg p-value correction for multiple hypothesis testing.
    """
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

def get_optimal_threshold_roc(
        foreground_scores_motif:np.ndarray,
        background_scores_motif:np.ndarray) -> float:
    foreground_scores_motif = foreground_scores_motif[foreground_scores_motif > 0]
    background_scores_motif = background_scores_motif[background_scores_motif > 0]
    all_scores = np.zeros(
        foreground_scores_motif.shape[0] + background_scores_motif.shape[0],
        dtype = float)
    all_scores[0:foreground_scores_motif.shape[0]] = foreground_scores_motif
    all_scores[foreground_scores_motif.shape[0]:] = background_scores_motif
    labels = np.zeros_like(all_scores, dtype = int)
    labels[0:foreground_scores_motif.shape[0]] = 1
    labels[foreground_scores_motif.shape[0]:] = 0
    fpr, tpr, thresholds = roc_curve(labels, all_scores)
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = thresholds[optimal_idx]
    return optimal_threshold

class DEMDatabase(FeatherRankingDatabase):
    def __init__(self, fname, fraction_overlap: float = 0.4):
        super().__init__(fname = fname, name = "DEM")
        self.db_regions = pr.PyRanges(region_names_to_coordinates(list(self.genes)))
        self.fraction_overlap = fraction_overlap
    def get_scores(
            self,
            regions: pr.PyRanges) -> Tuple[pd.DataFrame, pd.DataFrame]:
        target_to_db = target_to_query(
            target = regions,
            query = self.db_regions,
            fraction_overlap = self.fraction_overlap)
        db_regions_to_load = list(set(target_to_db["Query"]))
        db_scores = self.load(
            GeneSignature(
                name = "DEM",
                gene2weight=db_regions_to_load))
        return target_to_db, db_scores

class DEM(MotifEnrichmentResult):
    def __init__(
        self,
        foreground_regions: pr.PyRanges,
        background_regions: pr.PyRanges,
        name: str,
        species: Literal[
                "homo_sapiens", "mus_musculus", "drosophila_melanogaster"],
        adjpval_thr: float = 0.05,
        log2fc_thr: float = 1.0,
        mean_fg_thr: float = 0.0,
        motif_hit_thr: Optional[float] = None,
        path_to_motif_annotations: Optional[str] = None,
        annotation_version: str = 'v10nr_clust',
        annotation_to_use: list = ['Direct_annot', 'Motif_similarity_annot', 'Orthology_annot', 'Motif_similarity_and_Orthology_annot'],
        motif_similarity_fdr: float = 0.001,
        orthologous_identity_threshold: float = 0.0,
        motifs_to_use: Optional[list] = None):
        self.foreground_regions = foreground_regions
        self.background_regions = background_regions
        self.adjpval_thr = adjpval_thr
        self.log2fc_thr = log2fc_thr
        self.mean_fg_thr = mean_fg_thr
        self.motif_hit_thr = motif_hit_thr
        super().__init__(
            name = name,
            species = species,
            path_to_motif_annotations = path_to_motif_annotations,
            annotation_version = annotation_version,
            annotation_to_use = annotation_to_use,
            motif_similarity_fdr = motif_similarity_fdr,
            orthologous_identity_threshold = orthologous_identity_threshold,
            motifs_to_use = motifs_to_use)
    
    def run(self, dem_db: DEMDatabase):
        # Create logger
        level    = logging.INFO
        format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
        handlers = [logging.StreamHandler(stream=sys.stdout)]
        logging.basicConfig(level = level, format = format, handlers = handlers)
        log = logging.getLogger('DEM')

        log.info(f"Running DEM for {self.name}")
        foreground_region_to_db, foreground_scores = dem_db.get_scores(
            self.foreground_regions)
        background_region_to_db, background_scores = dem_db.get_scores(
            self.background_regions)
        self.regions_to_db = pd.concat(
            [foreground_region_to_db, background_region_to_db]) \
            .drop_duplicates() \
            .reset_index(drop = True)
        motif_names = foreground_scores.index
        background_scores = background_scores.loc[motif_names]
        foreground_scores_arr = np.array(foreground_scores)
        background_scores_arr = np.array(background_scores)
        
        # Perform Wilcoxon rank sum test
        stats, pvalues = ranksums_numba_multiple(
            X = foreground_scores_arr,
            Y = background_scores_arr)
        
        # Calculate log2FC
        logFC = get_log2_fc(
            fg_mat = foreground_scores_arr,
            bg_mat = background_scores_arr
        )
        # pvalue correction
        pvalues_adj = p_adjust_bh(pvalues)

        # Create result dataframe
        result = pd.DataFrame(
            data = {
                "Log2FC": logFC,
                "Adjusted_pval": pvalues_adj,
                "Mean_fg": mean_axis1(foreground_scores_arr),
                "Mean_bg": mean_axis1(background_scores_arr)},
            index = motif_names)
        
        # Threshold dataframe
        result = result.loc[
            np.logical_and.reduce(
                (
                    result["Adjusted_pval"] <= self.adjpval_thr,
                    result["Log2FC"]        >= self.log2fc_thr,
                    result["Mean_fg"]       >= self.mean_fg_thr
                )
            )
        ]

        result = result.sort_values([
            "Log2FC", "Adjusted_pval"], 
            ascending = [False, True])
        
        self.motif_enrichment = result
        log.info("Adding motif-to-TF annotation")
        self.add_motif_annotation()

        # Get motif hits
        result["Motif_hit_thr"] = None
        result["Motif_hits"] = None
        significant_motifs = result.index
        motif_hits = {}
        for motif in significant_motifs:
            if self.motif_hit_thr is None:
                thr = get_optimal_threshold_roc(
                    foreground_scores.loc[motif].values,
                    background_scores.loc[motif].values,)
            else:
                thr = self.motif_hit_thr 
            motif_hits[motif] = foreground_scores.loc[motif].loc[
                    foreground_scores.loc[motif] > thr].index.tolist()
            result.loc[motif, "Motif_hit_thr"] = thr
            result.loc[motif, "Motif_hits"] = len(motif_hits[motif])
        
        rs_motif_hits = {
            motif: list(set(
                self.regions_to_db.loc[
                    self.regions_to_db["Query"].isin(motif_hits[motif]),
                    "Target"].tolist()))
            for motif in motif_hits.keys()}
        self.motif_hits = {'Database': motif_hits, 'Region_set': rs_motif_hits}
        self.get_cistromes()

def get_foreground_and_background_regions(
        foreground_region_sets: List[pr.PyRanges],
        background_region_sets: List[pr.PyRanges],
        max_bg_regions: Optional[int] = None,
        genome_annotation: Optional[pd.DataFrame] = None,
        balance_number_of_promoters: bool = True,
        promoter_space: int = 1_000,
        seed: int = 555
) -> Tuple[pr.PyRanges, pr.PyRanges]:
    # Aggregate foreground and background regions into single dataframe and convert to region names
    foreground_region_names = list(set(coord_to_region_names(
        pd.concat(map(lambda ranges: ranges.df, foreground_region_sets)))))
    background_region_names = list(set(coord_to_region_names(
        pd.concat(map(lambda ranges: ranges.df, background_region_sets)))))
    
    background_region_names = list(
        set(background_region_names) - set(foreground_region_names))
    if max_bg_regions is not None:
        max_bg_regions = min(max_bg_regions, len(background_region_names))

    if balance_number_of_promoters:
        if genome_annotation is None:
            raise ValueError("Please provide genome_annotation if you wish to balance the number of promoters")
        tmp_pr_foreground = pr.PyRanges(region_names_to_coordinates(foreground_region_names))
        tmp_pr_background = pr.PyRanges(region_names_to_coordinates(background_region_names))
        
        # Define promoter locations
        promoter_annotation = genome_annotation.copy()
        promoter_annotation["End"]   = promoter_annotation["Start"] + promoter_space
        promoter_annotation["Start"] = promoter_annotation["Start"] - promoter_space

        # Calculate fraction of promoters in foreground
        pr_promoter_annotation = pr.PyRanges(
            promoter_annotation[["Chromosome", "Start", "End"]].drop_duplicates())
        nr_promoter_in_fg = len(
            tmp_pr_foreground \
                .count_overlaps(pr_promoter_annotation) \
                .df.query("NumberOverlaps > 0")
        )
        fraction_promoter_fg = nr_promoter_in_fg / len(foreground_region_names)
        
        # Calculate number of allowed promoter regions in background
        if max_bg_regions is not None:
            nr_allowed_promoters_background = int(max_bg_regions * fraction_promoter_fg)
        else:
            nr_allowed_promoters_background = int(len(background_region_names) * fraction_promoter_fg)

        # Seperate promoter and non promoter regions in background
        tmp_background_promoter_overlap = tmp_pr_background \
            .count_overlaps(pr_promoter_annotation).df
        background_region_names_no_promoter = coord_to_region_names(
            tmp_background_promoter_overlap.query("NumberOverlaps == 0")
        )
        background_region_names_promoter    = coord_to_region_names(
            tmp_background_promoter_overlap.query("NumberOverlaps > 0")
        )
        
        # Sample promoter regions from background
        nr_allowed_promoters_background = min(
            len(background_region_names_promoter), nr_allowed_promoters_background)
        random.Random(seed).shuffle(background_region_names_promoter)
        background_region_names_promoter = background_region_names_promoter[
            0:nr_allowed_promoters_background]
        
        # Sample non promoter regions
        if max_bg_regions is not None:
            nr_non_promoter_regions_background = max_bg_regions - nr_allowed_promoters_background
        else:
            nr_non_promoter_regions_background = len(background_region_names) - nr_allowed_promoters_background
        
        # Sample non promoter regions from background
        random.Random(seed).shuffle(background_region_names_no_promoter)
        background_region_names_no_promoter = background_region_names_no_promoter[
            0:nr_non_promoter_regions_background]
        
        background_region_names = [
            *background_region_names_promoter, *background_region_names_no_promoter]
        
    else:
        if max_bg_regions is not None:
            random.Random(seed).shuffle(background_region_names)
            background_region_names[0:max_bg_regions]
    foreground_regions = pr.PyRanges(region_names_to_coordinates(foreground_region_names))
    background_regions = pr.PyRanges(region_names_to_coordinates(background_region_names))
    return foreground_regions, background_regions