"""Motif enrichment analysis with the cisTarget algorithm (recovery of regions in motif rankings)."""

from __future__ import annotations

import logging
from os import PathLike

import numpy as np
import pandas as pd
import pyranges as pr

from pycistarget.database import FeatherCisTargetDatabase
from pycistarget.motif_enrichment_result import MotifEnrichmentResult, Species
from pycistarget.recovery import calculate_aucs, derive_rank_cutoff, leading_edge, recovery_curves
from pycistarget.utils import coord_to_region_names, region_set_motif_hits, target_to_query

log = logging.getLogger(__name__)


class cisTargetDatabase:
    """
    cisTarget rankings database.

    Contains a dataframe with motifs as rows, regions as columns and ranks as values,
    for the regions of the database that overlap with the region sets of interest.
    For more information on how to generate databases, please visit:
    https://github.com/aertslab/create_cisTarget_databases

    Attributes
    ----------
    regions_to_db: pd.DataFrame or dict of pd.DataFrame
        The mapping between query regions ("Target") and regions in the database ("Query").
        A dictionary with one mapping per region set when ``region_sets`` is a dictionary.
    db_rankings: pd.DataFrame
        A dataframe with motifs as rows, regions as columns and rank as values.
    total_regions: int
        Total number of regions in the database
    """

    def __init__(
        self,
        fname: str | PathLike,
        region_sets: dict[str, pr.PyRanges] | pr.PyRanges | None = None,
        name: str | None = None,
        fraction_overlap: float = 0.4,
    ):
        """
        Parameters
        ----------
        fname
            Path to feather file containing the cisTarget rankings database (regions_vs_motifs).
        region_sets
            Dictionary of PyRanges or a single PyRanges with the regions that are going to be analyzed.
            When None, the complete database is loaded.
        name
            Unused, kept for backwards compatibility.
        fraction_overlap
            Minimal overlap between query regions and regions in the database for the mapping.
        """
        self.regions_to_db, self.db_rankings, self.total_regions = self.load_db(
            fname, region_sets, name, fraction_overlap
        )

    def load_db(
        self,
        fname: str | PathLike,
        region_sets: dict[str, pr.PyRanges] | pr.PyRanges | None = None,
        name: str | None = None,
        fraction_overlap: float = 0.4,
    ) -> tuple[dict[str, pd.DataFrame] | pd.DataFrame | None, pd.DataFrame, int]:
        log.info("Reading cisTarget database")
        db = FeatherCisTargetDatabase(fname)
        db_regions = list(db.region_ids)
        # Region IDs can be prefixed, e.g. "hg38__chr1:100-600".
        prefix = db_regions[0].split("__")[0] if "__" in db_regions[0] else None
        if prefix is not None:
            db_regions = [x.split("__")[1] for x in db_regions]

        regions_to_db: dict[str, pd.DataFrame] | pd.DataFrame | None
        if region_sets is None:
            log.warning(
                "Loading complete cistarget database, this can take a long time and consumes a lot of memory!"
            )
            regions_to_db = None
            db_rankings = db.load_full()
        else:
            if isinstance(region_sets, dict):
                regions_to_db = {
                    key: target_to_query(regions, db_regions, fraction_overlap=fraction_overlap)
                    for key, regions in region_sets.items()
                }
                regions_in_db = set().union(*(m["Query"] for m in regions_to_db.values()))
            elif isinstance(region_sets, pr.PyRanges):
                regions_to_db = target_to_query(
                    region_sets, db_regions, fraction_overlap=fraction_overlap
                )
                regions_to_db.index = regions_to_db["Target"]
                regions_in_db = set(regions_to_db["Query"])
            else:
                raise ValueError(
                    "region_sets should be either a dict of PyRanges objects or a single PyRanges object, "
                    f"not {type(region_sets)}"
                )
            if prefix is not None:
                regions_in_db = {f"{prefix}__{x}" for x in regions_in_db}
            db_rankings = db.load(regions_in_db)
        if prefix is not None:
            db_rankings.columns = [x.split("__")[1] for x in db_rankings.columns]
        return regions_to_db, db_rankings, db.total_regions


class cisTarget(MotifEnrichmentResult):
    """
    Motif enrichment analysis on a set of regions with the cisTarget algorithm.

    Attributes
    ----------
    regions_to_db: pd.DataFrame
        A dataframe containing the mapping between query regions and regions in the database.
    region_set: pr.PyRanges
        A PyRanges containing region coordinates for the regions to be analyzed.
    name: str
        Analysis name
    species: str
        Species from which genomic coordinates come from
    auc_threshold: float
        The fraction of the ranked genome to take into account for the calculation of the
        Area Under the recovery Curve.
    nes_threshold: float
        The Normalized Enrichment Score (NES) threshold to select enriched features.
    rank_threshold: float
        The fraction of the ranked genome to take into account when creating a recovery curve.
    path_to_motif_annotations: str
        Path to motif annotations. If not provided, they will be downloaded from
        https://resources.aertslab.org based on the species name provided (only possible for mus_musculus,
        homo_sapiens and drosophila_melanogaster).
    annotation_version: str
        Motif collection version.
    annotation_to_use: List
        Annotation to use for forming cistromes. It can be 'Direct_annot' (direct evidence that the motif is
        linked to that TF), 'Motif_similarity_annot' (based on tomtom motif similarity), 'Orthology_annot'
        (based on orthology with a TF that is directly linked to that motif) or 'Motif_similarity_and_Orthology_annot'.
    motif_similarity_fdr: float
        Minimal motif similarity value to consider two motifs similar.
    orthologous_identity_threshold: float
        Minimal orthology value for considering two TFs orthologous.
    motifs_to_use: List
        A subset of motifs to use for the analysis (None: all).
    force: bool
        If the motif-to-TF annotations can not be loaded (e.g. no internet connection to download
        them), continue without annotations and cistromes (with a warning) instead of failing.
    motif_enrichment: pd.DataFrame
        A dataframe containing motif enrichment results
    motif_hits: Dict
        A dictionary containing regions that are considered enriched for each motif.
    cistromes: Dict
        A dictionary containing TF cistromes. Cistromes with no extension contain regions linked to directly
        annotated motifs, while '_extended' cistromes can contain regions linked to motifs annotated by
        similarity or orthology.

    References
    ----------
    Van de Sande B., Flerin C., et al. A scalable SCENIC workflow for single-cell gene regulatory network analysis.
    Nat Protoc. June 2020:1-30. doi:10.1038/s41596-020-0336-2
    """

    _analysis_type = "cisTarget"

    def __init__(
        self,
        region_set: pr.PyRanges,
        name: str,
        species: Species,
        auc_threshold: float = 0.005,
        nes_threshold: float = 3.0,
        rank_threshold: float = 0.05,
        path_to_motif_annotations: str | None = None,
        annotation_version: str = "v10nr_clust",
        annotation_to_use: list[str] | None = None,
        motif_similarity_fdr: float = 0.001,
        orthologous_identity_threshold: float = 0.0,
        motifs_to_use: list[str] | None = None,
        force: bool = False,
    ):
        self.region_set = region_set
        self.auc_threshold = auc_threshold
        self.nes_threshold = nes_threshold
        self.rank_threshold = rank_threshold
        super().__init__(
            name=name,
            species=species,
            path_to_motif_annotations=path_to_motif_annotations,
            annotation_version=annotation_version,
            annotation_to_use=(
                annotation_to_use
                if annotation_to_use is not None
                else [
                    "Direct_annot",
                    "Motif_similarity_annot",
                    "Orthology_annot",
                    "Motif_similarity_and_Orthology_annot",
                ]
            ),
            motif_similarity_fdr=motif_similarity_fdr,
            orthologous_identity_threshold=orthologous_identity_threshold,
            motifs_to_use=motifs_to_use,
            force=force,
        )

    def _set_empty_result(self):
        self.motif_enrichment = pd.DataFrame(
            data={
                "Logo": [],
                "Region_set": [],
                "Direct_annot": [],
                "Motif_similarity_annot": [],
                "Orthology_annot": [],
                "Motif_similarity_and_Orthology_annot": [],
                "NES": [],
                "AUC": [],
                "Rank_at_max": [],
            }
        )
        self.motif_hits = {"Database": {}, "Region_set": {}}
        self.cistromes = {"Database": {}, "Region_set": {}}

    def run_ctx(self, ctx_db: cisTargetDatabase):
        """
        Run cisTarget.

        Parameters
        ----------
        ctx_db
            A cistarget database object, loaded for (at least) the regions in ``region_set``.
        """
        if ctx_db.regions_to_db is None:
            raise ValueError(
                "The cisTarget database was loaded without region sets, "
                "load it with the region set(s) to analyze."
            )
        if isinstance(ctx_db.regions_to_db, dict):
            self.regions_to_db = ctx_db.regions_to_db[self.name]
        else:
            region_names = set(coord_to_region_names(self.region_set))
            self.regions_to_db = ctx_db.regions_to_db.loc[
                sorted(region_names & set(ctx_db.regions_to_db["Target"]))
            ]
        regions = np.array(list(dict.fromkeys(self.regions_to_db["Query"])), dtype=object)
        log.info(f"Running cisTarget for {self.name} which has {len(regions)} regions")

        db_rankings = ctx_db.db_rankings
        if self.motifs_to_use is not None:
            log.info("Using only user provided motifs")
            motifs_not_in_db = set(self.motifs_to_use) - set(db_rankings.index)
            if len(motifs_not_in_db) > 0:
                log.info(
                    "Some motifs provided by the parameter <motifs_to_use> are not in the "
                    f"rankings database: {motifs_not_in_db}"
                )
            db_rankings = db_rankings[db_rankings.index.isin(self.motifs_to_use)]
        motifs = db_rankings.index.to_numpy()
        rankings = db_rankings[regions].to_numpy()

        total_regions = ctx_db.total_regions
        rank_threshold = int(self.rank_threshold * total_regions)
        if derive_rank_cutoff(self.auc_threshold, total_regions) > rank_threshold:
            raise ValueError("auc_threshold must be smaller than or equal to rank_threshold.")

        # AUC and NES per motif.
        aucs = calculate_aucs(rankings, total_regions, self.auc_threshold)
        ness = (aucs - aucs.mean()) / aucs.std()
        enriched = ness >= self.nes_threshold
        if enriched.sum() == 0:
            log.info(f"No enriched motifs found for {self.name}")
            self._set_empty_result()
            return

        # Leading edge of the recovery curves of enriched motifs, compared to the
        # average + 2 standard deviations of the recovery curves of all motifs.
        rccs = recovery_curves(rankings, rank_threshold)
        avg2stdrcc = rccs.mean(axis=0) + 2.0 * rccs.std(axis=0)
        db_motif_hits: dict[str, list[str]] = {}
        rank_at_max: list[int] = []
        for motif, rcc, ranking in zip(
            motifs[enriched], rccs[enriched], rankings[enriched], strict=True
        ):
            db_motif_hits[motif], motif_rank_at_max = leading_edge(
                rcc, avg2stdrcc, ranking, regions
            )
            rank_at_max.append(motif_rank_at_max)

        motif_enrichment = pd.DataFrame(
            {
                "Region_set": self.name,
                "NES": ness[enriched],
                "AUC": aucs[enriched],
                "Rank_at_max": rank_at_max,
            },
            index=pd.Index(motifs[enriched], dtype=object),
        )
        self.motif_enrichment = motif_enrichment.sort_values("NES", ascending=False)
        log.info(f"Annotating motifs for {self.name}")
        self.add_motif_annotation()

        # Motif hits, in database coordinates and in region set coordinates.
        rs_motif_hits = region_set_motif_hits(self.regions_to_db, db_motif_hits)
        self.motif_hits = {"Database": db_motif_hits, "Region_set": rs_motif_hits}
        self.motif_enrichment["Motif_hits"] = [
            len(db_motif_hits[motif]) for motif in self.motif_enrichment.index
        ]
        log.info(f"Getting cistromes for {self.name}")
        self.get_cistromes()
