from pycistarget.utils import load_motif_annotations, get_cistromes_per_region_set
from typing import Optional, List, Dict, Literal
import pandas as pd
import ssl
from IPython.display import HTML
ssl._create_default_https_context = ssl._create_unverified_context
pd.set_option('display.max_colwidth', None)

class MotifEnrichmentResult:
    def __init__(
            self,
            name: str,
            species: Literal[
                "homo_sapiens", "mus_musculus", "drosophila_melanogaster"],
            path_to_motif_annotations: Optional[str] = None,
            annotation_version: str = "v10nr_clust",
            annotation_to_use: List[str] = [
                'Direct_annot',  'Orthology_annot'],
            motif_similarity_fdr: float = 0.001,
            orthologous_identity_threshold: float = 0.0,
            motifs_to_use: Optional[List[str]] = None):
        self.name = name
        self.species = species
        self.path_to_motif_annotations = path_to_motif_annotations
        self.annotation_version = annotation_version
        self.annotation_to_use = annotation_to_use
        self.motif_similarity_fdr = motif_similarity_fdr
        self.orthologous_identity_threshold = orthologous_identity_threshold
        self.motifs_to_use = motifs_to_use
        self.logo_url = f"https://motifcollections.aertslab.org/{self.annotation_version}/logos/"
        self.motif_enrichment = pd.DataFrame()
        self.motif_hits: Dict[str, Dict[str, List[str]]] = {
            "Database": dict(), "Region_set": dict()}
        self.cistromes: Dict[str, Dict[str, List[str]]] = {
            "Database": dict(), "Region_set": dict()}

    def add_motif_annotation(
            self,
            add_logo: bool = True):
        try:
            annot_df = load_motif_annotations(
                self.species,
                version = self.annotation_version,
                fname=self.path_to_motif_annotations,
                motif_similarity_fdr = self.motif_similarity_fdr,
                orthologous_identity_threshold = self.orthologous_identity_threshold)
            annot_df = annot_df[self.annotation_to_use]
            motif_enrichment_w_annot = self.motif_enrichment.merge(
                annot_df, how = "left", left_index = True, right_index = True)
        except:
            Warning("Unable to load motif-to-TF annotations!")
            motif_enrichment_w_annot = self.motif_enrichment
        if add_logo:
            motif_enrichment_w_annot['Logo']=[
                '<img src="' + self.logo_url + motif_name + '.png' + '" width="200" >'
                for motif_name in motif_enrichment_w_annot.index]
            # Put the logo at the first position
            motif_enrichment_w_annot = motif_enrichment_w_annot[
                [
                    *motif_enrichment_w_annot.columns[-1:],
                    *motif_enrichment_w_annot.columns[:-1]
                ]
            ]
        self.motif_enrichment = motif_enrichment_w_annot
    
    def get_cistromes(self):
        cistromes_db = get_cistromes_per_region_set(
            self.motif_enrichment, 
            self.motif_hits['Database'], 
            self.annotation_to_use)
        cistromes_rs = get_cistromes_per_region_set(
            self.motif_enrichment,
            self.motif_hits['Region_set'],
            self.annotation_to_use)
        self.cistromes = {"Database": cistromes_db, "Region_set": cistromes_rs}
    
    def show_result(self):
        return HTML(self.motif_enrichment.to_html(escape = False, col_space = 80))