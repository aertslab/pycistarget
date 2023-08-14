from pycistarget.motif_enrichment_cistarget import (
    cisTarget,
    cisTargetDatabase,
)
from pycistarget.input_output import write_cistarget
import pyranges as pr
from typing import List, Optional, Literal
import os


def run_cistarget_command(
    cistarget_db_fname: str,
    bed_fname: str,
    output_folder: str,
    fraction_overlap_w_cistarget_database: float,
    auc_threshold: float,
    nes_threshold: float,
    rank_threshold: float,
    path_to_motif_annotations: str,
    annotation_version: str,
    motif_similarity_fdr: float,
    orthologous_identity_threshold: float,
    species: str,
    annotations_to_use: List[str],
    name: Optional[str] = None,
    output_mode: Literal["tsv", "hdf5", "hdf5+"] = "hdf5",
    write_html: bool = True,
):
    """
    Run pycisTarget on region set

    Parameters
    ----------
    cistarget_db_fname: str
        Path to the cisTarget database.
    bed_fname: str
        Path to a bed file on which to run motif enrichment analysis.
    output_folder: str
        Path to the fodler where to store the results.
    fraction_overlap_w_cistarget_database: float
        Fraction of nucleotides, of regions in the bed file,
        that should overlap with regions in the cistarget database
        in order for them to be included in the analysis.
    auc_threshold: float
        Threshold on the AUC value for calling significant motifs.
    nes_threshold: float
        Threshold on the NES value for calling significant motifs.
        NES - Normalised Enrichment Score - is defined as (AUC - Avg(AUC)) / sd(AUC).
    rank_threshold: float
        The total number of ranked regions to take into account when creating a recovery curves.
    path_to_motif_annotations: str
        Path to the motif-to-TF annotations. Set to None to download these from the
        cisTarget webservers.
    annotation_version: str
        Version of the motif-to-TF annotation to use. This parameter is used
        to download the correct motif-to-TF data from the cisTarget webservers.
    motig_similarity_fdr: float
        Threshold on motif similarity scores for calling similar motifs.
    orthologous_identity_threshold: float
        Threshold on the protein-protein orthology score for calling orthologous motifs.
    species: str
        Species used for the analysis. This parameter is used to download the correct
        motif-to-TF annotations from the cisTarget webservers.
    annotations_to_use: List[str]
        Which annotations to use for annotation motifs to TFs.
    name: str, optional
        Name of this analysis. This name is appended to the output file name.
    output_mode: Literal["tsv", "hdf5", "hdf5+"]
        Specifies how the results will be saved to disk.
            tsv:    tab-seperated text file containing which motifs are enriched.
            hdf5:   a new hdf5 file will be created containing the full results.
            hdf5+:  an existing hdf5 file will be appended with the full motifs
                    enrichment results.
        Defaults to 'hdf5'
    write_html: bool
        Wether or not to save the results as an html file.
        Defaults to True
    """
    if name is None:
        # Get name from file name
        name = os.path.basename(bed_fname).replace(".bed", "").replace(".gz", "")
    regions: pr.PyRanges = pr.read_bed(bed_fname, as_df=False)
    # Read cisTarget database
    ctx_db = cisTargetDatabase(
        fname=cistarget_db_fname,
        region_sets=regions,
        name=name,
        fraction_overlap=fraction_overlap_w_cistarget_database,
    )
    # Set up cisTarget analysis
    cistarget_result = cisTarget(
        region_set=regions,
        name=name,
        specie=species,
        auc_threshold=auc_threshold,
        nes_threshold=nes_threshold,
        rank_threshold=rank_threshold,
        path_to_motif_annotations=path_to_motif_annotations,
        annotation_version=annotation_version,
        annotation=annotations_to_use,
        motif_similarity_fdr=motif_similarity_fdr,
        orthologous_identity_threshold=orthologous_identity_threshold,
    )
    # Run analysis
    cistarget_result.run_ctx(ctx_db)

    # Write results to file
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Write html report
    if write_html:
        cistarget_result.motif_enrichment.to_html(
            buf=os.path.join(output_folder, f"motif_enrichment_cistarget_{name}.html"),
            escape=False,
            col_space=80,
        )

    if output_mode == "tsv":
        cistarget_result.motif_enrichment.to_csv(
            path_or_buf=os.path.join(
                output_folder, f"motif_enrichment_cistarget_{name}.tsv"
            ),
            sep="\t",
            header=True,
            index=True,
        )
    elif output_mode == "hdf5":
        write_cistarget(
            cistarget=cistarget_result,
            path=os.path.join(
                output_folder, f"motif_enrichment_cistarget_{name}.hdf5"
            ),
            mode="w",
        )
    elif output_mode == "hdf5+":
        write_cistarget(
            cistarget=cistarget_result,
            path=os.path.join(
                output_folder, f"motif_enrichment_cistarget_{name}.hdf5"
            ),
            mode="a",
        )
    else:
        raise ValueError(f"Output mode: {output_mode} is not supported!")
