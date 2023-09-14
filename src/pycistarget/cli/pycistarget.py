import sys
import argparse
import pycistarget

LOGO = """

██████╗ ██╗   ██╗ ██████╗██╗███████╗████████╗ █████╗ ██████╗  ██████╗ ███████╗████████╗
██╔══██╗╚██╗ ██╔╝██╔════╝██║██╔════╝╚══██╔══╝██╔══██╗██╔══██╗██╔════╝ ██╔════╝╚══██╔══╝
██████╔╝ ╚████╔╝ ██║     ██║███████╗   ██║   ███████║██████╔╝██║  ███╗█████╗     ██║   
██╔═══╝   ╚██╔╝  ██║     ██║╚════██║   ██║   ██╔══██║██╔══██╗██║   ██║██╔══╝     ██║   
██║        ██║   ╚██████╗██║███████║   ██║   ██║  ██║██║  ██║╚██████╔╝███████╗   ██║   
╚═╝        ╚═╝    ╚═════╝╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝   ╚═╝   
                                                                                       
"""

# General constants for argument parser

VERSION = pycistarget.__version__
DESCRIPTION = "Motif enrichment analysis."

# Constants for cistarget argument parser

CISTARGET_CALL_COMMAND = "cistarget"
CISTARGET_COMMAND_DESCRIPTION = """
Run motif enrichment analysis using the cisTarget algorithm.
"""
DEM_CALL_COMMAND = "dem"
DEM_COMMAND_DESCRIPTION = """
Run motif enrichment analysis using the Differentially Enriched Motif (DEM) algorithm
"""

# cistarget default values
CISTARGET_DEFAULTS = dict(
    fraction_overlap_w_cistarget_database=0.4,
    auc_threshold=0.005,
    nes_threshold=3.0,
    rank_threshold=0.05,
    annotation_version="v10nr_clust",
    annotations_to_use=["Direct_annot", "Motif_similarity_annot", "Orthology_annot"],
    motif_similarity_fdr=0.001,
    orthologous_identity_threshold=0.0,
)

DEM_DEFAULTS = dict(
    fraction_overlap_w_dem_database= 0.4,
    promoter_space=1_000,
    adjpval_thr=0.05,
    log2fc_thr=1.0,
    mean_fg_thr=0.0,
    annotation_version='v10nr_clust',
    annotations_to_use=["Direct_annot", "Motif_similarity_annot", "Orthology_annot"],
    motif_similarity_fdr=0.001,
    orthologous_identity_threshold=0.0,
    seed=555)


def add_parser_for_motif_enrichment_cistarget(subparser: argparse._SubParsersAction):
    # Set up new argument parser
    parser: argparse.ArgumentParser = subparser.add_parser(
        name=CISTARGET_CALL_COMMAND,
        add_help=True,
        description=CISTARGET_COMMAND_DESCRIPTION,
    )

    # Define function to be run by this parser
    def cistarget(arg):
        from pycistarget.cli.commands import run_cistarget_command

        run_cistarget_command(
            cistarget_db_fname=arg.cistarget_db_fname,
            bed_fname=arg.bed_fname,
            output_folder=arg.output_folder,
            fraction_overlap_w_cistarget_database=arg.fraction_overlap_w_cistarget_database,
            auc_threshold=arg.auc_threshold,
            nes_threshold=arg.nes_threshold,
            rank_threshold=arg.rank_threshold,
            path_to_motif_annotations=arg.path_to_motif_annotations,
            annotation_version=arg.annotation_version,
            motif_similarity_fdr=arg.motif_similarity_fdr,
            orthologous_identity_threshold=arg.orthologous_identity_threshold,
            species=arg.species,
            annotations_to_use=arg.annotations_to_use,
            name=arg.name,
            output_mode=arg.output_mode,
            write_html=arg.write_html,
        )

    # Register the function with the parser
    parser.set_defaults(func=cistarget)
    # Define required Arguments
    parser.add_argument(
        "--cistarget_db_fname",
        dest="cistarget_db_fname",
        action="store",
        type=str,
        required=True,
        help="Path to the cisTarget rankings database (.regions_vs_motifs.rankings.feather).",
    )
    parser.add_argument(
        "--bed_fname",
        dest="bed_fname",
        action="store",
        type=str,
        required=True,
        help="Path to bed file on which to run motif enrichment analysis.",
    )
    parser.add_argument(
        "--output_folder",
        dest="output_folder",
        action="store",
        type=str,
        required=True,
        help="Path to the folder in which to write results.",
    )
    parser.add_argument(
        "--species",
        dest="species",
        action="store",
        type=str,
        required=True,
        help="""
            Species used for the analysis. This parameter is used to download the correct
            motif-to-TF annotations from the cisTarget webservers.""",
    )

    # Define optional arguments
    parser.add_argument(
        "--fr_overlap_w_ctx_db",
        dest="fraction_overlap_w_cistarget_database",
        action="store",
        type=float,
        required=False,
        help="""
            Fraction of nucleotides, of regions in the bed file,
            that should overlap with regions in the cistarget database
            in order for them to be included in the analysis.
            Defaults to: """ + str(CISTARGET_DEFAULTS["fraction_overlap_w_cistarget_database"]),
        default=CISTARGET_DEFAULTS["fraction_overlap_w_cistarget_database"],
    )
    parser.add_argument(
        "--auc_threshold",
        dest="auc_threshold",
        action="store",
        type=float,
        required=False,
        help="""
            Threshold on the AUC value for calling significant motifs.
            Defaults to: """+ str(CISTARGET_DEFAULTS["auc_threshold"]),
        default=CISTARGET_DEFAULTS["auc_threshold"],
    )
    parser.add_argument(
        "--nes_threshold",
        dest="nes_threshold",
        action="store",
        type=float,
        required=False,
        help="""
            Threshold on the NES value for calling significant motifs.
            NES - Normalised Enrichment Score - is defined as (AUC - Avg(AUC)) / sd(AUC).
            Defaults to: """ + str(CISTARGET_DEFAULTS["nes_threshold"]),
        default=CISTARGET_DEFAULTS["nes_threshold"],
    )
    parser.add_argument(
        "--rank_threshold",
        dest="rank_threshold",
        action="store",
        type=float,
        required=False,
        help="""
            The total number of ranked regions to take into account when creating a recovery curves.
            Defaults to: """ + str(CISTARGET_DEFAULTS["rank_threshold"]),
        default=CISTARGET_DEFAULTS["rank_threshold"],
    )
    parser.add_argument(
        "--path_to_motif_annotations",
        dest="path_to_motif_annotations",
        action="store",
        type=str,
        required=False,
        help="""
            Path to the motif-to-TF annotations.
            By default this will be downloaded from the cisTarget webservers.""",
        default=None,
    )
    parser.add_argument(
        "--annotation_version",
        dest="annotation_version",
        action="store",
        type=str,
        required=False,
        help="""
            Version of the motif-to-TF annotation to use. This parameter is used
            to download the correct motif-to-TF data from the cisTarget webservers.
            Defaults to: """ + str(CISTARGET_DEFAULTS["annotation_version"]),
        default=CISTARGET_DEFAULTS["annotation_version"],
    )
    parser.add_argument(
        "--motif_similarity_fdr",
        dest="motif_similarity_fdr",
        action="store",
        type=float,
        required=False,
        help=""""
            Threshold on motif similarity scores for calling similar motifs.
            Defaults to: """ + str(CISTARGET_DEFAULTS["motif_similarity_fdr"]),
        default=CISTARGET_DEFAULTS["motif_similarity_fdr"],
    )
    parser.add_argument(
        "--orthologous_identity_threshold",
        dest="orthologous_identity_threshold",
        action="store",
        type=float,
        required=False,
        help="""
            Threshold on the protein-protein orthology score for calling orthologous motifs.
            Defaults to: """ + str(CISTARGET_DEFAULTS["orthologous_identity_threshold"]),
        default=CISTARGET_DEFAULTS["orthologous_identity_threshold"],
    )
    parser.add_argument(
        "--annotations_to_use",
        dest="annotations_to_use",
        action="store",
        type=str,
        required=False,
        nargs="*",
        help="""
            Which annotations to use for annotation motifs to TFs.
            Defaults to: """ + ' '.join(CISTARGET_DEFAULTS["annotations_to_use"]),
        default=CISTARGET_DEFAULTS["annotations_to_use"],
    )
    parser.add_argument(
        "--name",
        dest="name",
        action="store",
        type=str,
        required=False,
        help="""
            Name of this analysis. This name is appended to the output file name.
            By default the file name of the bed file is used.""",
        default=None,
    )
    parser.add_argument(
        "--output_mode",
        dest="output_mode",
        action="store",
        choices=["tsv", "hdf5", "hdf5+"],
        required=False,
        help="""
            Specifies how the results will be saved to disk.
                tsv:    tab-seperated text file containing which motifs are enriched.
                hdf5:   a new hdf5 file will be created containing the full results.
                hdf5+:  an existing hdf5 file will be appended with the full motifs
                        enrichment results.
            Defaults to 'hdf5'""",
        default="hdf5"
    )
    parser.add_argument(
        "--write_html",
        dest="write_html",
        action="store_true",
        help="Wether or not to save the results as an html file.",
    )

def add_parser_for_motif_enrichment_dem(subparser: argparse._SubParsersAction):
    # Set up new argument parser
    parser: argparse.ArgumentParser = subparser.add_parser(
        name=DEM_CALL_COMMAND,
        add_help=True,
        description=DEM_CALL_COMMAND,
    )

    # Define function to be run by this parser
    def dem(arg):
        from pycistarget.cli.commands import run_dem_command
        run_dem_command(
            dem_db_fname=arg.dem_db_fname,
            paths_to_foreground_bed=arg.foreground_beds,
            paths_to_background_bed=arg.background_beds,
            output_folder=arg.output_folder,
            species=arg.species,
            fraction_overlap_w_dem_database=arg.fraction_overlap_w_dem_database,
            max_bg_regions=arg.max_bg_regions,
            path_to_genome_annotation=arg.genome_annotation,
            balance_number_of_promoters=arg.balance_number_of_promoters,
            promoter_space=arg.promoter_space,
            adjpval_thr=arg.adjpval_thr,
            log2fc_thr=arg.log2fc_thr,
            mean_fg_thr=arg.mean_fg_thr,
            motif_hit_thr=arg.motif_hit_thr,
            path_to_motif_annotations=arg.path_to_motif_annotations,
            annotation_version=arg.annotation_version,
            annotations_to_use=arg.annotations_to_use,
            motif_similarity_fdr=arg.motif_similarity_fdr,
            orthologous_identity_threshold=arg.orthologous_identity_threshold,
            seed=arg.seed,
            name=arg.name,
            output_mode=arg.output_mode,
            write_html=arg.write_html
        )

    # Register the function with the parser
    parser.set_defaults(func=dem)
    # Define required Arguments
    parser.add_argument(
        "--dem_db_fname",
        dest="dem_db_fname",
        action="store",
        type=str,
        required=True,
        help="Path to the DEM score database (.regions_vs_motifs.scores.feather).",
    )
    parser.add_argument(
        "--foreground_beds",
        dest="foreground_beds",
        action="store",
        nargs="+",
        type=str,
        required=True,
        help="Path(s) to bed file(s) to use as foreground regions."
    )
    parser.add_argument(
        "--background_beds",
        dest="background_beds",
        action="store",
        nargs="+",
        type=str,
        required=True,
        help="Path(s) to bed file(s) to use as background regions."
    )
    parser.add_argument(
        "--output_folder",
        dest="output_folder",
        action="store",
        type=str,
        required=True,
        help="Path to the folder in which to write results.",
    )
    parser.add_argument(
        "--species",
        dest="species",
        action="store",
        type=str,
        required=True,
        help="""
            Species used for the analysis. This parameter is used to download the correct
            motif-to-TF annotations from the cisTarget webservers.""",
    )

    # Define optional arguments
    parser.add_argument(
        "--fraction_overlap_w_dem_database",
        dest="fraction_overlap_w_dem_database",
        action="store",
        type=float,
        required=False,
        help="""
            Fraction of nucleotides, of regions in the bed file,
            that should overlap with regions in the dem database
            in order for them to be included in the analysis.
            Defaults to: """ + str(DEM_DEFAULTS["fraction_overlap_w_dem_database"]),
        default=DEM_DEFAULTS["fraction_overlap_w_dem_database"],
    )
    parser.add_argument(
        "--max_bg_regions",
        dest="max_bg_regions",
        action="store",
        type=int,
        required=False,
        help="""
            Maximum number of regions to use as background.
            Defaults to None (i.e. use all regions)""",
        default=None,
    )
    parser.add_argument(
        "--genome_annotation",
        dest="genome_annotation",
        action="store",
        type=str,
        required=False,
        help="""
        Path to genome annotation. 
        This parameter is required whe balance_number_of_promoters is set.
        Defaults to None.
        """,
        default=None
    )
    parser.add_argument(
        "--balance_number_of_promoters",
        dest="balance_number_of_promoters",
        action="store_true",
        help="""
        Set this flag to balance the number of promoter regions in fore- and background.
        When this is set a genome annotation must be provided using the 
        --genome_annotation parameter.
        """
    )
    parser.add_argument(
        "--promoter_space",
        dest="promoter_space",
        action="store",
        type=int,
        required=False,
        help="""
        Number of basepairs up- and downstream of the TSS that are considered as being
        the promoter for that gene.
        Defaults to: """ + str(DEM_DEFAULTS["promoter_space"]),
        default=DEM_DEFAULTS["promoter_space"]
    )
    parser.add_argument(
        "--adjpval_thr",
        dest="adjpval_thr",
        action="store",
        type=float,
        required=False,
        help="""
        Threshold on the Benjamini-Hochberg adjusted p-value
        from the Wilcoxon test performed on the motif score of foreground
        vs background regions for a motif to be considered as enriched.
        Defaults to: """ + str(DEM_DEFAULTS["adjpval_thr"]),
        default=DEM_DEFAULTS["adjpval_thr"]
    )
    parser.add_argument(
        "--log2fc_thr",
        dest="log2fc_thr",
        action="store",
        type=float,
        required=False,
        help="""
        Threshold on the log2 fold change of the motif score of foreground
        vs background regions for a motif to be considered as enriched.
        Defaults to: """ + str(DEM_DEFAULTS["log2fc_thr"]),
        default=DEM_DEFAULTS["log2fc_thr"]
    )
    parser.add_argument(
        "--mean_fg_thr",
        dest="mean_fg_thr",
        action="store",
        type=float,
        required=False,
        help="""
        Minimul mean signal in the foreground to consider a motif enriched.
        Defaults to: """ + str(DEM_DEFAULTS["mean_fg_thr"]),
        default=DEM_DEFAULTS["mean_fg_thr"]
    )
    parser.add_argument(
        "--motif_hit_thr",
        dest="motif_hit_thr",
        action="store",
        type=float,
        required=False,
        help="""
        Minimal CRM score to consider a region enriched for a motif. 
        Default: None (It will be automatically calculated based on precision-recall).
        """,
        default=None
    )
    parser.add_argument(
        "--path_to_motif_annotations",
        dest="path_to_motif_annotations",
        action="store",
        type=str,
        required=False,
        help="""
            Path to the motif-to-TF annotations.
            By default this will be downloaded from the cisTarget webservers.""",
        default=None,
    )
    parser.add_argument(
        "--annotation_version",
        dest="annotation_version",
        action="store",
        type=str,
        required=False,
        help="""
            Version of the motif-to-TF annotation to use. This parameter is used
            to download the correct motif-to-TF data from the cisTarget webservers.
            Defaults to: """ + str(DEM_DEFAULTS["annotation_version"]),
        default=DEM_DEFAULTS["annotation_version"],
    )
    parser.add_argument(
        "--motif_similarity_fdr",
        dest="motif_similarity_fdr",
        action="store",
        type=float,
        required=False,
        help=""""
            Threshold on motif similarity scores for calling similar motifs.
            Defaults to: """ + str(DEM_DEFAULTS["motif_similarity_fdr"]),
        default=DEM_DEFAULTS["motif_similarity_fdr"],
    )
    parser.add_argument(
        "--orthologous_identity_threshold",
        dest="orthologous_identity_threshold",
        action="store",
        type=float,
        required=False,
        help="""
            Threshold on the protein-protein orthology score for calling orthologous motifs.
            Defaults to: """ + str(DEM_DEFAULTS["orthologous_identity_threshold"]),
        default=DEM_DEFAULTS["orthologous_identity_threshold"],
    )
    parser.add_argument(
        "--annotations_to_use",
        dest="annotations_to_use",
        action="store",
        type=str,
        required=False,
        nargs="*",
        help="""
            Which annotations to use for annotation motifs to TFs.
            Defaults to: """ + ' '.join(DEM_DEFAULTS["annotations_to_use"]),
        default=DEM_DEFAULTS["annotations_to_use"],
    )
    parser.add_argument(
        "--name",
        dest="name",
        action="store",
        type=str,
        required=False,
        help="""
            Name of this analysis. This name is appended to the output file name.
            By default the file name of the bed file is used.""",
        default=None,
    )
    parser.add_argument(
        "--output_mode",
        dest="output_mode",
        action="store",
        choices=["tsv", "hdf5", "hdf5+"],
        required=False,
        help="""
            Specifies how the results will be saved to disk.
                tsv:    tab-seperated text file containing which motifs are enriched.
                hdf5:   a new hdf5 file will be created containing the full results.
                hdf5+:  an existing hdf5 file will be appended with the full motifs
                        enrichment results.
            Defaults to 'hdf5'""",
        default="hdf5"
    )
    parser.add_argument(
        "--write_html",
        dest="write_html",
        action="store_true",
        help="Wether or not to save the results as an html file.",
    )
    parser.add_argument(
        "--seed",
        dest="seed",
        action="store",
        type=int,
        required=False,
        help="""
        Random seed use for sampling background regions (if max_bg_regions is not None)
        Defaults to: """ + str(DEM_DEFAULTS["seed"]),
        default=DEM_DEFAULTS["seed"]
    )


def create_argument_parser() -> argparse.ArgumentParser:
    # Set up argument parser
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    subparsers = parser.add_subparsers()

    # Create parser for motif enrichment cistarget
    add_parser_for_motif_enrichment_cistarget(subparsers)

    # Create parser for motif enrichment dem
    add_parser_for_motif_enrichment_dem(subparsers)
    return parser


def main(argv=None) -> int:
    # Parse command line arguments
    parser = create_argument_parser()
    args = parser.parse_args(args=argv)

    if not hasattr(args, "func"):
        # pycistarget has been called without any command.
        # In this case the help is printed
        print(LOGO)
        print(f"pycistarget version: {VERSION}")
        parser.print_help()
    else:
        # pycistarget has been called with a command.
        # In this case run the command.
        args.func(args)
        
    return 0

if __name__ == "__main__":
    sys.exit(main())
