from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Callable

import pycistarget

LOGO = """

██████╗ ██╗   ██╗ ██████╗██╗███████╗████████╗ █████╗ ██████╗  ██████╗ ███████╗████████╗
██╔══██╗╚██╗ ██╔╝██╔════╝██║██╔════╝╚══██╔══╝██╔══██╗██╔══██╗██╔════╝ ██╔════╝╚══██╔══╝
██████╔╝ ╚████╔╝ ██║     ██║███████╗   ██║   ███████║██████╔╝██║  ███╗█████╗     ██║
██╔═══╝   ╚██╔╝  ██║     ██║╚════██║   ██║   ██╔══██║██╔══██╗██║   ██║██╔══╝     ██║
██║        ██║   ╚██████╗██║███████║   ██║   ██║  ██║██║  ██║╚██████╔╝███████╗   ██║
╚═╝        ╚═╝    ╚═════╝╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝   ╚═╝

"""

VERSION = pycistarget.__version__
DESCRIPTION = "Motif enrichment analysis."

CISTARGET_CALL_COMMAND = "cistarget"
CISTARGET_COMMAND_DESCRIPTION = "Run motif enrichment analysis using the cisTarget algorithm."
DEM_CALL_COMMAND = "dem"
DEM_COMMAND_DESCRIPTION = (
    "Run motif enrichment analysis using the Differentially Enriched Motif (DEM) algorithm."
)

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
    fraction_overlap_w_dem_database=0.4,
    promoter_space=1_000,
    adjpval_thr=0.05,
    log2fc_thr=1.0,
    mean_fg_thr=0.0,
    annotation_version="v10nr_clust",
    annotations_to_use=["Direct_annot", "Motif_similarity_annot", "Orthology_annot"],
    motif_similarity_fdr=0.001,
    orthologous_identity_threshold=0.0,
    seed=555,
)


def _add_common_arguments(parser: argparse.ArgumentParser, defaults: dict, analysis: str):
    """Arguments shared by the cistarget and dem commands."""
    parser.add_argument(
        "--output_folder",
        dest="output_folder",
        type=str,
        required=True,
        help="Path to the folder in which to write results.",
    )
    parser.add_argument(
        "--species",
        dest="species",
        type=str,
        required=True,
        help="Species used for the analysis. This parameter is used to download the correct "
        "motif-to-TF annotations from the cisTarget webservers.",
    )
    parser.add_argument(
        "--path_to_motif_annotations",
        dest="path_to_motif_annotations",
        type=str,
        default=None,
        help="Path to the motif-to-TF annotations. "
        "By default this will be downloaded from the cisTarget webservers.",
    )
    parser.add_argument(
        "--annotation_version",
        dest="annotation_version",
        type=str,
        default=defaults["annotation_version"],
        help="Version of the motif-to-TF annotation to use. This parameter is used to download "
        "the correct motif-to-TF data from the cisTarget webservers. Defaults to: %(default)s",
    )
    parser.add_argument(
        "--motif_similarity_fdr",
        dest="motif_similarity_fdr",
        type=float,
        default=defaults["motif_similarity_fdr"],
        help="Threshold on motif similarity scores for calling similar motifs. "
        "Defaults to: %(default)s",
    )
    parser.add_argument(
        "--orthologous_identity_threshold",
        dest="orthologous_identity_threshold",
        type=float,
        default=defaults["orthologous_identity_threshold"],
        help="Threshold on the protein-protein orthology score for calling orthologous motifs. "
        "Defaults to: %(default)s",
    )
    parser.add_argument(
        "--annotations_to_use",
        dest="annotations_to_use",
        type=str,
        nargs="*",
        default=defaults["annotations_to_use"],
        help="Which annotations to use for annotation motifs to TFs. Defaults to: "
        + " ".join(defaults["annotations_to_use"]),
    )
    parser.add_argument(
        "--name",
        dest="name",
        type=str,
        default=None,
        help="Name of this analysis. This name is appended to the output file name. "
        f"By default the file name of the bed file{'s' if analysis == 'dem' else ''} is used.",
    )
    parser.add_argument(
        "--output_mode",
        dest="output_mode",
        choices=["tsv", "hdf5", "hdf5+"],
        default="hdf5",
        help="Specifies how the results will be saved to disk. "
        "tsv: tab-separated text file containing which motifs are enriched. "
        "hdf5: a new hdf5 file will be created containing the full results. "
        "hdf5+: an existing hdf5 file will be appended with the full motifs enrichment results. "
        "Defaults to 'hdf5'",
    )
    parser.add_argument(
        "--write_html",
        dest="write_html",
        action="store_true",
        help="Whether or not to save the results as an html file.",
    )


def add_parser_for_motif_enrichment_cistarget(subparser: argparse._SubParsersAction):
    parser: argparse.ArgumentParser = subparser.add_parser(
        name=CISTARGET_CALL_COMMAND,
        add_help=True,
        description=CISTARGET_COMMAND_DESCRIPTION,
    )

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

    parser.set_defaults(func=cistarget)
    parser.add_argument(
        "--cistarget_db_fname",
        dest="cistarget_db_fname",
        type=str,
        required=True,
        help="Path to the cisTarget rankings database (.regions_vs_motifs.rankings.feather).",
    )
    parser.add_argument(
        "--bed_fname",
        dest="bed_fname",
        type=str,
        nargs="+",
        required=True,
        help="Path to bed file on which to run motif enrichment analysis. Several bed files can "
        "be given: each one is analysed separately (one output file per bed file).",
    )
    parser.add_argument(
        "--fr_overlap_w_ctx_db",
        dest="fraction_overlap_w_cistarget_database",
        type=float,
        default=CISTARGET_DEFAULTS["fraction_overlap_w_cistarget_database"],
        help="Fraction of nucleotides, of regions in the bed file, that should overlap with "
        "regions in the cistarget database in order for them to be included in the analysis. "
        "Defaults to: %(default)s",
    )
    parser.add_argument(
        "--auc_threshold",
        dest="auc_threshold",
        type=float,
        default=CISTARGET_DEFAULTS["auc_threshold"],
        help="Threshold on the AUC value for calling significant motifs. Defaults to: %(default)s",
    )
    parser.add_argument(
        "--nes_threshold",
        dest="nes_threshold",
        type=float,
        default=CISTARGET_DEFAULTS["nes_threshold"],
        help="Threshold on the NES value for calling significant motifs. "
        "NES - Normalised Enrichment Score - is defined as (AUC - Avg(AUC)) / sd(AUC). "
        "Defaults to: %(default)s",
    )
    parser.add_argument(
        "--rank_threshold",
        dest="rank_threshold",
        type=float,
        default=CISTARGET_DEFAULTS["rank_threshold"],
        help="The fraction of ranked regions to take into account when creating recovery curves. "
        "Defaults to: %(default)s",
    )
    _add_common_arguments(parser, CISTARGET_DEFAULTS, analysis="cistarget")


def add_parser_for_motif_enrichment_dem(subparser: argparse._SubParsersAction):
    parser: argparse.ArgumentParser = subparser.add_parser(
        name=DEM_CALL_COMMAND,
        add_help=True,
        description=DEM_COMMAND_DESCRIPTION,
    )

    def dem(arg):
        from pycistarget.cli.commands import run_dem_command, run_dem_one_vs_all_command

        if arg.one_vs_all_beds is not None:
            if arg.foreground_beds is not None or arg.background_beds is not None:
                parser.error(
                    "--one_vs_all_beds can not be combined with --foreground_beds or --background_beds"
                )
            if arg.name is not None:
                parser.error("--name can not be used with --one_vs_all_beds")
            if len(arg.one_vs_all_beds) < 2:
                parser.error("--one_vs_all_beds needs at least two bed files")
            dem_inputs = dict(paths_to_bed=arg.one_vs_all_beds)
            run: Callable[..., None] = run_dem_one_vs_all_command
        else:
            if arg.foreground_beds is None or arg.background_beds is None:
                parser.error(
                    "either --foreground_beds and --background_beds, or --one_vs_all_beds is required"
                )
            dem_inputs = dict(
                paths_to_foreground_bed=arg.foreground_beds,
                paths_to_background_bed=arg.background_beds,
                name=arg.name,
            )
            run = run_dem_command

        run(
            dem_db_fname=arg.dem_db_fname,
            **dem_inputs,
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
            output_mode=arg.output_mode,
            write_html=arg.write_html,
        )

    parser.set_defaults(func=dem)
    parser.add_argument(
        "--dem_db_fname",
        dest="dem_db_fname",
        type=str,
        required=True,
        help="Path to the DEM score database (.regions_vs_motifs.scores.feather).",
    )
    parser.add_argument(
        "--foreground_beds",
        dest="foreground_beds",
        nargs="+",
        type=str,
        default=None,
        help="Path(s) to bed file(s) to use as foreground regions.",
    )
    parser.add_argument(
        "--background_beds",
        dest="background_beds",
        nargs="+",
        type=str,
        default=None,
        help="Path(s) to bed file(s) to use as background regions.",
    )
    parser.add_argument(
        "--one_vs_all_beds",
        dest="one_vs_all_beds",
        nargs="+",
        type=str,
        default=None,
        help="Instead of --foreground_beds and --background_beds: run DEM for each of these bed "
        "files (foreground) against all other bed files (background), reading the database once. "
        "The analysis of x.bed is named x_vs_all.",
    )
    parser.add_argument(
        "--fraction_overlap_w_dem_database",
        dest="fraction_overlap_w_dem_database",
        type=float,
        default=DEM_DEFAULTS["fraction_overlap_w_dem_database"],
        help="Fraction of nucleotides, of regions in the bed file, that should overlap with "
        "regions in the dem database in order for them to be included in the analysis. "
        "Defaults to: %(default)s",
    )
    parser.add_argument(
        "--max_bg_regions",
        dest="max_bg_regions",
        type=int,
        default=None,
        help="Maximum number of regions to use as background. Defaults to None (i.e. use all regions)",
    )
    parser.add_argument(
        "--genome_annotation",
        dest="genome_annotation",
        type=str,
        default=None,
        help="Path to genome annotation. This parameter is required when "
        "--balance_number_of_promoters is set. Defaults to None.",
    )
    parser.add_argument(
        "--balance_number_of_promoters",
        dest="balance_number_of_promoters",
        action="store_true",
        help="Set this flag to balance the number of promoter regions in fore- and background. "
        "When this is set a genome annotation must be provided using the --genome_annotation parameter.",
    )
    parser.add_argument(
        "--promoter_space",
        dest="promoter_space",
        type=int,
        default=DEM_DEFAULTS["promoter_space"],
        help="Number of basepairs up- and downstream of the TSS that are considered as being "
        "the promoter for that gene. Defaults to: %(default)s",
    )
    parser.add_argument(
        "--adjpval_thr",
        dest="adjpval_thr",
        type=float,
        default=DEM_DEFAULTS["adjpval_thr"],
        help="Threshold on the Benjamini-Hochberg adjusted p-value from the Wilcoxon test performed "
        "on the motif score of foreground vs background regions for a motif to be considered as "
        "enriched. Defaults to: %(default)s",
    )
    parser.add_argument(
        "--log2fc_thr",
        dest="log2fc_thr",
        type=float,
        default=DEM_DEFAULTS["log2fc_thr"],
        help="Threshold on the log2 fold change of the motif score of foreground vs background "
        "regions for a motif to be considered as enriched. Defaults to: %(default)s",
    )
    parser.add_argument(
        "--mean_fg_thr",
        dest="mean_fg_thr",
        type=float,
        default=DEM_DEFAULTS["mean_fg_thr"],
        help="Minimal mean signal in the foreground to consider a motif enriched. "
        "Defaults to: %(default)s",
    )
    parser.add_argument(
        "--motif_hit_thr",
        dest="motif_hit_thr",
        type=float,
        default=None,
        help="Minimal CRM score to consider a region enriched for a motif. "
        "Default: None (It will be automatically calculated based on precision-recall).",
    )
    parser.add_argument(
        "--seed",
        dest="seed",
        type=int,
        default=DEM_DEFAULTS["seed"],
        help="Random seed used for sampling background regions (if --max_bg_regions is set or "
        "--balance_number_of_promoters is used). Defaults to: %(default)s",
    )
    _add_common_arguments(parser, DEM_DEFAULTS, analysis="dem")


def create_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("--version", action="version", version=f"pycistarget {VERSION}")
    subparsers = parser.add_subparsers()
    add_parser_for_motif_enrichment_cistarget(subparsers)
    add_parser_for_motif_enrichment_dem(subparsers)
    return parser


def main(argv=None) -> int:
    parser = create_argument_parser()
    args = parser.parse_args(args=argv)

    if not hasattr(args, "func"):
        # pycistarget has been called without any command: print help.
        print(LOGO)
        print(f"pycistarget version: {VERSION}")
        parser.print_help()
    else:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
            handlers=[logging.StreamHandler(stream=sys.stdout)],
        )
        args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
