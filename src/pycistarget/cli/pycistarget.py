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
CREATE_DB_CALL_COMMAND = "create_db"
CREATE_DB_COMMAND_DESCRIPTION = (
    "Create cisTarget motif databases (scores and rankings) from a FASTA file of regions and "
    "Cluster-Buster motif files. Gives the same databases as create_cistarget_motif_databases.py "
    "from create_cisTarget_databases, but faster."
)
PLOT_CALL_COMMAND = "plot"
PLOT_COMMAND_DESCRIPTION = 'Plot motif logos of enriched motifs with their scores (needs: pip install "pycistarget[plot]").'

# Species for which motif-to-TF annotations can be downloaded (see pycistarget.utils).
SPECIES = ["homo_sapiens", "mus_musculus", "drosophila_melanogaster"]

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
        metavar="SPECIES",
        help="Species used for the analysis. This parameter is used to download the correct "
        f"motif-to-TF annotations from the cisTarget webservers. Available: {', '.join(SPECIES)}. "
        "For other species, provide the annotations with --path_to_motif_annotations.",
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
    parser.add_argument(
        "--force",
        dest="force",
        action="store_true",
        help="Run even when the motif-to-TF annotations can not be loaded (e.g. no connection to "
        "the cisTarget webservers to download them, or no annotations for the species). The "
        "results then have no TF annotations and no cistromes.",
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
            force=arg.force,
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
            force=arg.force,
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


def add_parser_for_create_db(subparser: argparse._SubParsersAction):
    parser: argparse.ArgumentParser = subparser.add_parser(
        name=CREATE_DB_CALL_COMMAND,
        add_help=True,
        description=CREATE_DB_COMMAND_DESCRIPTION,
    )

    def create_db(arg):
        from pycistarget.create_database import create_cistarget_motif_databases

        create_cistarget_motif_databases(
            fasta=arg.fasta,
            motifs_dir=arg.motifs_dir,
            motifs_list=arg.motifs,
            db_prefix=arg.output_prefix,
            bg_padding=arg.bg_padding,
            mask=arg.mask,
            seed=arg.seed,
            n_threads=arg.threads,
        )

    parser.set_defaults(func=create_db)
    parser.add_argument(
        "--fasta",
        dest="fasta",
        type=str,
        required=True,
        help="FASTA file with the regions to score (e.g. created with "
        "create_fasta_with_padded_bg_from_bed.sh). The sequence names are the region IDs.",
    )
    parser.add_argument(
        "--motifs_dir",
        dest="motifs_dir",
        type=str,
        required=True,
        help="Directory with Cluster-Buster motif files (<motif ID>.cb).",
    )
    parser.add_argument(
        "--motifs",
        dest="motifs",
        type=str,
        required=True,
        help="File with the IDs of the motifs to score, one per line.",
    )
    parser.add_argument(
        "--output_prefix",
        dest="output_prefix",
        type=str,
        required=True,
        help="Prefix of the database files: <prefix>.motifs_vs_regions.scores.feather, "
        "<prefix>.regions_vs_motifs.scores.feather and "
        "<prefix>.regions_vs_motifs.rankings.feather.",
    )
    parser.add_argument(
        "--bg_padding",
        dest="bg_padding",
        type=int,
        default=0,
        help="Number of bases at the start and end of each sequence which are only used for "
        "calculating the background nucleotide frequencies, not for scoring motifs. "
        "Defaults to: %(default)s",
    )
    parser.add_argument(
        "--mask",
        dest="mask",
        action="store_true",
        help="Score masked (lowercase) bases as N.",
    )
    parser.add_argument(
        "--seed",
        dest="seed",
        type=int,
        default=None,
        help="Random seed for breaking ties when creating the rankings. With the same seed and "
        "input, the same rankings are created. By default a random seed is used.",
    )
    parser.add_argument(
        "--threads",
        dest="threads",
        type=int,
        default=0,
        help="Number of threads. Defaults to all cores.",
    )


def _add_plot_arguments(parser: argparse.ArgumentParser):
    """Arguments shared by the plot subcommands."""
    parser.add_argument(
        "--motifs_dir",
        dest="motifs_dir",
        type=str,
        required=True,
        help="Directory with Cluster-Buster motif files (<motif ID>.cb), used to draw the logos. "
        "For metaclusters, the matrix with the highest information content is drawn.",
    )
    parser.add_argument(
        "--output",
        dest="output",
        type=str,
        required=True,
        help="Output file. The format is set by the extension (e.g. .pdf, .png, .svg).",
    )
    parser.add_argument(
        "--threshold",
        dest="threshold",
        type=float,
        default=None,
        help="Only plot motifs with a score (NES for cisTarget, Log2FC for DEM) of at least this "
        "value. By default all enriched motifs are plotted.",
    )


def _save_plot(plot, output: str):
    import matplotlib.pyplot as plt

    plot.save(output)
    plt.close(plot.figure)
    logging.getLogger(__name__).info(f"Saved plot to {output}")


def add_parser_for_plot(subparser: argparse._SubParsersAction):
    parser: argparse.ArgumentParser = subparser.add_parser(
        name=PLOT_CALL_COMMAND,
        add_help=True,
        description=PLOT_COMMAND_DESCRIPTION,
    )
    parser.set_defaults(func=lambda arg: parser.print_help())
    plot_subparsers = parser.add_subparsers()

    barplot_parser = plot_subparsers.add_parser(
        name="barplot",
        add_help=True,
        description="Plot the logos of the enriched motifs of one region set, with a bar of their "
        "score (NES for cisTarget, Log2FC for DEM).",
    )

    def barplot(arg):
        from pycistarget.plotting import plot_logo_barplot

        _save_plot(
            plot_logo_barplot(
                arg.result, arg.motifs_dir, name=arg.name, threshold=arg.threshold, top=arg.top
            ),
            arg.output,
        )

    barplot_parser.set_defaults(func=barplot)
    barplot_parser.add_argument(
        "--result",
        dest="result",
        type=str,
        required=True,
        help="cisTarget or DEM result (hdf5 or tsv file written by pycistarget cistarget or dem).",
    )
    barplot_parser.add_argument(
        "--name",
        dest="name",
        type=str,
        default=None,
        help="Name of the result to plot, when the hdf5 file contains several results.",
    )
    barplot_parser.add_argument(
        "--top",
        dest="top",
        type=int,
        default=None,
        help="Only plot this number of motifs with the highest scores.",
    )
    _add_plot_arguments(barplot_parser)

    heatmap_parser = plot_subparsers.add_parser(
        name="heatmap",
        add_help=True,
        description="Plot the logos of the enriched motifs across region sets, with a heatmap of "
        "their score (NES for cisTarget, Log2FC for DEM) in each region set. Motifs are grouped by "
        "the region set with their highest score. Motifs that are not enriched in a region set "
        "have a grey cell.",
    )

    def heatmap(arg):
        from pycistarget.plotting import plot_logo_heatmap

        _save_plot(
            plot_logo_heatmap(arg.results, arg.motifs_dir, threshold=arg.threshold, top=arg.top),
            arg.output,
        )

    heatmap_parser.set_defaults(func=heatmap)
    heatmap_parser.add_argument(
        "--results",
        dest="results",
        type=str,
        nargs="+",
        required=True,
        help="cisTarget or DEM results (hdf5 or tsv files written by pycistarget cistarget or "
        "dem). Each result is a column of the heatmap.",
    )
    heatmap_parser.add_argument(
        "--top",
        dest="top",
        type=int,
        default=None,
        help="Only plot motifs that are among this number of motifs with the highest scores in "
        "at least one region set.",
    )
    _add_plot_arguments(heatmap_parser)


def create_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("--version", action="version", version=f"pycistarget {VERSION}")
    subparsers = parser.add_subparsers()
    add_parser_for_motif_enrichment_cistarget(subparsers)
    add_parser_for_motif_enrichment_dem(subparsers)
    add_parser_for_create_db(subparsers)
    add_parser_for_plot(subparsers)
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
