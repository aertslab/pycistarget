"""
Plot motif logos of enriched motifs, next to their enrichment scores.

Needs logomaker and matplotlib: ``pip install "pycistarget[plot]"``.

Logos are drawn from the Cluster-Buster motif files (``<motif ID>.cb``) of the motif collection. A
metacluster motif file contains the matrices of all motifs in the cluster: the matrix with the
highest information content is drawn. The score is the NES for cisTarget results and the Log2FC
for DEM results.

Example::

    from pycistarget.input_output import read_hdf5
    from pycistarget.plotting import plot_logo_barplot, plot_logo_heatmap

    results = read_hdf5("motif_enrichment_cistarget_Topic_1.hdf5")
    plot = plot_logo_barplot(results["Topic_1"], "v10nr_clust_public/singletons", threshold=5)
    plot.ax.set_xlim(0, 30)  # customise with matplotlib
    plot.save("Topic_1.pdf")

    plot = plot_logo_heatmap(
        ["motif_enrichment_dem_Topic_1_vs_all.hdf5", "motif_enrichment_dem_Topic_2_vs_all.hdf5"],
        "v10nr_clust_public/singletons",
        top=5,
    )
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import TYPE_CHECKING, Any, TypeAlias

import numpy as np
import pandas as pd

from pycistarget.motif_enrichment_result import MotifEnrichmentResult

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.colorbar import Colorbar
    from matplotlib.colors import Colormap
    from matplotlib.figure import Figure

log = logging.getLogger(__name__)

__all__ = [
    "LogoPlot",
    "draw_logo",
    "enrichment_scores",
    "information_matrix",
    "load_logo",
    "plot_logo_barplot",
    "plot_logo_heatmap",
    "read_cluster_buster_matrices",
    "read_enrichment_scores",
]

#: Score column of the motif enrichment table, per analysis type.
SCORE_COLUMNS = {"cisTarget": "NES", "DEM": "Log2FC"}
_TSV_PREFIXES = ("motif_enrichment_cistarget_", "motif_enrichment_dem_")

#: Default bar color and (single hue, light to dark) heatmap colors.
BAR_COLOR = "#2a78d6"
HEATMAP_COLORS = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
MISSING_COLOR = "#f0efec"
_INK = "#52514e"

# Margins (inches).
_MARGIN_LEFT = 0.2
_MARGIN_TOP = 0.5

#: A single result: a hdf5/tsv file, a cisTarget/DEM result, or scores per motif.
Result: TypeAlias = "str | PathLike[str] | MotifEnrichmentResult | pd.Series"
#: Several results: one of the above, a list of them, a dictionary of named results, or a table
#: of scores (motifs x region sets).
Results: TypeAlias = (
    "Result | Sequence[Result] | Mapping[str, MotifEnrichmentResult | pd.Series] | pd.DataFrame"
)


@dataclass
class LogoPlot:
    """A logo plot. Use the axes to customise it with matplotlib."""

    #: The matplotlib figure.
    figure: Figure
    #: Axes with the scores (bars or heatmap).
    ax: Axes
    #: Axes with the logos, one per motif (top to bottom).
    logo_axes: list[Axes]
    #: Motif IDs, top to bottom.
    motifs: list[str]
    #: Color bar of the heatmap (None for a barplot).
    colorbar: Colorbar | None = None

    def save(self, path: str | PathLike, **kwargs: Any):
        """Save the plot (format from the file extension). ``kwargs`` are passed to ``savefig``."""
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        self.figure.savefig(path, **{"dpi": 200, "facecolor": "white", **kwargs})


def _import_plotting():
    try:
        import logomaker
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError(
            'logomaker and matplotlib are required for plotting, install them with: pip install "pycistarget[plot]"'
        ) from e
    return logomaker, plt


# Logos.


def read_cluster_buster_matrices(path: str | PathLike) -> list[np.ndarray]:
    """Read the matrices (positions x A, C, G, T counts or frequencies) of a Cluster-Buster motif file."""
    matrices: list[list[list[float]]] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                matrices.append([])
            elif line and not line.startswith("#"):
                if not matrices:
                    raise ValueError(f"{path}: matrix rows before the first '>' header")
                matrices[-1].append([float(value) for value in line.split()[:4]])
    result = [np.array(rows, dtype=np.float64) for rows in matrices if rows]
    if not result or any(m.ndim != 2 or m.shape[1] != 4 for m in result):
        raise ValueError(f"{path}: not a Cluster-Buster motif file")
    return result


def information_matrix(matrix: np.ndarray) -> pd.DataFrame:
    """Information content (bits) of each base (columns A, C, G, T) at each position (rows)."""
    probabilities = matrix / matrix.sum(axis=1, keepdims=True).clip(min=np.finfo(float).tiny)
    with np.errstate(divide="ignore", invalid="ignore"):
        entropy = -np.nansum(probabilities * np.log2(probabilities), axis=1)
    return pd.DataFrame(probabilities * (2.0 - entropy)[:, None], columns=list("ACGT"))


def load_logo(motifs_dir: str | PathLike, motif: str) -> pd.DataFrame | None:
    """
    Logo (information matrix, see :func:`information_matrix`) of a motif from ``<motifs_dir>/<motif>.cb``.

    For a file with several matrices (metaclusters), the matrix with the highest total information
    content is returned. Returns None, with a warning, if the motif file is missing or invalid.
    """
    path = Path(motifs_dir) / f"{motif}.cb"
    if not path.is_file():
        log.warning(f"No logo for motif {motif}: {path.name} not found in {motifs_dir}")
        return None
    try:
        matrices = read_cluster_buster_matrices(path)
    except (OSError, ValueError) as e:
        log.warning(f"No logo for motif {motif}: {e}")
        return None
    logos = [information_matrix(m) for m in matrices]
    return max(logos, key=lambda logo: logo.to_numpy().sum())


def draw_logo(
    ax: Axes,
    logo: pd.DataFrame | None,
    width: int | None = None,
    color_scheme: str | dict = "classic",
):
    """
    Draw a logo (from :func:`load_logo`) on ``ax``, without axes.

    ``width`` sets the number of positions that fit on the axes (default: the length of the logo),
    so logos of different lengths can be drawn on the same scale. Nothing is drawn for a missing
    logo (None). ``color_scheme`` is a logomaker color scheme.
    """
    logomaker, _ = _import_plotting()
    if logo is not None:
        logomaker.Logo(logo, ax=ax, color_scheme=color_scheme, show_spines=False, baseline_width=0)
    ax.set_axis_off()
    ax.set_xlim(-0.5, (width if width is not None else len(logo) if logo is not None else 1) - 0.5)
    # 2 bits is the maximum; leave some space between rows of logos.
    ax.set_ylim(-0.15, 2.3)


# Scores.


def _result_name_from_tsv(path: Path) -> str:
    name = path.name.removesuffix(".tsv")
    for prefix in _TSV_PREFIXES:
        if name.startswith(prefix):
            return name[len(prefix) :]
    return name


def _scores_from_table(table: pd.DataFrame, source: object) -> tuple[str | None, pd.Series]:
    for analysis_type, column in SCORE_COLUMNS.items():
        if column in table.columns:
            return analysis_type, table[column].astype(np.float64)
    if len(table) == 0:  # no enriched motifs
        return None, pd.Series(dtype=np.float64)
    raise ValueError(f"{source}: the motif enrichment table has no NES or Log2FC column")


def _result_scores(result: MotifEnrichmentResult) -> tuple[str, pd.Series]:
    analysis_type, scores = _scores_from_table(result.motif_enrichment, result.name)
    return analysis_type or result._analysis_type, scores


def read_enrichment_scores(path: str | PathLike) -> dict[str, tuple[str | None, pd.Series]]:
    """
    Read the motif enrichment scores of the results in a hdf5 or tsv file.

    Returns a dictionary with the analysis name as key and as value the analysis type
    (``"cisTarget"``, ``"DEM"`` or None for a tsv file without enriched motifs) and the scores
    (NES or Log2FC) of the enriched motifs (index: motif IDs).
    """
    path = Path(path)
    if path.suffix == ".tsv":
        table = pd.read_csv(path, sep="\t", index_col=0)
        return {_result_name_from_tsv(path): _scores_from_table(table, path)}

    from pycistarget.input_output import read_hdf5

    return {name: _result_scores(result) for name, result in read_hdf5(path).items()}


def _collect(results: Any, name: str | None = None) -> list[tuple[str, str | None, pd.Series]]:
    """(name, analysis type, scores) of each result."""
    if isinstance(results, (str, PathLike)):
        return [(n, t, s) for n, (t, s) in read_enrichment_scores(results).items()]
    if isinstance(results, MotifEnrichmentResult):
        return [(results.name, *_result_scores(results))]
    if isinstance(results, pd.Series):
        # A column of a motif enrichment table (e.g. ``table["NES"]``) is named after the score.
        analysis_type = next((t for t, c in SCORE_COLUMNS.items() if c == results.name), None)
        if name is None:
            name = "" if results.name is None or analysis_type is not None else str(results.name)
        return [(name, analysis_type, results)]
    if isinstance(results, pd.DataFrame):
        return [(str(c), None, results[c]) for c in results.columns]
    if isinstance(results, Mapping):
        return [item for n, r in results.items() for item in _collect(r, name=str(n))]
    if isinstance(results, Sequence):
        return [item for r in results for item in _collect(r)]
    raise TypeError(f"Can not plot results of type {type(results).__name__}")


def enrichment_scores(results: Results) -> tuple[pd.DataFrame, str]:
    """
    Scores of the enriched motifs of one or more results, as plotted.

    Returns a table of scores (motifs x region sets, NaN for motifs that are not enriched in a region
    set) and the name of the score ("NES" for cisTarget, "Log2FC" for DEM results).

    ``results`` can be a hdf5 or tsv file written by ``pycistarget cistarget`` or ``dem``, a
    cisTarget or DEM result, a series of scores per motif, a list of these, a dictionary of named
    results (e.g. from :func:`pycistarget.input_output.read_hdf5`) or a table of scores (motifs x
    region sets).
    """
    collected = _collect(results)
    names = [n for n, _, _ in collected]
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        raise ValueError(f"Several results are named: {', '.join(duplicates)}")
    types = {t for _, t, _ in collected if t is not None}
    if len(types) > 1:
        raise ValueError("cisTarget and DEM results can not be plotted together.")
    score_name = SCORE_COLUMNS[types.pop()] if types else "Score"
    table = pd.DataFrame({n: s.astype(np.float64) for n, _, s in collected}, dtype=np.float64)
    return table, score_name


# Plots.


def _figure(
    n_rows: int,
    row_height: float,
    logo_width: float,
    data_width: float,
    right: float,
    bottom: float,
):
    """Figure with a column of logos and a column with the data (sizes and margins in inches)."""
    _, plt = _import_plotting()
    width = _MARGIN_LEFT + logo_width + data_width + right
    height = n_rows * row_height + _MARGIN_TOP + bottom
    fig = plt.figure(figsize=(width, height))
    gridspec = fig.add_gridspec(
        n_rows,
        2,
        width_ratios=[logo_width, data_width],
        hspace=0.0,
        wspace=0.05,
        left=_MARGIN_LEFT / width,
        right=1 - right / width,
        top=1 - _MARGIN_TOP / height,
        bottom=bottom / height,
    )
    return fig, gridspec


def _draw_logos(
    fig: Figure, gridspec, motifs: list[str], motifs_dir: str | PathLike, color_scheme
) -> list[Axes]:
    """Draw the logo of each motif in its own row of the first column of ``gridspec``."""
    logos = {motif: load_logo(motifs_dir, motif) for motif in motifs}
    # The same scale for all logos: one position has the same width in every row.
    width = max((len(logo) for logo in logos.values() if logo is not None), default=1)
    axes = []
    for i, motif in enumerate(motifs):
        ax = fig.add_subplot(gridspec[i, 0])
        draw_logo(ax, logos[motif], width=width, color_scheme=color_scheme)
        axes.append(ax)
    return axes


def _title(fig: Figure, title: str | None):
    if title:
        fig.suptitle(title, x=0.02, ha="left", y=1 - 0.1 / fig.get_figheight(), va="top")


def plot_logo_barplot(
    result: Results,
    motifs_dir: str | PathLike,
    *,
    name: str | None = None,
    threshold: float | None = None,
    top: int | None = None,
    title: str | None = None,
    color: str = BAR_COLOR,
    row_height: float = 0.4,
    logo_width: float = 2.4,
    bar_width: float = 3.0,
    color_scheme: str | dict = "classic",
) -> LogoPlot:
    """
    Plot the logos of the enriched motifs of one region set (rows), with a bar of their score.

    Parameters
    ----------
    result
        The result to plot: a hdf5 or tsv file written by ``pycistarget cistarget`` or ``dem``, a
        cisTarget or DEM result, or a series of scores per motif (see :func:`enrichment_scores`).
    motifs_dir
        Directory with the Cluster-Buster motif files (``<motif ID>.cb``) for the logos.
    name
        Name of the result to plot, if ``result`` contains several (e.g. a hdf5 file with several
        analyses or the dictionary returned by ``read_hdf5``).
    threshold
        Only plot motifs with a score of at least this value.
    top
        Only plot the ``top`` motifs with the highest scores.
    title
        Title of the plot. Defaults to the name of the result; use "" for no title.
    color
        Color of the bars.
    row_height, logo_width, bar_width
        Height of a row, width of the logos and width of the bars (inches).
    color_scheme
        logomaker color scheme of the logos.

    Returns
    -------
    The plot, with the figure and axes to customise it.
    """
    table, score_name = enrichment_scores(result)
    if name is not None:
        if name not in table.columns:
            raise ValueError(f"No result named {name}, available: {', '.join(table.columns)}")
        table = table[[name]]
    elif table.shape[1] != 1:
        raise ValueError(
            "Several results to plot, select one with name: " + ", ".join(table.columns)
        )
    result_name = str(table.columns[0])
    scores = table.iloc[:, 0].dropna().sort_values(ascending=False, kind="stable")
    if threshold is not None:
        scores = scores[scores >= threshold]
    if top is not None:
        scores = scores.iloc[:top]
    if len(scores) == 0:
        raise ValueError(
            f"No motifs to plot for {result_name} (no enriched motifs above the threshold)."
        )

    n = len(scores)
    motifs = [str(m) for m in scores.index]
    fig, gridspec = _figure(n, row_height, logo_width, bar_width, right=0.4, bottom=0.6)
    logo_axes = _draw_logos(fig, gridspec, motifs, motifs_dir, color_scheme)

    ax = fig.add_subplot(gridspec[:, 1])
    ax.barh(np.arange(n), scores.to_numpy(), height=0.6, color=color)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xlim(min(0.0, scores.min()), None)
    ax.set_xlabel(score_name)
    ax.grid(axis="x", color="#e4e3df", linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#8a8984")
    ax.tick_params(colors=_INK, left=False, labelleft=False)
    ax.xaxis.label.set_color(_INK)
    _title(fig, result_name if title is None else title)
    return LogoPlot(figure=fig, ax=ax, logo_axes=logo_axes, motifs=motifs)


def plot_logo_heatmap(
    results: Results,
    motifs_dir: str | PathLike,
    *,
    threshold: float | None = None,
    top: int | None = None,
    title: str | None = None,
    cmap: str | Colormap | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    missing_color: str = MISSING_COLOR,
    row_height: float = 0.4,
    logo_width: float = 2.4,
    column_width: float = 0.3,
    color_scheme: str | dict = "classic",
) -> LogoPlot:
    """
    Plot the logos of enriched motifs (rows) with a heatmap of their score in each region set.

    Motifs are grouped by the region set in which they have the highest score (in the order of the
    results) and sorted by that score. Motifs that are not enriched in a region set have no score
    (``missing_color``).

    Parameters
    ----------
    results
        The results to plot (one column each): hdf5 or tsv files written by ``pycistarget
        cistarget`` or ``dem``, cisTarget or DEM results, a dictionary of named results (e.g. from
        ``read_hdf5``) or a table of scores (motifs x region sets); see :func:`enrichment_scores`.
    motifs_dir
        Directory with the Cluster-Buster motif files (``<motif ID>.cb``) for the logos.
    threshold
        Only plot motifs with a score of at least this value in at least one region set.
    top
        Only plot motifs that are among the ``top`` motifs with the highest scores in at least one
        region set.
    title
        Title of the plot.
    cmap, vmin, vmax
        Colormap and range of the scores (default: light to dark blue, from the lowest to the
        highest plotted score).
    missing_color
        Color of the cells of motifs that are not enriched in a region set.
    row_height, logo_width, column_width
        Height of a row, width of the logos and width of a heatmap column (inches).
    color_scheme
        logomaker color scheme of the logos.

    Returns
    -------
    The plot, with the figure, axes and color bar to customise it.
    """
    _, plt = _import_plotting()
    from matplotlib.colors import LinearSegmentedColormap

    scores, score_name = enrichment_scores(results)
    if scores.shape[1] == 0:
        raise ValueError("No region sets to plot.")
    keep = scores.notna().any(axis=1)
    if threshold is not None:
        keep &= (scores >= threshold).any(axis=1)
    if top is not None:
        keep &= (scores.rank(axis=0, ascending=False, method="first") <= top).any(axis=1)
    scores = scores[keep]
    if len(scores) == 0:
        raise ValueError("No motifs to plot (no enriched motifs above the threshold).")

    best_region_set = scores.fillna(-np.inf).to_numpy().argmax(axis=1)
    scores = scores.iloc[np.lexsort((-scores.max(axis=1).to_numpy(), best_region_set))]

    n, m = scores.shape
    motifs = [str(motif) for motif in scores.index]
    # Room for the rotated region set names below and the color bar on the right.
    label_length = max(len(str(c)) for c in scores.columns)
    fig, gridspec = _figure(
        n, row_height, logo_width, m * column_width, right=1.0, bottom=0.2 + 0.09 * label_length
    )
    logo_axes = _draw_logos(fig, gridspec, motifs, motifs_dir, color_scheme)

    ax = fig.add_subplot(gridspec[:, 1])
    colormap = (
        plt.get_cmap(cmap)
        if cmap is not None
        else LinearSegmentedColormap.from_list("scores", HEATMAP_COLORS)
    )
    image = ax.imshow(
        np.ma.masked_invalid(scores.to_numpy()),
        aspect="auto",
        cmap=colormap.with_extremes(bad=missing_color),
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
        extent=(-0.5, m - 0.5, n - 0.5, -0.5),
    )
    # White lines between the cells.
    ax.set_xticks(np.arange(m + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(n + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="white", linewidth=1)
    ax.tick_params(which="minor", length=0)
    ax.set_xticks(np.arange(m), [str(c) for c in scores.columns], rotation=90)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(colors=_INK, left=False, labelleft=False, bottom=False)

    colorbar_height = min(1.5, n * row_height) / fig.get_figheight()
    bbox = ax.get_position()
    cax = fig.add_axes(
        (
            bbox.x1 + 0.15 / fig.get_figwidth(),
            bbox.y1 - colorbar_height,
            0.12 / fig.get_figwidth(),
            colorbar_height,
        )
    )
    colorbar = fig.colorbar(image, cax=cax)
    colorbar.set_label(score_name, color=_INK)
    colorbar.outline.set_visible(False)
    cax.tick_params(colors=_INK)
    _title(fig, title)
    return LogoPlot(figure=fig, ax=ax, logo_axes=logo_axes, motifs=motifs, colorbar=colorbar)
