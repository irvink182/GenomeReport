###############################################################################
# GenomeReport Style
###############################################################################
import math
import numpy as np
import pandas as pd

from matplotlib import (
	colormaps,
	cm,)

from matplotlib.patches import (
	Rectangle,
	Patch,
)

from matplotlib.transforms import Bbox

from .protein_features import (
	PROTEIN_FEATURES
)

###############################################################################
# Figure sizes
###############################################################################

FIGSIZE = (10, 5)
CONSOLE_TEXT_WIDTH = 45
AGE_BINS = 20
BOXPLOT_FIGSIZE = (6, 5)
HEATMAP_FIGSIZE = (10, 5)
LOLLIPOP_FIGSIZE = (14, 5)
MAP_FIGSIZE = (8, 10)

###############################################################################
# Typography
###############################################################################

TITLE_SIZE = 16
SUBTITLE_SIZE = 14
LABEL_SIZE = 13
TICK_SIZE = 10
LEGEND_SIZE = 9
LEGEND_TITLE_SIZE=10
ANNOTATION_SIZE = 8

PANEL_TITLE_SIZE = 14

HEATMAP_ANNOTATION_SIZE = 6

LOLLIPOP_LABEL_SIZE = 8
MAP_LABEL_SIZE = 10

###############################################################################
# Colours
###############################################################################

COLOR_A = "#2C7FB8"
COLOR_B = "#1B9E77"

SUBTYPE_COLORS = {

    "A": COLOR_A,

    "B": COLOR_B,

}

COLOR_SENTINEL = "#2E8B57"
COLOR_NON_SENTINEL = "#D95F02"

OTHER_COLOR = "#7d7d7d"
UNASSIGNED_COLOR = "#C7C7C7"
OTHER_LABEL = "Other"
UNASSIGNED_LABEL = "unassigned"

HEATMAP_TOTAL_COLOR = "D9D9D9"

LOLLIPOP_STEM_COLOUR = "#C7C7C7"

MAP_EDGECOLOR = "white"
MAP_MISSING_COLOR = "#C7C7C7"

CLUSTERED_COLOR = "#1F77B4"
SINGLETON_COLOR = "#C7C7C7"


###############################################################################
# Colormaps
###############################################################################

LINEAGE_CMAP_A = colormaps["Blues"]
LINEAGE_CMAP_B = colormaps["BuGn"]
GENOTYPE_CMAP_A = cm.Blues
GENOTYPE_CMAP_B = cm.Greens
GENOTYPE_COLORS = {

    "ON1": GENOTYPE_CMAP_A(0.90),
    "NA1": GENOTYPE_CMAP_A(0.60),
    "GA5": GENOTYPE_CMAP_A(0.30),

    "BA9": GENOTYPE_CMAP_B(0.90),
    "GB3": GENOTYPE_CMAP_B(0.45),

}

MAP_CMAP = "Blues"

###############################################################################
# Legend
###############################################################################

LEGEND_FRAME = True
LEGEND_ALPHA = 1.0
LEGEND_EDGECOLOR = "#C7C7C7"
LEGEND_LINEWIDTH = 0.8
LEGEND_TARGET_ROWS = 3
LEGEND_COLUMNS = 5

LEGEND_SPACING = 0.05
LEGEND_BBOX = (1.02, 1)
LEGEND_BORDER_PAD = 0.4
LEGEND_HANDLE_LENGTH = 1.3
LEGEND_HANDLE_HEIGHT = 0.8
LEGEND_HANDLE_TEXT_PAD = 0.4
LEGEND_LABEL_SPACING = 0.35
LEGEND_COLUMN_SPACING = 1.0
FIGURE_BOTTOM = 0.15

###############################################################################
# Layout
###############################################################################

BAR_WIDTH = 0.80
BAR_HEIGHT = 0.45
LABEL_OFFSET = 0.02
GRID_ALPHA = 0.30
GRID_WIDTH = 0.8
AXES_WIDTH = 1.2

HEATMAP_BAR_RATIO = [0.2, 9.5, 0.12]

LOLLIPOP_MARKER_MIN = 40
LOLLIPOP_MARKER_MAX = 250
LOLLIPOP_BACKBONE_WIDTH = 3
LOLLIPOP_LABEL_THRESHOLD = 1
LOLLIPOP_LABEL_STEP = 0.20

# Distance (pixels) below which labels are considered overlapping
LOLLIPOP_LABEL_CLUSTER_PIXELS = 25
LOLLIPOP_LABEL_ROTATION = 45
LOLLIPOP_LABEL_MAX_VERTICAL = 2
LOLLIPOP_LABEL_X_STEP = 3.0


LOLLIPOP_LABEL_ROTATION = 90
LOLLIPOP_LABEL_OFFSET = 0.25
LOLLIPOP_LABEL_COLUMN_PIXELS = 22
LOLLIPOP_LABEL_LEVEL_PIXELS = 18
LOLLIPOP_CONNECTOR_COLOR = "#808080"
LOLLIPOP_CONNECTOR_WIDTH = 0.8

MAP_LINEWIDTH = 1.2
MAP_LABEL_ALPHA = 0.85
MAP_LABEL_BOXSTYLE = "round,pad=0.2"

###############################################################################
# Protein track
###############################################################################

PROTEIN_TRACK_Y = -0.08
PROTEIN_DOMAIN_HEIGHT = 0.35
PROTEIN_FEATURE_HEIGHT = 0.20
PROTEIN_DOMAIN_LABEL_SIZE = 8
PROTEIN_FEATURE_LABEL_SIZE = 6
PROTEIN_DOMAIN_LINEWIDTH = 1.0
PROTEIN_FEATURE_LINEWIDTH = 0.3
PROTEIN_FEATURE_OFFSET = 0.03
PROTEIN_MIN_FEATURE_WIDTH = 4
PROTEIN_TRACK_ALPHA = 0.85
PROTEIN_TRACK_EDGECOLOUR = "black"

###############################################################################
# Other parameters
###############################################################################

MEAN_LINESTYLE = ":"
MEDIAN_LINESTYLE = "--"
REFERENCE_LINEWIDTH = 1.5
REFERENCE_COLOR = "black"

GROUP_LABELS = {

    "A": "RSV-A",
    "B": "RSV-B",

}

LOLLIPOP_MIN_FREQUENCY = 10

FEATURE_COLOURS = {
    "signal": "#D9D9D9",
    "transmembrane": "#636363",
    "hypervariable": "#F28E2B",
    "domain": "#4E79A7",
    "motif": "#B07AA1",
    "epitope": "#59A14F",
    "drug_target": "#E15759",
}

FEATURE_LABELS = {
    "signal": "Signal peptide",
    "transmembrane": "Transmembrane region",
    "hypervariable": "Hypervariable region",
    "domain": "Conserved domain",
    "motif": "Motif",
    "epitope": "Epitope",
    "drug_target": "Antiviral target",
}

CONSOLE_MUTATION_FREQUENCY = 20

###############################################################################
# Styling functions
###############################################################################

def style_axes(ax):
    """
    Apply GenomeReport style.
    """

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.grid(
        axis="y",
        alpha=GRID_ALPHA,
        linewidth=GRID_WIDTH,
    )

    ax.tick_params(
    	axis="both",
    	labelsize=TICK_SIZE,
    )

    ax.set_axisbelow(True)

def style_labels(
    ax,
    xlabel=None,
    ylabel=None,
):

    if xlabel is not None:

        ax.set_xlabel(
            xlabel,
            fontsize=LABEL_SIZE,
        )

    if ylabel is not None:

        ax.set_ylabel(
            ylabel,
            fontsize=LABEL_SIZE,
        )

def style_title(
    ax,
    title,
):
    ax.set_title(
        title,
        fontsize=TITLE_SIZE,
        fontweight="bold",
        loc="left",
    )

def style_legend(
    ax,
    handles=None,
    labels=None,
    title=None,
    loc="best",
    bbox_to_anchor=None,
    ncol=1,
    adjust_figure=False,
):
    """
    Apply the standard legend style.
    """

    legend = ax.legend(
    	handles=handles,
    	labels=labels,
        title=title,
        loc=loc,
        bbox_to_anchor=bbox_to_anchor,
        ncol=ncol,
        frameon=LEGEND_FRAME,
        fancybox=True,
        framealpha=LEGEND_ALPHA,
        fontsize=LEGEND_SIZE,
        title_fontsize=LABEL_SIZE,
    )

    legend.get_frame().set_edgecolor(
        LEGEND_EDGECOLOR,
    )

    legend.get_frame().set_linewidth(
        0.8,
    )

    if adjust_figure:

        fig = ax.figure

        fig.canvas.draw()

        bbox = legend.get_window_extent().transformed(
            fig.transFigure.inverted()
        )

        fig.subplots_adjust(
            bottom=bbox.height + 0.06,
        )

    return legend

def style_bubble_legend(
    ax,
    handles,
    labels,
):
    """
    Standard legend for mutation frequency.
    """

    return style_legend(
        ax=ax,
        handles=handles,
        labels=labels,
        title="Frequency\n(% genomes)",
        loc="upper left",
        bbox_to_anchor=(1.03, 0.88),
    )

def style_feature_legend(
    ax,
    subtype,
    protein,
    bbox_to_anchor=(1.03, 0.38),
):
    """
    Legend for protein annotations.

    Includes both domain- and feature-level annotations
    present in the selected protein.
    """

    features = PROTEIN_FEATURES[subtype][protein]

    used_types = []

    # Domains first, features second
    for level in ("domain", "feature"):

        for feature in features:

            if (
                not feature["display"]
                or feature["level"] != level
            ):
                continue

            if feature["type"] not in used_types:

                used_types.append(
                    feature["type"]
                )

    handles = [

        Patch(
            facecolor=FEATURE_COLOURS[t],
            edgecolor="black",
            label=FEATURE_LABELS[t],
        )

        for t in used_types
    ]

    return style_legend(
        ax=ax,
        handles=handles,
        title="Protein annotation",
        loc="upper left",
        bbox_to_anchor=bbox_to_anchor,
    )

def style_colorbar(
    cbar,
    label=None,
    ticks=None,
    formatter=None,
    tick_length=4,
):

    if label is not None:

        cbar.set_label(
            label,
            fontsize=LABEL_SIZE,
        )

    cbar.ax.tick_params(
        labelsize=TICK_SIZE,
        length=tick_length,
    )

    if ticks is not None:
        cbar.set_ticks(ticks)

    if formatter is not None:
        cbar.ax.yaxis.set_major_formatter(
            formatter
        )

def style_map(ax):
    """
    Standard styling for geographic maps.
    """

    ax.set_axis_off()

    ax.set_aspect(
        "equal",
        adjustable="box",
    )

    ax.set_xlim(-19.0, 5.0)
    ax.set_ylim(27.0, 44.8)

def style_heatmap_axes(
    ax,
    hide_ylabels=False,
):

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    ax.tick_params(
        axis="both",
        length=0,
        labelsize=TICK_SIZE,
    )

    if hide_ylabels:
    	ax.set_yticklabels([])

def format_year_axis(
    ax,
    years,
    fontsize=10,
    every=1,
):
    """
    Format the x-axis for yearly plots.
    """

    years = (
        pd.Series(years)
        .astype(int)
        .sort_values()
    )

    ticks = years[::every]

    ax.set_xticks(
        ticks,
    )

    ax.set_xticklabels(
        ticks.astype(str),
        rotation=45,
        ha="right",
        fontsize=fontsize,
    )

    ax.tick_params(
        axis="x",
        pad=2,
    )

def style_protein_axis(
    ax,
):
	ax.set_yticks([])
	ax.tick_params(
		left=False,
		labelleft=False,
	)

	ax.spines["left"].set_visible(False)
	ax.spines["right"].set_visible(False)
	ax.spines["top"].set_visible(False)

def draw_protein_backbone(
    ax,
    length,
):
	ax.hlines(
		y=0,
		xmin=1,
		xmax=length,
		linewidth=LOLLIPOP_BACKBONE_WIDTH,
		color="black",
	)

###############################################################################
# Annotation helpers
###############################################################################
def add_bar_labels(
    ax,
    df,
    category,
    value="Total",
    position=None,
    threshold=0.08,
    inside=True,
    prefix="",
):
    """
    Add value labels to a horizontal barplot.
    """

    if position is None:
        position = value

    if isinstance(position, str):

        max_position = df[position].max()

    else:

        max_position = position

    label_offset = max_position * LABEL_OFFSET

    for _, row in df.iterrows():

        # Label text
        label = row[value]

        if isinstance(label, float):

            if label.is_integer():
                label = f"{int(label):,}"
            else:
                label = f"{label:.2f}"

        else:

            label = f"{label:,}"

        label = f"{prefix}{label}"

        # Label position
        if isinstance(position, str):
            xpos = row[position]
            max_position = df[position].max()
        else:
            xpos = position
            max_position = position

            if inside:
                xpos -= max_position * LABEL_OFFSET

        # Draw label
        if inside and xpos >= max_position * threshold:

            ax.text(
                xpos - label_offset,
                row[category],
                label,
                ha="right",
                va="center",
                fontsize=TICK_SIZE,
                color="white",
                fontweight="bold",
            )

        else:

            ax.text(
                xpos + label_offset,
                row[category],
                label,
                ha="left",
                va="center",
                fontsize=TICK_SIZE,
            )

def add_grouped_bar_labels(
    ax,
    y,
    values,
    offset,
    xmax,
    threshold=0.08,
    inside=True,
):
    """
    Add labels to grouped horizontal barplots.

    """

    label_offset = xmax * LABEL_OFFSET

    for yy, xx in zip(y, values):

        if pd.isna(xx):
            continue

        label = (
            f"{xx:,.0f}"
            if float(xx).is_integer()
            else f"{xx:.2f}"
        )

        if inside and xx >= xmax * threshold:

            ax.text(
                xx - label_offset,
                yy + offset,
                label,
                ha="right",
                va="center",
                fontsize=ANNOTATION_SIZE,
                color="white",
                fontweight="bold",
            )

        else:

            ax.text(
                xx + label_offset,
                yy + offset,
                label,
                ha="left",
                va="center",
                fontsize=ANNOTATION_SIZE,
                color="black",
            )

def add_heatmap_bar_labels(
    ax,
    y,
    values,
):
    """
    Add labels to the sample-size bars of a heatmap.
    """

    xmax = max(values)

    offset = xmax * LABEL_OFFSET

    for yy, xx in zip(y, values):

        ax.text(
            xx + offset,
            yy,
            f"{xx:,}",
            ha="left",
            va="center",
            fontsize=ANNOTATION_SIZE,
        )

def add_heatmap_row_labels(
    ax,
    labels,
    totals,
    title=None,
):
    """
    Draw row labels for categorical heatmaps.
    """

    y = np.arange(len(labels))

    ax.set_xlim(0, 1)
    ax.set_ylim(len(labels) - 0.5, -0.5)
    ax.axis("off")

    for yy, (label, total) in enumerate(zip(labels, totals)):

        ax.text(
            0,
            yy,
            f"{label} (n={total:,})",
            ha="right",
            va="center",
            fontsize=8,
        )

    if title is None:
    	ax.text(
    		0,
    		1.02,
    		title,
    		transform=ax.transAxes,
    		ha="left",
    		va="bottom",
    		fontsize=LABEL_SIZE,
    		fontweight="bold",
    	)

def style_map_label(
    ax,
    x,
    y,
    label,
):
    """
    Draw a map label.
    """

    ax.text(
        x,
        y,
        label,
        ha="center",
        va="center",
        fontsize=MAP_LABEL_SIZE,
        fontweight="bold",
        color="black",
        bbox=dict(
            facecolor="white",
            edgecolor="none",
            alpha=MAP_LABEL_ALPHA,
            boxstyle=MAP_LABEL_BOXSTYLE,
        ),
        zorder=10,
    )

def draw_map_labels(
    ax,
    gdf,
    value_column="Total",
):
    """
    Draw sequence counts for each autonomous community.
    """

    for _, row in gdf.iterrows():

        if row[value_column] == 0:
            continue

        x = row["x_cvi"]
        y = row["y_cvi"]

        style_map_label(
            ax=ax,
            x=x,
            y=y,
            label=f"{int(row[value_column])}",
        )

def get_domains(features):

    return [
        f
        for f in features
        if (
            f["level"] == "domain"
            and f["display"]
        )
    ]

def get_features(features):

    return [
        f
        for f in features
        if (
            f["level"] == "feature"
            and f["display"]
        )
    ]

def draw_domains(
    ax,
    features,
):
    """
    Draw protein domains.
    """

    for feature in get_domains(features):

        if (
            feature["level"] != "domain"
            or not feature["display"]
        ):
            continue

        center = (
            feature["start"]
            + feature["end"]
        ) / 2

        rect = Rectangle(
            (
                feature["start"],
                PROTEIN_TRACK_Y,
            ),
            width=(
                feature["end"]
                - feature["start"]
                + 1
            ),
            height=PROTEIN_DOMAIN_HEIGHT,
            facecolor=FEATURE_COLOURS[
                feature["type"]
            ],
            edgecolor="white",
            linewidth=PROTEIN_DOMAIN_LINEWIDTH,
            alpha=PROTEIN_TRACK_ALPHA,
            clip_on=False,
            zorder=2,
        )

        ax.add_patch(rect)

        ax.text(
            center,
            PROTEIN_TRACK_Y + PROTEIN_DOMAIN_HEIGHT / 2,
            feature["label"],
            ha="center",
            va="center",
            fontsize=PROTEIN_DOMAIN_LABEL_SIZE,
            color="black",
            fontweight="bold",
            clip_on=False,
            zorder=3,
        )

def draw_features(
    ax,
    features,
):
    """
    Draw motifs, epitopes and drug targets.
    """

    feature_y = (
        PROTEIN_TRACK_Y
        + PROTEIN_DOMAIN_HEIGHT
        + PROTEIN_FEATURE_OFFSET
    )

    for feature in get_features(features):

        if (
            feature["level"] != "feature"
            or not feature["display"]
        ):
            continue

        center = (
            feature["start"]
            + feature["end"]
        ) / 2

        ax.plot(
            [center, center],
            [
                PROTEIN_TRACK_Y
                + PROTEIN_DOMAIN_HEIGHT,
                feature_y,
            ],
            color="black",
            linewidth=0.8,
            zorder=4,
        )

        width = max(
            feature["end"] - feature["start"] + 1,
            PROTEIN_MIN_FEATURE_WIDTH,
        )

        rect = Rectangle(
            (
                feature["start"],
                feature_y,
            ),
            width=width,
            height=PROTEIN_FEATURE_HEIGHT,
            facecolor=FEATURE_COLOURS[
                feature["type"]
            ],
            edgecolor=PROTEIN_TRACK_EDGECOLOUR,
            linewidth=PROTEIN_FEATURE_LINEWIDTH,
            alpha=PROTEIN_TRACK_ALPHA,
            clip_on=False,
            zorder=5,
        )

        ax.add_patch(rect)

        # Feature label
        label = feature["label"]

        target = feature.get(
            "target"
        )

        if (
            target is not None
            and target != label
        ):
            label = (
                f"{label}\n"
                f"({target})"
            )

        ax.text(
            center,
            feature_y
            + PROTEIN_FEATURE_HEIGHT
            + PROTEIN_FEATURE_OFFSET
            + 0.01,
            label,
            ha="left",
            va="bottom",
            rotation=0,
            fontsize=PROTEIN_FEATURE_LABEL_SIZE + 1.5,
            linespacing=1.2,
            clip_on=False,
            zorder=6,
        )

def draw_annotation_levels(ax):
    """
    Draw labels indicating the two protein annotation levels.
    """

    ax.text(
        -5,
        0.28,
        "Feature",
        ha="right",
        va="center",
        fontsize=9,
        fontweight="normal",
    )

    ax.text(
        -5,
        -0.28,
        "Domain",
        ha="right",
        va="center",
        fontsize=9,
        fontweight="normal",
    )

def protein_track(
    ax,
    subtype,
    protein,
):
    """
    Draw annotated protein track.
    """

    features = PROTEIN_FEATURES[subtype][protein]
    
    draw_domains(
        ax=ax,
        features=features,
    )

    draw_features(
        ax=ax,
        features=features,
    )

    draw_annotation_levels(
    	ax=ax,
    )

def arrange_lollipop_labels(
    ax,
    labels,
    top_fraction=0.96,
    bottom_fraction=0.88,
    rotation=35,
    min_dx=10,
):
    """
    Arrange lollipop labels in two horizontal rows.
    """

    if len(labels) == 0:
        return

    # Axis limits
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    y_top = ymax * top_fraction
    y_bottom = ymax * bottom_fraction

    # Sort labels by x position
    labels = sorted(
        labels,
        key=lambda d: d["x"],
    )

    # Split into two rows
    row_top = labels[::2]
    row_bottom = labels[1::2]

    # Spread labels horizontally
    def spread(row):

        if len(row) == 0:
            return

        xs = [d["x"] for d in row]

        for i in range(1, len(xs)):

            if xs[i] - xs[i-1] < min_dx:

                xs[i] = xs[i-1] + min_dx

        if xs[-1] > xmax:

            shift = xs[-1] - xmax

            xs = [
                x - shift
                for x in xs
            ]

        for i in range(len(xs)-2, -1, -1):

            if xs[i+1] - xs[i] < min_dx:

                xs[i] = xs[i+1] - min_dx

        return xs

    xs_top = spread(row_top)
    xs_bottom = spread(row_bottom)

    # Draw top row
    for xpos, label in zip(xs_top, row_top):

        ax.plot(
            [label["x"], xpos],
            [label["y"], y_top],
            color="0.55",
            lw=0.8,
            zorder=2,
        )

        ax.text(
            xpos,
            y_top,
            label["text"],
            fontsize=LOLLIPOP_LABEL_SIZE,
            rotation=rotation,
            ha="left",
            va="bottom",
            zorder=5,
        )

    # Draw bottom row
    for xpos, label in zip(xs_bottom, row_bottom):

        ax.plot(
            [label["x"], xpos],
            [label["y"], y_bottom],
            color="0.55",
            lw=0.8,
            zorder=2,
        )

        ax.text(
            xpos,
            y_bottom,
            label["text"],
            fontsize=LOLLIPOP_LABEL_SIZE,
            rotation=rotation,
            ha="left",
            va="bottom",
            zorder=5,
        )