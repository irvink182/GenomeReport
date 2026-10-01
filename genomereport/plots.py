"""
Plots.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MaxNLocator
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm
from shapely.affinity import translate
from matplotlib.patches import Rectangle
from shapely.ops import unary_union
from adjustText import adjust_text
import geopandas as gpd
import pandas as pd
import numpy as np
import math
import warnings

from .config import (
    Columns, 
    MutationColumns,
    Metadata,
    SUBTYPES,
    PROTEIN_ORDER,
    PROTEIN_LENGTHS,
)

from .style import (
    FIGSIZE,
    BOXPLOT_FIGSIZE,
    COLOR_A,
    COLOR_B,
    SUBTYPE_COLORS,
    COLOR_SENTINEL,
    COLOR_NON_SENTINEL,
    OTHER_COLOR,
    UNASSIGNED_COLOR,
    OTHER_LABEL,
    UNASSIGNED_LABEL,
    HEATMAP_TOTAL_COLOR,
    LINEAGE_CMAP_A,
    LINEAGE_CMAP_B,
    GENOTYPE_CMAP_A,
    GENOTYPE_CMAP_B,
    GENOTYPE_COLORS,
    BAR_WIDTH,
    BAR_HEIGHT,
    MEAN_LINESTYLE,
    MEDIAN_LINESTYLE,
    REFERENCE_LINEWIDTH,
    REFERENCE_COLOR,
    LABEL_SIZE,
    TICK_SIZE,
    LABEL_OFFSET,
    LEGEND_TARGET_ROWS,
    LEGEND_COLUMNS,
    LEGEND_BBOX,
    LEGEND_SPACING,
    FIGURE_BOTTOM,
    style_axes,
    style_labels,
    style_title,
    style_colorbar,
    style_heatmap_axes,
    style_legend,
    format_year_axis,
    add_bar_labels,
    add_grouped_bar_labels,
    add_heatmap_bar_labels,
    add_heatmap_row_labels,
    GROUP_LABELS,
    HEATMAP_ANNOTATION_SIZE,
    HEATMAP_FIGSIZE,
    HEATMAP_BAR_RATIO,
    LOLLIPOP_FIGSIZE,
    LOLLIPOP_STEM_COLOUR,
    LOLLIPOP_LABEL_SIZE,
    LOLLIPOP_MARKER_MIN,
    LOLLIPOP_MARKER_MAX,
    LOLLIPOP_BACKBONE_WIDTH,
    LOLLIPOP_LABEL_THRESHOLD,
    style_protein_axis,
    draw_protein_backbone,
    LOLLIPOP_MIN_FREQUENCY,
    protein_track,
    style_bubble_legend,
    style_feature_legend,
    arrange_lollipop_labels,
    LOLLIPOP_LABEL_OFFSET,
    draw_map_labels,
    style_map_label,
    MAP_FIGSIZE,
    MAP_CMAP,
    MAP_LINEWIDTH,
    MAP_EDGECOLOR,
    MAP_MISSING_COLOR,
    MAP_LABEL_SIZE,
    MAP_LABEL_ALPHA,
    MAP_LABEL_BOXSTYLE,
    style_map,
    TITLE_SIZE,
    CLUSTERED_COLOR,
    SINGLETON_COLOR,
)

from .virus_dynamics import (
    pivot_category,
    complete_time_series,
    filter_category_matrix,
)

from .geography import (
    SPAIN_CCAA_GEOJSON,
    AUTONOMOUS_COMMUNITY_MAP,
    MAP_LABEL_OFFSETS,
)

from .tables import(
    lineage_cluster_heatmap_table,
)

from .clusters import(
    prepare_cluster_plot_data,
)

# Open figures at the same time
mpl.rcParams["figure.max_open_warning"] = 100

#############################
# Colors
#############################

def get_category_subtype(category):
    """
    Return RSV subtype from a lineage/genotype name.
    """

    category = str(category).strip().upper()

    if category.startswith("A"):
        return "A"

    if category.startswith("B"):
        return "B"

    return None

def get_category_colors(
    categories,
    cmap_a,
    cmap_b,
    subtype=None,
):
    """
    Return colours for RSV categories.
    """

    colors = {}

    # Real categories
    real_categories = [
        category
        for category in categories
        if category != OTHER_LABEL
        and str(category).lower() != UNASSIGNED_LABEL.lower()
    ]

    # RSV-A
    if subtype in ("A", None):

        categories_a = [
            category
            for category in real_categories
            if get_category_subtype(category) == "A"
        ]

        if categories_a:

            levels = np.linspace(
                0.90,
                0.15,
                len(categories_a),
            )

            for category, level in zip(categories_a, levels):

                colors[category] = cmap_a(level)

    # RSV-B
    if subtype in ("B", None):

        categories_b = [
            category
            for category in real_categories
            if get_category_subtype(category) == "B"
        ]

        if categories_b:

            levels = np.linspace(
                0.90,
                0.15,
                len(categories_b),
            )

            for category, level in zip(categories_b, levels):

                colors[category] = cmap_b(level)

    # Special categories
    colors[OTHER_LABEL] = OTHER_COLOR
    colors[UNASSIGNED_LABEL] = UNASSIGNED_COLOR

    # Preserve original order
    return [
        colors.get(category, OTHER_COLOR)
        for category in categories
    ]

def get_lineage_colors(
    categories,
    subtype=None,
):
    """
    Return colours for lineage plots.
    """

    return get_category_colors(
        categories=categories,
        cmap_a=LINEAGE_CMAP_A,
        cmap_b=LINEAGE_CMAP_B,
        subtype=subtype,
    )

def get_genotype_colors(categories):
    """
    Return colours for genotype plots.
    """

    colors = []

    for category in categories:

        if category == OTHER_LABEL:

            colors.append(OTHER_COLOR)

        elif str(category).lower() == UNASSIGNED_LABEL.lower():

            colors.append(UNASSIGNED_COLOR)

        else:

            colors.append(
                GENOTYPE_COLORS.get(
                    category,
                    OTHER_COLOR,
                )
            )

    return colors

#############################
# Generic Plots
#############################
def style_layout(fig):
    """
    Apply a consistent layout
    while suppressing harmless
    Matplotlib layout warnings.
    """

    with warnings.catch_warnings():

        warnings.filterwarnings(
            "ignore",
            category=UserWarning,
            message=".*tight_layout.*",
        )

        fig.tight_layout()

def temporal_stacked_barplot(
    table,
    x,
    columns,
    colors,
    labels,
    ylabel,
    percentage=False,
    show_totals=False,
    legend_title=None,
    legend_large=False,
    title=None,
):
    """
    Generic stacked barplot for temporal distributions.
    """

    df = table.copy()

    # Normalize to percentages
    if percentage:

        for column in columns:

            df[column] = (
                df[column]
                / df["Total"]
                * 100
            )

    # Figure
    fig, ax = plt.subplots(
        figsize=FIGSIZE,
    )

    bottom = np.zeros(
        len(df),
    )

    for column, color, label in zip(
        columns,
        colors,
        labels,
    ):

        ax.bar(
            df[x],
            df[column],
            width=BAR_WIDTH,
            bottom=bottom,
            color=color,
            label=label,
        )

        bottom += df[column]

    # Labels
    style_labels(
        ax,
        xlabel="Collection year",
        ylabel=ylabel,
    )

    # Sample size
    if show_totals:

        ymax = 100 if percentage else bottom.max()

        for year, total in zip(
            df[x],
            df["Total"],
        ):

            ax.text(
                year,
                ymax + 2,
                f"{int(total)}",
                ha="center",
                va="bottom",
                fontsize=TICK_SIZE,
                fontweight="bold",
            )

    # Axis
    if percentage:

        ax.set_ylim(
            0,
            108,
        )

        ax.set_yticks(
            [
                0,
                25,
                50,
                75,
                100,
            ]
        )

    format_year_axis(
        ax,
        df[x],
        every=1,
    )

    # Style
    if title is not None:

        style_title(
            ax,
            title,
        )

    style_axes(
        ax,
    )

    if legend_large:
        style_legend(
            ax,
            title=legend_title,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.16),
            ncol=5,
        )

        fig.subplots_adjust(
            bottom=0.34,
        )
    else:

        style_legend(
            ax,
            title=legend_title,
        )

    style_layout(
        fig,
    )

    return fig

def categorical_stacked_barplot(
    table,
    category,
    columns=None,
    colors=None,
    labels=None,
    xlabel="Number of sequences",
    x_margin=1.05,
    normalize=False,
    show_totals=True,
    title=None,
):
    """
    Horizontal stacked barplot for categorical variables.
    """

    if columns is None:
        columns = ["A", "B"]

    if colors is None:
        colors = [COLOR_A, COLOR_B]

    if labels is None:
        labels = columns

    df = table.copy()

    for column in columns:
        if column not in df.columns:
            df[column] = 0

    df = df.sort_values(
        "Total",
        ascending=True,
    )

    fig_height = max(
        5,
        len(df) * 0.45,
    )

    fig, ax = plt.subplots(
        figsize=(10, fig_height),
    )

    left = np.zeros(len(df))

    for i, column in enumerate(columns):

        if normalize:

            values = (
                df[column]
                / df["Total"]
                * 100
            )

        else:

            values = df[column]

        ax.barh(
            df[category],
            values,
            left=left,
            color=colors[i],
            label=labels[i],
        )

        left += values

    # Labels
    style_labels(
        ax,
        xlabel=xlabel,
    )

    # Total labels
    if show_totals:

        if normalize:

            # Draw the absolute number inside the complete 100% bar
            add_bar_labels(
                ax,
                df,
                category,
                value="Total",
                position=100,
                inside=True,
            )

        else:

            add_bar_labels(
                ax,
                df,
                category,
                value="Total",
            )

    # Axes
    if normalize:

        ax.set_xlim(
            0,
            100,
        )

    else:

        ax.set_xlim(
            0,
            df["Total"].max() * x_margin,
        )

    # Style
    if title is not None:

        style_title(
            ax,
            title,
        )

    style_axes(
        ax,
    )

    style_legend(
        ax,
    )

    style_layout(
        fig,
    )

    return fig

def categorical_barplot(
    table,
    category,
    value="Total",
    xlabel="Number of sequences",
    color=COLOR_A,
    top_n=None,
    x_margin=1.05,
    sort=True,
    title=None,
):
    """
    Horizontal barplot for a single-category table.
    """

    df = table.copy()

    # Keep only the Top N categories if requested
    if top_n is not None:
        df = (
            df.sort_values(
                value,
                ascending=False,
            )
            .head(top_n)
        )

    # Sort so that the largest bar appears at the top
    if sort:
        df = df.sort_values(
            value,
            ascending=True,
        )

    fig_height = max(
        5,
        len(df) * 0.45,
    )

    fig, ax = plt.subplots(
        figsize=(10, fig_height),
    )

    ax.barh(
        df[category],
        df[value],
        color=color,
    )

    style_labels(
        ax,
        xlabel=xlabel,
    )

    add_bar_labels(
        ax,
        df,
        category,
        value=value,
    )

    ax.set_xlim(
        0,
        df[value].max() * x_margin,
    )

    if title is not None:
        style_title(
            ax,
            title,
        )
    
    style_axes(ax)

    style_layout(fig)

    return fig

def virus_dynamics_areaplot(
    table,
    category,
    subtype,
    relative=False,
    top_n=15,
    min_count=None,
    min_percentage=None,
    title=None,
):
    """
    Generic virus dynamics area plot.
    """

    matrix = pivot_category(
        table,
        index="Year",
        columns=category,
        values="Count",
    )

    matrix = complete_time_series(
        matrix
    )

    matrix = filter_category_matrix(
        matrix,
        index="Year",
        top_n=top_n,
        min_count=min_count,
        min_percentage=min_percentage,
        subtype_order=(
            subtype is None
            and category == Columns.LINEAGE
        ),
    )

    if relative:
        counts = matrix.drop(columns="Year")

        percentages = (
            counts.div(
                counts.sum(axis=1),
                axis=0,
            ) *100
        )

        matrix = pd.concat(
            [
                matrix[["Year"]],
                percentages,
            ],
            axis=1,
        )

    columns = matrix.columns[1:].tolist()

    if category == Columns.LINEAGE:

        colors = get_lineage_colors(
            columns,
            subtype=subtype,
        )

        legend_title = "Lineage"

    else:

        colors = get_genotype_colors(
            columns,
        )

        legend_title = "Genotype"

    single_year = (
        matrix["Year"].nunique() == 1
    )

    if single_year:

        return temporal_stacked_barplot(
            table=matrix,
            x="Year",
            columns=columns,
            colors=colors,
            labels=columns,
            ylabel=(
                "Relative abundance (%)"
                if relative
                else "Number of sequences"
            ),
            percentage=False,
            show_totals=False,
            legend_title=legend_title,
            legend_large=(
                category == Columns.LINEAGE
            ),
            title=title,
        )

    if category == Columns.LINEAGE:

        return stacked_area_plot_large_legend(
            matrix,
            colors=colors,
            ylabel=(
                "Relative abundance (%)"
                if relative
                else "Number of sequences"
            ),
            percentage=relative,
            legend_title=legend_title,
            title=title,
        )
        
    return stacked_area_plot(
        matrix,
        colors=colors,
        ylabel=(
            "Relative abundance (%)"
            if relative
            else "Number of sequences"
        ),
        percentage=relative,
        legend_title=legend_title,
        title=title,
    )

def stacked_area_plot(
    table,
    colors,
    xlabel="Collection year",
    ylabel="Number of sequences",
    legend_title=None,
    percentage=False,
    figsize=FIGSIZE,
    title=None,
):
    """
    Generic stacked area plot for datasets with a small legend
    (e.g. genotypes).
    """

    fig, ax = plt.subplots(
        figsize=figsize,
    )

    df = table.copy()

    x = df.iloc[:, 0]
    y = df.iloc[:, 1:]

    ax.stackplot(
        x,
        y.T,
        labels=y.columns,
        colors=colors,
        alpha=0.90,
    )

    style_labels(
        ax,
        xlabel,
        ylabel,
    )

    if percentage:
        ax.set_ylim(
            0,
            100,
        )

    format_year_axis(
        ax,
        x,
    )

    style_axes(
        ax,
    )

    if title is not None:
        style_title(
            ax,
            title,
        )

    handles, labels = ax.get_legend_handles_labels()

    if labels:

        style_legend(
            ax,
            title=legend_title,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.10),
            ncol=min(
                LEGEND_COLUMNS,
                len(labels),
            ),
        )

        fig.subplots_adjust(
            bottom=0.20,
        )

    return fig

def stacked_area_plot_large_legend(
    table,
    colors,
    xlabel="Collection year",
    ylabel="Number of sequences",
    legend_title=None,
    percentage=False,
    figsize=(12, 7.5),
    title=None,
):
    """
    Generic stacked area plot for datasets with many categories
    (e.g. viral lineages).
    """

    fig, ax = plt.subplots(
        figsize=figsize,
    )

    df = table.copy()

    x = df.iloc[:, 0]
    y = df.iloc[:, 1:]

    ax.stackplot(
        x,
        y.T,
        labels=y.columns,
        colors=colors,
        alpha=0.90,
    )

    style_labels(
        ax,
        xlabel,
        ylabel,
    )

    if percentage:
        ax.set_ylim(
            0,
            100,
        )

    format_year_axis(
        ax,
        x,
    )

    style_axes(
        ax,
    )

    if title is not None:
        style_title(
            ax,
            title,
        )

    handles, labels = ax.get_legend_handles_labels()

    if labels:

        style_legend(
            ax,
            title=legend_title,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.16),
            ncol=5,
        )

        fig.subplots_adjust(
            bottom=0.34,
        )

    return fig

def histogram_plot(
    data,
    bins=50,
    color=COLOR_A,
    xlabel=None,
    ylabel="Number of sequences",
    figsize=FIGSIZE,
    show_mean=False,
    show_median=False,
    title=None,
):
    """
    Generic histogram.
    """

    fig, ax = plt.subplots(
        figsize=figsize,
    )

    ax.hist(
        data,
        bins=bins,
        color=color,
        edgecolor="white",
        linewidth=0.8,
        alpha=0.90,
    )

    if show_mean:
        mean = np.mean(data)

        ax.axvline(
            mean,
            color=REFERENCE_COLOR,
            linestyle=MEAN_LINESTYLE,
            linewidth=REFERENCE_LINEWIDTH,
            label=f"Mean = {mean:.1f}",
        )

    if show_median:
        median = np.median(data)

        ax.axvline(
            median,
            color=REFERENCE_COLOR,
            linestyle=MEDIAN_LINESTYLE,
            linewidth=REFERENCE_LINEWIDTH,
            label=f"Median = {median:.1f}",
        )

    if show_mean or show_median:
        style_legend(
            ax,
        )

    if title is not None:
        style_title(
            ax,
            title,
        )

    style_labels(
        ax,
        xlabel,
        ylabel,
    )
    
    style_axes(ax)

    style_layout(fig)

    return fig

def grouped_histogram_plot(
    dataframe,
    value_column,
    group_column,
    groups,
    colours,
    bins=50,
    xlabel=None,
    ylabel="Number of sequences",
    figsize=(12,8),
    show_mean=False,
    show_median=False,
    title=None,
):

    fig, axes = plt.subplots(
        len(groups),
        1,
        figsize=figsize,
        sharex=True,
    )

    if len(groups) == 1:
        axes = [axes]

    for ax, group, colour in zip(
        axes,
        groups,
        colours,
    ):

        values = (
            pd.to_numeric(
                dataframe.loc[
                    dataframe[group_column] == group,
                    value_column,
                ],
                errors="coerce",
            )
            .dropna()
        )

        ax.hist(
            values,
            bins=bins,
            color=colour,
            edgecolor="white",
            linewidth=0.8,
            alpha=0.90,
        )

        if show_mean:
            mean = values.mean()

            ax.axvline(
                mean,
                color=REFERENCE_COLOR,
                linestyle=MEAN_LINESTYLE,
                linewidth=REFERENCE_LINEWIDTH,
                label=f"Mean = {mean:.1f}",
            )

        if show_median:
            median = values.median()

            ax.axvline(
                median,
                color=REFERENCE_COLOR,
                linestyle=MEDIAN_LINESTYLE,
                linewidth=REFERENCE_LINEWIDTH,
                label=f"Median = {median:.1f}",
            )

        style_title(
            ax,
            f"RSV-{group}",
        )

        style_labels(
            ax,
            None,
            ylabel,
        )

        style_axes(ax)

        if show_mean or show_median:
            style_legend(ax)

    axes[-1].set_xlabel(xlabel)

    if title is not None:
        fig.suptitle(
            title,
            fontsize=22,
            fontweight="bold",
        )

    style_layout(fig)

    return fig

def boxplot_plot(
    table,
    group,
    value,
    colors,
    xlabel=None,
    ylabel=None,
    figsize=BOXPLOT_FIGSIZE,
    title=None,
):
    """
    Generic boxplot.
    """

    fig, ax = plt.subplots(figsize=figsize)

    groups = [
        g for g in table[group].dropna().unique()
    ]

    data = [
        table.loc[
            table[group] == g,
            value,
        ]
        for g in groups
    ]

    bp = ax.boxplot(
        data,
        tick_labels=groups,
        patch_artist=True,
        showfliers=True,
        flierprops={
            "marker": "o",
            "markersize": 3,
            "markerfacecolor": OTHER_COLOR,
            "markeredgecolor": OTHER_COLOR,
            "alpha": 0.5,
        }
    )

    # Style
    for patch, group in zip(bp["boxes"], groups):
        patch.set_facecolor(
            colors.get(
                group,
                OTHER_COLOR,
            )
        )

        patch.set_alpha(0.9)
        patch.set_linewidth(1.2)
    
    for median in bp["medians"]:

        median.set_color("black")
        median.set_linewidth(1.5)

    for whisker in bp["whiskers"]:

        whisker.set_linewidth(1.2)

    for cap in bp["caps"]:

        cap.set_linewidth(1.2)

    if title is not None:
        style_title(
            ax,
            title,
        )

    style_labels(
        ax,
        xlabel,
        ylabel,
    )
    
    style_axes(ax)

    style_layout(fig)

    return fig

def grouped_categorical_barplot(
    table,
    category,
    hue,
    value,
    colors,
    xlabel,
    order=PROTEIN_ORDER,
    x_margin=1.05,
    bar_height=BAR_HEIGHT,
    ax=None,
    show_legend=True,
    title=None,
):
    """
    Grouped horizontal barplot.
    """

    created_figure = ax is None

    pivot = (
        table
        .pivot(
            index=category,
            columns=hue,
            values=value,
        )
    )

    if order is not None:
        pivot = (
            pivot
            .reindex(order)
            .dropna(how="all")
        )

    pivot = pivot.iloc[::-1]

    if ax is None:

        fig_height = max(
            5,
            len(pivot) * 0.45,
        )

        fig, ax = plt.subplots(
            figsize=(10, fig_height),
        )

    else:

        fig = ax.figure

    y = np.arange(len(pivot))

    groups = list(
        pivot.columns
    )

    offsets = np.linspace(
        -bar_height / 2,
        bar_height / 2,
        len(groups),
    )

    xmax = pivot.max().max()

    for group, offset in zip(groups, offsets):

        values = pivot[group]

        ax.barh(
            y + offset,
            values,
            height=bar_height,
            color=colors[group],
            label= GROUP_LABELS.get(
                group,
                group,
            )
        )

        add_grouped_bar_labels(
            ax=ax,
            y=y,
            values=values,
            offset=offset,
            xmax=xmax,
        )
    
    ax.set_yticks(y)
    ax.set_yticklabels(
        pivot.index,
    )

    style_labels(
        ax,
        xlabel=xlabel,
    )

    ax.set_xlim(
        0,
        pivot.max().max() * x_margin,
    )

    if title is not None:
        style_title(
            ax,
            title,
        )

    style_axes(ax)

    if show_legend:
        style_legend(
            ax,
        )

    if created_figure:
        style_layout(fig)

    return fig

def categorical_heatmap(
    table,
    cmap,
    percentages=None,
    row_totals=None,
    xlabel=None,
    ylabel=None,
    row_total_label=None,
    colorbar_label=None,
    title=None,
    figsize=HEATMAP_FIGSIZE,
    annotate=True,
    show_row_totals=True,
    column_labels=None,
    vmin=None,
    vmax=None,
    center=None,
    annotate_zero=False,
):
    """
    Generic categorical heatmap.

    If 'percentages' is provided, it is used for the heatmap colours,
    while 'table' provides the absolute counts shown as annotations.
    Otherwise, row-wise percentages are calculated from 'table'.
    """

    # Data
    index = table.columns[0]

    counts = (
        table
        .set_index(index)
    )

    if row_totals is None:
        row_totals = counts.sum(axis=1)

    if percentages is None:

        percentages = (
            counts.div(
                row_totals,
                axis=0,
            ) * 100
        )

    else:

        percentages = (
            percentages
            .set_index(index)
            .reindex(
                index=counts.index,
                columns=counts.columns,
            )
        )

    # Figure
    fig = plt.figure(
        figsize=figsize,
    )

    if show_row_totals:

        gs = GridSpec(
            1,
            3,
            width_ratios=HEATMAP_BAR_RATIO,
            wspace=0.05,
        )

        ax_labels = fig.add_subplot(gs[0])
        ax = fig.add_subplot(gs[1], sharey=ax_labels)
        cax = fig.add_subplot(gs[2])

    else:

        gs = GridSpec(
            1,
            2,
            width_ratios=[30, 1],
            wspace=0.05,
        )

        ax = fig.add_subplot(gs[0])
        cax = fig.add_subplot(gs[1])

    # Row labels
    if show_row_totals:

        add_heatmap_row_labels(
            ax_labels,
            counts.index,
            row_totals,
            title=(
                row_total_label
                if row_total_label is not None
                else ylabel
            ),
        )

    # Heatmap
    if vmin is None:
        vmin = np.nanmin(percentages.values)
    if vmax is None:
        vmax = np.nanmax(percentages.values)
    if center is not None:
        norm = TwoSlopeNorm(
            vmin=vmin,
            vcenter=center,
            vmax=vmax,
        )

        image = ax.imshow(
            percentages.values,
            cmap=cmap,
            aspect="auto",
            interpolation="nearest",
            norm=norm,
        )
    else:
        image = ax.imshow(
            percentages.values,
            cmap=cmap,
            aspect="auto",
            interpolation="nearest",
            vmin=vmin,
            vmax=vmax,
        )


    # Axes
    ax.set_yticks(
        np.arange(
            counts.shape[0]
        )
    )

    if show_row_totals:

        ax.tick_params(
            axis="y",
            left=False,
            labelleft=False,
        )

    else:

        ax.set_yticklabels(
            counts.index,
        )

        ax.tick_params(
            axis="y",
            left=False,
            labelleft=True,
            pad=8,
        )

    ax.set_xticks(
        np.arange(
            counts.shape[1]
        )
    )

    labels = (
        column_labels
        if column_labels is not None
        else counts.columns
    )

    ax.set_xticklabels(
        labels,
        rotation=45,
        ha="right",
    )

    ax.set_xticks(
        np.arange(
            counts.shape[1] + 1
        ) - 0.5,
        minor=True,
    )

    ax.set_yticks(
        np.arange(
            counts.shape[0] + 1
        ) - 0.5,
        minor=True,
    )

    ax.grid(
        which="minor",
        color="white",
        linewidth=0.8,
    )

    ax.tick_params(
        which="minor",
        bottom=False,
        left=False,
    )

    if xlabel is not None:

        style_labels(
            ax,
            xlabel=xlabel,
        )

    # Cell annotations
    if annotate:

        for i in range(counts.shape[0]):

            for j in range(counts.shape[1]):

                value = counts.iloc[i, j]

                if float(value).is_integer():
                    label = f"{int(value)}"
                else:
                    label = f"{value:.2f}".rstrip("0").rstrip(".")

                if value == 0 and not annotate_zero:
                    continue

                colour = (
                    "white"
                    if image.norm(
                        percentages.iloc[i, j]
                    ) >= 0.5
                    else "black"
                )

                ax.text(
                    j,
                    i,
                    label,
                    ha="center",
                    va="center",
                    fontsize=HEATMAP_ANNOTATION_SIZE,
                    color=colour,
                    fontweight="bold",
                )

    # Colorbar
    cbar = fig.colorbar(
        image,
        cax=cax,
    )

    style_colorbar(
        cbar,
        label=colorbar_label,
    )
    
    # Title
    if title is not None:

        style_title(
            ax,
            title,
        )

    style_layout(
        fig,
    )

    return fig

def lollipop_plot(
    table,
    x,
    y,
    bubble_size,
    protein_length,
    label=None,
    colour="#4E79A7",
    xlabel=None,
    ylabel=None,
    title=None,
    figsize=LOLLIPOP_FIGSIZE,
    label_threshold=LOLLIPOP_LABEL_THRESHOLD,
    track=None,
    annotation_legend=None,
):
    """
    Generic lollipop plot.
    """
    # Data
    df = (
        table
        .sort_values(x)
        .reset_index(drop=True)
    )

    subtype = table[Columns.SUBTYPE].iat[0]

    # Bubble sizes
    sizes = np.interp(
        df[bubble_size],
        (
            df[bubble_size].min(),
            df[bubble_size].max(),
        ),
        (
            LOLLIPOP_MARKER_MIN,
            LOLLIPOP_MARKER_MAX,
        ),
    )

    # Figure
    fig = plt.figure(
        figsize=figsize,
    )

    gs = GridSpec(
        2,
        1,
        height_ratios=[8, 1.6],
        hspace=0.03,
    )

    ax = fig.add_subplot(gs[0])

    ax_protein = fig.add_subplot(
        gs[1],
        sharex=ax,
    )

    # Lollipop stems
    ax.vlines(
        x=df[x],
        ymin=0,
        ymax=df[y],
        color=LOLLIPOP_STEM_COLOUR,
        linewidth=1.2,
        zorder=1,
    )

    # Lollipop heads
    bubble_color = SUBTYPE_COLORS[subtype]
    
    scatter = ax.scatter(
        df[x],
        df[y],
        s=sizes,
        color=bubble_color,
        edgecolor="black",
        linewidth=0.7,
        zorder=3,
    )

    # Bubble size and protein annotation legends
    legend_values = [10, 25, 50, 75, 100]

    legend_sizes = np.interp(
        legend_values,
        (
            df[bubble_size].min(),
            df[bubble_size].max(),
        ),
        (
            LOLLIPOP_MARKER_MIN,
            LOLLIPOP_MARKER_MAX,
        ),
    )

    handles = [
        plt.scatter(
            [],
            [],
            s=size,
            color=colour,
            edgecolor="black",
        )
        for size in legend_sizes
    ]

    bubble_legend = style_bubble_legend(
        ax,
        handles=handles,
        labels=[f"{v}%" for v in legend_values],
    )

    ax.add_artist(bubble_legend)

    if annotation_legend is not None:
        annotation_legend(ax)
    
    # Labels
    if label is not None:

        selected = (
            df[
                df[y] >= label_threshold
            ]
            .sort_values(x)
            .reset_index(drop=True)
        )

        labels = []

        for _, row in selected.iterrows():

            labels.append(
                {
                    "x": row[x],
                    "y": row[y],
                    "text": row[label],
                }
            )

        arrange_lollipop_labels(
            ax=ax,
            labels=labels,
        )

    # Main axis
    style_axes(ax)

    style_labels(
        ax,
        ylabel=ylabel,
    )

    ax.yaxis.set_major_locator(
        MaxNLocator(integer=True)
    )

    ax.set_xlim(
        1,
        protein_length + 15,
    )

    ax.set_ylim(
        0,
        df[y].max() + 1,
    )

    #X-axis ticks accorfing to proteing length
    if protein_length <= 350:
        tick_step = 25
    elif protein_length <= 700:
        tick_step = 50
    elif protein_length <= 1500:
        tick_step = 100
    else:
        tick_step = 250

    ax.set_xticks(
        np.arange(
            tick_step,
            protein_length + 1,
            tick_step,
        )
    )

    ax.grid(
        axis="x",
        alpha=0.15,
        linewidth=0.6,
    )

    # Protein track
    if track is not None:
        track(ax_protein)

    # Protein axis
    ax_protein.set_xlim(
        1,
        protein_length + 5,
    )

    ax_protein.set_ylim(
        -1,
        1,
    )

    ax_protein.set_yticks([])

    style_labels(
        ax_protein,
        xlabel=xlabel,
    )

    ax_protein.spines["left"].set_visible(False)
    ax_protein.spines["right"].set_visible(False)
    ax_protein.spines["top"].set_visible(False)

    ax_protein.tick_params(
        left=False,
        labelleft=False,
    )

    # Title
    if title is not None:

        style_title(
            ax,
            title,
        )

    style_layout(fig)

    return fig

#############################
# Geography
#############################
def draw_geography(
    ax,
    gdf,
    value_column="Total",
):
    """
    Draw the geographic distribution of genomes across Spain.
    """
    # Map
    gdf.plot(
        column=value_column,
        cmap=MAP_CMAP,
        norm=LogNorm(
            vmin=1,
            vmax=max(
                1,
                gdf[value_column].max(),
            ),
        ),
        linewidth=MAP_LINEWIDTH,
        edgecolor=MAP_EDGECOLOR,
        legend=False,
        ax=ax,
        missing_kwds={
            "color": MAP_MISSING_COLOR,
        },
    )
    
    # Country outline
    country = unary_union(
        gdf.geometry
    )

    gpd.GeoSeries(
        [country]
    ).boundary.plot(
        ax=ax,
        color="#555555",
        linewidth=1.2,
        zorder=5,
    )

    # Labels
    draw_map_labels(
        ax=ax,
        gdf=gdf,
        value_column=value_column,
    )

    style_map(
        ax,
    )

def geography_map(
    table,
    title=None,
):

    # Read map
    gdf = gpd.read_file(
        SPAIN_CCAA_GEOJSON,
    )

    # Community names
    table = table.copy()

    table[
        Columns.COMMUNITY
    ] = (
        table[
            Columns.COMMUNITY
        ].replace(
            AUTONOMOUS_COMMUNITY_MAP
        )
    )

    # Merge metadata
    value_column = "Total"

    gdf = gdf.merge(
        table,
        left_on="rotulo",
        right_on=Columns.COMMUNITY,
        how="left",
    )

    gdf[value_column] = (
        gdf[value_column]
        .fillna(0)
    )

    # Figure
    fig, ax = plt.subplots(
        figsize=MAP_FIGSIZE,
    )

    #Draw maps
    draw_geography(
        ax=ax,
        gdf=gdf,
        value_column=value_column,
    )

    #Title
    style_title(
        ax,
        title,
    )

    return fig

def _available_rsv_subtypes(
    table,
):
    """
    Return plotting configuration for RSV subtypes
    available in the table.
    """

    subtype_config = {
        "A": (COLOR_A, "RSV-A"),
        "B": (COLOR_B, "RSV-B"),
    }

    subtypes = [
        subtype
        for subtype in subtype_config
        if subtype in table.columns
    ]

    colors = [
        subtype_config[subtype][0]
        for subtype in subtypes
    ]

    labels = [
        subtype_config[subtype][1]
        for subtype in subtypes
    ]

    return subtypes, colors, labels

def geography_community_barplot(table):
    """
    Distribution of sequences by autonomous community.
    """

    return categorical_stacked_barplot(
        table,
        category=Columns.COMMUNITY,
        columns=["A", "B"],
        colors=[COLOR_A, COLOR_B],
        labels=["RSV-A", "RSV-B"],
        title="RSV distribution by Community"
    )

def sampling_strategy_barplot(table):
    """
    Distribution sampling strategy by autonomous community.
    """

    table = table[
        ~table[
            Columns.COMMUNITY
        ].isin(
            [
                "Unknown",
                "NA",
            ]
        )
    ]

    return categorical_stacked_barplot(
        table,
        category=Columns.COMMUNITY,
        columns=["Sentinel", "Non-sentinel"],
        colors=[COLOR_SENTINEL, COLOR_NON_SENTINEL],
        labels=["Sentinel", "Non-sentinel"],
        xlabel="Number of Sequences",
        title="Sampling strategy by Community",
    )

def sequence_length_histogram(dataset):

    return grouped_histogram_plot(
        dataframe=dataset.df,
        value_column=Columns.SEQUENCE_LENGTH,
        group_column=Columns.SUBTYPE,
        groups=["A", "B"],
        colours=[COLOR_A, COLOR_B],
        bins=50,
        xlabel="Sequence length (bp)",
        show_median=True,
        title="Sequence length distribution",
    )

#############################
# Laboratory
#############################

def laboratory_originating_barplot(table):

    return categorical_stacked_barplot(
        table,
        category=Columns.ORIGINATING_LAB,
        columns=["A", "B"],
        colors=[COLOR_A, COLOR_B],
        labels=["RSV-A", "RSV-B"],
        title="RSV distribution by Originating Lab"
    )

def laboratory_institution_barplot(table):

    return categorical_stacked_barplot(
        table,
        category=Columns.INSTITUTION,
        columns=["A", "B"],
        colors=[COLOR_A, COLOR_B],
        labels=["RSV-A", "RSV-B"],
        title="RSV distribution by Institution"
    )


def laboratory_submitter_barplot(table):

    return categorical_stacked_barplot(
        table,
        category=Columns.SUBMITTER_LAB,
        columns=["A", "B"],
        colors=[COLOR_A, COLOR_B],
        labels=["RSV-A", "RSV-B"],
        title="RSV distribution by Submitter Institution"
    )


#############################
# Patient
#############################

def patient_status_barplot(table):
    """
    Patient status distribution.
    """

    return categorical_stacked_barplot(
        table,
        category=Columns.STATUS_NORMALIZED,
        xlabel="Number of sequences",
        columns=["A", "B"],
        colors=[COLOR_A, COLOR_B],
        labels=["RSV-A", "RSV-B"],
        title="RSV distribution by Patient Status"
    )

def patient_gender_barplot(table):

    return categorical_stacked_barplot(
        table,
        category=Columns.GENDER_NORMALIZED,
        xlabel="Number of sequences",
        columns=["A", "B"],
        colors=[COLOR_A, COLOR_B],
        labels=["RSV-A", "RSV-B"],
        title="RSV distribution by Gender"
    )

def patient_age_groups_barplot(table):
    """
    Age distribution by age groups.
    """
    
    return categorical_stacked_barplot(
        table,
        category="Age group",
        xlabel="Number of sequences",
        columns=["A", "B"],
        colors=[COLOR_A, COLOR_B],
        labels=["RSV-A", "RSV-B"],
        title="RSV distribution by Age"
    )

#############################
# Temporal
#############################

def temporal_year_barplot(
    table,
):
    """
    Annual number of sequences.
    """

    subtypes, colors, labels = (
        _available_rsv_subtypes(table)
    )

    return temporal_stacked_barplot(
        table,
        x="Year",
        columns=subtypes,
        colors=colors,
        labels=labels,
        ylabel="Number of sequences",
        title="RSV distribution through time",
    )

def temporal_year_percentage(
    table,
):
    """
    Relative annual distribution of RSV-A and RSV-B.
    """

    subtypes, colors, labels = (
        _available_rsv_subtypes(table)
    )

    return temporal_stacked_barplot(
        table,
        x="Year",
        columns=subtypes,
        colors=colors,
        labels=labels,
        ylabel="Relative abundance (%)",
        percentage=True,
        show_totals=True,
        title="RSV abundance distribution through time",
    )

#######################
# Sequencing
#######################

def sequencing_technology_barplot(table):

    return categorical_stacked_barplot(
        table,
        category=Columns.TECHNOLOGY_NORMALIZED,
        xlabel="Number of sequences",
        columns=["A", "B"],
        colors=[COLOR_A, COLOR_B],
        labels=["RSV-A", "RSV-B"],
        title="RSV distribution by Sequencing technology"
    )


def sequencing_assembly_barplot(table):

    return categorical_stacked_barplot(
        table,
        category=Columns.ASSEMBLY_NORMALIZED,
        xlabel="Number of sequences",
        columns=["A", "B"],
        colors=[COLOR_A, COLOR_B],
        labels=["RSV-A", "RSV-B"],
        title="RSV distribution by Assembly method"
    )

#############################
# Virus lineage & genotype overview
#############################

def virus_lineage_barplot(
    table,
    subtype,
):
    """
    Distribution of viral lineages.
    """

    color = (
        COLOR_A
        if subtype == "A"
        else COLOR_B
    )

    title = (
        "RSV-A lineage distribution"
        if subtype == "A"
        else "RSV-B lineage distribution"
    )

    return categorical_barplot(
        table,
        category=Columns.LINEAGE,
        color=color,
        xlabel="Number of sequences",
        top_n=15,
        x_margin=1.03,
        title=title,
    )

def virus_genotype_barplot(
    table,
    subtype,
):
    """
    Distribution of viral genotypes.
    """

    color = (
        COLOR_A
        if subtype == "A"
        else COLOR_B
    )

    title = (
        "RSV-A genotype distribution"
        if subtype == "A"
        else "RSV-B genotype distribution"
    )


    return categorical_barplot(
        table,
        category=Columns.GENOTYPE,
        color=color,
        xlabel="Number of sequences",
        x_margin=1.03,
        title=title,
    )

#############################
# Virus lineage & genotype dynamics
#############################

def lineage_year_areaplot(
    table,
    subtype,
    relative=False,
):

    if subtype is None:
        virus = "RSV"
    else:
        virus = f"RSV-{subtype}"

    title = f"{virus} lineage dynamics"

    if relative:
        title += " (%)"

    return virus_dynamics_areaplot(
        table=table,
        category=Columns.LINEAGE,
        subtype=subtype,
        relative=relative,
        title=title,
    )

def genotype_year_areaplot(
    table,
    subtype,
    relative=False,
):

    if subtype is None:
        virus = "RSV"
    else:
        virus = f"RSV-{subtype}"

    title = f"{virus} genotype dynamics"

    if relative:
        title += " (%)"

    return virus_dynamics_areaplot(
        table=table,
        category=Columns.GENOTYPE,
        subtype=subtype,
        relative=relative,
        title=title,
    )

def community_lineage_heatmap(
    table,
    subtype,
):
    """
    Distribution of major viral lineages across autonomous communities.
    """

    return categorical_heatmap(
        table=table,
        cmap=(
            LINEAGE_CMAP_A
            if subtype == "A"
            else LINEAGE_CMAP_B
        ),
        colorbar_label="Relative abundance (%)",
        title=f"RSV-{subtype} major lineages distribution",
    )

def transmission_lineage_heatmap(
    table,
    percentages,
    lineage_totals,
    subtype,
    threshold,
):
    """
    Transmission clustering of major viral lineages through time.
    """

    return categorical_heatmap(
        table,
        cmap=(
            LINEAGE_CMAP_A
            if subtype == "A"
            else LINEAGE_CMAP_B
        ),
        percentages=percentages,
        row_totals=lineage_totals,
        xlabel="Collection year",
        row_total_label="Total\nsequences",
        colorbar_label="Yearly clustered sequences (%)",
        title=(
            f"Transmission clustering by lineage "
            f"(RSV-{subtype}, {threshold} SNP threshold)"
        ),
    )

#############################
# Mutations overview
#############################

def mutation_histogram(
    mutation_counts,
    subtype=None,
):
    """
    Histogram of amino acid substitutions per genome.
    """

    df = mutation_counts.copy()

    if subtype is not None:
        df = df[
            df[Columns.SUBTYPE] == subtype
        ]

    data = df[
        MutationColumns.COUNT
    ]

    if subtype == "A":
        color = COLOR_A

    elif subtype == "B":
        color = COLOR_B

    else:
        color = OTHER_COLOR

    return histogram_plot(
        data=data,
        color=color,
        bins=50,
        xlabel="Number of amino acid substitutions",
        ylabel="Number of sequences",
        show_mean=True,
        show_median=True,
        title="AA substitutions distribution"
    )

def mutation_boxplot(
    mutation_counts,
):
    """
    Distribution of amino acid substitutions by subtype.
    """

    return boxplot_plot(
        table=mutation_counts,
        group=Columns.SUBTYPE,
        value=MutationColumns.COUNT,
        colors=SUBTYPE_COLORS,
        xlabel="Subtype",
        ylabel="Number of amino acid substitutions",
        title="AA substitutions distribution by RSV subtype"
    )

def protein_mutation_events(
    protein_summary,
    ax=None,
):
    return grouped_categorical_barplot(
        table=protein_summary,
        category=MutationColumns.PROTEIN,
        hue=Columns.SUBTYPE,
        value=MutationColumns.TOTAL_EVENTS,
        colors=SUBTYPE_COLORS,
        xlabel="Total mutation events",
        order=PROTEIN_ORDER,
        ax=ax,
        show_legend=False,
    )

def protein_distinct_substitutions(
    protein_summary,
    ax=None,
):
    return grouped_categorical_barplot(
        table=protein_summary,
        category=MutationColumns.PROTEIN,
        hue=Columns.SUBTYPE,
        value=MutationColumns.DISTINCT_SUBSTITUTIONS,
        colors=SUBTYPE_COLORS,
        xlabel="Number of distinct AA substitutions",
        order=PROTEIN_ORDER,
        ax=ax,
        show_legend=True,
    )

def protein_mutation_summary_plot(
    protein_summary,
):

    fig, axes = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(10, 10),
        sharey=True,
        constrained_layout=True,
    )

    protein_mutation_events(
        protein_summary,
        ax=axes[0],
    )

    protein_distinct_substitutions(
        protein_summary,
        ax=axes[1],
    )

    style_title(
        axes[0],
        "A. Total mutation events",
    )

    style_title(
        axes[1],
        "B. Number of distinct AA substitutions",
    )

    return fig

def protein_lineage_heatmap(
    table,
    percentages,
    summary,
    subtype,
    title,
):
    """
    Protein-by-lineage heatmap.
    """

    # Column labels
    lineage_order = table.columns[1:]

    genomes = (
        summary.loc[
            summary[Columns.SUBTYPE] == subtype,
            [
                Columns.LINEAGE,
                MutationColumns.GENOMES,
            ],
        ]
        .drop_duplicates()
        .set_index(Columns.LINEAGE)
        .loc[lineage_order]
    )

    column_labels = [
        f"{lineage}\n({n:,})"
        for lineage, n in zip(
            genomes.index,
            genomes[MutationColumns.GENOMES],
        )
    ]

    # Plot
    return categorical_heatmap(
        table=table,
        percentages=percentages,
        cmap=(
            LINEAGE_CMAP_A
            if subtype == "A"
            else LINEAGE_CMAP_B
        ),
        xlabel=None,
        ylabel="Protein",
        colorbar_label="Relative abundance (%)",
        show_row_totals=False,
        column_labels=column_labels,
        title=title,
    )

def protein_lineage_comparative_heatmap(
    table,
    summary,
    subtype,
    title,
):
    """
    Protein-by-lineage comparative heatmap.

    Colours represent robust z-scores of the relative
    distribution of distinct amino acid substitutions
    across lineages, calculated independently per protein.
    """

    # Column labels
    lineage_order = table.columns[1:]

    genomes = (
        summary.loc[
            summary[Columns.SUBTYPE] == subtype,
            [
                Columns.LINEAGE,
                MutationColumns.GENOMES,
            ],
        ]
        .drop_duplicates()
        .set_index(Columns.LINEAGE)
        .loc[lineage_order]
    )

    column_labels = [
        f"{lineage}\n({n:,})"
        for lineage, n in zip(
            genomes.index,
            genomes[MutationColumns.GENOMES],
        )
    ]

    return categorical_heatmap(
        table=table,
        percentages=table,
        cmap="RdBu_r",
        xlabel=None,
        ylabel="Protein",
        colorbar_label=(
            "Deviation from typical\n"
            "protein distribution\n"
            "(robust z-score)"
        ),
        show_row_totals=False,
        column_labels=column_labels,
        title=title,
        vmin=-3,
        vmax=3,
        center=0,
        annotate_zero=True,
    )

def protein_lollipop(
    table,
    min_frequency=LOLLIPOP_MIN_FREQUENCY,
):
    """
    Protein lollipop plot.
    """

    # Keep only frequent substitutions
    table = (
        table[
            table[MutationColumns.MOST_COMMON_FREQUENCY]
            >= min_frequency
        ]
        .copy()
    )
    
    lineage = (
        table[Columns.LINEAGE]
        .iat[0]
    )

    subtype = (
        table[Columns.SUBTYPE]
        .iat[0]
    )

    protein = (
        table[MutationColumns.PROTEIN]
        .iat[0]
    )

    # Protein track function
    def track(ax):
        protein_track(
            ax=ax,
            subtype=subtype,
            protein=protein,
        )
    
    # Protein legend
    def annotation_legend(ax):

        return style_feature_legend(
            ax=ax,
            subtype=subtype,
            protein=protein,
        )

    return lollipop_plot(
        table=table,
        x=MutationColumns.POSITION,
        y=MutationColumns.DISTINCT_SUBSTITUTIONS,
        bubble_size=MutationColumns.MOST_COMMON_FREQUENCY,
        protein_length=PROTEIN_LENGTHS[subtype][protein],
        label=MutationColumns.LABEL,
        colour=(
            LINEAGE_CMAP_A(0.7)
            if subtype == "A"
            else LINEAGE_CMAP_B(0.7)
        ),
        xlabel="Amino acid position",
        ylabel="Number of distinct substitutions",
        title=(
            f"RSV-{subtype} {lineage} "
            f"{protein} protein mutation hotspots"
        ),
        track=track,
        annotation_legend=annotation_legend,
    )

#############################
# Clustering plots
#############################
def cluster_year_barplot(
    table,
    subtype,
    threshold,
):
    """
    Relative distribution of clustered and singleton sequences
    through time.
    """

    df = (
        table[
            table["Subtype"] == subtype
        ]
        .copy()
    )

    color = (
        COLOR_A
        if subtype == "A"
        else COLOR_B
    )

    return temporal_stacked_barplot(
        df,
        x=Columns.COLLECTION_YEAR,
        columns=[
            "Clustered",
            "Singleton",
        ],
        colors=[
            color,
            SINGLETON_COLOR,
        ],
        labels=[
            "Clustered",
            "Singleton",
        ],
        ylabel="Sequences (%)",
        percentage=True,
        show_totals=True,
        title=(
            f"Transmission clustering through time "
            f"(RSV-{subtype}, {threshold} SNP threshold)"
        ),
    )

def cluster_community_barplot(
    table,
    subtype,
    threshold,
):
    """
    Relative distribution of clustered and singleton sequences
    by autonomous community.
    """

    df = (
        table[
            table["Subtype"] == subtype
        ]
        .copy()
    )

    clustered_color = (
        COLOR_A
        if subtype == "A"
        else COLOR_B
    )

    return categorical_stacked_barplot(
        df,
        category=Columns.COMMUNITY,
        columns=[
            "Clustered",
            "Singleton",
        ],
        colors=[
            clustered_color,
            SINGLETON_COLOR,
        ],
        labels=[
            "Clustered",
            "Singleton",
        ],
        xlabel="Sequences (%)",
        normalize=True,
        show_totals=True,
        title=f"Transmission clustering by Community (RSV-{subtype}, {threshold} SNP threshold)"
    )