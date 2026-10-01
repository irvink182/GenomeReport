"""
Cluster statistics and summary tables.
"""

import pandas as pd

from .config import (
    Columns,
    ClusterMetrics,
    ClusterColumns,
)

SINGLETON = "Singleton"

###############################################################################
# Utilities
###############################################################################

def cluster_columns(dataset):
    """
    Return all cluster columns present in the dataset.
    """

    columns = [
        column
        for column in dataset.df.columns
        if column.startswith(ClusterColumns.PREFIX)
    ]

    return sorted(
        columns,
        key=lambda x: int(
            x.replace(
                ClusterColumns.PREFIX,
                "",
            )
        ),
    )


def cluster_thresholds(dataset):
    """
    Return all SNP thresholds available in the dataset.
    """

    return [
        int(
            column.replace(
                ClusterColumns.PREFIX,
                "",
            )
        )
        for column in cluster_columns(dataset)
    ]


def cluster_column(threshold):
    """
    Return the cluster column name for a SNP threshold.
    """

    return (
        f"{ClusterColumns.PREFIX}{threshold}"
    )

###############################################################################
# Summary
###############################################################################

def has_clustering_data(
    dataset,
):
    """
    Check whether clustering information is available.
    """

    thresholds = cluster_thresholds(
        dataset
    )

    if not thresholds:
        return False

    return any(
        dataset.df[
            cluster_column(threshold)
        ].notna().any()
        for threshold in thresholds
    )

def cluster_summary(
    dataset,
    groupby=None,
):
    """
    Summarise transmission clusters for each SNP threshold.
    """

    df = dataset.df.copy()

    if df.empty:
        return pd.DataFrame()

    # Groups
    if groupby is None:

        groups = [None]

    else:

        groups = (
            df[groupby]
            .dropna()
            .sort_values()
            .unique()
        )

    # Summary
    summary_tables = []

    for group in groups:

        if group is None:

            subset = df.copy()

        else:

            subset = df[
                df[groupby] == group
            ].copy()

        summary = {}

        for threshold in cluster_thresholds(dataset):

            column = cluster_column(threshold)

            # Split transmission clusters / singletons
            transmission = subset[
                subset[column] != SINGLETON
            ]

            singletons = subset[
                subset[column] == SINGLETON
            ]

            # Cluster sizes
            cluster_sizes = (
                transmission
                .groupby(column)
                .size()
            )

            clustered_sequences = len(
                transmission
            )

            singleton_sequences = len(
                singletons
            )

            clustered_sequences_pct = round(
                100 * clustered_sequences / len(subset),
                2,
            )

            singleton_sequences_pct = round(
                100 * singleton_sequences / len(subset),
                2,
            )

            # Cluster statistics
            if len(cluster_sizes) > 0:

                largest = int(
                    cluster_sizes.max()
                )

                mean = round(
                    cluster_sizes.mean(),
                    2,
                )

                median = round(
                    cluster_sizes.median(),
                    2,
                )

            else:

                largest = 0
                mean = 0
                median = 0

            # Store summary
            summary[f"{threshold} SNP"] = {

                ClusterMetrics.SEQUENCES:
                    len(subset),

                ClusterMetrics.TRANSMISSION_CLUSTERS:
                    transmission[column].nunique(),

                ClusterMetrics.CLUSTERED_SEQUENCES:
                    clustered_sequences,

                ClusterMetrics.CLUSTERED_SEQUENCES_PERCENT:
                    clustered_sequences_pct,

                ClusterMetrics.SINGLETON_SEQUENCES:
                    singleton_sequences,

                ClusterMetrics.SINGLETON_SEQUENCES_PERCENT:
                    singleton_sequences_pct,

                ClusterMetrics.LARGEST_CLUSTER_SIZE:
                    largest,

                ClusterMetrics.MEAN_CLUSTER_SIZE:
                    mean,

                ClusterMetrics.MEDIAN_CLUSTER_SIZE:
                    median,

            }

        summary = (
            pd.DataFrame(summary)
            .reset_index()
            .rename(
                columns={
                    "index": "Metric",
                }
            )
        )

        if groupby is not None:

            summary.insert(
                0,
                groupby,
                group,
            )

        summary_tables.append(
            summary
        )

    return (
        pd.concat(
            summary_tables,
            ignore_index=True,
        )
    )

###############################################################################
# Generic tables
###############################################################################

def _cluster_by_variable(
    dataset,
    variable,
    groupby=None,
):
    """
    Generate a long-format table describing the distribution of
    transmission clusters across a metadata variable.
    """

    df = dataset.df.copy()

    if df.empty:
        return pd.DataFrame()

    # Groups
    if groupby is None:

        groups = [None]

    else:

        groups = (
            df[groupby]
            .dropna()
            .sort_values()
            .unique()
        )

    # Generate tables
    results = []

    for threshold in cluster_thresholds(dataset):

        cluster = cluster_column(threshold)

        for group in groups:

            if group is None:

                subset = df.copy()

            else:

                subset = df[
                    df[groupby] == group
                ].copy()

            # Group by cluster and metadata
            table = (
                subset
                .groupby(
                    [
                        cluster,
                        variable,
                    ],
                    dropna=False,
                )
                .size()
                .reset_index(
                    name="Sequences",
                )
            )

            table = table.rename(
                columns={
                    cluster: "Cluster",
                }
            )

            table.insert(
                0,
                "Threshold",
                threshold,
            )

            if groupby is not None:

                table.insert(
                    1,
                    groupby,
                    group,
                )

            results.append(
                table
            )

    if not results:

        return pd.DataFrame()

    return (
        pd.concat(
            results,
            ignore_index=True,
        )
    )

###############################################################################
# Generic summary tables
###############################################################################

def summarize_cluster_table(
    table,
    sort_by=None,
):
    """
    Summarise transmission clusters by metadata variable.
    """

    if table.empty:
        return table

    # Metadata column
    metadata_column = next(
        column
        for column in table.columns
        if column not in (
            "Threshold",
            "Subtype",
            "Cluster",
            "Sequences",
        )
    )

    # Split transmission clusters / singletons
    transmission = table[
        table["Cluster"] != SINGLETON
    ].copy()

    singletons = table[
        table["Cluster"] == SINGLETON
    ].copy()

    # Group columns
    group_columns = ["Threshold"]

    if "Subtype" in table.columns:
        group_columns.append(
            "Subtype"
        )

    group_columns.append(
        metadata_column
    )

    # Transmission summary
    transmission_summary = (
        transmission
        .groupby(
            group_columns,
            dropna=False,
        )
        .agg(
            {
                "Cluster": "nunique",
                "Sequences": "sum",
            }
        )
        .reset_index()
        .rename(
            columns={
                "Cluster": "Transmission clusters",
                "Sequences": "Clustered sequences",
            }
        )
    )

    # Singleton summary
    singleton_summary = (
        singletons
        .groupby(
            group_columns,
            dropna=False,
        )
        .agg(
            {
                "Sequences": "sum",
            }
        )
        .reset_index()
        .rename(
            columns={
                "Sequences": "Singleton sequences",
            }
        )
    )

    # Merge summaries
    summary = (
        transmission_summary
        .merge(
            singleton_summary,
            on=group_columns,
            how="outer",
        )
        .fillna(0)
    )

    # Total sequences and percentages
    summary["Sequences"] = (
        summary["Clustered sequences"]
        +
        summary["Singleton sequences"]
    )

    summary["Clustered sequences (%)"] = (
        100
        * summary["Clustered sequences"]
        / summary["Sequences"]
    ).round(2)

    summary["Singleton sequences (%)"] = (
        100
        * summary["Singleton sequences"]
        / summary["Sequences"]
    ).round(2)

    # Convert to wide format
    index_columns = [metadata_column]

    if "Subtype" in summary.columns:
        index_columns.append(
            "Subtype"
        )

    metrics = [
        "Transmission clusters",
        "Clustered sequences",
        "Clustered sequences (%)",
        "Singleton sequences",
        "Singleton sequences (%)",
    ]

    wide_tables = []

    for metric in metrics:

        wide = (
            summary
            .pivot_table(
                index=index_columns,
                columns="Threshold",
                values=metric,
                fill_value=0,
            )
            .reset_index()
        )

        wide.columns.name = None

        wide = wide.rename(
            columns={
                threshold: f"{int(threshold)} SNP"
                for threshold in wide.columns
                if isinstance(
                    threshold,
                    (int, float),
                )
            }
        )

        wide.insert(
            len(index_columns),
            "Metric",
            metric,
        )

        wide_tables.append(
            wide
        )

    summary = pd.concat(
        wide_tables,
        ignore_index=True,
    )

    # Sort
    metric_order = {
        "Transmission clusters": 0,
        "Clustered sequences": 1,
        "Clustered sequences (%)": 2,
        "Singleton sequences": 3,
        "Singleton sequences (%)": 4,
    }

    summary["Metric_order"] = (
        summary["Metric"]
        .map(metric_order)
    )

    if sort_by == Columns.COLLECTION_YEAR:

        summary = summary.sort_values(
            by=[
                metadata_column,
                "Subtype",
                "Metric_order",
            ],
            ascending=[
                False,
                True,
                True,
            ],
        )

    elif sort_by == Columns.LINEAGE:

        snp_columns = sorted(
            [
                c
                for c in summary.columns
                if c.endswith("SNP")
            ],
            key=lambda x: int(
                x.split()[0]
            ),
        )

        last_threshold = snp_columns[-1]

        order = (
            summary[
                summary["Metric"] == "Clustered sequences"
            ][
                [
                    metadata_column,
                    "Subtype",
                    last_threshold,
                ]
            ]
            .rename(
                columns={
                    last_threshold: "_order",
                }
            )
        )

        summary = (
            summary
            .merge(
                order,
                on=[
                    metadata_column,
                    "Subtype",
                ],
                how="left",
            )
            .sort_values(
                by=[
                    "Subtype",
                    "_order",
                    "Metric_order",
                ],
                ascending=[
                    True,
                    False,
                    True,
                ],
            )
            .drop(
                columns="_order",
            )
        )

    else:

        summary = summary.sort_values(
            by=[
                metadata_column,
                "Subtype",
                "Metric_order",
            ],
            ascending=[
                True,
                True,
                True,
            ],
        )

    summary = (
        summary
        .drop(
            columns="Metric_order",
        )
        .reset_index(
            drop=True,
        )
    )

    return summary

###############################################################################
# Tables
###############################################################################

def cluster_by_lineage(
    dataset,
    groupby=None,
):
    """
    Generate a long-format transmission cluster table by lineage.
    """

    return _cluster_by_variable(
        dataset,
        Columns.LINEAGE,
        groupby,
    )

def cluster_by_year(
    dataset,
    groupby=None,
):
    """
    Generate a long-format transmission cluster table by lineage.
    """

    return _cluster_by_variable(
        dataset,
        Columns.COLLECTION_YEAR,
        groupby,
    )

def cluster_by_community(
    dataset,
    groupby=None,
):
    """
    Generate a long-format transmission cluster table by lineage.
    """

    return _cluster_by_variable(
        dataset,
        Columns.COMMUNITY,
        groupby,
    )

def prepare_cluster_plot_data(
    table,
    threshold,
    subtype=None,
):
    """
    Prepare transmission cluster summary for plotting.
    Returns absolute counts. Percentages are computed by the plotting
    functions when normalize=True.
    """

    if table.empty:
        return table

    # Threshold column
    if isinstance(threshold, int):
        threshold = f"{threshold} SNP"

    # Metadata column
    metadata_column = next(
        column
        for column in table.columns
        if column not in (
            "Subtype",
            "Metric",
        )
        and not column.endswith("SNP")
    )

    # Clustered sequences (absolute)
    clustered = (
        table[
            table["Metric"] == "Clustered sequences"
        ][
            [
                metadata_column,
                "Subtype",
                threshold,
            ]
        ]
        .rename(
            columns={
                threshold: "Clustered",
            }
        )
    )

    # Singleton sequences (absolute)
    singleton = (
        table[
            table["Metric"] == "Singleton sequences"
        ][
            [
                metadata_column,
                "Subtype",
                threshold,
            ]
        ]
        .rename(
            columns={
                threshold: "Singleton",
            }
        )
    )

    # Merge
    df = clustered.merge(
        singleton,
        on=[
            metadata_column,
            "Subtype",
        ],
    )

    # Total
    df["Total"] = (
        df["Clustered"]
        + df["Singleton"]
    )

    # Optional subtype filter
    if subtype is not None:

        df = df[
            df["Subtype"] == subtype
        ].copy()

    if metadata_column == Columns.COLLECTION_YEAR:
    	df[metadata_column] = df[
    		metadata_column
    	].astype(int)

    	df = df.sort_values(
    		metadata_column,
    	)

    return df

def cluster_by_lineage_year(
    dataset,
    groupby=None,
):
    """
    Summarise clustered and total sequences by lineage and
    collection year.
    """

    df = dataset.df.copy()

    if df.empty:
        return pd.DataFrame()

    # Groups
    if groupby is None:

        groups = [None]

    else:

        groups = (
            df[groupby]
            .dropna()
            .sort_values()
            .unique()
        )

    # Generate tables
    results = []

    for threshold in cluster_thresholds(dataset):

        cluster = cluster_column(threshold)

        for group in groups:

            if group is None:

                subset = df.copy()

            else:

                subset = df[
                    df[groupby] == group
                ].copy()

            # Total sequences
            total = (
                subset
                .groupby(
                    [
                        Columns.LINEAGE,
                        Columns.COLLECTION_YEAR,
                    ],
                    dropna=False,
                )
                .size()
                .reset_index(
                    name="Total sequences",
                )
            )

            # Clustered sequences
            clustered = (
                subset[
                    subset[cluster] != SINGLETON
                ]
                .groupby(
                    [
                        Columns.LINEAGE,
                        Columns.COLLECTION_YEAR,
                    ],
                    dropna=False,
                )
                .size()
                .reset_index(
                    name="Clustered sequences",
                )
            )

            # Merge
            table = (
                total
                .merge(
                    clustered,
                    on=[
                        Columns.LINEAGE,
                        Columns.COLLECTION_YEAR,
                    ],
                    how="left",
                )
                .fillna(0)
            )

            table["Clustered sequences"] = (
                table["Clustered sequences"]
                .astype(int)
            )

            table.insert(
                0,
                "Threshold",
                threshold,
            )

            if groupby is not None:

                table.insert(
                    1,
                    groupby,
                    group,
                )

            results.append(
                table
            )

    if not results:

        return pd.DataFrame()

    return pd.concat(
        results,
        ignore_index=True,
    )