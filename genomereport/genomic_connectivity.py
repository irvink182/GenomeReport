import itertools

import pandas as pd

from .config import (
    Columns,
    TRANSMISSION_CLUSTER_THRESHOLDS,
)


SINGLETON_LABEL = "Singleton"


def _prepare_connectivity_data(
    dataset,
    threshold,
    group_variable,
):
    """
    Prepare clustered sequences for genomic connectivity analysis.
    """

    cluster_column = f"Cluster_{threshold}"

    required_columns = [
        Columns.SUBTYPE,
        Columns.ACCESSION,
        Columns.COLLECTION_DATE,
        group_variable,
        cluster_column,
    ]

    missing = [
        column
        for column in required_columns
        if column not in dataset.df.columns
    ]

    if missing:
        raise ValueError(
            "Missing columns required for genomic connectivity analysis: "
            + ", ".join(missing)
        )

    df = dataset.df.copy()

    # Keep clustered sequences only
    df = df[
        df[cluster_column].notna()
        &
        (df[cluster_column] != SINGLETON_LABEL)
    ].copy()

    # Remove sequences without grouping information
    df = df[
        df[group_variable].notna()
        &
        (df[group_variable] != "")
        &
        (df[group_variable] != "Unknown")
    ].copy()

    df[Columns.COLLECTION_DATE] = pd.to_datetime(
        df[Columns.COLLECTION_DATE],
        errors="coerce",
    )

    df["Threshold"] = threshold
    df["Cluster"] = df[cluster_column]

    return df

def cluster_connectivity_summary(
    dataset,
    group_variable,
    thresholds=TRANSMISSION_CLUSTER_THRESHOLDS,
):
    """
    Summarize genomic connectivity for each genomic cluster.

    A cluster is classified as:

        SINGLE_GROUP
            All sequences belong to the same group.

        MULTI_GROUP
            Sequences belong to two or more groups.

    Parameters
    ----------
    dataset
        GenomeReport Dataset object.

    group_variable : str
        Metadata column used to define groups
        (e.g. autonomous_community, institution).

    thresholds : iterable
        SNP clustering thresholds to analyse.

    """

    rows = []

    for threshold in thresholds:

        df = _prepare_connectivity_data(
            dataset=dataset,
            threshold=threshold,
            group_variable=group_variable,
        )

        if df.empty:
            continue

        for (
            subtype,
            cluster,
        ), cluster_df in df.groupby(
            [
                Columns.SUBTYPE,
                "Cluster",
            ],
            observed=False,
        ):

            groups = sorted(
                cluster_df[group_variable]
                .dropna()
                .astype(str)
                .unique()
            )

            lineages = []

            if Columns.LINEAGE in cluster_df.columns:
                lineages = sorted(
                    cluster_df[Columns.LINEAGE]
                    .dropna()
                    .astype(str)
                    .unique()
                )

            seasons = []

            if "epidemiological_season" in cluster_df.columns:
                seasons = sorted(
                    cluster_df["epidemiological_season"]
                    .replace("", pd.NA)
                    .dropna()
                    .astype(str)
                    .unique()
                )

            dates = (
                cluster_df[Columns.COLLECTION_DATE]
                .dropna()
            )

            if dates.empty:
                first_date = pd.NaT
                last_date = pd.NaT
                duration_days = pd.NA
            else:
                first_date = dates.min()
                last_date = dates.max()

                duration_days = (
                    last_date - first_date
                ).days

            n_groups = len(groups)

            connectivity = (
                "MULTI_GROUP"
                if n_groups > 1
                else "SINGLE_GROUP"
            )

            rows.append({
                Columns.SUBTYPE: subtype,
                "Threshold": threshold,
                "Cluster": cluster,
                "n_sequences": len(cluster_df),
                "n_groups": n_groups,
                "Groups": ";".join(groups),
                "n_lineages": len(lineages),
                "Lineages": ";".join(lineages),
                "n_seasons": len(seasons),
                "Seasons": ";".join(seasons),
                "First_date": first_date,
                "Last_date": last_date,
                "Duration_days": duration_days,
                "Connectivity": connectivity,
            })

    columns = [
        Columns.SUBTYPE,
        "Threshold",
        "Cluster",
        "n_sequences",
        "n_groups",
        "Groups",
        "n_lineages",
        "Lineages",
        "n_seasons",
        "Seasons",
        "First_date",
        "Last_date",
        "Duration_days",
        "Connectivity",
    ]

    if not rows:
        return pd.DataFrame(
            columns=columns,
        )

    table = pd.DataFrame(rows)

    return (
        table
        .sort_values(
            [
                Columns.SUBTYPE,
                "Threshold",
                "n_groups",
                "n_sequences",
                "Cluster",
            ],
            ascending=[
                True,
                True,
                False,
                False,
                True,
            ],
        )
        .reset_index(drop=True)
    )

def select_connectivity_group(dataset):
    """
    Select the grouping variable for genomic connectivity analysis.

    Autonomous community is preferred when multiple communities
    are represented. Institution is used when the dataset contains
    a single community but multiple institutions.
    """

    df = dataset.df

    communities = (
        df[Columns.COMMUNITY]
        .replace("", pd.NA)
        .replace("Unknown", pd.NA)
        .dropna()
        .nunique()
    )

    institutions = (
        df[Columns.INSTITUTION]
        .replace("", pd.NA)
        .replace("Unknown", pd.NA)
        .dropna()
        .nunique()
    )

    if communities > 1:
        return Columns.COMMUNITY

    if institutions > 1:
        return Columns.INSTITUTION

    return None

def cluster_connectivity_samples(
    dataset,
    group_variable,
    thresholds=TRANSMISSION_CLUSTER_THRESHOLDS,
):
    """
    Return sequence-level information for genomic clusters.

    Each row represents one sequence belonging to a non-singleton
    genomic cluster. Connectivity is defined according to the number
    of distinct groups represented in the cluster:

        SINGLE_GROUP
            All sequences in the cluster belong to the same group.

        MULTI_GROUP
            Sequences in the cluster belong to two or more groups.

    Parameters
    ----------
    dataset
        GenomeReport Dataset object.

    group_variable : str
        Metadata column used to define groups
        (e.g. autonomous_community, institution).

    thresholds : iterable
        SNP clustering thresholds to analyse.

    """

    tables = []

    for threshold in thresholds:

        df = _prepare_connectivity_data(
            dataset=dataset,
            threshold=threshold,
            group_variable=group_variable,
        )

        if df.empty:
            continue

        # Number of groups represented in each cluster.
        cluster_n_groups = (
            df.groupby(
                [
                    Columns.SUBTYPE,
                    "Cluster",
                ],
                observed=False,
            )[group_variable]
            .nunique()
        )

        df["Connectivity"] = [
            (
                "MULTI_GROUP"
                if cluster_n_groups.loc[(subtype, cluster)] > 1
                else "SINGLE_GROUP"
            )
            for subtype, cluster in zip(
                df[Columns.SUBTYPE],
                df["Cluster"],
            )
        ]

        # Columns to export.
        output_columns = [
            Columns.SUBTYPE,
            "Threshold",
            "Cluster",
            Columns.ACCESSION,
            Columns.VIRUS_NAME,
            group_variable,
            Columns.LINEAGE,
            Columns.COLLECTION_DATE,
            "epidemiological_season",
            "gisaid_epiweek",
            "season_week",
            "season_epiweek_custom",
            "Connectivity",
        ]

        # Keep only columns available in the dataset.
        output_columns = [
            column
            for column in output_columns
            if column in df.columns
        ]

        table = df[
            output_columns
        ].copy()

        # Rename the selected grouping variable to a generic name.
        table = table.rename(
            columns={
                group_variable: "Group",
            }
        )

        tables.append(table)

    columns = [
        Columns.SUBTYPE,
        "Threshold",
        "Cluster",
        Columns.ACCESSION,
        Columns.VIRUS_NAME,
        "Group",
        Columns.LINEAGE,
        Columns.COLLECTION_DATE,
        "epidemiological_season",
        "gisaid_epiweek",
        "season_week",
        "season_epiweek_custom",
        "Connectivity",
    ]

    if not tables:
        return pd.DataFrame(
            columns=columns,
        )

    table = pd.concat(
        tables,
        ignore_index=True,
    )

    return (
        table
        .sort_values(
            [
                Columns.SUBTYPE,
                "Threshold",
                "Cluster",
                Columns.COLLECTION_DATE,
                "Group",
                Columns.ACCESSION,
            ]
        )
        .reset_index(drop=True)
    )


def group_connectivity(
    dataset,
    group_variable,
    thresholds=TRANSMISSION_CLUSTER_THRESHOLDS,
):
    """
    Summarize genomic connectivity between pairs of groups.

    For each subtype and SNP threshold, every cluster containing
    sequences from two or more groups contributes one shared cluster
    to each pairwise combination of those groups.

    Parameters
    ----------
    dataset
        GenomeReport Dataset object.

    group_variable : str
        Metadata column used to define groups
        (e.g. autonomous_community, institution).

    thresholds : iterable
        SNP clustering thresholds to analyse.
    """

    rows = []

    for threshold in thresholds:

        df = _prepare_connectivity_data(
            dataset=dataset,
            threshold=threshold,
            group_variable=group_variable,
        )

        if df.empty:
            continue

        # Collect shared clusters for each pair of groups.
        pair_clusters = {}

        for (
            subtype,
            cluster,
        ), cluster_df in df.groupby(
            [
                Columns.SUBTYPE,
                "Cluster",
            ],
            observed=False,
        ):

            groups = sorted(
                cluster_df[group_variable]
                .dropna()
                .astype(str)
                .unique()
            )

            # Single-group clusters do not contribute
            # to between-group connectivity.
            if len(groups) < 2:
                continue

            for group_1, group_2 in itertools.combinations(
                groups,
                2,
            ):
                key = (
                    subtype,
                    group_1,
                    group_2,
                )

                pair_clusters.setdefault(
                    key,
                    set(),
                ).add(cluster)

        for (
            subtype,
            group_1,
            group_2,
        ), clusters in pair_clusters.items():

            clusters = sorted(clusters)

            rows.append({
                Columns.SUBTYPE: subtype,
                "Threshold": threshold,
                "Group_1": group_1,
                "Group_2": group_2,
                "n_shared_clusters": len(clusters),
                "Shared_clusters": ";".join(clusters),
            })

    columns = [
        Columns.SUBTYPE,
        "Threshold",
        "Group_1",
        "Group_2",
        "n_shared_clusters",
        "Shared_clusters",
    ]

    if not rows:
        return pd.DataFrame(
            columns=columns,
        )

    table = pd.DataFrame(rows)

    return (
        table
        .sort_values(
            [
                Columns.SUBTYPE,
                "Threshold",
                "n_shared_clusters",
                "Group_1",
                "Group_2",
            ],
            ascending=[
                True,
                True,
                False,
                True,
                True,
            ],
        )
        .reset_index(drop=True)
    )