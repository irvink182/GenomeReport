"""
Virus dynamics.

Utilities to aggregate viral metadata through time and space.
"""

import pandas as pd

from .config import Columns

###############################################################################
# Generic aggregation
###############################################################################

def aggregate_category(
    dataset,
    category,
    by,
    subtype=None,
):
    """
    Aggregate a categorical variable by one or more metadata fields.

    """

    df = dataset.df.copy()

    # Filter subtype
    if subtype is not None:

        df = df[
            df[Columns.SUBTYPE] == subtype
        ]

    # Extract collection year
    if by == "Year" or (
        isinstance(by, list)
        and "Year" in by
    ):

        df["Year"] = (
            pd.to_datetime(
                df[Columns.COLLECTION_DATE],
                errors="coerce",
            )
            .dt.year
        )

    # Accept one or multiple grouping variables
    if isinstance(by, str):
        by = [by]

    grouping = by + [category]

    # Count observations
    table = (
        df.groupby(grouping)
        .size()
        .reset_index(name="Count")
    )

    # Remove empty categories
    table = table[
        table[category]
        .astype(str)
        .str.strip()
        != ""
    ]

    table = table.sort_values(
        by=grouping
    ).reset_index(drop=True)

    return table

###############################################################################
# Matrix utilities
###############################################################################

def pivot_category(
    table,
    index,
    columns,
    values="Count",
    fill_value=0,
):
    """
    Convert a long-format table into a matrix.
    """

    matrix = (
        table.pivot_table(
            index=index,
            columns=columns,
            values=values,
            fill_value=fill_value,
            aggfunc="sum",
        )
        .reset_index()
    )

    matrix.columns.name = None

    return matrix


def complete_time_series(
    matrix,
    year_column="Year",
    fill_value=0,
):
    """
    Complete missing years in a temporal matrix.
    """

    years = pd.DataFrame({
        year_column: range(
            int(matrix[year_column].min()),
            int(matrix[year_column].max()) + 1,
        )
    })

    matrix = years.merge(
        matrix,
        on=year_column,
        how="left",
    )

    return matrix.fillna(fill_value)


def sort_category_columns(
    matrix,
    index,
    other_label="Other",
    subtype_order=False,
):
    """
    Sort category columns by total abundance.

    'Other' and 'unassigned' are always placed last.
    """

    df = matrix.copy()

    totals = (
        df
        .drop(columns=index)
        .sum(axis=0)
        .sort_values(
            ascending=False,
        )
    )

    ordered = totals.index.tolist()

    if subtype_order:
        a = [
            c
            for c in ordered
            if str(c).startswith("A")
        ]
        b = [
            c
            for c in ordered
            if str(c).startswith("B")
        ]
        other = [
            c
            for c in ordered
            if not str(c).startswith(("A", "B"))
        ]

        ordered = a + b + other

    if other_label in ordered:
        ordered.remove(other_label)
        ordered.append(other_label)

    unassigned = [
        c
        for c in ordered
        if str(c).lower() == "unassigned"
    ]

    for column in unassigned:
        ordered.remove(column)
        ordered.append(column)

    return df[
        [index] + ordered
    ]


def filter_category_matrix(
    matrix,
    index,
    top_n=None,
    min_count=None,
    min_percentage=None,
    keep_unassigned=True,
    other_label="Other",
    subtype_order=False,
):
    """
    Filter a category matrix prior to plotting.

    Categories can be filtered according to:

        • Top N
        • Minimum count
        • Minimum percentage

    Remaining categories are merged into "Other".
    """

    df = matrix.copy()

    values = df.drop(columns=index)

    # Total abundance
    totals = values.sum(axis=0)

    # Categories to retain
    keep = set()

    if top_n is not None:

        keep.update(
            totals
            .sort_values(ascending=False)
            .head(top_n)
            .index
        )

    if min_count is not None:

        keep.update(
            totals[
                totals >= min_count
            ].index
        )

    if min_percentage is not None:

        percentages = (
            totals / totals.sum() * 100
        )

        keep.update(
            percentages[
                percentages >= min_percentage
            ].index
        )

    if (
        top_n is None
        and min_count is None
        and min_percentage is None
    ):
        keep = set(totals.index)

    # Preserve unassigned
    if keep_unassigned:

        keep.update([
            c
            for c in values.columns
            if str(c).lower() == "unassigned"
        ])

    # Build filtered matrix
    result = pd.DataFrame()

    result[index] = df[index]

    other = []

    for column in values.columns:

        if column in keep:

            result[column] = values[column]

        else:

            other.append(column)

    # Merge remaining categories
    if other:

        result[other_label] = (
            values[other]
            .sum(axis=1)
        )

    # Sort columns
    result = sort_category_columns(
        result,
        index=index,
        other_label=other_label,
        subtype_order=subtype_order,
    )

    return result

###############################################################################
# Generic wrappers
###############################################################################

def category_summary(
    dataset,
    category,
    by,
    subtype=None,
):
    """
    Aggregate a category by one or more metadata variables.
    """

    return aggregate_category(
        dataset,
        category=category,
        by=by,
        subtype=subtype,
    )

###############################################################################
# Virus-specific wrappers
###############################################################################

def lineage_by_year(
    dataset,
    subtype=None,
):
    """
    Lineage distribution by year.
    """

    return category_summary(
        dataset,
        category=Columns.LINEAGE,
        by="Year",
        subtype=subtype,
    )


def genotype_by_year(
    dataset,
    subtype=None,
):
    """
    Genotype distribution by year.
    """

    return category_summary(
        dataset,
        category=Columns.GENOTYPE,
        by="Year",
        subtype=subtype,
    )


def lineage_by_community(
    dataset,
    subtype=None,
):
    """
    Lineage distribution by Autonomous Community.
    """

    return category_summary(
        dataset,
        category=Columns.LINEAGE,
        by=Columns.COMMUNITY,
        subtype=subtype,
    )


def genotype_by_community(
    dataset,
    subtype=None,
):
    """
    Genotype distribution by Autonomous Community.
    """

    return category_summary(
        dataset,
        category=Columns.GENOTYPE,
        by=Columns.COMMUNITY,
        subtype=subtype,
    )


def lineage_by_year_community(
    dataset,
    subtype=None,
):
    """
    Lineage distribution by year and Autonomous Community.
    """

    return category_summary(
        dataset,
        category=Columns.LINEAGE,
        by=[
            "Year",
            Columns.COMMUNITY,
        ],
        subtype=subtype,
    )


def genotype_by_year_community(
    dataset,
    subtype=None,
):
    """
    Genotype distribution by year and Autonomous Community.
    """

    return category_summary(
        dataset,
        category=Columns.GENOTYPE,
        by=[
            "Year",
            Columns.COMMUNITY,
        ],
        subtype=subtype,
    )