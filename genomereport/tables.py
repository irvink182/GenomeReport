"""
Summary tables.
"""

from pathlib import Path

import pandas as pd
import numpy as np

from .config import (
    Columns, 
    Metadata,
    SUBTYPES,
    MutationColumns,
    PROTEIN_ORDER,
    PROTEIN_LENGTHS,
    LOLLIPOP_PROTEINS,
)

from .style import (
    OTHER_LABEL,
    LOLLIPOP_LABEL_THRESHOLD,
    CONSOLE_MUTATION_FREQUENCY,
)

from .normalize import (
    is_missing,
    normalize_sampling_strategy,
)

from .virus_dynamics import(
    lineage_by_community,)

from .mutations import (
    explode_mutations,)

from .protein_features import (
    PROTEIN_FEATURES,
)

from .clusters import(
    cluster_by_lineage_year,
)

#######################
# Generic
#######################

def count_unique(df, column):
    """
    Count unique non-empty values.
    """

    if column not in df.columns:
        return None

    return (
        df[column]
        .replace("", pd.NA)
        .dropna()
        .nunique()
    )

def count_duplicates(df, column):
    """
    Count duplicated values.
    """

    if column not in df.columns:
        return None

    return len(df) - df[column].nunique()

def count_duplicate_patients(df):
    """
    Count duplicated patient records based on
    the 'Additional host information' field.
    """

    if Columns.ADDITIONAL_HOST_INFORMATION not in df.columns:
        return None

    return (
        df[Columns.ADDITIONAL_HOST_INFORMATION]
        .fillna("")
        .str.contains(
            "refers to same patient",
            case=False,
            regex=False,
        )
        .sum()
    )

def sequence_length_summary(df):
    """
    Genome length summary statistics.
    """

    if Columns.SEQUENCE_LENGTH not in df.columns:
        return None

    length = (
        pd.to_numeric(
            df[Columns.SEQUENCE_LENGTH],
            errors="coerce",
        )
        .dropna()
    )

    if length.empty:
        return None

    return (
        f"{int(length.min()):,} - "
        f"{int(length.max()):,} bp "
        f"(median: {int(length.median()):,})"
    )

def collection_period(df):
    """
    Return collection date range using only
    complete collection dates (YYYY-MM-DD).
    """

    if Columns.COLLECTION_DATE not in df.columns:
        return None

    collection_date = (
        df[Columns.COLLECTION_DATE]
        .astype("string")
        .str.strip()
    )

    # Keep only complete dates (YYYY-MM-DD)
    complete_dates = collection_date.str.fullmatch(
        r"\d{4}-\d{2}-\d{2}",
        na=False,
    )

    dates = pd.to_datetime(
        collection_date[complete_dates],
        format="%Y-%m-%d",
        errors="coerce",
    ).dropna()

    if dates.empty:
        return None

    return (
        f"{dates.min().date()} - {dates.max().date()}"
    )

def complete_years(table, column="Year"):
    """
    Fill missing years with zero counts.
    """

    years = pd.DataFrame({
        column: range(
            table[column].min(),
            table[column].max() + 1,
        )
    })

    table = (
        years
        .merge(table, on=column, how="left")
        .fillna(0)
    )

    for col in table.columns:
        if col != column:
            table[col] = table[col].astype(int)

    return table

def count_by(
    df,
    column,
    subtype_column=Columns.SUBTYPE,
    total_column="Total",
    sort_by=None,
):
    """
    Count sequences by a categorical variable and subtype.
    """

    subtypes = (
        df[subtype_column]
        .dropna()
        .unique()
        .tolist()
    )

    df = df[
        ~df[column].apply(is_missing)
    ].copy()

    table = (
        df.groupby(
            [column, subtype_column],
            observed=False,
        )
        .size()
        .unstack(fill_value=0)
        .reindex(
            columns=subtypes,
            fill_value=0,
        )
        .rename_axis(None, axis=1)
        .reset_index()
    )

    table[total_column] = (
        table[subtypes]
        .sum(axis=1)
    )

    table["Percentage"] = (
        table[total_column]
        / table[total_column].sum()
        * 100
    ).round(2)

    if sort_by is None:
        sort_by = column

    return (
        table
        .sort_values(sort_by)
        .reset_index(drop=True)
    )

def count_by_bins(
    df,
    column,
    bins,
    labels,
    category_name=None,
    unknown_label="Unknown",
):
    """
    Count sequences after binning a continuous variable.
    """

    df = df.copy()

    values = pd.to_numeric(
        df[column],
        errors="coerce",
    )

    if category_name is None:
        category_name = f"{column}_group"

    df[category_name] = pd.cut(
        values,
        bins=bins,
        labels=labels,
        right=False,
    )

    df[category_name] = (
        df[category_name]
        .cat.add_categories([unknown_label])
        .fillna(unknown_label)
    )

    return count_by(
        df,
        category_name,
    )

def virus_category(
    dataset,
    category,
    subtype,
):
    """
    Count viral categories (lineages, genotypes, etc.)
    for a given RSV subtype.
    """

    return count_by_subtype(
        dataset.df,
        category=category,
        subtype=subtype,
    )

def count_by_subtype(
    df,
    category,
    subtype,
):
    """
    Count values of a variable within one RSV subtype.
    """

    df = df[
        df[Columns.SUBTYPE] == subtype
    ].copy()

    table = (
        df.groupby(category)
        .size()
        .reset_index(name="Total")
        .sort_values(
            "Total",
            ascending=False,
        )
    )

    table["Percentage"] = (
        table["Total"]
        / table["Total"].sum()
        * 100
    ).round(2)

    return table

#######################
# Dataset
#######################

def dataset_summary(dataset):
    """
    Build dataset overview table.
    """

    df = dataset.df

    summary = {

        "Dataset": dataset.name,

        "Total sequences": len(df),

        "Unique accession IDs":
            count_unique(df, Columns.ACCESSION),

        "Duplicated accession IDs":
            count_duplicates(df, Columns.ACCESSION),

        "Duplicate patients":
            count_duplicate_patients(df),

        "RSV-A":
            (df[Columns.SUBTYPE] == "A").sum(),

        "RSV-B":
            (df[Columns.SUBTYPE] == "B").sum(),

        "Collection period":
            collection_period(df),

        "Sequence length":
            sequence_length_summary(df),

        "Institutions":
            count_unique(df, Columns.INSTITUTION),

        "Originating laboratories":
            count_unique(df, Columns.ORIGINATING_LAB),

        "Submitter laboratories":
            count_unique(df, Columns.SUBMITTER_LAB),

        "Communities":
            count_unique(df, Columns.COMMUNITY),

        "Cities":
            count_unique(df, Columns.CITY),

        "Lineages":
            count_unique(df, Columns.LINEAGE),

        "Genotypes":
            count_unique(df, Columns.GENOTYPE),

        "Sequencing technologies":
            count_unique(df, Columns.TECHNOLOGY),

        "Assembly methods":
            count_unique(df, Columns.ASSEMBLY),

    }

    return pd.DataFrame(
        summary.items(),
        columns=[
            "Metric",
            "Value",
        ],
    )

def metadata_completeness(dataset, columns=None):
    """
    Summarize metadata completeness for selected columns.
    """

    df = dataset.df

    if columns is None:
        columns = Metadata.CORE

    rows = []

    for column in columns:

        missing = df[column].apply(is_missing).sum()

        complete = len(df) - missing

        rows.append({
            "Variable": column,
            "Complete": complete,
            "Missing": missing,
            "Missing (%)": round(
                missing / len(df) * 100,
                2,
            ),
        })

    return (
        pd.DataFrame(rows)
        .sort_values("Missing (%)")
        .reset_index(drop=True)
    )

#######################
# Geography
#######################

def geography_community(dataset):
    """
    Number of sequences by autonomous community.
    """

    return count_by(
        dataset.df,
        Columns.COMMUNITY,
    )

    table["Percentage"] = (
        table["Total"] /
        table["Total"].sum() * 100
    ).round(2)

    return table.sort_values(
        "Total",
        ascending=False,
    )

def sampling_strategy_community(
    dataset,
):
    """
    Number of sentinel and non-sentinel sequences
    by autonomous community.
    """

    table = dataset.df.copy()

    table["Sampling group"] = (
        table[
            Columns.SAMPLING_STRATEGY
        ].apply(
            normalize_sampling_strategy
        )
    )

    table = (
        table
        .groupby(
            [
                Columns.COMMUNITY,
                "Sampling group",
            ]
        )
        .size()
        .unstack(
            fill_value=0,
        )
    )

    for column in [
        "Sentinel",
        "Non-sentinel",
    ]:
        if column not in table.columns:
            table[column] = 0

    table["Total"] = (
        table[
            [
                "Sentinel",
                "Non-sentinel",
            ]
        ].sum(
            axis=1,
        )
    )

    # Total
    table["Total"] = (
        table["Sentinel"] +
        table["Non-sentinel"]
    )

    # Percentages
    table["Sentinel (%)"] = (
        table["Sentinel"]
        / table["Total"]
        * 100
    ).round(2)

    table["Non-sentinel (%)"] = (
        table["Non-sentinel"]
        / table["Total"]
        * 100
    ).round(2)

    table = table.sort_values(
        "Total",
        ascending=False,
    )

    #print (table.index.name)
    table = table.reset_index()
    #print(table.columns.tolist())

    return table[
        [
            Columns.COMMUNITY,
            "Sentinel",
            "Non-sentinel",
            "Total",
            "Sentinel (%)",
            "Non-sentinel (%)",
        ]
    ]

def summarize_sampling_strategy(table):

    return (
        table[
            [
                "Sentinel",
                "Non-sentinel",
            ]
        ]
        .sum()
        .rename_axis(
            "Sampling strategy",
        )
        .reset_index(
            name="Total",
        )
    )

#######################
# Laboratory
#######################

def laboratory_institution(dataset):

    return count_by(
        dataset.df,
        Columns.INSTITUTION,
    )

    table["Percentage"] = (
        table["Total"] /
        table["Total"].sum() * 100
    ).round(2)

    return table.sort_values(
        "Total",
        ascending=False,
    )

def laboratory_originating(dataset):

    return count_by(
        dataset.df,
        Columns.ORIGINATING_LAB,
    )
    
    table["Percentage"] = (
        table["Total"] /
        table["Total"].sum() * 100
    ).round(2)

    return table.sort_values(
        "Total",
        ascending=False,
    )

def laboratory_submitter(dataset):

    return count_by(
        dataset.df,
        Columns.SUBMITTER_LAB,
    )

    table["Percentage"] = (
        table["Total"] /
        table["Total"].sum() * 100
    ).round(2)

    return table.sort_values(
        "Total",
        ascending=False,
    )

def temporal_year(dataset):
    """
    Number of sequences by collection year.
    """

    df = dataset.df.copy()

    # Use normalized collection year
    df["Year"] = pd.to_numeric(
        df[Columns.COLLECTION_YEAR],
        errors="coerce",
    )

    # Remove records without year information
    df = df.dropna(
        subset=["Year"]
    )

    # Convert to integer
    df["Year"] = (
        df["Year"]
        .astype(int)
    )

    # Count sequences
    table = count_by(
        df,
        "Year",
    )

    return complete_years(table)

#######################
# Sequencing
#######################

def sequencing_technology(dataset):
    """
    Number of sequences by sequencing technology.
    """

    return count_by(
        dataset.df,
        Columns.TECHNOLOGY_NORMALIZED,
    )

def sequencing_assembly(dataset):
    """
    Number of sequences by assembly method.
    """

    return count_by(
        dataset.df,
        Columns.ASSEMBLY_NORMALIZED,
    )

#######################
# Virus lineage & genotype
#######################

def virus_lineage(
    dataset,
    subtype,
):
    """
    Number of sequences by lineage.
    """

    return count_by_subtype(
        dataset.df,
        category=Columns.LINEAGE,
        subtype=subtype,
    )

def virus_genotype(
    dataset,
    subtype,
):
    """
    Number of sequences by genotype.
    """

    return count_by_subtype(
        dataset.df,
        category=Columns.GENOTYPE,
        subtype=subtype,
    )

def major_lineages(
    lineages,
    coverage=80,
):
    """
    Select the major circulating lineages covering a given
    cumulative percentage of genomes.
    """

    rows = []

    for subtype, table in lineages.items():

        df = table.copy()

        df["Cumulative (%)"] = (
            df["Percentage"]
            .cumsum()
            .round(2)
        )

        df["Selected"] = False

        # Lineages below the threshold
        df.loc[
            df["Cumulative (%)"] <= coverage,
            "Selected",
        ] = True

        # Include the first lineage exceeding the threshold
        above = df[
            df["Cumulative (%)"] > coverage
        ]

        if not above.empty:

            df.loc[
                above.index[0],
                "Selected",
            ] = True

        df.insert(
            0,
            Columns.SUBTYPE,
            subtype,
        )

        rows.append(df)

    return pd.concat(
        rows,
        ignore_index=True,
    )

def get_selected_lineage(
    major_lineages,
    subtype,
):
    """
    Return the most abundant selected lineage for a subtype.
    """

    return (
        major_lineages
        .loc[
            (major_lineages["Subtype"] == subtype)
            &
            (major_lineages["Selected"]),
        ]
        .sort_values(
            "Total",
            ascending=False,
        )
        .iloc[0]["Lineage"]
    )

def community_lineage_heatmap_table(
    lineage_table,
    major_lineages,
    subtype,
):
    """
    Prepare the count matrix for the community lineage heatmap.
    """

    # Copy table
    df = lineage_table.copy()

    df = df[
        df[Columns.COMMUNITY] != "Unknown"
    ]

    # Major lineages
    selected = (
        major_lineages.loc[
            (major_lineages[Columns.SUBTYPE] == subtype)
            &
            (major_lineages["Selected"])
        ][Columns.LINEAGE]
        .tolist()
    )

    # Collapse minor lineages
    df[Columns.LINEAGE] = (
        df[Columns.LINEAGE]
        .where(
            df[Columns.LINEAGE].isin(selected),
            OTHER_LABEL,
        )
    )

    # Merge "Other"
    df = (
        df
        .groupby(
            [
                Columns.COMMUNITY,
                Columns.LINEAGE,
            ],
            as_index=False,
        )["Count"]
        .sum()
    )

    # Count matrix
    counts = (
        df
        .pivot(
            index=Columns.COMMUNITY,
            columns=Columns.LINEAGE,
            values="Count",
        )
        .fillna(0)
        .astype(int)
    )

    # Column order
    column_order = selected + [OTHER_LABEL]

    counts = counts.reindex(
        columns=column_order,
        fill_value=0,
    )

    # Sort communities by total number of sequences
    counts = counts.loc[
        counts.sum(axis=1)
        .sort_values(ascending=False)
        .index
    ]

    counts = counts.reset_index()
    
    return counts

def lineage_cluster_heatmap_table(
    lineage_clusters,
    major_lineages,
    subtype,
    threshold,
):
    """
    Prepare the lineage transmission heatmap table.
    """
    # Filter
    df = lineage_clusters.copy()

    df = df[
        (df["Threshold"] == threshold)
        &
        (df[Columns.SUBTYPE] == subtype)
    ]

    # Major lineages
    selected = (
        major_lineages.loc[
            (major_lineages[Columns.SUBTYPE] == subtype)
            &
            (major_lineages["Selected"])
        ][Columns.LINEAGE]
        .tolist()
    )

    df = df[
        df[Columns.LINEAGE].isin(selected)
    ]

    # Clustering percentage within lineage-year
    year_totals = (
        df.groupby(Columns.COLLECTION_YEAR)["Clustered sequences"]
            .sum()
            .rename("Year total")
    )

    df = df.merge(
        year_totals,
        on=Columns.COLLECTION_YEAR,
    )

    df["Clustered (%)"] = (
        df["Clustered sequences"]
        / df["Year total"]
        * 100
    ).fillna(0)

    # Pivot tables
    counts = (
        df
        .pivot_table(
            index=Columns.LINEAGE,
            columns=Columns.COLLECTION_YEAR,
            values="Clustered sequences",
            aggfunc="sum",
            fill_value=0,
        )
        .astype(int)
    )

    percentages = (
        df
        .pivot_table(
            index=Columns.LINEAGE,
            columns=Columns.COLLECTION_YEAR,
            values="Clustered (%)",
            aggfunc="mean",
            fill_value=0,
        )
    )

    # Total sequences per lineage (for row labels)
    lineage_totals = (
        df
        .groupby(Columns.LINEAGE)["Total sequences"]
        .sum()
        .astype(int)
    )

    # Order years
    years = sorted(
        counts.columns,
        key=int,
    )

    counts = counts.reindex(
        columns=years,
        fill_value=0,
    )

    percentages = percentages.reindex(
        columns=years,
        fill_value=0,
    )

    # Order lineages by total number of sequences
    order = (
        lineage_totals
        .sort_values(
            ascending=False,
        )
        .index
    )

    counts = counts.loc[order]
    percentages = percentages.loc[order]
    lineage_totals = lineage_totals.loc[order]

    return (
        counts.reset_index(),
        percentages.reset_index(),
        lineage_totals,
    )

#######################
# Patient
#######################

def patient_host(dataset):
    """
    Number of sequences by host.
    """

    return count_by(
        dataset.df,
        Columns.HOST,
    )

def patient_gender(dataset):
    """
    Number of sequences by gender.
    """

    return count_by(
        dataset.df,
        Columns.GENDER_NORMALIZED,
    )

def patient_status(dataset):
    """
    Number of sequences by patient status.
    """

    return count_by(
        dataset.df,
        Columns.STATUS_NORMALIZED,
    )

def patient_age_summary(dataset):
    """
    Summary statistics of patient age.
    """

    df = dataset.df.copy()

    age = (
        pd.to_numeric(
            df[Columns.AGE_NORMALIZED],
            errors="coerce",
        )
        .dropna()
    )

    rows = [

        {
            "Metric": "Patients with age",
            "Value": len(age),
        },

        {
            "Metric": "Missing",
            "Value": dataset.n_rows - len(age),
        },

        {
            "Metric": "Mean age",
            "Value": round(age.mean(), 1),
        },

        {
            "Metric": "Median age",
            "Value": round(age.median(), 1),
        },

        {
            "Metric": "Median (IQR)",
            "Value": (
                f"{age.median():.1f} "
                f"({age.quantile(0.25):.1f}–"
                f"{age.quantile(0.75):.1f})"
            ),
        },

        {
            "Metric": "Minimum",
            "Value": age.min(),
        },

        {
            "Metric": "Maximum",
            "Value": age.max(),
        },

    ]

    return pd.DataFrame(rows)

def patient_age_groups(dataset):
    """
    Number of sequences by age group.
    """

    bins = [
        -1,
        1,
        5,
        18,
        50,
        65,
        150,
    ]

    labels = [
        "<1 year",
        "1–4 years",
        "5–17 years",
        "18–49 years",
        "50–64 years",
        "≥65 years",
    ]

    return count_by_bins(
        dataset.df,
        Columns.AGE_NORMALIZED,
        bins,
        labels,
        category_name="Age group",
    )

#######################
# Mutations
#######################

def protein_lineage_summary(
    mutations,
    major_lineages,
):
    """
    Summarise amino acid substitutions by lineage and protein.

    Only the selected major lineages are retained.
    """

    # Copy table
    df = mutations.copy()

    # Number of genomes per lineage
    genomes = (
        df
        .groupby(
            [
                Columns.SUBTYPE,
                Columns.LINEAGE,
            ],
            as_index=False,
        )
        .agg(
            **{
                MutationColumns.GENOMES: (
                    Columns.ACCESSION,
                    "nunique",
                )
            }
        )
    )

    # Keep genomes with mutations
    df = df[
        df[MutationColumns.HAS_MUTATIONS]
    ]

    # Keep only major lineages
    selected = (
        major_lineages
        .loc[
            major_lineages["Selected"],
            [
                Columns.SUBTYPE,
                Columns.LINEAGE,
            ],
        ]
        .drop_duplicates()
    )

    df = df.merge(
        selected,
        on=[
            Columns.SUBTYPE,
            Columns.LINEAGE,
        ],
        how="inner",
    )

    # Protein summary
    summary = (
        df
        .groupby(
            [
                Columns.SUBTYPE,
                Columns.LINEAGE,
                MutationColumns.PROTEIN,
            ],
            as_index=False,
        )

        .agg(
            **{
                MutationColumns.TOTAL_EVENTS: (
                    MutationColumns.AA_SUBSTITUTION,
                    "size",
                ),
                MutationColumns.DISTINCT_SUBSTITUTIONS: (
                    MutationColumns.AA_SUBSTITUTION,
                    "nunique",
                ),
            }
        )
    )

    # Merge genome counts
    summary = summary.merge(
        genomes,
        on=[
            Columns.SUBTYPE,
            Columns.LINEAGE,
        ],
        how="left",
    )

    # Protein length
    summary[MutationColumns.PROTEIN_LENGTH] = summary.apply(
        lambda row: PROTEIN_LENGTHS[
            row[Columns.SUBTYPE]
        ][
            row[MutationColumns.PROTEIN]
        ],
        axis=1,
    )

    # Normalised metrics by genome
    summary[MutationColumns.EVENTS_PER_GENOME] = (
        summary[MutationColumns.TOTAL_EVENTS]
        / summary[MutationColumns.GENOMES]
    ).round(2)

    # Normalised metrics by protein length
    summary[MutationColumns.EVENTS_PER_100_AA] = (
        summary[MutationColumns.TOTAL_EVENTS]
        / summary[MutationColumns.PROTEIN_LENGTH]
        * 100
    ).round(2)

    summary[MutationColumns.DISTINCT_PER_100_AA] = (
        summary[MutationColumns.DISTINCT_SUBSTITUTIONS]
        * 100
        / summary[MutationColumns.PROTEIN_LENGTH]
    ).round(2)

    # Protein order
    summary[MutationColumns.PROTEIN] = pd.Categorical(
        summary[MutationColumns.PROTEIN],
        categories=PROTEIN_ORDER,
        ordered=True,
    )

    summary = (
        summary
        .sort_values(
            [
                Columns.SUBTYPE,
                Columns.LINEAGE,
                MutationColumns.PROTEIN,
            ]
        )

        .reset_index(
            drop=True,
        )
    )
    return summary

def protein_lineage_comparative_heatmap_table(
    percentages,
):
    """
    Calculate robust z-scores of the relative distribution of
    distinct amino acid substitutions across lineages.

    Robust z-scores are calculated independently for each protein
    using the median and median absolute deviation (MAD).
    """

    # Protein column
    protein_column = percentages.columns[0]

    # Keep protein names
    proteins = percentages[protein_column]

    # Numerical matrix
    values = percentages.iloc[:, 1:].copy()

    # Median across lineages for each protein
    median = values.median(
        axis=1,
    )

    # Median absolute deviation (MAD)
    mad = (
        values
        .sub(
            median,
            axis=0,
        )
        .abs()
        .median(
            axis=1,
        )
    )

    # Robust z-score
    robust_z = (
        values
        .sub(
            median,
            axis=0,
        )
        .div(
            mad.replace(0, np.nan),
            axis=0,
        )
        * 0.6745
    )

    # Proteins without variability across lineages
    robust_z = robust_z.fillna(0)

    # Restore protein column
    robust_z.insert(
        0,
        protein_column,
        proteins.values,
    )

    return robust_z

def protein_position_summary(
    mutations,
    major_lineages,
):
    """
    Summarise amino acid substitutions by lineage, protein and position.

    One row corresponds to one amino acid position within one protein and one
    lineage.
    """

    # Copy table
    df = mutations.copy()

    # Keep only major lineages
    selected = (
        major_lineages
        .loc[
            major_lineages["Selected"],
            [
                Columns.SUBTYPE,
                Columns.LINEAGE,
            ],
        ]
        .drop_duplicates()
    )

    df = df.merge(
        selected,
        on=[
            Columns.SUBTYPE,
            Columns.LINEAGE,
        ],
        how="inner",
    )

    # Number of genomes per lineage
    genomes = (
        df
        .groupby(
            [
                Columns.SUBTYPE,
                Columns.LINEAGE,
            ],
            as_index=False,
        )
        .agg(
            **{
                MutationColumns.GENOMES: (
                    Columns.ACCESSION,
                    "nunique",
                )
            }
        )
    )

    # Keep only genomes with mutations
    df = df[
        df[MutationColumns.HAS_MUTATIONS]
    ]

    # Position summary
    summary = (
        df
        .groupby(
            [
                Columns.SUBTYPE,
                Columns.LINEAGE,
                MutationColumns.PROTEIN,
                MutationColumns.POSITION,
                MutationColumns.REFERENCE_AA,
            ],
            as_index=False,
        )
        .agg(
            **{
                MutationColumns.TOTAL_EVENTS: (
                    MutationColumns.AA_SUBSTITUTION,
                    "size",
                ),
                MutationColumns.DISTINCT_SUBSTITUTIONS: (
                    MutationColumns.AA_SUBSTITUTION,
                    "nunique",
                ),
                MutationColumns.GENOMES_MUTATED: (
                    Columns.ACCESSION,
                    "nunique",
                ),
                MutationColumns.ALTERNATIVE_AA: (
                    MutationColumns.ALTERNATIVE_AA,
                    lambda x: "/".join(
                        sorted(x.dropna().unique())
                    ),
                ),
                MutationColumns.MOST_COMMON_SUBSTITUTION: (
                    MutationColumns.AA_SUBSTITUTION,
                    lambda x: x.value_counts().idxmax(),
                ),
                MutationColumns.MOST_COMMON_COUNT: (
                    MutationColumns.AA_SUBSTITUTION,
                    lambda x: x.value_counts().iloc[0],
                ),
            }
        )
    )

    # Genome counts
    summary = summary.merge(
        genomes,
        on=[
            Columns.SUBTYPE,
            Columns.LINEAGE,
        ],
        how="left",
    )

    # Mutation frequency
    summary[MutationColumns.FREQUENCY] = (
        summary[MutationColumns.GENOMES_MUTATED]
        / summary[MutationColumns.GENOMES]
        * 100
    ).round(2)

    #Most common count
    summary[MutationColumns.MOST_COMMON_FREQUENCY] = (
        summary[MutationColumns.MOST_COMMON_COUNT]
        / summary[MutationColumns.GENOMES]
        * 100
    ).round(2)

    #Dominance
    summary[MutationColumns.DOMINANCE] = (
        summary[MutationColumns.MOST_COMMON_COUNT]
        / summary[MutationColumns.GENOMES_MUTATED]
        * 100
    ).round(1)

    #Label protein
    summary[MutationColumns.LABEL] = (
        summary[MutationColumns.MOST_COMMON_SUBSTITUTION]
        .str.split("_")
        .str[-1]
    )

    # Protein order
    summary[MutationColumns.PROTEIN] = pd.Categorical(
        summary[MutationColumns.PROTEIN],
        categories=PROTEIN_ORDER,
        ordered=True,
    )

    # Sort
    summary = (
        summary
        .sort_values(
            [
                Columns.SUBTYPE,
                Columns.LINEAGE,
                MutationColumns.PROTEIN,
                MutationColumns.POSITION,
            ]
        )
        .reset_index(
            drop=True,
        )
    )

    # Column order
    summary = summary[
        [
            Columns.SUBTYPE,
            Columns.LINEAGE,
            MutationColumns.PROTEIN,
            MutationColumns.POSITION,
            MutationColumns.REFERENCE_AA,
            MutationColumns.ALTERNATIVE_AA,
            MutationColumns.MOST_COMMON_SUBSTITUTION,
            MutationColumns.GENOMES,
            MutationColumns.GENOMES_MUTATED,
            MutationColumns.FREQUENCY,
            MutationColumns.TOTAL_EVENTS,
            MutationColumns.DISTINCT_SUBSTITUTIONS,    
            MutationColumns.MOST_COMMON_COUNT,        
            MutationColumns.MOST_COMMON_FREQUENCY,
            MutationColumns.DOMINANCE,
            MutationColumns.LABEL,
        ]
    ]

    return summary

def protein_annotation_summary(
    summary,
    label_threshold=LOLLIPOP_LABEL_THRESHOLD,
):
    """
    Summarise hotspot mutations by annotated protein regions.

    One row corresponds to one annotated protein region (domain or feature)
    for one subtype, lineage and protein.
    """
    
    # Output rows
    rows = []

    # Iterate over subtype / lineage / protein
    for (subtype, lineage, protein), df in (
        summary.groupby(
            [
                Columns.SUBTYPE,
                Columns.LINEAGE,
                MutationColumns.PROTEIN,
            ],
            observed=False,
        )
    ):

        # Skip proteins without annotation
        if (
            subtype not in PROTEIN_FEATURES
            or protein not in PROTEIN_FEATURES[subtype]
        ):
            continue

        # Iterate over annotated regions
        for feature in PROTEIN_FEATURES[subtype][protein]:

            # Mutations within the region
            mutations = df[
                (
                    df[MutationColumns.POSITION]
                    >= feature["start"]
                )
                &
                (
                    df[MutationColumns.POSITION]
                    <= feature["end"]
                )
            ]

            major_mutations = (
                mutations[
                    mutations[
                        MutationColumns.FREQUENCY
                    ] >= CONSOLE_MUTATION_FREQUENCY
                ]
                .sort_values(
                    MutationColumns.POSITION,
                )
            )

            # Store row
            rows.append(
                {
                    Columns.SUBTYPE:
                        subtype,

                    Columns.LINEAGE:
                        lineage,

                    MutationColumns.PROTEIN:
                        protein,

                    "feature":
                        feature["label"],

                    "level":
                        feature["level"],

                    "type":
                        feature["type"],

                    "function":
                        feature["function"],

                    "start":
                        feature["start"],

                    "end":
                        feature["end"],

                    "mutations":
                        (
                            "; ".join(
                                f"{m} ({f:.0f}%)"
                                for m, f in zip(
                                    mutations[MutationColumns.LABEL],
                                    mutations[MutationColumns.FREQUENCY],
                                )
                            )
                            if not mutations.empty
                            else "—"
                        ),

                    "n_mutations":
                        mutations.shape[0],

                    "major_mutations": 
                        (
                            "; ".join(
                                f"{label} ({freq:.0f}%)"
                                for label, freq in zip(
                                    major_mutations[
                                        MutationColumns.LABEL
                                    ],
                                    major_mutations[
                                        MutationColumns.FREQUENCY
                                    ],
                                )
                            )
                            if not major_mutations.empty
                            else "No substitutions above threshold"
                        ),

                    "max_frequency":
                        (
                            mutations[
                                MutationColumns.FREQUENCY
                            ].max()
                            if len(mutations)
                            else 0
                        ),

                    "max_dominance":
                        (
                            mutations[
                                MutationColumns.DOMINANCE
                            ].max()
                            if len(mutations)
                            else 0
                        ),
                }
            )

    # Output table
    table = pd.DataFrame(rows)

    table = table.sort_values(
        [
            Columns.SUBTYPE,
            Columns.LINEAGE,
            MutationColumns.PROTEIN,
            "start",
        ]
    ).reset_index(
        drop=True,
    )

    return table

def protein_summary_report(
    annotation_summary,
    major_lineages,
):
    """
    Prepare a console summary of hotspot mutations for the selected
    lineage of each subtype.
    """

    rows = []

    for subtype in SUBTYPES:

        lineage = get_selected_lineage(
            major_lineages,
            subtype,
        )

        table = annotation_summary[
            (annotation_summary[Columns.SUBTYPE] == subtype)
            &
            (annotation_summary[Columns.LINEAGE] == lineage)
            &
            (
                annotation_summary[
                    MutationColumns.PROTEIN
                ].isin(LOLLIPOP_PROTEINS)
            )
        ].copy()

        table["Hotspot mutations"] = table[
            "hotspot_mutations"
        ].replace(
            {
                "—": "No substitutions above threshold",
            }
        )

        rows.append(table)

    return (
        pd.concat(rows)
        .reset_index(drop=True)
    )

def protein_lineage_heatmap_table(
    summary,
    subtype,
    value,
    normalize=None,
):
    """
    Prepare the protein-by-lineage heatmap table.
    """

    # Filter subtype
    df = (
        summary[
            summary[Columns.SUBTYPE] == subtype
        ]
        .copy()
    )

    # Protein × lineage matrix (values displayed in the cells)
    counts = (
        df
        .pivot(
            index=MutationColumns.PROTEIN,
            columns=Columns.LINEAGE,
            values=value,
        )
    )

    # Protein order
    counts = counts.reindex(
        PROTEIN_ORDER
    )

    # Lineage order
    lineage_order = (
        summary.loc[
            summary[Columns.SUBTYPE] == subtype,
            Columns.LINEAGE,
        ]
        .drop_duplicates()
        .tolist()
    )

    counts = counts.reindex(
        columns=lineage_order,
    )

    counts = counts.fillna(0)

    # Colours
    if normalize == "lineage":

        percentages = (
            counts
            .div(
                counts.sum(axis=0),
                axis=1,
            )
            .fillna(0)
            * 100
        )

    else:

        percentages = counts.copy()

    return (
        counts.reset_index(),
        percentages.reset_index(),
    )

def protein_lollipop_table(
    summary,
    subtype,
    lineage,
    protein,
):
    """
    Prepare the protein position table for a lollipop plot.
    """

    table = (
        summary[
            (summary[Columns.SUBTYPE] == subtype)
            &
            (summary[Columns.LINEAGE] == lineage)
            &
            (summary[MutationColumns.PROTEIN] == protein)
        ]
        .copy()
        .sort_values(
            MutationColumns.POSITION
        )
        .reset_index(drop=True)
    )

    return table

def lineage_labels(
    summary,
    subtype,
):
    """
    Create lineage labels including the number of genomes.
    -------
    
    """
    # Genome counts    
    genomes = (
        summary.loc[
            summary[Columns.SUBTYPE] == subtype,
            [
                Columns.LINEAGE,
                MutationColumns.GENOMES,
            ],
        ]
        .drop_duplicates()
    )

    # Labels
    labels = [
        f"{lineage}\n(n={genomes:,})"
        for lineage, genomes in zip(
            genomes[Columns.LINEAGE],
            genomes[MutationColumns.GENOMES],
        )
    ]
    return labels