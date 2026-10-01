"""
Mutations.

Utilities to aggregate viral mutations metadata.
"""

import pandas as pd
import re

from .config import (
	Columns,
	MutationColumns,
	MutationMetrics,
)

###############################################################################
# Parsing
###############################################################################

def parse_aa_substitutions(substitutions):
    """
    Parse the GISAID AA Substitutions field.

    """

    if pd.isna(substitutions):

        return []

    substitutions = substitutions.strip()

    if substitutions in ("", "()"):

        return []

    substitutions = substitutions.strip("()")

    mutations = []

    for item in substitutions.split(","):

        item = item.strip()

        match = re.match(
            r"(.+?)_([A-Z])(\d+)([A-Z])",
            item,
        )

        if match is None:
            continue

        protein, ref, position, alt = match.groups()

        mutations.append({
        	MutationColumns.PROTEIN: protein,
        	MutationColumns.MUTATION: f"{ref}{position}{alt}",
        	MutationColumns.AA_SUBSTITUTION: item,
        	MutationColumns.REFERENCE_AA: ref,
        	MutationColumns.POSITION: int(position),
        	MutationColumns.ALTERNATIVE_AA: alt,}
        )

    return mutations

def explode_mutations(
    dataset,
):
    """
    Expand the GISAID AA substitutions field into one mutation per row.
    """

    rows = []

    df = dataset.df.copy()

    # Collection year
    df[Columns.YEAR] = (
        pd.to_datetime(
            df[Columns.COLLECTION_DATE],
            errors="coerce",
        )
        .dt.year
    )

    # Iterate over sequences
    for _, record in df.iterrows():

        substitutions = parse_aa_substitutions(
            record[Columns.AA_SUBSTITUTIONS]
        )

        # Genome without amino acid substitutions
        if not substitutions:

            rows.append({
                Columns.VIRUS_NAME: record[Columns.VIRUS_NAME],
                Columns.ACCESSION: record[Columns.ACCESSION],
                Columns.SUBTYPE: record[Columns.SUBTYPE],
                Columns.GENOTYPE: record[Columns.GENOTYPE],
                Columns.LINEAGE: record[Columns.LINEAGE],
                Columns.COMMUNITY: record[Columns.COMMUNITY],
                Columns.YEAR: record[Columns.YEAR],
                MutationColumns.PROTEIN: None,
                MutationColumns.MUTATION: None,
                MutationColumns.AA_SUBSTITUTION: None,
                MutationColumns.REFERENCE_AA: None,
                MutationColumns.POSITION: None,
                MutationColumns.ALTERNATIVE_AA: None,
                MutationColumns.HAS_MUTATIONS: False,
            })

            continue

        # Genome with substitutions
        for mutation in substitutions:

            rows.append({
                Columns.VIRUS_NAME: record[Columns.VIRUS_NAME],
                Columns.ACCESSION: record[Columns.ACCESSION],
                Columns.SUBTYPE: record[Columns.SUBTYPE],
                Columns.GENOTYPE: record[Columns.GENOTYPE],
                Columns.LINEAGE: record[Columns.LINEAGE],
                Columns.COMMUNITY: record[Columns.COMMUNITY],
                Columns.YEAR: record[Columns.YEAR],
                MutationColumns.HAS_MUTATIONS: True,
                MutationColumns.PROTEIN: mutation[
                    MutationColumns.PROTEIN
                ],
                MutationColumns.MUTATION: mutation[
                    MutationColumns.MUTATION
                ],
                MutationColumns.AA_SUBSTITUTION: mutation[
                    MutationColumns.AA_SUBSTITUTION
                ],
                MutationColumns.REFERENCE_AA: mutation[
                    MutationColumns.REFERENCE_AA
                ],
                MutationColumns.POSITION: mutation[
                    MutationColumns.POSITION
                ],
                MutationColumns.ALTERNATIVE_AA: mutation[
                    MutationColumns.ALTERNATIVE_AA
                ],
            })

    mutations = pd.DataFrame(rows)

    mutations[MutationColumns.POSITION] = (
    	mutations[MutationColumns.POSITION]
    	.astype("Int64")
    )

    mutations["Year"] = (
    	mutations["Year"]
    	.astype("Int64")
    )

    mutations[MutationColumns.HAS_MUTATIONS] = (
    	mutations[MutationColumns.HAS_MUTATIONS]
    	.astype(bool)
    )

    return mutations

###############################################################################
# Summary statistics
###############################################################################

def mutation_counts(
    mutations,
    subtype=None,
):
    """
    Count the number of amino acid substitutions per genome.

    """

    df = mutations.copy()

    # Filter subtype
    if subtype is not None:

        df = df[
            df[Columns.SUBTYPE] == subtype
        ]

    # Count mutations
    counts = (

        df.groupby(

            [
                Columns.VIRUS_NAME,
                Columns.ACCESSION,
                Columns.SUBTYPE,
                Columns.GENOTYPE,
                Columns.LINEAGE,
                Columns.COMMUNITY,
                Columns.YEAR,
            ],

            dropna=False,

        )[MutationColumns.HAS_MUTATIONS]

        .sum()

        .reset_index(

            name=MutationColumns.COUNT

        )

    )

    return counts

###############################################################################
# Generic aggregation
###############################################################################

def mutation_summary(
    mutation_counts,
):
    """
    Generate summary statistics for amino acid substitutions.
    """

    def summarize(df):

    	sequences = len(df)

    	with_mutations = (
    		df[MutationColumns.COUNT] > 0
    	).sum()

    	without_mutations = (
        	df[MutationColumns.COUNT] == 0
    	).sum()

    	return {

        	MutationMetrics.SEQUENCES:
            	sequences,

        	MutationMetrics.WITH_MUTATIONS:
            	with_mutations,

        	MutationMetrics.WITHOUT_MUTATIONS:
            	without_mutations,

        	MutationMetrics.PERCENT_WITH_MUTATIONS:
        		(
            		round(
                		with_mutations / sequences * 100,
                		2,
            		)
            		if sequences > 0
            		else 0
            	),

        	MutationMetrics.MEAN:
            	round(
                	df[MutationColumns.COUNT].mean(),
                	2,
            	) if sequences else 0,

        	MutationMetrics.MEDIAN:
            	round(
                	df[MutationColumns.COUNT].median(),
                	2,
            	) if sequences else 0,

        	MutationMetrics.MAX:
            	int(
                	df[MutationColumns.COUNT].max()
            	) if sequences else 0,

    	}

    # Summary
    all_summary = summarize(
    	mutation_counts,
    )

    a_summary = summarize(
    	mutation_counts[
    		mutation_counts[Columns.SUBTYPE] == "A"
    	]
    )

    b_summary = summarize(
    	mutation_counts[
    		mutation_counts[Columns.SUBTYPE] == "B"
    	]
    )

    summary = pd.DataFrame(
    	{
    		"Metric": all_summary.keys(),
			"All": all_summary.values(),
			"RSV-A": a_summary.values(),
			"RSV-B": b_summary.values(),
		}
	)

    return summary

def mutation_frequencies(
    mutations,
    subtype=None,
):
    """
    Calculate amino acid substitution frequencies.
    """

    df = mutations.copy()

    if subtype is not None:

    	df = df[
    		df[Columns.SUBTYPE] == subtype
    	]

    genomes = (
        df[
            [
                Columns.SUBTYPE,
                Columns.ACCESSION,
            ]
        ]
        .drop_duplicates()
        .groupby(
            Columns.SUBTYPE,
        )
        .size()
        .rename(
            MutationColumns.GENOMES,
        )
    )

    frequencies = (
        df[
            df[
                MutationColumns.HAS_MUTATIONS
            ]
        ]
        .groupby(
            [
                Columns.SUBTYPE,
                MutationColumns.PROTEIN,
                MutationColumns.AA_SUBSTITUTION,
                MutationColumns.MUTATION,
                MutationColumns.POSITION,
            ],
            as_index=False,
        )
        .agg(

            **{
                MutationColumns.COUNT: (
                    Columns.ACCESSION,
                    "nunique",
                )
            }

        )
    )

    frequencies[
        MutationColumns.GENOMES
    ] = (
        frequencies[
            Columns.SUBTYPE
        ].map(genomes)
    )

    frequencies[
        MutationColumns.FREQUENCY
    ] = (
        frequencies[
            MutationColumns.COUNT
        ]
        / frequencies[
            MutationColumns.GENOMES
        ]
        * 100
    ).round(2)

    frequencies = (
        frequencies
        .sort_values(
            [
                Columns.SUBTYPE,
                MutationColumns.COUNT,
                MutationColumns.PROTEIN,
                MutationColumns.POSITION,
            ],
            ascending=[
                True,
                False,
                True,
                True,
            ],
        )
        .reset_index(
            drop=True,
        )
    )

    return frequencies

def protein_summary(
    mutations,
    subtype=None,
):
    """
    Summarise amino acid substitutions by protein.
    """

    df = mutations.copy()

    if subtype is not None:

    	df = df[
    		df[Columns.SUBTYPE] == subtype
    	]

    genomes = (
        df[
            [
                Columns.SUBTYPE,
                Columns.ACCESSION,
            ]
        ]
        .drop_duplicates()
        .groupby(
            Columns.SUBTYPE,
        )
        .size()
    )

    summary = (
        df[
            df[
                MutationColumns.HAS_MUTATIONS
            ]
        ]
        .groupby(
            [
                Columns.SUBTYPE,
                MutationColumns.PROTEIN,
            ],
            as_index=False,
        )
        .agg(

            **{

                MutationColumns.GENOMES_WITH_MUTATIONS: (
                    Columns.ACCESSION,
                    "nunique",
                ),

                MutationColumns.TOTAL_EVENTS: (
                    MutationColumns.AA_SUBSTITUTION,
                    "count",
                ),

                MutationColumns.DISTINCT_SUBSTITUTIONS: (
                    MutationColumns.AA_SUBSTITUTION,
                    "nunique",
                ),

                MutationColumns.MAX_POSITION: (
                    MutationColumns.POSITION,
                    "max",
                ),

            }

        )
    )

    summary[
        MutationColumns.GENOMES
    ] = (
        summary[
            Columns.SUBTYPE
        ].map(genomes)
    )

    summary[
        MutationColumns.MEAN_PER_GENOME
    ] = (
        summary[
            MutationColumns.TOTAL_EVENTS
        ]
        / summary[
            MutationColumns.GENOMES
        ]
    ).round(2)

    summary[
        MutationColumns.PERCENT_GENOMES_WITH_MUTATIONS
    ] = (
        summary[
            MutationColumns.GENOMES_WITH_MUTATIONS
        ]
        / summary[
            MutationColumns.GENOMES
        ]
        * 100
    ).round(2)

    summary[
        MutationColumns.MEAN_PER_MUTATED_GENOME
    ] = (
        summary[
            MutationColumns.TOTAL_EVENTS
        ]
        / summary[
            MutationColumns.GENOMES_WITH_MUTATIONS
        ]
    ).round(2)

    summary = (
        summary
        .sort_values(
            [
                Columns.SUBTYPE,
                MutationColumns.TOTAL_EVENTS,
            ],
            ascending=[
                True,
                False,
            ],
        )
        .reset_index(
            drop=True,
        )
    )

    summary = summary[
    	[
    		Columns.SUBTYPE,
    		MutationColumns.PROTEIN,
    		MutationColumns.GENOMES,
    		MutationColumns.GENOMES_WITH_MUTATIONS,
    		MutationColumns.PERCENT_GENOMES_WITH_MUTATIONS,
    		MutationColumns.TOTAL_EVENTS,
    		MutationColumns.DISTINCT_SUBSTITUTIONS,
    		MutationColumns.MEAN_PER_GENOME,
    		MutationColumns.MEAN_PER_MUTATED_GENOME,
    		MutationColumns.MAX_POSITION,
    	]
    ]

    return summary

def filter_mutations(
    mutations,
    subtype=None,
    lineage=None,
):
    """
    Filter mutation table by subtype and/or lineage.
    """

    df = mutations.copy()

    if subtype is not None:

        df = df[
            df[Columns.SUBTYPE] == subtype
        ]

    if lineage is not None:

        df = df[
            df[Columns.LINEAGE] == lineage
        ]

    return df

###############################################################################
# Wrappers
###############################################################################

def protein_summary_report(
    protein_summary,
):
    """
    Format protein summary for console output.
    """

    df = protein_summary.copy()

    report = (
        df.pivot(
            index=MutationColumns.PROTEIN,
            columns=Columns.SUBTYPE,
            values=[
                MutationColumns.PERCENT_GENOMES_WITH_MUTATIONS,
                MutationColumns.MEAN_PER_GENOME,
                MutationColumns.DISTINCT_SUBSTITUTIONS,
            ],
        )
    )

    metric_labels = {
        MutationColumns.PERCENT_GENOMES_WITH_MUTATIONS:
            "Seq. mut. (%)",
        MutationColumns.MEAN_PER_GENOME:
            "Mean",
        MutationColumns.DISTINCT_SUBSTITUTIONS:
            "Unique AA subs",
    }

    report.columns = [
        f"{metric_labels[metric]} {subtype}"
        for metric, subtype in report.columns
    ]

    return (
        report
        .reset_index()
    )
