"""
GenomeReport sections.

Defines the organization of the report, including the
table_ids and figure_ids that belong to each section.
"""

TRANSMISSION_CLUSTER_THRESHOLDS = (
    0,
    5,
)

cluster_figures = []

for threshold in TRANSMISSION_CLUSTER_THRESHOLDS:

    cluster_figures.extend(

        [
            f"transmission_cluster_by_year_{threshold}_snps_A",
            f"transmission_cluster_by_year_{threshold}_snps_B",
            f"transmission_cluster_by_community_{threshold}_snps_A",
            f"transmission_cluster_by_community_{threshold}_snps_B",
            f"transmission_lineage_heatmap_A_{threshold}_snps",
            f"transmission_lineage_heatmap_B_{threshold}_snps",
        ]

    )

REPORT_SECTIONS = {

    "overview": {

        "name": "Dataset overview",

        "table_ids": (

            "dataset_summary",
            "metadata_completeness",

        ),

        "figure_ids": (

        ),

    },

    "geography": {

        "name": "Geographical distribution",

        "table_ids": (

            "temporal_year",
            "geography_community",
            "sampling_strategy_community",

        ),

        "figure_ids": (

            "temporal_year",
            "temporal_year_percentage",

            "geography_map",
            "geography_community_barplot",

            "sampling_strategy_barplot",

        ),

    },

    "laboratory": {

        "name": "Laboratory metadata",

        "table_ids": (

            "laboratory_originating",
            "laboratory_submitter",
            "laboratory_institution",

        ),

        "figure_ids": (

            "laboratory_originating_barplot",
            "laboratory_submitter_barplot",
            "laboratory_institution_barplot",

        ),

    },

    "sequencing": {

        "name": "Sequencing metadata",

        "table_ids": (

            "sequencing_technology",
            "sequencing_assembly",

        ),

        "figure_ids": (

            "sequence_length_histogram",

            "sequencing_technology_barplot",
            "sequencing_assembly_barplot",

        ),

    },

    "patient": {

        "name": "Patient metadata",

        "table_ids": (

            "patient_host",
            "patient_gender",
            "patient_status",
            "patient_age_summary",
            "patient_age_groups",

        ),

        "figure_ids": (

            "patient_gender_barplot",
            "patient_status_barplot",
            "patient_age_groups_barplot",

        ),

    },

    "virus": {

        "name": "Virus diversity",

        "table_ids": (

            "virus_lineage_A",
            "virus_lineage_B",

            "virus_genotype_A",
            "virus_genotype_B",

            "major_lineages",

            "lineage_all_per_year",
            "lineage_A_per_year",
            "lineage_B_per_year",

            "genotype_A_per_year",
            "genotype_B_per_year",

            "lineage_A_by_community",
            "lineage_B_by_community",

            "genotype_A_by_community",
            "genotype_B_by_community",

            "lineage_A_per_year_by_community",
            "lineage_B_per_year_by_community",

            "major_lineages_community_heatmap_A",
            "major_lineages_community_heatmap_B",

        ),

        "figure_ids": (

            "virus_lineage_A_barplot",
            "virus_lineage_B_barplot",

            "virus_genotype_A_barplot",
            "virus_genotype_B_barplot",

            "lineage_all_counts_per_year",
            "lineage_all_percentage_per_year",

            "lineage_A_counts_per_year",
            "lineage_A_percentage_per_year",

            "lineage_B_counts_per_year",
            "lineage_B_percentage_per_year",

            "genotype_A_counts_per_year",
            "genotype_A_percentage_per_year",

            "genotype_B_counts_per_year",
            "genotype_B_percentage_per_year",

            "community_lineage_heatmap_A",
            "community_lineage_heatmap_B",

        ),

    },

    "mutations": {

        "name": "Protein variation",

        "table_ids": (

            "mutations",
            "mutation_counts",
            "mutation_summary",
            "mutation_frequencies",

            "protein_summary",
            "protein_summary_report",

            "protein_position_summary",
            "protein_annotation_summary",

            "protein_lineage_summary",

            "protein_lineage_distinct_substitutions_heatmap_A",
            "protein_lineage_distinct_substitutions_heatmap_A_percentages",
            "protein_lineage_distinct_substitutions_heatmap_B",
            "protein_lineage_distinct_substitutions_heatmap_B_percentages",

            "protein_lineage_comparative_heatmap_A",
            "protein_lineage_comparative_heatmap_B",

            "lollipop_A_F",
            "lollipop_A_G",
            "lollipop_A_L",

            "lollipop_B_F",
            "lollipop_B_G",
            "lollipop_B_L",

        ),

        "figure_ids": (

            "mutations_histogram",
            "mutations_boxplot",

            "protein_summary_substitutions",

            "protein_lineage_distinct_substitutions_heatmap_A",
            "protein_lineage_distinct_substitutions_heatmap_B",

            "protein_lineage_comparative_heatmap_A",
            "protein_lineage_comparative_heatmap_B",

            "lollipop_A_F",
            "lollipop_A_G",
            "lollipop_A_L",

            "lollipop_B_F",
            "lollipop_B_G",
            "lollipop_B_L",

        ),

    },

    "clusters": {

        "name": "Transmission clusters",

        "table_ids": (
            "transmission_cluster_summary",
            "transmission_cluster_by_lineage_long",
            "transmission_cluster_by_year_long",
            "transmission_cluster_by_community_long",
            "transmission_cluster_by_lineage_summary",
            "transmission_cluster_by_year_summary",
            "transmission_cluster_by_community_summary",

        ),

        "figure_ids": tuple(
            cluster_figures
        ),

    },

}

SECTION_ORDER = (
    "overview",
    "geography",
    "laboratory",
    "sequencing",
    "patient",
    "virus",
    "mutations",
    "clusters",
)