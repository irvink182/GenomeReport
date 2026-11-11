from pathlib import Path
from datetime import datetime
import shutil
import os
import matplotlib.pyplot as plt

from .sections import (
    REPORT_SECTIONS,
)

from .html import (
    HTMLReport,
)

from .display import (
    ReportLogger,
    print_header,
    print_dataset_summary,
    print_table,
    print_footer,
    info,
    warning,
    success,
    print_protein_summary,
    summarize_transmission_clusters,
    summarize_transmission_clusters_by_year,
    summarize_transmission_clusters_by_community,
)

from .io import (
    write_table,
    write_figure,
)

from .config import (
    Columns, 
    MutationColumns,
    CATEGORY_NAMES,
    LOLLIPOP_PROTEINS,
    TRANSMISSION_CLUSTER_THRESHOLDS,
)


from .tables import (
    dataset_summary,
    metadata_completeness,
    temporal_year,
    geography_community,
    sampling_strategy_community,
    summarize_sampling_strategy,
    laboratory_institution,
    laboratory_originating,
    laboratory_submitter,
    sequencing_technology,
    sequencing_assembly,
    patient_host,
    patient_gender,
    patient_status,
    patient_age_summary,
    patient_age_groups,
    virus_lineage,
    virus_genotype,
    major_lineages,
    community_lineage_heatmap_table,
    protein_lineage_summary,
    lineage_labels,
    protein_lineage_heatmap_table,
    protein_position_summary,
    protein_lollipop_table,
    get_selected_lineage,
    protein_annotation_summary,
    lineage_cluster_heatmap_table,
    protein_lineage_comparative_heatmap_table,
)

from .plots import(
    temporal_year_barplot,
    temporal_year_percentage,
    geography_community_barplot,
    geography_map,
    laboratory_originating_barplot,
    laboratory_institution_barplot,
    laboratory_submitter_barplot,
    sampling_strategy_barplot,
    sequence_length_histogram,
    sequencing_technology_barplot,
    sequencing_assembly_barplot,
    patient_status_barplot,
    patient_gender_barplot,
    patient_age_groups_barplot,
    virus_lineage_barplot,
    virus_genotype_barplot,
    lineage_year_areaplot,
    genotype_year_areaplot,
    community_lineage_heatmap,
    mutation_histogram,
    mutation_boxplot,
    protein_mutation_events,
    protein_distinct_substitutions,
    protein_mutation_summary_plot,
    protein_lineage_heatmap,
    protein_lollipop,
    cluster_year_barplot,
    cluster_community_barplot,
    community_lineage_heatmap,
    transmission_lineage_heatmap,
    protein_lineage_comparative_heatmap,
    genomic_cluster_spread_scatter,
    genomic_cluster_composition_barplot,
    genomic_cluster_timeline,
)

from .virus_dynamics import (
    lineage_by_year,
    genotype_by_year,
    lineage_by_community,
    genotype_by_community,
    lineage_by_year_community,
    genotype_by_year_community,

)

from .mutations import(
    explode_mutations,
    mutation_counts,
    mutation_summary,
    mutation_frequencies,
    protein_summary,
    protein_summary_report,
)

from .clusters import(
    has_clustering_data,
    cluster_summary,
    summarize_cluster_table,
    cluster_by_lineage,
    cluster_by_year,
    cluster_by_community,
    prepare_cluster_plot_data,
    cluster_by_lineage_year,
    )

from .genomic_connectivity import (
    cluster_connectivity_summary,
    cluster_connectivity_samples,
    group_connectivity,
    select_connectivity_group,
)

class DatasetReport:

    def __init__(
        self,
        dataset,
        output="output",
        phylogenies=None,
        input_file=None,
    ):

        self.dataset = dataset

        self.output = Path(output)

        self.generated_at = (
            datetime.now().astimezone()
        )

        self.input_file = (
            Path(input_file)
            if input_file is not None
            else None
        )

        self.phylogenies = [
            Path(path)
            for path in (
                phylogenies or []
            )
        ]

        self.tables = {}

        self.figures = {}

        self.statistics = {}

    def save_phylogenies(
        self,
    ):
        """
        Copy user-provided phylogenetic figures
        to the report output directory.
        """

        if not self.phylogenies:
            return

        output_dir = (
            self.output
            / "figures"
            / "phylogeny"
        )

        output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        for path in self.phylogenies:

            if not path.exists():
                raise FileNotFoundError(
                    f"Phylogenetic figure not found: {path}"
                )
            if path.suffix.lower() != ".png":
                raise ValueError(
                    "Phylogenetic figures must be PNG files: "
                    f"{path}"
                )

            destination = (
                output_dir
                / path.name
            )

            shutil.copy2(
                path,
                destination,
            )

    def build(self):

        self.output.mkdir(
            parents=True,
            exist_ok=True,
        )

        logger = ReportLogger(
            os.path.join(
                self.output,
                "report_summary.txt",
            )
        )

        clustering_available = has_clustering_data(
            self.dataset
        )
        
        print_header(
            self.dataset,
            self.output,
        )

        ###########################################################################
        # Tables
        ###########################################################################

        self.tables["dataset_summary"] = dataset_summary(
            self.dataset
        )

        self.tables["metadata_completeness"] = metadata_completeness(
            self.dataset
        )

        self.tables["temporal_year"] = (temporal_year(
            self.dataset)
        )

        self.tables["geography_community"] = (geography_community(
            self.dataset)
        )

        self.tables["sampling_strategy_community"] = sampling_strategy_community(
            self.dataset)

        self.tables["laboratory_institution"] = (laboratory_institution(
            self.dataset)
        )

        self.tables["laboratory_originating"] = (laboratory_originating(
            self.dataset)
        )

        self.tables["laboratory_submitter"] = (laboratory_submitter(
            self.dataset)
        )        

        self.tables["sequencing_technology"] = (sequencing_technology(
            self.dataset)
        )    

        self.tables["sequencing_assembly"] = (sequencing_assembly(
            self.dataset)
        )

        self.tables["patient_host"] = (
            patient_host(self.dataset)
        )

        self.tables["patient_gender"] = (
            patient_gender(self.dataset)
        )

        self.tables["patient_status"] = (
            patient_status(self.dataset)
        )

        self.tables["patient_age_summary"] = (
            patient_age_summary(self.dataset)
        )

        self.tables["patient_age_groups"] = (
            patient_age_groups(self.dataset)
        )

        # Virus
        self.tables["virus_lineage_A"] = (
            virus_lineage(
                self.dataset,
                subtype="A",
            )
        )

        self.tables["virus_lineage_B"] = (
            virus_lineage(
                self.dataset,
                subtype="B",
            )
        )

        self.tables["virus_genotype_A"] = (
            virus_genotype(
                self.dataset,
                subtype="A",
            )
        )

        self.tables["virus_genotype_B"] = (
            virus_genotype(
                self.dataset,
                subtype="B",
            )
        )

        self.tables["major_lineages"] = (
            major_lineages(
                lineages={
                    "A": self.tables["virus_lineage_A"],
                    "B": self.tables["virus_lineage_B"],
                },
                coverage=80,
            )
        )

        # Virus dynamics
        self.tables["lineage_A_per_year"] = (
            lineage_by_year(
                self.dataset,
                subtype="A",
            )
        )

        self.tables["lineage_B_per_year"] = (
            lineage_by_year(
                self.dataset,
                subtype="B",
            )
        )

        self.tables["lineage_all_per_year"] = (
            lineage_by_year(
                self.dataset,
                subtype=None,
            )
        )

        self.tables["genotype_A_per_year"] = (
            genotype_by_year(
                self.dataset,
                subtype="A",
            )
        )

        self.tables["genotype_B_per_year"] = (
            genotype_by_year(
                self.dataset,
                subtype="B",
            )
        )

        self.tables["lineage_A_by_community"] = (
            lineage_by_community(
                self.dataset,
                subtype="A",
            )
        )

        self.tables["lineage_B_by_community"] = (
            lineage_by_community(
                self.dataset,
                subtype="B",
            )
        )

        self.tables["genotype_A_by_community"] = (
            genotype_by_community(
                self.dataset,
                subtype="A",
            )
        )

        self.tables["genotype_B_by_community"] = (
            genotype_by_community(
                self.dataset,
                subtype="B",
            )
        )

        self.tables["lineage_A_per_year_by_community"] = (
            lineage_by_year_community(
                self.dataset,
                subtype="A",
            )
        )

        self.tables["lineage_B_per_year_by_community"] = (
            lineage_by_year_community(
                self.dataset,
                subtype="B",
            )
        )

        self.tables["major_lineages_community_heatmap_A"] = (
            community_lineage_heatmap_table(
                lineage_table=self.tables["lineage_A_by_community"],
                major_lineages=self.tables["major_lineages"],
                subtype="A",
                )
        )

        self.tables["major_lineages_community_heatmap_B"] = (
            community_lineage_heatmap_table(
                lineage_table=self.tables["lineage_B_by_community"],
                major_lineages=self.tables["major_lineages"],
                subtype="B",
                )
        )

        # Mutations/proteins tables
        self.tables["mutations"] = explode_mutations(
            self.dataset,
        )
        
        self.tables["mutation_counts"] = mutation_counts(
            self.tables["mutations"],
        )

        self.tables["mutation_summary"] = mutation_summary(
            self.tables["mutation_counts"],
        )

        self.tables["mutation_frequencies"] = mutation_frequencies(
            self.tables["mutations"],
        )
        
        self.tables["protein_summary"] = protein_summary(
            self.tables["mutations"],
        )
        
        self.tables["protein_summary_report"] = protein_summary_report(
            self.tables["protein_summary"],
        )

        self.tables["protein_lineage_summary"] = (
            protein_lineage_summary(
                mutations=self.tables["mutations"],
                major_lineages=self.tables["major_lineages"],
            )
        )

        # Mutation Heatmaps exploratory tables
        (
            self.tables["protein_lineage_distinct_substitutions_heatmap_A"],
            self.tables["protein_lineage_distinct_substitutions_heatmap_A_percentages"],
        ) = protein_lineage_heatmap_table(
            summary=self.tables["protein_lineage_summary"],
            subtype="A",
            value=MutationColumns.DISTINCT_SUBSTITUTIONS,
            normalize="lineage",
        )

        (
            self.tables["protein_lineage_distinct_substitutions_heatmap_B"],
            self.tables["protein_lineage_distinct_substitutions_heatmap_B_percentages"],
        ) = protein_lineage_heatmap_table(
            summary=self.tables["protein_lineage_summary"],
            subtype="B",
            value=MutationColumns.DISTINCT_SUBSTITUTIONS,
            normalize="lineage",
        )

        # Mutation Heatmaps comparative tables
        self.tables[
            "protein_lineage_comparative_heatmap_A"
        ] = protein_lineage_comparative_heatmap_table(
            percentages=self.tables[
                "protein_lineage_distinct_substitutions_heatmap_A_percentages"
            ],
        )

        self.tables[
            "protein_lineage_comparative_heatmap_B"
        ] = protein_lineage_comparative_heatmap_table(
            percentages=self.tables[
                "protein_lineage_distinct_substitutions_heatmap_B_percentages"
            ],
        )


        self.tables["protein_position_summary"] = (
            protein_position_summary(
                mutations=self.tables["mutations"],
                major_lineages=self.tables["major_lineages"],
            )
        )

        self.tables["protein_annotation_summary"] = (
            protein_annotation_summary(
                summary=self.tables["protein_position_summary"],
            )
        )

        summary = self.tables["protein_position_summary"]

        subtypes = (
            self.tables["major_lineages"]
            .loc[
                self.tables["major_lineages"]["Selected"],
                "Subtype",
            ]
            .dropna()
            .unique()
        )

        for subtype in subtypes:
            lineage = get_selected_lineage(
                self.tables["major_lineages"],
                subtype
            )

            for protein in LOLLIPOP_PROTEINS:
                self.tables[
                    f"lollipop_{subtype}_{protein}"
                ] = (
                    protein_lollipop_table(
                        summary=summary,
                        subtype=subtype,
                        lineage=lineage,
                        protein=protein,
                    )
                )

        if has_clustering_data(
            self.dataset
        ):
            summary = cluster_summary(
                self.dataset,
                groupby=Columns.SUBTYPE,
            )

            self.tables["transmission_cluster_summary"] = summary

            lineage = cluster_by_lineage(
                self.dataset,
                groupby=Columns.SUBTYPE,
            )

            year = cluster_by_year(
                self.dataset,
                groupby=Columns.SUBTYPE,
            )

            community = cluster_by_community(
                self.dataset,
                groupby=Columns.SUBTYPE,
            )

            lineage_year = cluster_by_lineage_year(
                self.dataset,
                groupby=Columns.SUBTYPE,
            )

            self.tables["transmission_cluster_by_lineage_long"] = lineage
            self.tables["transmission_cluster_by_year_long"] = year
            self.tables["transmission_cluster_by_community_long"] = community
            self.tables["transmission_cluster_by_lineage_year"] = lineage_year

            self.tables["transmission_cluster_by_lineage_summary"] = (
                summarize_cluster_table(lineage)
            )

            self.tables["transmission_cluster_by_year_summary"] = (
                summarize_cluster_table(year)
            )

            self.tables["transmission_cluster_by_community_summary"] = (
                summarize_cluster_table(community)
            )

            #Genomic connectivity tables
            self.connectivity_group_variable = (
                select_connectivity_group(
                    self.dataset
                )
            )

            group_variable = (
                self.connectivity_group_variable
            )

            if group_variable is not None:
                self.tables[
                    "cluster_connectivity_summary"
                ] = cluster_connectivity_summary(
                    self.dataset,
                    group_variable=group_variable,
                )

                self.tables[
                    "cluster_connectivity_samples"
                ] = cluster_connectivity_samples(
                    self.dataset,
                    group_variable=group_variable,
                )

                self.tables[
                    "group_connectivity"
                ] = group_connectivity(
                    self.dataset,
                    group_variable=group_variable,
                )


        ###########################################################################
        # Figures
        ###########################################################################

        # Overviwew
        self.figures["temporal_year"] = temporal_year_barplot(
            self.tables["temporal_year"]
        )

        self.figures["temporal_year_percentage"] = (
            temporal_year_percentage(
                self.tables["temporal_year"]
            )
        )

        self.figures["sequence_length_histogram"] = (
            sequence_length_histogram(
                self.dataset
            )
        )

        self.figures["geography_community_barplot"] = geography_community_barplot(
            self.tables["geography_community"]
        )

        self.figures["laboratory_originating_barplot"] = laboratory_originating_barplot(
            self.tables["laboratory_originating"]
        )

        self.figures["laboratory_institution_barplot"] = laboratory_institution_barplot(
            self.tables["laboratory_institution"]
        )

        self.figures["laboratory_submitter_barplot"] = laboratory_submitter_barplot(
            self.tables["laboratory_submitter"]
        )

        self.figures["sampling_strategy_barplot"] = sampling_strategy_barplot(
            self.tables["sampling_strategy_community"]
        )

        self.figures["sequencing_technology_barplot"] = sequencing_technology_barplot(
            self.tables["sequencing_technology"]
        )

        self.figures["sequencing_assembly_barplot"] = sequencing_assembly_barplot(
            self.tables["sequencing_assembly"]
        )

        if not self.tables["patient_gender"].empty:
            self.figures["patient_gender_barplot"] = (
                patient_gender_barplot(
                    self.tables["patient_gender"]
                )
            )

        if not self.tables["patient_status"].empty:
            self.figures["patient_status_barplot"] = (
                patient_status_barplot(
                    self.tables["patient_status"]
                )
            )

        if not self.tables["patient_age_groups"].empty:
            self.figures["patient_age_groups_barplot"] = (
                patient_age_groups_barplot(
                    self.tables["patient_age_groups"]
                )
            )

        if not self.tables["virus_lineage_A"].empty:
            self.figures["virus_lineage_A_barplot"] = (
                virus_lineage_barplot(
                    self.tables["virus_lineage_A"],
                    subtype="A",
                )
            )

        if not self.tables["virus_lineage_B"].empty:
            self.figures["virus_lineage_B_barplot"] = (
                virus_lineage_barplot(
                    self.tables["virus_lineage_B"],
                    subtype="B",
                )
            )

        if not self.tables["virus_genotype_A"].empty:
            self.figures["virus_genotype_A_barplot"] = (
                virus_genotype_barplot(
                    self.tables["virus_genotype_A"],
                    subtype="A",
                )
            )

        if not self.tables["virus_genotype_B"].empty:
            self.figures["virus_genotype_B_barplot"] = (
                virus_genotype_barplot(
                    self.tables["virus_genotype_B"],
                    subtype="B",
                )
            )

        self.figures["geography_map"] = (
            geography_map(
                table=self.tables["geography_community"],
                title="Geographic distribution of RSV sequences across Spain",
            )
        )

        # Virus dynamics plots
        if not self.tables["lineage_A_per_year"].empty:
            self.figures["lineage_A_counts_per_year"] = (
                lineage_year_areaplot(
                    self.tables["lineage_A_per_year"],
                    subtype="A",
                    relative=False,
                )
            )

            self.figures["lineage_A_percentage_per_year"] = (
                lineage_year_areaplot(
                    self.tables["lineage_A_per_year"],
                    subtype="A",
                    relative=True,
                )
            )

        if not self.tables["lineage_B_per_year"].empty:
            self.figures["lineage_B_counts_per_year"] = (
                lineage_year_areaplot(
                    self.tables["lineage_B_per_year"],
                    subtype="B",
                    relative=False,
                )
            )

            self.figures["lineage_B_percentage_per_year"] = (
                lineage_year_areaplot(
                    self.tables["lineage_B_per_year"],
                    subtype="B",
                    relative=True,
                )
            )

        self.figures["lineage_all_counts_per_year"] = (
            lineage_year_areaplot(
                self.tables["lineage_all_per_year"],
                subtype=None,
                relative=False,
            )
        )

        self.figures["lineage_all_percentage_per_year"] = (
            lineage_year_areaplot(
                self.tables["lineage_all_per_year"],
                subtype=None,
                relative=True,
            )
        )


        if not self.tables["genotype_A_per_year"].empty:
            self.figures["genotype_A_counts_per_year"] = (
                genotype_year_areaplot(
                    self.tables["genotype_A_per_year"],
                    subtype="A",
                    relative=False,
                )
            )


            self.figures["genotype_A_percentage_per_year"] = (
                genotype_year_areaplot(
                    self.tables["genotype_A_per_year"],
                    subtype="A",
                    relative=True,
                )
            )
        
        if not self.tables["genotype_B_per_year"].empty:
            self.figures["genotype_B_counts_per_year"] = (
                genotype_year_areaplot(
                    self.tables["genotype_B_per_year"],
                    subtype="B",
                    relative=False,
                )
            )

            self.figures["genotype_B_percentage_per_year"] = (
                genotype_year_areaplot(
                    self.tables["genotype_B_per_year"],
                    subtype="B",
                    relative=True,
                )
            )

        if not self.tables[
            "major_lineages_community_heatmap_A"
        ].empty:
            self.figures["community_lineage_heatmap_A"] = (
                community_lineage_heatmap(
                    table=self.tables[
                        "major_lineages_community_heatmap_A"
                    ],
                    subtype="A",
                )
            )

        if not self.tables[
            "major_lineages_community_heatmap_B"
        ].empty:
            self.figures["community_lineage_heatmap_B"] = (
                community_lineage_heatmap(
                    table=self.tables[
                        "major_lineages_community_heatmap_B"
                    ],
                    subtype="B",
                )
            )

        # Mutations/Proteins plots
        self.figures["mutations_histogram"] = (
            mutation_histogram(
                self.tables["mutation_counts"],
            )
        )

        self.figures["mutations_boxplot"] = (
            mutation_boxplot(
                self.tables["mutation_counts"],
            )
        )

        self.figures["protein_summary_substitutions"] = (
            protein_mutation_summary_plot(
                self.tables["protein_summary"],
            )
        )

        if (
            self.tables[
                "protein_lineage_distinct_substitutions_heatmap_A"
            ].shape[1] > 1
            and self.tables[
                "protein_lineage_distinct_substitutions_heatmap_A_percentages"
            ].shape[1] > 1
        ):
            self.figures[
                "protein_lineage_distinct_substitutions_heatmap_A"
            ] = protein_lineage_heatmap(
                table=self.tables[
                    "protein_lineage_distinct_substitutions_heatmap_A"
                ],
                percentages=self.tables[
                    "protein_lineage_distinct_substitutions_heatmap_A_percentages"
                ],
                summary=self.tables[
                    "protein_lineage_summary"
                ],
                subtype="A",
                title=(
                    "RSV-A distribution of distinct "
                    "amino acid substitutions"
                ),
            )

        if (
            self.tables[
                "protein_lineage_distinct_substitutions_heatmap_B"
            ].shape[1] > 1
            and self.tables[
                "protein_lineage_distinct_substitutions_heatmap_B_percentages"
            ].shape[1] > 1
        ):
            self.figures[
                "protein_lineage_distinct_substitutions_heatmap_B"
            ] = protein_lineage_heatmap(
                table=self.tables[
                    "protein_lineage_distinct_substitutions_heatmap_B"
                ],
                percentages=self.tables[
                    "protein_lineage_distinct_substitutions_heatmap_B_percentages"
                ],
                summary=self.tables[
                    "protein_lineage_summary"
                ],
                subtype="B",
                title=(
                    "RSV-B distribution of distinct "
                    "amino acid substitutions"
                ),
            )


        if (
            self.tables[
                "protein_lineage_comparative_heatmap_A"
            ].shape[1] > 1
        ):
            self.figures[
                "protein_lineage_comparative_heatmap_A"
            ] = protein_lineage_comparative_heatmap(
                table=self.tables[
                    "protein_lineage_comparative_heatmap_A"
                ],
                summary=self.tables[
                    "protein_lineage_summary"
                ],
                subtype="A",
                title=(
                    "RSV-A relative amino acid substitution patters "
                    "across lineages"
                ),
            )

        if (
            self.tables[
                "protein_lineage_comparative_heatmap_B"
            ].shape[1] > 1
        ):
            self.figures[
                "protein_lineage_comparative_heatmap_B"
            ] = protein_lineage_comparative_heatmap(
                table=self.tables[
                    "protein_lineage_comparative_heatmap_B"
                ],
                summary=self.tables[
                    "protein_lineage_summary"
                ],
                subtype="B",
                title=(
                    "RSV-B relative amino acid substitution patters "
                    "across lineages"
                ),
            )

        # Lollipop plots
        for subtype in ["A", "B"]:
            for protein in LOLLIPOP_PROTEINS:
                key = f"lollipop_{subtype}_{protein}"

                if (
                    key in self.tables
                    and not self.tables[key].empty
                ):
                    self.figures[key] = (
                        protein_lollipop(
                            table=self.tables[key],
                        )
                    )
        
        # Clustering plots
        if has_clustering_data(
            self.dataset
        ):
            cluster_subtypes = (
                self.tables["major_lineages"][Columns.SUBTYPE]
                .dropna()
                .unique()
            )
            
            for threshold in TRANSMISSION_CLUSTER_THRESHOLDS:
                #BY year
                cluster_year = prepare_cluster_plot_data(
                    self.tables["transmission_cluster_by_year_summary"],
                    threshold=threshold,
                )

                for subtype in cluster_subtypes:
                    self.figures[
                        f"transmission_cluster_by_year_{threshold}_snps_{subtype}"
                    ] = cluster_year_barplot(
                        cluster_year,
                        subtype=subtype,
                        threshold=threshold,
                    )

                #BY Community
                cluster_community = prepare_cluster_plot_data(
                    self.tables["transmission_cluster_by_community_summary"],
                    threshold=threshold,
                )

                for subtype in cluster_subtypes:
                    self.figures[
                        f"transmission_cluster_by_community_{threshold}_snps_{subtype}"
                    ] = cluster_community_barplot(
                        cluster_community,
                        subtype=subtype,
                        threshold=threshold,
                    )

            # Transmission lineage heatmaps
            for threshold in TRANSMISSION_CLUSTER_THRESHOLDS:
                for subtype in cluster_subtypes:
                    counts, percentages, lineage_totals = lineage_cluster_heatmap_table(
                        self.tables["transmission_cluster_by_lineage_year"],
                        self.tables["major_lineages"],
                        subtype=subtype,
                        threshold=threshold,
                    )
                    self.figures[
                        f"transmission_lineage_heatmap_{subtype}_{threshold}_snps"
                    ] = transmission_lineage_heatmap(
                        table=counts,
                        percentages=percentages,
                        lineage_totals=lineage_totals,
                        subtype=subtype,
                        threshold=threshold,
                    )

            #Genomic connectivity plots
            if group_variable is not None:
                for threshold in TRANSMISSION_CLUSTER_THRESHOLDS:
                    
                    #Cluster spread
                    self.figures[
                        f"genomic_cluster_spread_{threshold}snp"
                    ] = genomic_cluster_spread_scatter(
                        self.tables["cluster_connectivity_summary"],
                        threshold=threshold,
                        group_variable=group_variable,
                    )

                    #Cluster composition
                    for subtype in cluster_subtypes:

                        self.figures[
                            f"genomic_cluster_composition_{threshold}snp_{subtype}"
                        ] = genomic_cluster_composition_barplot(
                            self.tables["cluster_connectivity_samples"],
                            subtype=subtype,
                            threshold=threshold,
                            group_variable=group_variable,
                        )

                        #Cluster timeline
                        self.figures[
                            f"genomic_cluster_timeline_{threshold}snp_{subtype}"
                        ] = genomic_cluster_timeline(
                            self.tables["cluster_connectivity_samples"],
                            subtype=subtype,
                            threshold=threshold,
                            group_variable=group_variable,
                        )




        ###########################################################################
        # Console
        ###########################################################################

        print_dataset_summary(
            self.tables["dataset_summary"],
            logger=logger,
        )

        print_table(
            "Metadata completeness",
            self.tables["metadata_completeness"],
            logger=logger,
        )

        print_table(
            "Temporal distribution by year",
            self.tables["temporal_year"],
            logger=logger,
        )

        print_table(
            "Geographic distribution by autonomous community",
            self.tables["geography_community"],
            logger=logger,
        )

        print_table(
            "Originating laboratories",
            self.tables["laboratory_originating"],
            max_rows=8,
            sort_by="Total",
            logger=logger,
        )

        print_table(
            "Submitter institutions",
            self.tables["laboratory_submitter"],
            max_rows=8,
            sort_by="Total",
            logger=logger,
        )

        print_table(
            "Sampling strategy",
            summarize_sampling_strategy(
                self.tables[
                    "sampling_strategy_community"
                ]
            ),
            logger=logger,
        )

        print_table(
            "Sequecing technology",
            self.tables["sequencing_technology"],
            sort_by="Total",
            logger=logger,
        )

        print_table(
            "Assembly methods",
            self.tables["sequencing_assembly"],
            max_rows=8,
            sort_by="Total",
            logger=logger,
        )

        print_table(
            "Host",
            self.tables["patient_host"],
            logger=logger,
        )

        print_table(
            "Gender",
            self.tables["patient_gender"],
            sort_by="Total",
            logger=logger,
        )

        print_table(
            "Patient status",
            self.tables["patient_status"],
            sort_by="Total",
            logger=logger,
        )

        print_table(
            "Age summary",
            self.tables["patient_age_summary"],
            logger=logger,
        )

        print_table(
            "Age groups",
            self.tables["patient_age_groups"],
            logger=logger,
        )

        print_table(
            "Virus lineages (RSV-A)",
            self.tables["virus_lineage_A"],
            max_rows=6,
            logger=logger,
        )

        print_table(
            "Virus lineages (RSV-B)",
            self.tables["virus_lineage_B"],
            max_rows=6,
            logger=logger,
        )

        print_table(
            "Mutation summary",
            self.tables["mutation_summary"],
            logger=logger,
        )
        
        print_table(
            "Protein mutation summary",
            self.tables["protein_summary_report"],
            logger=logger,
        )

        print_protein_summary(
            annotation_summary=self.tables[
                "protein_annotation_summary"
            ],
            major_lineages=self.tables[
                "major_lineages"
            ],
            logger=logger,
        )

        if clustering_available:

            print_table(
                "Transmission clustering",
                summarize_transmission_clusters(
                    self.tables["transmission_cluster_summary"]
                ),
                logger=logger,
            )

            print_table(
                "Transmission clusters by autonomous community",
                summarize_transmission_clusters_by_community(
                    self.tables["transmission_cluster_by_community_summary"]
                    ),
                max_rows=10,
                logger=logger,
            )

            print_table(
                "Transmission clusters by year",
                summarize_transmission_clusters_by_community(
                    self.tables["transmission_cluster_by_year_summary"]
                ),
                logger=logger,
            )

        logger.close()


        ###########################################################################
        # Save
        ###########################################################################

        self.save()

        self.save_phylogenies()

        HTMLReport(
            self,
        ).build()

        print_footer()

        return self

    def save(self):
        """
        Save all report tables and figures.
        """

        # Create output directories
        tables_dir = (
            self.output /
            "tables"
        )

        figures_dir = (
            self.output /
            "figures"
        )

        tables_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        figures_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        for section in REPORT_SECTIONS:

            (
                tables_dir /
                section
            ).mkdir(
                exist_ok=True,
            )

            (
                figures_dir /
                section
            ).mkdir(
                exist_ok=True,
            )

        # Export tables
        for section, content in REPORT_SECTIONS.items():

            for table_id in content["table_ids"]:

                if table_id not in self.tables:
                    continue

                write_table(
                    self.tables[
                        table_id
                    ],
                    (
                        tables_dir /
                        section /
                        f"{table_id}.tsv"
                    ),
                )

        # Export figures
        for section, content in REPORT_SECTIONS.items():

            for figure_id in content["figure_ids"]:

                if figure_id not in self.figures:
                    continue

                write_figure(
                    self.figures[
                        figure_id
                    ],
                    (
                        figures_dir /
                        section /
                        figure_id
                    ),
                )

                plt.close(
                    self.figures[figure_id]
                )

        success(
            "Report successfully exported."
        )