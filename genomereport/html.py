"""
HTML report generation.
"""

from pathlib import Path
from html import escape
from .config import Columns

import pandas as pd

ENGLISH_MONTHS = (
    "January",
    "February",
    "March",
    "April",
    "May",
    "June",
    "July",
    "August",
    "September",
    "October",
    "November",
    "December",
)

class HTMLReport:
    """
    GenomeReport HTML report.
    """

    def __init__(
        self,
        report,
    ):

        self.report = report

    def format_date_english(
        self,
        date,
        include_time=False,
    ):
        """
        Format dates using English month names
        independently of the system locale.
        """

        month = ENGLISH_MONTHS[
            date.month - 1
        ]

        formatted = (
            f"{date.day:02d} "
            f"{month} "
            f"{date.year}"
        )

        if include_time:
            formatted += (
                f", {date:%H:%M %Z}"
            )

        return formatted

    def format_month_year_english(
        self,
        date,
    ):
        """
        Format month and year using English
        month names independently of the system locale.
        """

        return (
            f"{ENGLISH_MONTHS[date.month - 1]} "
            f"{date.year}"
        )

    def build(self):
        """
        Build HTML report.
        """

        package_dir = Path(__file__).parent

        # Read HTML template
        template = (
            package_dir
            / "templates"
            / "report.html"
        )

        html = template.read_text(
            encoding="utf-8",
        )

        # Read stylesheet
        stylesheet = (
            package_dir
            / "static"
            / "style.css"
        )

        style = stylesheet.read_text(
            encoding="utf-8",
        )

        # Read JavaScript
        script_file = (
            package_dir
            / "static"
            / "report.js"
        )

        script = script_file.read_text(
            encoding="utf-8",
        )

        # Build report content
        content = self.build_content()

        # Replace template placeholder
        html = html.replace(
            "{{ STYLE }}",
            style,
        )
        
        html = html.replace(
            "{{ CONTENT }}",
            content,
        )

        html = html.replace(
            "{{ SCRIPT }}",
            script,
        )

        # Save report
        output = (
            self.report.output /
            "report.html"
        )

        output.write_text(
            html,
            encoding="utf-8",
        )

    def report_metadata(
        self,
    ):
        """
        Render report metadata.
        """

        generated = self.format_date_english(
            self.report.generated_at
        )
        
        return f"""
    <div class="report-metadata">

        <div class="report-metadata-item">
            <span class="report-metadata-label">
                Author
            </span>
            <span class="report-metadata-value">
                Irving Cancino-Muñoz
            </span>
        </div>

        <div class="report-metadata-item">
            <span class="report-metadata-label">
                Institution
            </span>
            <span class="report-metadata-value">
                UV-I2SySBIO & FISABIO Public Health
            </span>
        </div>

        <div class="report-metadata-item">
            <span class="report-metadata-label">
                Generated
            </span>
            <span class="report-metadata-value">
                {escape(generated)}
            </span>
        </div>
    </div>
    """

    def report_text(
        self,
        text,
    ):
        """
        Render descriptive report text.
        """

        return f"""
    <div class="report-text">
        <p>
            {text}
        </p>
    </div>
    """

    def table(
        self,
        table_id,
        interactive=False,
        max_rows=None,
    ):
        """
        Render a report table.

        Parameters
        ----------
        table_id
            Identifier of the table stored in report.tables.

        interactive
            Add search, sorting and pagination controls.

        max_rows
            Optional maximum number of rows included in the HTML.
            If None, all rows are included.
        """

        if table_id not in self.report.tables:
            return ""

        # Get table
        table = self.report.tables[
            table_id
        ].copy()

        if max_rows is not None:
            table = table.head(
                max_rows
            )

        def format_value(
            value,
        ):
            """
            Format table returns for HTML display.
            """

            if pd.isna(value):
                return "NA"

            return escape(
                str(value)
            )

        # CSS classes
        classes = [
            "report-table",
        ]

        if interactive:
            classes.append(
                "interactive-table"
            )

        class_name = " ".join(
            classes
        )

        # Table header
        header = "".join(
            f"<th>{escape(str(column))}</th>"
            for column in table.columns
        )

        # Table body
        rows = []

        for row in table.itertuples(
            index=False,
            name=None,
        ):
            cells = "".join(
                f"<td>{format_value(value)}</td>"
                for value in row
            )

            rows.append(
                f"<tr>{cells}</tr>"
            )

        body = "\n".join(
            rows
        )

        return f"""
    <div class="table-container">

        <table class="{class_name}">

            <thead>
                <tr>
                    {header}
                </tr>
            </thead>

            <tbody>
                {body}
            </tbody>

        </table>

    </div>
    """

    def figure(
        self,
        figure_id,
        section,
        caption=None,
        alt=None,
        max_width=None,
    ):
        """
        Render a report figure.
        """

        figure_path = (
            Path("figures")
            / section
            / f"{figure_id}.png"
        )
        if alt is None:
            alt = figure_id.replace(
                "_",
                " ",
            )

        caption_html = ""

        if caption is not None:
            
            caption_html = f"""
            <figcaption class="figure-caption">
                {escape(caption)}
            </figcaption>
            """

        image_style = ""

        if max_width is not None:
            image_style= (
                f' style="max-width: {max_width}px;"'
            )

        return f"""
    <figure class="figure-container">

        <img
            src="{figure_path.as_posix()}"
            alt="{escape(alt)}"
            loading="lazy"
            {image_style}
        >

        {caption_html}

    </figure>
    """

    def has_clustering(
        self,
    ):
        """
        Check whether clustering results are available.
        """

        return (
            "transmission_cluster_summary"
            in self.report.tables
        )

    def connectivity_group_label(
        self,
        plural=False,
    ):
        """
        Return the display label used for genomic
        connectivity grouping.
        """

        group_variable = getattr(
            self.report,
            "connectivity_group_variable",
            None,
        )

        if group_variable == Columns.COMMUNITY:
            return (
                "Communities"
                if plural
                else "Community"
            )

        if group_variable == Columns.INSTITUTION:
            return (
                "Institutions"
                if plural
                else "Institution"
            )

        return (
            "Groups"
            if plural
            else "Group"
        )

    def build_content(
        self,
    ):
        """
        Build HTML report content.
        """

        html = ""
        html += self.header()
        html += self.summary()
        html += self.sample_characteristics()
        html += self.laboratory_contribution()
        html += self.virus_diversity()
        html += self.protein_summary()

        #Clustering section
        if self.has_clustering():
            html += self.clustering()

        #Phylogenies
        if self.report.phylogenies:
            html += self.phylogeny()

        #Analysis information footer
        html += self.analysis_information()

        return html

    def header(
        self,
    ):
        """
        Report header.
        """

        return f"""
    <header class="report-header">

        <h1>
            Respiratory Syncytial Virus (RSV) Sequence Data Report 
        </h1>
        <p class="report-subtitle">
            Descriptive genomic surveillance report for
            RSV sequences from Spain
        </p>

        {self.report_metadata()}

    </header>
    """

    def table_metric(
        self,
        table_id,
        metric,
    ):
        """
        Get a value from a Metric/Value table.
        """

        table = self.report.tables[
            table_id
        ]

        values = table.loc[
            table["Metric"] == metric,
            "Value",
        ]

        if values.empty:
            return None

        return values.iloc[0]

    def dataset_metric(
        self,
        metric,
    ):
        """
        Get a value from the dataset summary table.
        """

        return self.table_metric(
            table_id="dataset_summary",
            metric=metric,
        )

    def metric_card(
        self,
        value,
        label,
        css_class=None,
    ):
        """
        Render a summary metric card.
        """

        classes = [
            "metric-card",
        ]

        if css_class is not None:
            classes.append(
                css_class
            )

        class_name = " ".join(
            classes
        )

        return f"""
    <div class="{class_name}">

        <span class="metric-value">
            {escape(str(value))}
        </span>

        <span class="metric-label">
            {escape(label)}
        </span>

    </div>
    """

    def secondary_metric(
        self,
        value,
        label,
    ):
        """
        Render a compact secondary dataset metric.
        """

        return f"""
    <div class="secondary-metric">

        <span class="secondary-metric-value">
            {escape(str(value))}
        </span>

        <span class="secondary-metric-label">
            {escape(label)}
        </span>

    </div>
    """

    def metadata_completeness(
        self,
    ):
        """
        Render metadata completeness as horizontal bars.
        """

        table = self.report.tables[
            "metadata_completeness"
        ]

        rows = []

        for row in table.itertuples(
            index=False,
        ):

            variable = row[0]
            missing_percentage = float(
                row[3]
            )

            completeness = (
                100.0
                - missing_percentage
            )

            rows.append(
                f"""
        <div class="completeness-row">

            <div class="completeness-variable">
                {escape(str(variable))}
            </div>

            <div class="completeness-bar">
                <div
                    class="completeness-fill"
                    style="width: {completeness:.2f}%"
                ></div>
            </div>

            <div class="completeness-value">
                {completeness:.1f}%
            </div>

        </div>
        """
            )
        return "\n".join(
            rows
        )

    def metadata_completeness_text(
        self,
        n=3,
    ):
        """
        Generate descriptive text for variables
        with the lowest metadata completeness.
        """

        table = self.report.tables[
            "metadata_completeness"
        ].copy()

        table["Completeness (%)"] = (
            100
            - table["Missing (%)"]
        )

        # Sampling strategy uses empty values to represent
        # sentinel-surveillance records and is therefore
        # not treated as missing metadata in this summary.
        table = table[
            table["Variable"]
            != "sampling_strategy"
        ]

        lowest = (
            table
            .sort_values(
                "Completeness (%)",
                ascending=True,
            )
            .head(n)
        )

        values = []
        for row in lowest.itertuples(
            index=False,
        ):
            variable = (
                str(row[0])
                .replace("_", " ")
            )

            completeness = row[-1]

            values.append(
                f"{variable} "
                f"({completeness:.2f}%)"
            )

            variables = (
            ", ".join(values[:-1])
            + " and "
            + values[-1]
        )

        return (
            f"The variables with the lowest metadata completeness "
            f"were {variables}. "
            f"For sampling strategy, empty values represent "
            f"sentinel-surveillance records by default, while "
            f"non-sentinel sampling strategies are explicitly recorded "
            f"in the sampling_strategy field."
        )

    def collapsible_table(
        self,
        table_id,
        label="View data",
        interactive=True,
    ):
        """
        Render a table inside a collapsible HTML block.
        """

        table_html = self.table(
            table_id=table_id,
            interactive=interactive,
        )

        if not table_html:
            return ""

        return f"""
    <details class="collapsible-table">
        <summary>
            {escape(label)}
        </summary>

        <div class="collapsible-table-content">
            {table_html}
        </div>
    </details>
    """

    def figure_grid(
        self,
        figures,
        columns=2,
    ):
        """
        Render multiple figures in a responsive grid.
        """

        content = "\n".join(
            figures
        )

        return f"""
    <div class="figure-grid figure-grid-{columns}">

        {content}

    </div>
    """

    def temporal_distribution_text(
        self,
    ):
        """
        Generate descriptive text for the temporal
        distribution of RSV sequences.
        """

        table = self.report.tables[
            "temporal_year"
        ].copy()

        first_year = int(
            table["Year"].min()
        )
        last_year = int(
            table["Year"].max()
        )

        peak = table.loc[
            table["Total"].idxmax()
        ]
        peak_year = int(
            peak["Year"]
        )
        peak_total = int(
            peak["Total"]
        )

        subtype_counts = []

        if "A" in table.columns:
            subtype_counts.append(
                f"{int(peak['A']):,} RSV-A"
            )

        if "B" in table.columns:
            subtype_counts.append(
                f"{int(peak['B']):,} RSV-B"
            )

        if len(subtype_counts) == 2:
            peak_description = (
                f" ({subtype_counts[0]} and "
                f"{subtype_counts[1]})"
            )
        else:
            peak_description = ""

        return (
            f"The temporal distribution includes RSV sequences "
            f"collected between {first_year} and {last_year}. "
            f"The largest number of sequences was recorded in "
            f"{peak_year}, with {peak_total:,} sequences"
            f"{peak_description}."
        )

    def geographic_distribution_text(
        self,
        n=3,
    ):
        """
        Generate descriptive text for the geographic
        distribution of RSV sequences.
        """

        table = self.report.tables[
            "geography_community"
        ].copy()

        table = table.sort_values(
            "Total",
            ascending=False,
        )

        n_communities = (
            table[
                "autonomous_community"
            ]
            .nunique()
        )

        top = table.head(n)

        values = []

        for row in top.itertuples(
            index=False,
        ):

            values.append(
                f"{row.autonomous_community} "
                f"(n = {int(row.Total):,}; "
                f"{row.Percentage:.2f}%)"
            )

        if len(values) == 1:
            communities = values[0]
        else:
            communities = (
                ", ".join(values[:-1])
                + " and "
                + values[-1]
            )

        most_represented = top.iloc[0]

        community = (
            most_represented[
                "autonomous_community"
            ]
        )

        subtype_counts = []

        if "A" in table.columns:
            subtype_counts.append(
                f"{int(most_represented['A']):,} RSV-A"
            )

        if "B" in table.columns:
            subtype_counts.append(
                f"{int(most_represented['B']):,} RSV-B"
            )

        if len(subtype_counts) == 2:
            subtype_text = (
                f" The dataset includes "
                f"{subtype_counts[0]} and "
                f"{subtype_counts[1]} sequences "
                f"from this Community."
            )
        else:
            subtype_text = ""
        
        if n_communities == 1:
            return (
                f"The dataset includes sequences from "
                f"1 Community: "
                f"{community} "
                f"(n = {int(most_represented['Total']):,}; "
                f"{most_represented['Percentage']:.2f}%). "
                f"{subtype_text}"
            )
        
        if len(subtype_counts) == 2:
            subtype_text = (
                f" {community}, the most represented Community, "
                f"included {subtype_counts[0]} and "
                f"{subtype_counts[1]} sequences."
            )
        else:
            subtype_text = ""

        return (
            f"The dataset includes sequences from "
            f"{n_communities} Communities. "
            f"The largest contributions were from "
            f"{communities}. "
            f"{subtype_text}"
        )

    def sampling_strategy_text(
        self,
    ):
        """
        Generate descriptive text for the sampling strategy.
        """
        table = self.report.tables[
            "sampling_strategy_community"
        ].copy()
        sentinel = int(
            table["Sentinel"].sum()
        )
        non_sentinel = int(
            table["Non-sentinel"].sum()
        )
        total = int(
            table["Total"].sum()
        )
        sentinel_percentage = (
            sentinel
            / total
            * 100
        )
        non_sentinel_percentage = (
            non_sentinel
            / total
            * 100
        )
        return (
            f"Of the {total:,} sequences included in the dataset, "
            f"{sentinel:,} ({sentinel_percentage:.2f}%) corresponded "
            f"to sentinel surveillance and "
            f"{non_sentinel:,} ({non_sentinel_percentage:.2f}%) "
            f"to non-sentinel sampling."
        )

    def laboratory_contribution_text(
        self,
        n=3,
    ):
        """
        Generate descriptive text for laboratory contribution.
        """

        originating = (
            self.report.tables[
                "laboratory_originating"
            ]
            .copy()
            .sort_values(
                "Total",
                ascending=False,
            )
        )
        submitter = (
            self.report.tables[
                "laboratory_submitter"
            ]
            .copy()
            .sort_values(
                "Total",
                ascending=False,
            )
        )

        def top_labs(
            table,
            column,
        ):
            """
            Format the most represented laboratories.
            """
            values = []

            for row in table.head(n).itertuples(
                index=False,
            ):

                lab = getattr(
                    row,
                    column,
                )
                values.append(
                    f"{lab} "
                    f"(n = {int(row.Total):,}; "
                    f"{row.Percentage:.2f}%)"
                )

            if len(values) == 1:
                return values[0]
            if len(values) == 2:
                return (
                    values[0]
                    + " and "
                    + values[1]
                )

            return (
                ", ".join(values[:-1])
                + " and "
                + values[-1]
            )

        originating_labs = top_labs(
            originating,
            "originating_lab",
        )
        submitter_labs = top_labs(
            submitter,
            "submitter_lab_name",
        )
        n_originating = len(
            originating
        )
        n_submitter = len(
            submitter
        )

        return (
            f"The dataset includes sequences from "
            f"{n_originating} originating laboratories and "
            f"{n_submitter} submitter laboratories. "
            f"The largest contributions by originating laboratory "
            f"were {originating_labs}. "
            f"The most represented submitter laboratories were "
            f"{submitter_labs}."
        )

    def sequencing_characteristics_text(
        self,
        n=3,
    ):
        """
        Generate descriptive text for sequencing
        technologies and assembly methods.
        """

        technology = (
            self.report.tables[
                "sequencing_technology"
            ]
            .copy()
            .sort_values(
                "Total",
                ascending=False,
            )
        )

        assembly = (
            self.report.tables[
                "sequencing_assembly"
            ]
            .copy()
            .sort_values(
                "Total",
                ascending=False,
            )
        )

        def top_categories(
            table,
            column,
        ):
            """
            Format the most represented categories.
            """

            values = []

            for row in table.head(n).itertuples(
                index=False,
            ):

                category = getattr(
                    row,
                    column,
                )

                values.append(
                    f"{category} "
                    f"(n = {int(row.Total):,}; "
                    f"{row.Percentage:.2f}%)"
                )

            if len(values) == 1:
                return values[0]

            if len(values) == 2:
                return (
                    values[0]
                    + " and "
                    + values[1]
                )

            return (
                ", ".join(values[:-1])
                + " and "
                + values[-1]
            )

        top_technologies = top_categories(
            technology,
            "sequencing_technology",
        )

        top_assemblies = top_categories(
            assembly,
            "assembly_method",
        )

        assembly_available = int(
            assembly["Total"].sum()
        )

        return (
            f"The most frequently reported sequencing technologies "
            f"were {top_technologies}. "
            f"Assembly method was available for "
            f"{assembly_available:,} sequences, with "
            f"{top_assemblies} as the most frequently "
            f"reported methods."
        )

    def virus_lineage_composition_text(
        self,
        n=3,
    ):
        """
        Generate descriptive text for RSV lineage composition.
        """

        paragraphs = []

        for subtype in ("A", "B"):
            key = f"virus_lineage_{subtype}"

            if key not in self.report.tables:
                continue

            table = (
                self.report.tables[key]
                .copy()
                .sort_values(
                    "Total",
                    ascending=False,
                )
            )

            if table.empty:
                continue

            values = []

            for row in table.head(n).itertuples(
                index=False,
            ):
                values.append(
                    f"{row.Lineage} "
                    f"(n = {int(row.Total):,}; "
                    f"{row.Percentage:.2f}%)"
                )

            if len(values) == 1:
                top_lineages = values[0]
            elif len(values) == 2:
                top_lineages = (
                    values[0]
                    + " and "
                    + values[1]
                )
            else:
                top_lineages = (
                    ", ".join(values[:-1])
                    + " and "
                    + values[-1]
                )

            n_lineages = len(table)

            paragraphs.append(
                f"A total of {n_lineages} RSV-{subtype} lineages "
                f"were represented in the dataset. "
                f"The most frequently represented RSV-{subtype} "
                f"lineages were {top_lineages}."
            )

        return " ".join(
            paragraphs
        )

    def mutation_burden_text(
        self,
    ):
        """
        Generate descriptive text for overall
        amino acid substitution burden.
        """

        table = self.report.tables[
            "mutation_summary"
        ]

        def metric(
            name,
            column="All",
        ):
            """
            Get a metric from the mutation summary table.
            """

            if column not in table.columns:
                return None

            values = table.loc[
                table["Metric"] == name,
                column,
            ]

            if values.empty:
                return None

            return values.iloc[0]

        sequences = metric(
            "Sequences"
        )
        with_mutations = metric(
            "Sequences with mutations"
        )
        genomes_percentage = metric(
            "Genomes with mutations (%)"
        )
        median_all = metric(
            "Median substitutions / genome"
        )

        subtype_medians = []

        for subtype in ("A", "B"):
            median = metric(
                "Median substitutions / genome",
                f"RSV-{subtype}",
            )

            if median is not None:
                subtype_medians.append(
                    f"{float(median):g} for RSV-{subtype}"
                )

        text = (
            f"Amino acid substitutions were identified in "
            f"{int(with_mutations):,} of {int(sequences):,} genomes "
            f"({float(genomes_percentage):.2f}%). "
            f"The median number of substitutions per genome was "
            f"{float(median_all):g}"
        )

        if len(subtype_medians) == 2:
            text += (
                f", with subtype-specific medians of "
                f"{subtype_medians[0]} and "
                f"{subtype_medians[1]}"
            )

        return text + "."

    def protein_level_variation_text(
        self,
    ):
        """
        Generate descriptive text for protein-level
        variation.
        """

        table = self.report.tables[
            "protein_summary"
        ].copy()

        paragraphs = []

        available_subtypes = (
            table["Subtype"]
            .dropna()
            .unique()
        )

        for subtype in available_subtypes:
            subtype_table = table[
                table["Subtype"] == subtype
            ]

            if subtype_table.empty:
                continue

            top_events = (
                subtype_table
                .sort_values(
                    "Total mutation events",
                    ascending=False,
                )
                .iloc[0]
            )   

            top_substitutions = (
                subtype_table
                .sort_values(
                    "Distinct AA substitutions",
                    ascending=False,
                )
                .iloc[0]
            )

            paragraphs.append(
                f"For RSV-{subtype}, the highest number of mutation "
                f"events was observed in the "
                f"{top_events['Protein']} protein "
                f"(n = {int(top_events['Total mutation events']):,}), "
                f"while the largest number of distinct amino acid "
                f"substitutions was observed in the "
                f"{top_substitutions['Protein']} protein "
                f"(n = "
                f"{int(top_substitutions['Distinct AA substitutions']):,})."
            )

        if (
            "A" in available_subtypes
            and "B" in available_subtypes
        ):
            comparison = (
                table.pivot(
                    index="Protein",
                    columns="Subtype",
                    values="Distinct AA substitutions",
                )
                .dropna(
                    subset=["A", "B"]
                )
            )

            if not comparison.empty:
                comparison["Absolute difference"] = (
                    comparison["B"]
                    - comparison["A"]
                ).abs()

                protein_difference = (
                    comparison[
                        "Absolute difference"
                    ]
                    .idxmax()
                )

                rsv_a_difference = int(
                    comparison.loc[
                        protein_difference,
                        "A",
                    ]
                )
                rsv_b_difference = int(
                    comparison.loc[
                        protein_difference,
                        "B",
                    ]
                )

                absolute_difference = int(
                    comparison.loc[
                        protein_difference,
                        "Absolute difference",
                    ]
                )

                paragraphs.append(
                    f"The largest absolute difference in distinct "
                    f"amino acid substitutions between subtypes was "
                    f"observed for the {protein_difference} protein, "
                    f"with {rsv_a_difference:,} substitutions in RSV-A "
                    f"and {rsv_b_difference:,} in RSV-B "
                    f"(absolute difference: "
                    f"{absolute_difference:,})."
                )

            return " ".join(
                paragraphs
            )

    def protein_lineage_comparative_text(
        self,
    ):
        """
        Generate descriptive text for comparative
        protein variation across RSV lineages.
        """

        table = self.report.tables[
            "protein_lineage_summary"
        ].copy()
        
        paragraphs = []

        for subtype in ("A", "B"):
            subtype_table = table[
                table["Subtype"].eq(subtype)
            ].copy()

            if subtype_table.empty:
                continue
            
            n_lineages = subtype_table[
                "Lineage"
            ].nunique()

            if n_lineages < 2:
                continue

            matrix = subtype_table.pivot(
                index="Protein",
                columns="Lineage",
                values="Distinct AA substitutions",
            )

            percentages = (
                matrix
                .div(
                    matrix.sum(axis=0),
                    axis=1,
                )
                * 100
            )
            median = percentages.median(
                axis=1
            )
            mad = (
                percentages
                .sub(
                    median,
                    axis=0,
                )
                .abs()
                .median(
                    axis=1
                )
            )
            robust_z = (
                percentages
                .sub(
                    median,
                    axis=0,
                )
                .div(
                    mad.replace(
                        0,
                        float("nan"),
                    ),
                    axis=0,
                )
                * 0.6745
            )
            robust_z = robust_z.fillna(0)
            values = (
                robust_z
                .rename_axis(
                    index="Protein",
                    columns="Lineage",
                )
                .stack()
                .rename("Robust z-score")
                .reset_index()
            )
            highest = values.loc[
                values[
                    "Robust z-score"
                ].idxmax()
            ]
            lowest = values.loc[
                values[
                    "Robust z-score"
                ].idxmin()
            ]

            paragraphs.append(
                f"Across {n_lineages} RSV-{subtype} lineages, "
                f"the largest positive deviation from the typical "
                f"protein distribution was observed for "
                f"{highest['Protein']} in lineage "
                f"{highest['Lineage']} "
                f"(robust z-score = "
                f"{highest['Robust z-score']:.2f}), while the largest "
                f"negative deviation was observed for "
                f"{lowest['Protein']} in lineage "
                f"{lowest['Lineage']} "
                f"(robust z-score = "
                f"{lowest['Robust z-score']:.2f})."
            )

        return " ".join(
            paragraphs
        )

    def clustering_summary_text(
        self,
    ):
        """
        Generate descriptive clustering summary text
        for the 0-SNP and 5-SNP thresholds.
        """

        table = self.report.tables[
            "transmission_cluster_summary"
        ].copy()

        sentences = []

        for threshold in ("0 SNP", "5 SNP"):

            results = []

            for subtype in ("A", "B"):

                subtype_table = table[
                    table["Subtype"].eq(
                        subtype
                    )
                ]

                if subtype_table.empty:
                    continue

                clusters = subtype_table.loc[
                    subtype_table["Metric"].eq(
                        "Transmission clusters"
                    ),
                    threshold,
                ]

                singletons = subtype_table.loc[
                    subtype_table["Metric"].eq(
                        "Singleton sequences"
                    ),
                    threshold,
                ]

                if (
                    clusters.empty
                    or singletons.empty
                ):
                    continue

                clusters = int(
                    clusters.iloc[0]
                )

                singletons = int(
                    singletons.iloc[0]
                )

                results.append(
                    f"RSV-{subtype} had "
                    f"{clusters:,} genomic clusters and "
                    f"{singletons:,} singleton sequences"
                )

            if results:

                sentences.append(
                    f"At the {threshold} threshold, "
                    + " and ".join(results)
                    + "."
                )

        return " ".join(
            sentences
        )


    def genomic_connectivity_text(
        self,
    ):
        """
        Generate descriptive text for genomic connectivity.
        """

        table = self.report.tables[
            "cluster_connectivity_summary"
        ].copy()

        group_label = self.connectivity_group_label()

        group_label_plural = (
            self.connectivity_group_label(
                plural=True
            )
        )

        paragraphs = [
            (
                f"Genomic connectivity was evaluated using "
                f"{group_label} as the grouping variable. "
                f"A multi-group cluster is defined as a genomic cluster "
                f"containing sequences from at least two different "
                f"{group_label_plural} groups, indicating genomic connectivity "
                f"between groups. A single-group cluster contains "
                f"sequences from only one {group_label} group."
            )
        ]

        for threshold in (0, 5):

            threshold_table = table[
                table["Threshold"].eq(
                    threshold
                )
            ].copy()


            if threshold_table.empty:
                continue

            for subtype in ("A", "B"):

                subtype_table = threshold_table[
                    threshold_table[
                        Columns.SUBTYPE
                    ].eq(subtype)
                ].copy()

                if subtype_table.empty:
                    continue

                multi_group = subtype_table[
                    subtype_table["Connectivity"].eq(
                        "MULTI_GROUP"
                    )
                ].copy()

                single_group = subtype_table[
                    subtype_table["Connectivity"].eq(
                        "SINGLE_GROUP"
                    )
                ].copy()

                n_multi = len(
                    multi_group
                )
                n_single = len(
                    single_group
                )

                sentence = (
                    f"At the {threshold}-SNP threshold, RSV-{subtype} "
                    f"had {n_multi:,} multi-group genomic clusters and "
                    f"{n_single:,} single-group genomic clusters."
                )

                if not multi_group.empty:

                    largest = multi_group.loc[
                        multi_group[
                            "n_sequences"
                        ].idxmax()
                    ]

                    broadest = multi_group.loc[
                        multi_group[
                            "n_groups"
                        ].idxmax()
                    ]

                    sentence += (
                        f" The largest RSV-{subtype} multi-group cluster "
                        f"was {largest['Cluster']}, comprising "
                        f"{int(largest['n_sequences']):,} sequences "
                        f"from {int(largest['n_groups']):,} "
                        f"{group_label_plural}."
                    )

                    sentence += (
                        f" The broadest RSV-{subtype} multi-group cluster "
                        f"was {broadest['Cluster']}, including sequences "
                        f"from {int(broadest['n_groups']):,} "
                        f"{group_label_plural} "
                        f"(n = {int(broadest['n_sequences']):,} sequences)."
                    )

                paragraphs.append(
                    sentence
                )

        return " ".join(
            paragraphs
        )


    def genomic_cluster_composition_text(
        self,
    ):
        """
        Generate descriptive text for genomic cluster composition.
        """

        group_label = self.connectivity_group_label()
        

        return (
            f"Composition of multi-group genomic clusters according "
            f"to {group_label}. Clusters are ordered by decreasing "
            f"number of sequences. An asterisk (*) indicates clusters "
            f"detected in two or more epidemiological seasons."
        )


    def genomic_cluster_timeline_text(
        self,
    ):
        """
        Generate descriptive text for genomic cluster timelines.
        """

        group_label = (
            self.connectivity_group_label()
        )

        group_label_plural = self.connectivity_group_label(
            plural=True
        )

        return (
            f"Temporal distribution of sequences belonging to "
            f"multi-group genomic clusters. Colors indicate the "
            f"{group_label} associated with each sequence. "
            f"Shaded areas indicate epidemiological seasons. "
            f"Cluster lineage is indicated next to the cluster "
            f"identifier, and an asterisk (*) indicates clusters "
            f"detected in two or more epidemiological seasons."
        )

    def genomic_connectivity_caption(
        self,
        plot_type,
        threshold,
        subtype=None,
    ):
        """
        Generate captions for genomic connectivity figures.
        """

        group_label = self.connectivity_group_label()

        group_label_plural = (
            self.connectivity_group_label(
                plural=True
            )
        )

        if plot_type == "spread":
            return (
                f"Genomic connectivity of multi-group clusters "
                f"across {group_label_plural} at the "
                f"{threshold}-SNP threshold."
            )

        if plot_type == "composition":
            return (
                f"{group_label} composition of multi-group "
                f"RSV-{subtype} genomic clusters at the "
                f"{threshold}-SNP threshold. An asterisk (*) "
                f"indicates clusters detected in two or more "
                f"epidemiological seasons."
            )

        if plot_type == "timeline":
            return (
                f"Temporal distribution of multi-group RSV-{subtype} "
                f"genomic clusters at the {threshold}-SNP threshold. "
                f"Colors indicate {group_label}, shaded areas indicate "
                f"epidemiological seasons, and an asterisk (*) indicates "
                f"clusters detected in two or more epidemiological seasons."
            )

        return None

    def summary(
        self,
    ):
        """
        Dataset summary.
        """

        # Main dataset metrics
        total = self.dataset_metric(
            "Total sequences"
        )
        rsv_a = self.dataset_metric(
            "RSV-A"
        )
        rsv_b = self.dataset_metric(
            "RSV-B"
        )
        collection_period = self.dataset_metric(
            "Collection period"
        )
        # Secondary dataset metrics
        communities = self.dataset_metric(
            "Communities"
        )
        lineages = self.dataset_metric(
            "Lineages"
        )
        institutions = self.dataset_metric(
            "Institutions"
        )
        submitter_labs = self.dataset_metric(
            "Submitter laboratories"
        )
        technologies = self.dataset_metric(
            "Sequencing technologies"
        )
        assembly_methods = self.dataset_metric(
            "Assembly methods"
        )
        sequence_length = self.dataset_metric(
            "Sequence length"
        )

        # Format sequence counts
        total = f"{int(total):,}"
        rsv_a = f"{int(rsv_a):,}"
        rsv_b = f"{int(rsv_b):,}"

        #Format collection period
        start_date, end_date = (
            collection_period
            .split(" - ")
        )
        start_date = pd.to_datetime(
            start_date
        )
        end_date = pd.to_datetime(
            end_date
        )
        collection_period = (
            f"{start_date:%m/%Y}"
            " - "
            f"{end_date:%m/%Y}"
        )

        collection_period_text = (
            f"{self.format_month_year_english(start_date)}"
            " and "
            f"{self.format_month_year_english(end_date)}"
        )

        return f"""
    <section class="report-section">

        <h2>Dataset overview</h2>

        {self.report_text(
            (
                f"This report includes GISAID metadata for "
                f"{total} RSV sequences collected between "
                f"{collection_period_text}. "
                f"The dataset comprises {rsv_a} RSV-A and "
                f"{rsv_b} RSV-B sequences from "
                f"{communities} Autonomous Communities and/or Cities in Spain."
            )
        )}

        {self.report_text(
            self.metadata_completeness_text()
        )}

        <p class="section-description">
            Summary of the genomic surveillance dataset included
            in this report.
        </p>

        <div class="metric-grid">
            
            {self.metric_card(
                value=total,
                label="Sequences",
            )}
            {self.metric_card(
                value=rsv_a,
                label="RSV-A",
                css_class="rsv_a",
            )}
            {self.metric_card(
                value=rsv_b,
                label="RSV-B",
                css_class="rsv_b",
            )}
            {self.metric_card(
                value=collection_period,
                label="Collection period",
            )}

        </div>

        <div class="dataset-characteristics">
            <h3>Dataset characteristics</h3>

            <div class="secondary-metric-grid">
                
                {self.secondary_metric(
                    value=communities,
                    label="Communities",
                )}
                {self.secondary_metric(
                    value=lineages,
                    label="Lineages",
                )}
                {self.secondary_metric(
                    value=institutions,
                    label="Laboratories",
                )}
                {self.secondary_metric(
                    value=submitter_labs,
                    label="Submitter Institutions",
                )}
                {self.secondary_metric(
                    value=technologies,
                    label="Sequencing technologies",
                )}
                {self.secondary_metric(
                    value=assembly_methods,
                    label="Assembly methods",
                )}
            
            </div>

            <div class="dataset-details">

                <span class="dataset-detail-label">
                    Sequence length
                </span>

                <span class="dataset-detail-value">
                    {escape(str(sequence_length))}
                </span>

            </div>

        </div>

        <div class="report-subsection">
            <h3>Metadata completeness</h3>

            <p class="subsection-description">
                Proportion of sequences with available information
                for each metadata variable.
            </p>

            <div class="completeness-container">

                {self.metadata_completeness()}
            </div>

            {self.collapsible_table(
                table_id="metadata_completeness",
                label="View data",
                interactive=True,
            )}

        </div>

        <div class="report-subsection">
            <h3>Temporal distribution</h3>

            <p class="subsection-description">
                Distribution of RSV sequences according to collection year
                and viral subtype.
            </p>

            {self.report_text(
                self.temporal_distribution_text()
            )}
            
            {self.figure(
                figure_id="temporal_year",
                section="geography",
                caption=(
                    "Annual distribution of RSV-A & RSV-B sequences "
                    "included in the dataset."
                ),
            )}

            {self.collapsible_table(
                table_id="temporal_year",
                label="View data",
            )}

        </div>

        <div class="report-subsection">
            <h3>Geographic distribution</h3>

            <p class="subsection-description">
                Distribution of RSV sequences across Communities.
            </p>

            {self.report_text(
                self.geographic_distribution_text()
            )}
            
            {self.figure_grid(
                figures=[
                    self.figure(
                        figure_id="geography_map",
                        section="geography",
                        caption=(
                            "Geographic distribution RSV-A & RSV-B sequences "
                            "across Communities."
                        ),
                    ),
                    self.figure(
                        figure_id="geography_community_barplot",
                        section="geography",
                        caption=(
                            "Number of RSV-A & RSV-B sequences by Community."
                        ),
                    ),
                ],
                columns=2,
            )}

        </div>

        <div class="report-subsection">
            <h3>Sampling strategy</h3>

            <p class="subsection-description">
                Distribution of RSV sequences according to 
                Sampling strategy.
            </p>

            {self.report_text(
                self.sampling_strategy_text()
            )}
            
            {self.figure(
                figure_id="sampling_strategy_barplot",
                section="geography",
                caption=(
                    "Distribution of RSV-A & RSV-B sequences "
                    "by sampling strategy."
                ),
                max_width=800,
            )}

            {self.collapsible_table(
                table_id="sampling_strategy_community",
                label="View data",
            )}

        </div>

    </section>
    """

    def sample_characteristics(
        self,
    ):
        """
        Patient characteristics.
        """

        patients_with_age = self.table_metric(
            "patient_age_summary",
            "Patients with age",
        )
        median_age = self.table_metric(
            "patient_age_summary",
            "Median age",
        )
        median_iqr = self.table_metric(
            "patient_age_summary",
            "Median (IQR)",
        )
        minimum_age = self.table_metric(
            "patient_age_summary",
            "Minimum",
        )
        maximum_age = self.table_metric(
            "patient_age_summary",
            "Maximum",
        )
        iqr = (
            str(median_iqr)
            .split("(", 1)[1]
            .rstrip(")")
        )
        age_range = (
            f"{float(minimum_age):g}"
            "-"
            f"{float(maximum_age):g}"
        )
        patients_with_age = (
            f"{int(patients_with_age):,}"
        )
        median_age = (
            f"{float(median_age):g}"
        )

        return f"""
    <section class="report-section">
        <h2>Patient characteristics</h2>

        <p class="section-description">
            Description of patient metadata
        </p>

        <div class="report-subsection">
            <h3>Patient characteristics</h3>

            <p class="subsection-description">
                Distribution of available patient demographic metadata
                in the dataset.
            </p>
            
            {self.figure_grid(
                figures=[
                    self.figure(
                        figure_id="patient_gender_barplot",
                        section="patient",
                        caption=(
                            "Distribution RSV-A & RSV-B sequences "
                            "by patient gender."
                        ),
                    ),
                    self.figure(
                        figure_id="patient_age_groups_barplot",
                        section="patient",
                        caption=(
                            "Distribution RSV-A & RSV-B sequences "
                            "by patient age groups"
                        ),
                    ),
                    self.figure(
                        figure_id="patient_status_barplot",
                        section="patient",
                        caption=(
                            "Distribution RSV-A & RSV-B sequences "
                            "by patient status"
                        ),
                    ),
                ],
                columns=3,
            )}

        </div>

                <div class="patient-age-summary">
            <h4>Age summary</h4>

            <div class="secondary-metric-grid age-summary-grid">
                {self.secondary_metric(
                    value=patients_with_age,
                    label="Patients with age",
                )}
                {self.secondary_metric(
                    value=f"{median_age} years",
                    label="Median age",
                )}
                {self.secondary_metric(
                    value=f"{iqr} years",
                    label="IQR",
                )}
                {self.secondary_metric(
                    value=f"{age_range} years",
                    label="Range",
                )}

            </div>

            {self.collapsible_table(
                table_id="patient_age_summary",
                label="View age summary",
                interactive=False,
            )}

        </div>

    </section>
    """

    def laboratory_contribution(
        self,
    ):
        """
        Patient characteristics.
        """

        return f"""
    <section class="report-section">
        <h2>Laboratory contribution and Sequencing characteristics</h2>

        <p class="section-description">
            Laboratory contribution and sequencing methods.
        </p>

        <div class="report-subsection">
            <h3>Laboratory contribution</h3>

            <p class="subsection-description">
                Distribution of RSV sequences according to 
                originating and subbmitting institutions.
            </p>

            {self.report_text(
                self.laboratory_contribution_text()
            )}
            
            {self.figure(
                figure_id="laboratory_originating_barplot",
                section="laboratory",
                caption=(
                    "Distribution of RSV-A & RSV-B sequences "
                    "by originating laboratory."
                ),
                max_width=700,
            )}

            {self.collapsible_table(
                table_id="laboratory_originating",
                label="View data",
                interactive=True,
            )}

            {self.figure(
                figure_id="laboratory_submitter_barplot",
                section="laboratory",
                caption=(
                    "Distribution of RSV-A & RSV-B sequences "
                    "by submitter institution."
                ),
                max_width=800,
            )}

        </div>

        <div class="report-subsection">
            <h3>Sequencing characteristics</h3>

            <p class="subsection-description">
                Distribution of sequencing technologies and 
                assembly methods in the dataset.
            </p>

            {self.report_text(
                self.sequencing_characteristics_text()
            )}
            
            {self.figure_grid(
                figures=[
                    self.figure(
                        figure_id="sequencing_technology_barplot",
                        section="sequencing",
                        caption=(
                            "Distribution RSV-A & RSV-B sequences "
                            "by sequencing technology."
                        ),
                    ),
                    self.figure(
                        figure_id="sequencing_assembly_barplot",
                        section="sequencing",
                        caption=(
                            "Distribution RSV-A & RSV-B sequences "
                            "by assembly method."
                        ),
                    ),
                ],
                columns=2,
            )}

            <div class="data-links">

                {self.collapsible_table(
                    table_id="sequencing_technology",
                    label="View sequencing technology data",
                    interactive=True,
                )}

                {self.collapsible_table(
                    table_id="sequencing_assembly",
                    label="View sequencing method data",
                    interactive=True,
                )}

            </div>

        </div>

    </section>
    """

    def virus_diversity(
        self,
    ):
        """
        Virus diversity and lineage distribution.
        """

        return f"""
    <section class="report-section">
        <h2>Virus diversity</h2>

        <div class="report-subsection">
            <h3>Lineage composition</h3>

            <p class="subsection-description">
                Distribution of RSV-A and RSV-B sequences
                according to viral lineage.
            </p>

            {self.report_text(
                self.virus_lineage_composition_text()
            )}
            
            {self.figure_grid(
                figures=[
                    self.figure(
                        figure_id="virus_lineage_A_barplot",
                        section="virus",
                        caption=(
                            "Distribution of RSV-A sequences "
                            "by lineage."
                        ),
                    ),
                    self.figure(
                        figure_id="virus_lineage_B_barplot",
                        section="virus",
                        caption=(
                            "Distribution of RSV-B sequences "
                            "by lineage."
                        ),
                    ),
                ],
                columns=2,
            )}
            
            <div class="data-links">
                {self.collapsible_table(
                    table_id="virus_lineage_A",
                    label="View RSV-A lineage data",
                    interactive=True,
                )}
                {self.collapsible_table(
                    table_id="virus_lineage_B",
                    label="View RSV-B lineage data",
                    interactive=True,
                )}
            </div>
        
        <div class="report-subsection">
            <h3>Temporal lineage distribution</h3>

            <p class="subsection-description">
                Annual relative distribution of RSV lineages.
                Up to 15 lineages are displayed individually, while
                the remaining lineages are grouped as Other.
            </p>

            {self.figure(
                figure_id="lineage_all_percentage_per_year",
                section="virus",
                caption=(
                    "Annual relative distribution of RSV lineages. "
                    "Lineages not displayed individually are grouped "
                    "as Other."
                ),
            )}

            {self.collapsible_table(
                table_id="lineage_all_per_year",
                label="View temporal lineage data",
                interactive=True,
            )}
        </div>

        <div class="report-subsection">
            <h3>Geographic lineage distribution</h3>

            <p class="subsection-description">
                Geographic distribution of major RSV lineages across
                Communities. Major lineages were defined independently
                for RSV-A & RSV-B as the most frequent lineages reaching
                at least 80% of cumulative representation.
                Remaining lineages were grouped as Other.
            </p>

            {self.figure(
                figure_id="community_lineage_heatmap_A",
                section="virus",
                caption=(
                    "Distribution of major RSV-A lineages across "
                    "Communities. Cell values indicate the number "
                    "of sequences and colour intensity "
                    "represents relative abundance within each Community."
                ),
            )}

            {self.figure(
                figure_id="community_lineage_heatmap_B",
                section="virus",
                caption=(
                    "Distribution of major RSV-B lineages across "
                    "Communities. Cell values indicate the number "
                    "of sequences and colour intensity "
                    "represents relative abundance within each Community."
                ),
            )}

            <div class="data-links">
                {self.collapsible_table(
                    table_id="major_lineages_community_heatmap_A",
                    label="View RSV-A geographic lineage data",
                    interactive=True,
                )}
                {self.collapsible_table(
                    table_id="major_lineages_community_heatmap_B",
                    label="View RSV-B geographic lineage data",
                    interactive=True,
                )}
            </div>

        </div>
    
    </section>
    """

    def protein_summary(
        self,
    ):
        """
        Virus diversity and lineage distribution.
        """

        return f"""
    <section class="report-section">
        <h2>Protein summary</h2>

        <div class="report-subsection">
            <h3>Overall mutation burden</h3>

            <p class="subsection-description">
                Distribution of amino acid substitution per sequence
                across the dataset and by RSV subtype.
            </p>

            {self.report_text(
                self.mutation_burden_text()
            )}
            
            {self.figure_grid(
                figures=[
                    self.figure(
                        figure_id="mutations_histogram",
                        section="mutations",
                        caption=(
                            "Distribution of the number of amino acid "
                            "substitutions per sequence."
                        ),
                    ),
                    self.figure(
                        figure_id="mutations_boxplot",
                        section="mutations",
                        caption=(
                            "Distribution of amino acid substitutions "
                            "per sequence by RSV subtype."
                        ),
                    ),
                ],
                columns=2,
            )}
            
            {self.collapsible_table(
                table_id="mutation_summary",
                label="View mutation summary",
                interactive=True,
            )}

        </div>

        <div class="report-subsection">
            <h3>Protein-level variation</h3>

            <p class="subsection-description">
                Distribution of amino acid substitutions across
                RSV-A & RSV-B proteins.
            </p>

            {self.report_text(
                self.protein_level_variation_text()
            )}

            {self.figure(
                figure_id="protein_summary_substitutions",
                section="mutations",
                caption=(
                    "Distribution of amino acid substitutions across "
                    "RSV proteins by subtype."
                ),
                max_width=800,
            )}

            {self.collapsible_table(
                table_id="protein_summary",
                label="View protein-level data",
                interactive=True,
            )}
        </div>

        <div class="report-subsection">
            <h3>Lineage-specific protein variation</h3>

            <p class="subsection-description">
                Distribution of distinct amino acid substitutions across
                proteins in major RSV lineages. Major lineages were defined
                independently for RSV-A and RSV-B as the most frequent
                lineages reaching at least 80% cumulative representation.
            </p>

            {self.figure(
                figure_id="protein_lineage_distinct_substitutions_heatmap_A",
                section="mutations",
                caption=(
                    "Distribution of distinct amino acid substitutions "
                    "across proteins in major RSV-A lineages. Values "
                    "within cells indicate the number of distinct amino "
                    "acid substitutions, while cell colour represents "
                    "the percentage of all distinct substitutions within "
                    "each lineage."
                ),
                max_width=800,
            )}

            {self.figure(
                figure_id="protein_lineage_distinct_substitutions_heatmap_B",
                section="mutations",
                caption=(
                    "Distribution of distinct amino acid substitutions "
                    "across proteins in major RSV-B lineages. Values "
                    "within cells indicate the number of distinct amino "
                    "acid substitutions, while cell colour represents "
                    "the percentage of all distinct substitutions within "
                    "each lineage."
                ),
                max_width=800,
            )}

            <div class="data-links">
                {self.collapsible_table(
                    table_id="protein_lineage_distinct_substitutions_heatmap_A",
                    label="View RSV-A heatmap data",
                    interactive=True,
                )}
                {self.collapsible_table(
                    table_id="protein_lineage_distinct_substitutions_heatmap_B",
                    label="View RSV-B heatmap data",
                    interactive=True,
                )}
            </div>

        <div class="report-subsection">
            <h3>Comparative protein variation across lineages</h3>

            <p class="subsection-description">
                Relative protein-level amino acid substitution patterns
                across RSV lineages. For each protein, robust z-scores
                describe deviations from its typical relative distribution
                across lineages.
            </p>

            {self.report_text(
                self.protein_lineage_comparative_text()
            )}

            {self.figure(
                figure_id="protein_lineage_comparative_heatmap_A",
                section="mutations",
                caption=(
                    "Relative protein-level amino acid substitution "
                    "patterns across RSV-A lineages. Positive robust "
                    "z-scores indicate a higher relative contribution "
                    "of a protein compared with its typical distribution "
                    "across lineages, while negative values indicate "
                    "a lower relative contribution."
                ),
                max_width=800,
            )}

            {self.figure(
                figure_id="protein_lineage_comparative_heatmap_B",
                section="mutations",
                caption=(
                    "Relative protein-level amino acid substitution "
                    "patterns across RSV-B lineages. Positive robust "
                    "z-scores indicate a higher relative contribution "
                    "of a protein compared with its typical distribution "
                    "across lineages, while negative values indicate "
                    "a lower relative contribution."
                ),
                max_width=800,
            )}
                {self.collapsible_table(
                    table_id="protein_lineage_summary",
                    label="View protein-lineage data",
                    interactive=True,
                )}
            </div>

        <div class="report-subsection">
            <h3>Positional protein variation</h3>

            <p class="subsection-description">
                Positional distribution of amino acid variation along the F, G,
                and L proteins for the most represented lineage of each RSV
                subtype. Lollipop height represents the number of distinct amino
                acid substitutions observed at each position, while marker size
                represents the frequency of the most common substitution whithin
                the lineage. Only positions where the most common substitution
                reach a frequency of at least 10% are shown. Annotated 
                protein features are displayed along the protein track.
            </p>

            <h4>RSV-A</h4>
            {self.figure(
                figure_id="lollipop_A_F",
                section="mutations",
                caption=(
                    "Amino acid variation across the F protein."
                ),
            )}
            {self.figure(
                figure_id="lollipop_A_G",
                section="mutations",
                caption=(
                    "Amino acid variation across the G protein."
                ),
            )}
            {self.figure(
                figure_id="lollipop_A_L",
                section="mutations",
                caption=(
                    "Amino acid variation across the L protein."
                ),
            )}

            <h4>RSV-B</h4>
            {self.figure(
                figure_id="lollipop_B_F",
                section="mutations",
                caption=(
                    "Amino acid variation across the F protein."
                ),
            )}
            {self.figure(
                figure_id="lollipop_B_G",
                section="mutations",
                caption=(
                    "Amino acid variation across the G protein."
                ),
            )}
            {self.figure(
                figure_id="lollipop_B_L",
                section="mutations",
                caption=(
                    "Amino acid variation across the L protein."
                ),
            )}
        </div>

    </section>
    """

    def clustering(
        self,
    ):
        """
        Build clustering analysis section.
        """

        html = f"""
    <section class="report-section">
        <h2>Clustering analysis</h2>

        {self.report_text(
            (
                "Genomic clustering was evaluated separately for "
                "RSV-A and RSV-B using pairwise SNP distances. "
                "Sequences not connected to any other "
                "sequence at a given threshold were classified as singletons. "
                "Results were evaluated across multiple "
                "SNP distance thresholds, with selected 0-SNP and "
                "5-SNP thresholds shown in the figures below."
            )
        )}

        <div class="report-subsection">
            <h3>Clustering summary</h3>

            <p class="subsection-description">
                Summary of genomic clustering across the evaluated
                SNP distance thresholds for RSV-A and RSV-B.
            </p>

            {self.report_text(
                self.clustering_summary_text()
            )}

            {self.collapsible_table(
                table_id="transmission_cluster_summary",
                label="View clustering summary",
                interactive=True,
            )}

        </div>

        <div class="report-subsection">
            <h3>Clustering geographic distribution</h3>

            <p class="subsection-description">
                Distribution of clustered and singleton sequences across
                Communities at 0-SNP and 5-SNP distance thresholds.
                Bars represent the percentage of sequences
                classified as clustered or singleton within each community,
                while numbers indicate the total number of sequences from
                that community.
            </p>

            <h4>0-SNP threshold</h4>

            {self.figure_grid(
                figures=[
                    self.figure(
                        figure_id="transmission_cluster_by_community_0_snps_A",
                        section="clusters",
                        caption=(
                            "Distribution of clustered and singleton "
                            "RSV-A sequences across Communities "
                            "at 0-SNP threshold."
                        ),
                    ),
                    self.figure(
                        figure_id="transmission_cluster_by_community_0_snps_B",
                        section="clusters",
                        caption=(
                            "Distribution of clustered and singleton "
                            "RSV-B sequences across Communities "
                            "at 0-SNP threshold."
                        ),
                    ),
                ],
                columns=2,
            )}

            <h4>5-SNP threshold</h4>

            {self.figure_grid(
                figures=[
                    self.figure(
                        figure_id="transmission_cluster_by_community_5_snps_A",
                        section="clusters",
                        caption=(
                            "Distribution of clustered and singleton "
                            "RSV-A sequences across Communities "
                            "at 5-SNP threshold."
                        ),
                    ),
                    self.figure(
                        figure_id="transmission_cluster_by_community_5_snps_B",
                        section="clusters",
                        caption=(
                            "Distribution of clustered and singleton "
                            "RSV-B sequences across Communities "
                            "at 5-SNP threshold."
                        ),
                    ),
                ],
                columns=2,
            )}

            {self.collapsible_table(
                table_id="transmission_cluster_by_community_summary",
                label="View geographic clustering data",
                interactive=True,
            )}

        </div>

        <div class="report-subsection">
            <h3>Clustering temporal distribution</h3>

            <p class="subsection-description">
                Distribution of clustered and singleton sequences by
                Year at 0-SNP and 5-SNP distance thresholds.
                Bars represent the percentage of sequences
                classified as clustered or singleton within each year,
                while numbers above the bars indicate the total number of sequences 
                collected that year.
            </p>

            <h4>0-SNP threshold</h4>

            {self.figure_grid(
                figures=[
                    self.figure(
                        figure_id="transmission_cluster_by_year_0_snps_A",
                        section="clusters",
                        caption=(
                            "Temporal distribution of clustered and singleton "
                            "RSV-A sequences at 0-SNP threshold."
                        ),
                    ),
                    self.figure(
                        figure_id="transmission_cluster_by_year_0_snps_B",
                        section="clusters",
                        caption=(
                            "Temporal distribution of clustered and singleton "
                            "RSV-B sequences at 0-SNP threshold."
                        ),
                    ),
                ],
                columns=2,
            )}

            <h4>5-SNP threshold</h4>

            {self.figure_grid(
                figures=[
                    self.figure(
                        figure_id="transmission_cluster_by_year_5_snps_A",
                        section="clusters",
                        caption=(
                            "Temporal distribution of clustered and singleton "
                            "RSV-A sequences at 5-SNP threshold."
                        ),
                    ),
                    self.figure(
                        figure_id="transmission_cluster_by_year_5_snps_B",
                        section="clusters",
                        caption=(
                            "Temporal distribution of clustered and singleton "
                            "RSV-B sequences at 5-SNP threshold."
                        ),
                    ),
                ],
                columns=2,
            )}

            {self.collapsible_table(
                table_id="transmission_cluster_by_year_summary",
                label="View temporal clustering data",
                interactive=True,
            )}

        </div>

        <div class="report-subsection">
            <h3>Clustering lineage distribution</h3>

            <p class="subsection-description">
                Distribution of clustered sequences across RSV major lineages
                and collection years at the 0-SNP and 5-SNP distance
                thresholds. Values within cells indicate the number of
                clustered sequences for each lineage and year, while cell
                colour represents the percentage of all clustered sequences
                from that year assigned to the corresponding lineage.
                The total number of sequences for each lineage, including
                both clustered and singleton sequences, is shown in
                parentheses next to the lineage name.
            </p>

            <h4>0-SNP threshold</h4>

            {self.figure_grid(
                figures=[
                    self.figure(
                        figure_id="transmission_lineage_heatmap_A_0_snps",
                        section="clusters",
                        caption=(
                            "Distribution of clustered RSV-A sequences "
                            "across lineages and collection years at "
                            "0-SNP threshold."
                        ),
                    ),
                    self.figure(
                        figure_id="transmission_lineage_heatmap_B_0_snps",
                        section="clusters",
                        caption=(
                            "Distribution of clustered RSV-B sequences "
                            "across lineages and collection years at "
                            "0-SNP threshold."
                        ),
                    ),
                ],
                columns=2,
            )}

            <h4>5-SNP threshold</h4>

            {self.figure_grid(
                figures=[
                    self.figure(
                        figure_id="transmission_lineage_heatmap_A_5_snps",
                        section="clusters",
                        caption=(
                            "Distribution of clustered RSV-A sequences "
                            "across lineages and collection years at "
                            "5-SNP threshold."
                        ),
                    ),
                    self.figure(
                        figure_id="transmission_lineage_heatmap_B_5_snps",
                        section="clusters",
                        caption=(
                            "Distribution of clustered RSV-B sequences "
                            "across lineages and collection years at "
                            "5-SNP threshold."
                        ),
                    ),
                ],
                columns=2,
            )}

            {self.collapsible_table(
                table_id="transmission_cluster_by_lineage_summary",
                label="View lineage clustering data",
                interactive=True,
            )}

        </div>

        <div class="report-subsection">
            <h3>Genomic connectivity</h3>

            <p class="subsection-description">
            {self.report_text(
                self.genomic_connectivity_text()
            )}
            </p>

            <h4>Cluster connectivity</h4>

            <p class="subsection-description">
                Relationship between genomic cluster size and the number of
                represented groups. Only clusters detected in two or more 
                groups are shown.

            </p>

            {self.figure_grid(
                figures=[
                    self.figure(
                        figure_id="genomic_cluster_spread_0snp",
                        section="genomic_connectivity",
                        caption=self.genomic_connectivity_caption(
                            plot_type="spread",
                            threshold=0,
                        ),
                    ),
                    self.figure(
                        figure_id="genomic_cluster_spread_5snp",
                        section="genomic_connectivity",
                        caption=self.genomic_connectivity_caption(
                            plot_type="spread",
                            threshold=5,
                        ),
                    ),
                ],
                columns=2,
            )}

            {self.collapsible_table(
                table_id="cluster_connectivity_summary",
                label="View clustering connectivity summary",
                interactive=True,
            )}

            <h4>Cluster composition</h4>

            <p class="subsection-description">
            {self.report_text(
                self.genomic_cluster_composition_text()
            )}
            </p>

            <h5>0-SNP threshold</h5>

            {self.figure_grid(
                figures=[
                    self.figure(
                        figure_id="genomic_cluster_composition_0snp_A",
                        section="genomic_connectivity",
                        caption=self.genomic_connectivity_caption(
                            plot_type="composition",
                            threshold=0,
                            subtype="A",
                        ),
                    ),
                    self.figure(
                        figure_id="genomic_cluster_composition_0snp_B",
                        section="genomic_connectivity",
                        caption=self.genomic_connectivity_caption(
                            plot_type="composition",
                            threshold=0,
                            subtype="B",
                        ),
                    ),
                ],
                columns=2,
            )}

            <h5>5-SNP threshold</h5>

            {self.figure_grid(
                figures=[
                    self.figure(
                        figure_id="genomic_cluster_composition_5snp_A",
                        section="genomic_connectivity",
                        caption=self.genomic_connectivity_caption(
                            plot_type="composition",
                            threshold=5,
                            subtype="A",
                        ),
                    ),
                    self.figure(
                        figure_id="genomic_cluster_composition_5snp_B",
                        section="genomic_connectivity",
                        caption=self.genomic_connectivity_caption(
                            plot_type="composition",
                            threshold=5,
                            subtype="B",
                        ),
                    ),
                ],
                columns=2,
            )}

            <h4>Temporal distribution of connected clusters</h4>

            <p class="subsection-description">
            {self.report_text(
                self.genomic_cluster_timeline_text()
            )}
            </p>

            <h5>0-SNP threshold</h5>

            {self.figure_grid(
                figures=[
                    self.figure(
                        figure_id="genomic_cluster_timeline_0snp_A",
                        section="genomic_connectivity",
                        caption=self.genomic_connectivity_caption(
                            plot_type="timeline",
                            threshold=0,
                            subtype="A",
                        ),
                    ),
                    self.figure(
                        figure_id="genomic_cluster_timeline_0snp_B",
                        section="genomic_connectivity",
                        caption=self.genomic_connectivity_caption(
                            plot_type="timeline",
                            threshold=0,
                            subtype="B",
                        ),
                    ),
                ],
                columns=2,
            )}

            <h5>5-SNP threshold</h5>

            {self.figure_grid(
                figures=[
                    self.figure(
                        figure_id="genomic_cluster_timeline_5snp_A",
                        section="genomic_connectivity",
                        caption=self.genomic_connectivity_caption(
                            plot_type="timeline",
                            threshold=5,
                            subtype="A",
                        ),
                    ),
                    self.figure(
                        figure_id="genomic_cluster_timeline_5snp_B",
                        section="genomic_connectivity",
                        caption=self.genomic_connectivity_caption(
                            plot_type="timeline",
                            threshold=5,
                            subtype="B",
                        ),
                    ),
                ],
                columns=2,
            )}

            <div class="data-links">

                {self.collapsible_table(
                    table_id="cluster_connectivity_samples",
                    label="View cluster connectivity sample data",
                    interactive=True,
                )}


                {self.collapsible_table(
                    table_id="group_connectivity",
                    label="View group connectivity data",
                    interactive=True,
                )}

            </div>

        </div>

    </section>
    """
        return html

    def phylogeny(
        self,
    ):
        """
        Build phylogenetic analysis section.
        """

        figures = []

        for path in self.report.phylogenies:

            figure_path = (
                Path("figures")
                / "phylogeny"
                / path.name
            )

            figures.append(
                f"""
        
    <figure class="figure-container">
        <img
            src="{figure_path.as_posix()}"
            alt="{escape(path.stem)}"
            loading="lazy"
        >
        <figcaption class="figure-caption">
            Phylogenetic representation of the sequences
            included in the analysis.
        </figcaption>

    </figure>
    """
            )

        figures_html = "\n".join(
            figures
        )

        return f"""
    <section class="report-section">
        <h2>Phylogenetic analysis</h2>

        <p class="section-description">
            Phylogenetic representations provided for the sequences
            included in the analysis.
        </p>

        {figures_html}

    </section>
    """

    def analysis_information(
        self,
    ):
        """
        Render report generation and analysis information.
        """

        generated = self.format_date_english(
            self.report.generated_at,
            include_time=True,
        )
        input_file = (
            self.report.input_file.name
            if self.report.input_file is not None
            else "NA"
        )
        version = "1.1.0"

        return f"""
    <footer class="report-footer">
        <div class="analysis-information">
            
            <h2>Analysis information</h2>
            <div class="analysis-columns">

                <div class="analysis-column">

                    <h3>Report details</h3>

                    <div class="analysis-info-grid">

                        <div class="analysis-info-label">
                            Input dataset:
                        </div>
                        <div>
                            {escape(input_file)}
                        </div>

                        <div class="analysis-info-label">
                            Report date generated:
                        </div>
                        <div>
                            {escape(generated)}
                        </div>

                        <div class="analysis-info-label">
                            GenomeReport version:
                        </div>
                        <div>
                            {escape(version)}
                        </div>
                    </div>
                </div>

                <div class="analysis-column">

                    <h3>Analysis parameters</h3>

                    <div class="analysis-info-grid">

                        <div class="analysis-info-label">
                            Major lineages:
                        </div>
                        <div>
                            ≥80% cumulative representation
                        </div>

                        <div class="analysis-info-label">
                            Mutations filter (lollipop plots):
                        </div>
                        <div>
                            ≥10% frequency
                        </div>

                        <div class="analysis-info-label">
                            Clustering thresholds:
                        </div>
                        <div>
                            0, 5, 10, 15, 20 SNPs
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <div class="report-attribution">

            <span>
                GenomeReport
            </span>

            <span class="footer-separator">
                ·
            </span>

            <span>
                Developed by Irving Cancino-Muñoz
            </span>

            <span class="footer-separator">
                ·
            </span>

            <a
                href="TU_GITHUB"
                target="_blank"
                rel="noopener noreferrer"
            >
                https://github.com/irvink182/GenomeReport
            </a>

        </div>

    </footer>
    """
