"""
Console output utilities.
"""
import numbers
from pathlib import Path
from .style import(
    CONSOLE_TEXT_WIDTH,
    CONSOLE_MUTATION_FREQUENCY,
)
from .config import (
    SUBTYPES,
    LOLLIPOP_PROTEINS,
    Columns,
    MutationColumns,
)

from .tables import (
    get_selected_lineage,
)

import pandas as pd
import math

###########
# Constants
###########

LINE = "=" * 70
SUBLINE = "-" * 70

###########
# Messages
###########

def info(message):
    print(f"[INFO] {message}")

def warning(message):
    print(f"[WARNING] {message}")

def success(message):
    print(f"[OK] {message}")

###########
# Report
###########

def print_header(dataset, output):

    print()
    print(LINE)
    print("GenomeReport")
    print(LINE)

    print(f"Input : {dataset.filepath}")
    print(f"Output: {Path(output)}")

def print_footer():

    print()
    print(LINE)
    print("Report successfully generated.")
    print(LINE)
    print()

###########
# Classes
###########
class ReportLogger:
    """
    Print to stdout and optionally write to a text file.
    """

    def __init__(self, filename=None):
        self.file = None

        if filename is not None:
            self.file = open(
                filename,
                "w",
                encoding="utf-8",
            )

    def write(
        self,
        text="",
    ):

        print(text)

        if self.file is not None:
            self.file.write(text)
            self.file.write("\n")

    def close(self):
        if self.file is not None:
            self.file.close()

###########
# Tables
###########

def truncate(text, width=45):
    """
    Truncate long strings for console display.
    """

    text = str(text)

    if len(text) <= width:
        return text

    return text[:width - 3] + "..."


def format_value(value, column):
    """
    Format values for console output.
    """

    if pd.isna(value):
        return ""

    # Text
    if isinstance(value, str):
        return value

    # Year
    if column == "Year":
        return str(int(value))

    # Percentages
    if column.endswith("(%)") or column == "Percentage":
        return f"{float(value):.2f}"

    # Integers
    if isinstance(value, numbers.Integral):
        return f"{value:,}"

    # Floats that are actually integers (common with iterrows)
    if isinstance(value, numbers.Real):

        if float(value).is_integer():
            return f"{int(value):,}"

        return f"{value:.2f}"

    return truncate(value)

def print_dataset_summary(
    summary,
    logger=None,
):

    lines = []

    lines.append("")
    lines.append("Dataset summary")
    lines.append("-" * 70)

    width = summary["Metric"].str.len().max() + 2

    for _, row in summary.iterrows():

        value = row["Value"]

        if isinstance(value, numbers.Integral):
            value = f"{value:,}"

        lines.append(
            f"{row['Metric']:.<{width}} {str(value):>12}"
        )

    text = "\n".join(lines)

    if logger is None:
        print(text)
    else:
        logger.write(text)

def print_protein_summary(
    annotation_summary,
    major_lineages,
    logger=None,
):
    """
    Print a summary of hotspot mutations by protein region for the
    predominant lineage of each RSV subtype.
    """

    lines = []

    lines.append("")
    lines.append("Protein hotspot summary")
    lines.append("=" * 80)

    lines.append(
        f"Only substitutions with a frequency ≥ "
        f"{CONSOLE_MUTATION_FREQUENCY}%"
    )

    subtypes = (
        major_lineages
        .loc[
            major_lineages["Selected"],
            Columns.SUBTYPE,
        ]
        .dropna()
        .unique()
    )

    for subtype in subtypes:

        lineage = get_selected_lineage(
            major_lineages,
            subtype,
        )

        lines.append("")
        lines.append(
            f"RSV-{subtype} ({lineage})"
        )

        for protein in LOLLIPOP_PROTEINS:

            table = (
                annotation_summary[
                    (
                        annotation_summary[
                            Columns.SUBTYPE
                        ] == subtype
                    )
                    &
                    (
                        annotation_summary[
                            Columns.LINEAGE
                        ] == lineage
                    )
                    &
                    (
                        annotation_summary[
                            MutationColumns.PROTEIN
                        ] == protein
                    )
                ]
                .sort_values("start")
            )

            if table.empty:
                continue

            lines.append("")
            lines.append(
                f"  {protein} protein"
            )

            for _, row in table.iterrows():

                mutation = row["major_mutations"]

                if mutation in ("", "—"):
                    mutation = (
                        "No substitutions above threshold"
                    )

                lines.append(
                    f"    {row['feature']:<15}{mutation}"
                )

    lines.append(SUBLINE)

    text = "\n".join(lines)

    if logger is None:
        print(text)
    else:
        logger.write(text)

def summarize_transmission_clusters(
    table,
):
    """
    Prepare a compact transmission clustering summary for console output.
    """

    metrics = [
        "Transmission clusters",
        "Clustered sequences",
        "Clustered sequences (%)",
        "Singleton sequences",
        "Singleton sequences (%)",
        "Largest cluster size",
    ]

    summary = (
        table[
            table["Metric"].isin(metrics)
        ]
        .copy()
    )

    summary = summary.rename(
        columns={
            "Subtype": "RSV",
        }
    )

    return summary

def summarize_transmission_clusters_by_year(
    table,
):
    """
    Prepare transmission cluster counts by collection year.
    """

    summary = (
        table[
            table["Metric"] == "Transmission clusters"
        ]
        .copy()
    )

    summary = summary.rename(
        columns={
            Columns.COLLECTION_YEAR: "Year",
            Columns.SUBTYPE: "RSV",
        }
    )

    return summary.sort_values("Year")

def summarize_transmission_clusters_by_community(
    table,
):
    """
    Prepare transmission cluster counts by autonomous community.
    """

    summary = (
        table[
            table["Metric"] == "Transmission clusters"
        ]
        .copy()
    )

    summary = summary.rename(
        columns={
            Columns.COMMUNITY: "Community",
            Columns.SUBTYPE: "RSV",
        }
    )

    return summary

def print_table(
    title,
    table,
    max_rows=20,
    sort_by=None,
    ascending=False,
    logger=None,
):
    """
    Print a pandas DataFrame.
    """

    table = table.copy()

    if sort_by is not None:

        table = table.sort_values(
            sort_by,
            ascending=ascending,
        )

    lines = []

    lines.append("")
    lines.append(title)
    lines.append(SUBLINE)

    # Empty table
    if table.empty:

        lines.append("No data available.")
        lines.append(SUBLINE)

    else:

        # Column widths
        widths = []

        for column in table.columns:

            values = [
                format_value(value, column)
                for value in table[column]
            ]

            width = max(
                len(str(column)),
                max(len(v) for v in values),
            )

            widths.append(width)

        # Header
        lines.append(
            "  ".join(
                f"{column:<{width}}"
                for column, width in zip(
                    table.columns,
                    widths,
                )
            )
        )

        lines.append(
            "  ".join(
                "-" * width
                for width in widths
            )
        )

        # Rows
        for _, row in table.head(max_rows).iterrows():

            values = [

                format_value(value, column)

                for column, value in zip(
                    table.columns,
                    row,
                )
            ]

            lines.append(
                "  ".join(
                    f"{value:<{width}}"
                    for value, width in zip(
                        values,
                        widths,
                    )
                )
            )

        # Truncated table
        if len(table) > max_rows:

            lines.append("")
            lines.append(
                f"... ({len(table) - max_rows} more rows)"
            )

        lines.append(SUBLINE)

    text = "\n".join(lines)

    if logger is None:
        print(text)
    else:
        logger.write(text)