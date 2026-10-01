"""
Command-line interface.
"""

import argparse

from pathlib import Path
from importlib.metadata import version
from ..dataset import Dataset
from ..report import DatasetReport


def main():

    parser = argparse.ArgumentParser(
        prog="genomereport",
        description="Generate genomic dataset reports.",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"GenomeReport {version('genomereport')}",
    )

    parser.add_argument(
        "input",
        help="Input metadata table.",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="results",
        help="Output directory.",
    )

    parser.add_argument(
        "--phylogeny",
        nargs="+",
        type=Path,
        default=None,
        metavar="PNG",
        help=(
            "One or more phylogenetic PNG figures "
            "to include in the report."
        ),
    )

    args = parser.parse_args()

    dataset = Dataset(args.input)

    report = DatasetReport(
        dataset,
        output=args.output,
        phylogenies=args.phylogeny,
        input_file=args.input,
    )

    report.build()


if __name__ == "__main__":
    main()