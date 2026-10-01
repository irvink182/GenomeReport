"""
Dataset report CLI.
"""

from pathlib import Path

from ..dataset import Dataset
from ..report import DatasetReport


def run(args):

    dataset = Dataset(args.input)

    report = DatasetReport(
        dataset,
        output=args.output,
        phylogenies=ars.phylogeny,
    )

    report.build()


def add_parser(subparsers):

    parser = subparsers.add_parser(
        "report",
        help="Generate a dataset report.",
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

    parser.set_defaults(func=run)