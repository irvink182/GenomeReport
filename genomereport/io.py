"""
Input/output utilities.
"""

from pathlib import Path
import matplotlib.pyplot as plt

import pandas as pd

from .config import Report

SUPPORTED_EXTENSIONS = {
    ".tsv",
    ".csv",
    ".txt",
    ".xlsx",
}

##############
# Tables
##############
def read_table(filepath):
    """
    Read a tabular file into a pandas DataFrame.

    Parameters
    ----------
    filepath : str or Path

    Returns
    -------
    pandas.DataFrame
    """

    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(
            f"File not found: {filepath}"
        )

    suffix = filepath.suffix.lower()

    if suffix not in SUPPORTED_EXTENSIONS:
        raise ValueError(
            f"Unsupported format: {suffix}"
        )

    if suffix in {".tsv", ".txt"}:

        return pd.read_csv(
            filepath,
            sep="\t",
            keep_default_na=False,
            dtype=str,
        )

    if suffix == ".csv":

        return pd.read_csv(
            filepath,
            keep_default_na=False,
            dtype=str,
        )

    if suffix == ".xlsx":

        return pd.read_excel(
            filepath,
            keep_default_na=False,
            dtype=str,
        )


def write_table(df, filepath):
    """
    Write a DataFrame.
    """

    filepath = Path(filepath)

    suffix = filepath.suffix.lower()

    if suffix in {".tsv", ".txt"}:

        df.to_csv(
            filepath,
            sep="\t",
            index=False,
        )

    elif suffix == ".csv":

        df.to_csv(
            filepath,
            index=False,
        )

    elif suffix == ".xlsx":

        df.to_excel(
            filepath,
            index=False,
        )

    else:

        raise ValueError(
            f"Unsupported format: {suffix}"
        )

##############
# Figures
##############

def write_figure(fig, filepath):
    """
    Save matplotlib figure.
    """

    filepath = Path(filepath)

    for fmt in Report.FIGURE_FORMATS:

        kwargs = {
            "bbox_inches": "tight",
        }

        if fmt == "png":
            kwargs["dpi"] = Report.FIGURE_DPI

        fig.savefig(
            filepath.with_suffix(f".{fmt}"),
            **kwargs,
        )

    plt.close(fig)

##############
# HTML
##############