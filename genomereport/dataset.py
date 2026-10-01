"""
Dataset object.
"""

from pathlib import Path

from .io import read_table

from .config import Columns

from .normalize import (
    normalize_patient_age,
    normalize_sequencing_technology,
    normalize_assembly_method,
    normalize_patient_status,
    normalize_gender,
)


class Dataset:

    def __init__(self, filepath):

        self.filepath = Path(filepath)

        self.name = self.filepath.stem

        self.extension = self.filepath.suffix

        self.df = read_table(filepath)

        self.prepare()

        self.tables = {}

        self.figures = {}

        self.statistics = {}

    def prepare(self):
        """
        Prepare dataset for downstream analyses.
        """

        self.normalize()

        return self

    def normalize(self):
        """
        Normalize metadata.
        """

        ###########################################################################
        # Patient age
        ###########################################################################

        self.df[Columns.AGE_NORMALIZED] = (
            self.df[Columns.AGE]
            .apply(normalize_patient_age)
        )

        ###########################################################################
        # Sequencing
        ###########################################################################

        self.df[Columns.TECHNOLOGY_NORMALIZED] = (
            self.df[Columns.TECHNOLOGY]
            .apply(normalize_sequencing_technology)
        )

        self.df[Columns.ASSEMBLY_NORMALIZED] = (
            self.df[Columns.ASSEMBLY]
            .apply(normalize_assembly_method)
        )

        ###########################################################################
        # Patient Status
        ###########################################################################

        self.df[Columns.STATUS_NORMALIZED] = (
            self.df[Columns.STATUS]
            .apply(normalize_patient_status)
        )

        ###########################################################################
        # Patient Gender
        ###########################################################################

        self.df[Columns.GENDER_NORMALIZED] = (
            self.df[Columns.GENDER]
            .apply(normalize_gender)
        )

    @property
    def n_rows(self):

        return len(self.df)

    @property
    def n_columns(self):

        return len(self.df.columns)

    @property
    def columns(self):

        return list(self.df.columns)

    @property
    def n_communities(self):

        return self.df[
            Columns.COMMUNITY
        ].nunique()