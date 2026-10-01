"""
Normalization utilities.
"""

import pandas as pd

from .config import (
    MISSING_VALUES,
    YEAR_UNITS,
    MONTH_UNITS,
    WEEK_UNITS,
    DAY_UNITS,
    SAMPLING_STRATEGY_MAPPING,
    TECHNOLOGY_MAPPING,
    ASSEMBLY_MAPPING,
    STATUS_MAPPING,
    GENDER_MAPPING,
)

def is_missing(value):
    """
    Return True if a value should be considered missing.
    """

    if pd.isna(value):
        return True

    value = str(value).strip().upper()

    return value in MISSING_VALUES

def normalize_patient_age(value):
    """
    Normalize patient age to years.

    Parameters
    ----------
    value : str

    Returns
    -------
    float or None
        Age expressed in years.
        Returns None for missing or invalid values.
    """

    if is_missing(value):
        return None

    value = str(value).strip()

    # Already numeric
    try:
        return round(float(value), 3)

    except ValueError:
        pass

    # Parse "<number> <unit>"
    tokens = value.lower().split()

    if len(tokens) != 2:
        return None

    try:
        number = float(tokens[0])

    except ValueError:
        return None

    unit = tokens[1]

    # Years
    if unit in YEAR_UNITS:
        return round(number, 3)

    # Months
    if unit in MONTH_UNITS:
        return round(number / 12, 3)

    # Weeks
    if unit in WEEK_UNITS:
        return round(number *7 / 365.25, 3)

    # Days
    if unit in DAY_UNITS:
        return round(number / 365.25, 3)

    return None

def is_valid_patient_age(age):
    """
    Return True if age is biologically plausible.
    """

    if age is None:
        return False

    return 0 <= age <= 120

def normalize_mapping(value, mapping):
    """
    Normalize a value using a keyword mapping.
    """

    if is_missing(value):
        return "NA"

    value = " ".join(
        str(value).lower().split()
    )

    for category, keywords in mapping.items():

        for keyword in keywords:

            if keyword in value:
                return category

    return "Other"

def normalize_sequencing_technology(value):
    """
    Normalize sequencing technology into broad categories.
    """

    return normalize_mapping(
        value,
        TECHNOLOGY_MAPPING,
    )

def normalize_assembly_method(value):
    """
    Normalize assembly methods into broad categories.
    """

    return normalize_mapping(
        value,
        ASSEMBLY_MAPPING,
    )

def normalize_patient_status(value):
    """
    Normalize patient status into broad categories.
    """

    return normalize_mapping(
        value,
        STATUS_MAPPING,
    )

def normalize_gender(value):
    """
    Normalize patient gender.
    """

    return normalize_mapping(
        value,
        GENDER_MAPPING,
    )

def normalize_sampling_strategy(
    value,
):
    """
    Normalize sampling strategy into broad surveillance categories.
    """

    if is_missing(
        value,
    ):

        return "Sentinel"

    return normalize_mapping(
        value,
        SAMPLING_STRATEGY_MAPPING,
    )