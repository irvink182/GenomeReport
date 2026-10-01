###############################################################################
# Missing values
###############################################################################

MISSING_VALUES = {
    "",
    "NA",
    "N/A",
    "UNKNOWN",
    "NOT AVAILABLE",
    "NONE",
    "NULL",
}

###############################################################################
# Age units
###############################################################################

YEAR_UNITS = {
    "year",
    "years",
    "yr",
    "yrs",
    "y",
}

MONTH_UNITS = {
    "month",
    "months",
    "mes",
    "meses",
    "mo",
    "mos",
    "moth",
    "moths",
}

WEEK_UNITS = {
    "week",
    "weeks",
    "semana",
    "semanas",
}

DAY_UNITS = {
    "day",
    "days",
    "dia",
    "dias",
    "día",
    "días",
}

###############################################################################
# Other variables
###############################################################################

SUBTYPES = [
    "RSV-A",
    "RSV-B",
]

INTEGER_COLUMNS = {
    "Complete",
    "Missing",
    "Total",
    "A",
    "B",
}

FLOAT_COLUMNS = {
    "Missing (%)",
    "Percentage",
}

SAMPLING_STRATEGY_MAPPING = {

    "Sentinel": [
        "sentinel surveillance",
    ],

    "Non-sentinel": [
        "non-sentinel",
    ],

}

TECHNOLOGY_MAPPING = {
    "Illumina": [
        "illumina",
        "miseq",
        "nextseq",
        "novaseq",
        "iseq",
        "miniseq",
        "eseq",
    ],
    "Nanopore": [
        "nanopore",
        "minion",
        "gridion",
    ],
    "Sanger": [
        "sanger",
    ],
    "Ion Torrent": [
        "ion torrent",
    ],
    "DNBSEQ": [
        "dnbseq",
    ],
}

ASSEMBLY_MAPPING = {

    "SeqScape": [
        "seqscape",
        "secscape",
    ],

    "MICGENT": [
        "micgent",
        "ngs-mstb",
    ],

    "Geneious": [
        "geneious",
    ],

    "IRMA": [
        "irma",
    ],

    "Genome Detective": [
        "genome detective",
    ],

    "Galaxy": [
        "galaxy",
    ],

    "SPAdes": [
        "spades",
    ],

    "Ampliseq": [
        "ampliseq",
    ],

    "ViralRecon": [
        "viralrecon",
        "nf-core/viralrecon",
    ],

    "ARTIC": [
        "artic",
    ],

    "CLC": [
        "clc",
    ],

    "INSaFLU": [
        "insaflu",
    ],

    "iVar": [
        "ivar",
    ],

    "Other": [],
}

STATUS_MAPPING = {

    "Hospitalized": [
        "hospitalized",
    ],

    "Emergency": [
        "emergency room",
        "emergency ward",
        "emergency",
    ],

    "Outpatient": [
        "outpatient",
        "ambulatory",
    ],

    "Primary care": [
        "primary-care",
        "primary care",
    ],

    "Released": [
        "released",
    ],

    "Live": [
        "live",
    ],

    "Other": [],
}

GENDER_MAPPING = {

    "Female": [
        "female",
        "f",
    ],

    "Male": [
        "male",
        "m",
    ],

    "Unknown": [
        "u",
        "unknown",
        "U",
    ],

    "Other": [],
}

TRANSMISSION_CLUSTER_THRESHOLDS = (0, 5)

###############################################################################
# Protein lengths (GISAID reference sequences)
#
# RSV-A: hRSV/A/England/397/2017
# RSV-B: hRSV/B/Australia/VIC-RCH056/2019
###############################################################################
PROTEIN_ORDER = [

    "NS1",
    "NS2",
    "N",
    "P",
    "M",
    "SH",
    "G",
    "F",
    "M2-2",
    "M2-1",
    "L",

]

PROTEIN_LENGTHS = {
    
    "A": {
        "NS1": 139,
        "NS2": 124,
        "N": 391,
        "P": 241,
        "M": 256,
        "SH": 64,
        "G": 321,
        "F": 574,
        "M2-2": 90,
        "M2-1": 194,
        "L": 2165,
    },

    "B": {
        "NS1": 139,
        "NS2": 84,
        "N": 391,
        "P": 241,
        "M": 256,
        "SH": 65,
        "G": 312,
        "F": 574,
        "M2-2": 56,
        "M2-1": 194,
        "L": 2165,
    }
}

LOLLIPOP_PROTEINS = [
    "G",
    "F",
    "L",
]


###############################################################################
# Defining Classes
###############################################################################

class Columns:

    # General
    ACCESSION = "Accession ID"
    VIRUS_NAME = "Virus name"
    SUBTYPE = "Subtype"
    SAMPLING_STRATEGY = "Sampling strategy"
    ADDITIONAL_HOST_INFORMATION = "Additional host information"

    # Temporal
    COLLECTION_DATE = "Collection date"
    YEAR = "Year"
    COLLECTION_YEAR = "gisaid_collection_year"

    # Geography
    LOCATION = "Location"
    COMMUNITY = "autonomous_community"
    CITY = "city"

    # Laboratory
    INSTITUTION = "institution"
    ORIGINATING_LAB = "originating_lab"
    SUBMITTER_LAB = "submitter_lab_name"

    # Sequencing
    SEQUENCE_LENGTH = "Sequence Length"
    TECHNOLOGY = "Sequencing technology"
    ASSEMBLY = "Assembly method"

    # Virus
    LINEAGE = "Lineage"
    GENOTYPE = "Genotype"

    # Patient
    HOST = "Host"
    GENDER = "Gender"
    AGE = "Patient age"
    STATUS = "Patient status"

    # Mutations
    AA_SUBSTITUTIONS = "AA Substitutions"

    # Normalized columns
    AGE_NORMALIZED = "patient_age"
    TECHNOLOGY_NORMALIZED = "sequencing_technology"
    ASSEMBLY_NORMALIZED = "assembly_method"
    STATUS_NORMALIZED = "patient_status"
    GENDER_NORMALIZED = "gender"

class MutationColumns:

    PROTEIN = "Protein"
    MUTATION = "Mutation"
    AA_SUBSTITUTION = "AA substitution"
    POSITION = "Position"
    REFERENCE_AA = "Reference AA"
    ALTERNATIVE_AA = "Alternative AA"
    HAS_MUTATIONS = "Has mutations"

    COUNT = "Count"
    GENOMES = "Genomes"
    FREQUENCY = "Frequency (%)"
    GENOMES_WITH_MUTATIONS = "Genomes with mutations"
    TOTAL_EVENTS = "Total mutation events"
    DISTINCT_SUBSTITUTIONS = "Distinct AA substitutions"
    MEAN_PER_GENOME = "Mean substitutions / genome"
    MAX_POSITION = "Max position"
    PERCENT_GENOMES_WITH_MUTATIONS = "Genomes with mutations (%)"
    MEAN_PER_MUTATED_GENOME = (
    "Mean substitutions / mutated genome"
    )

    PROTEIN_LENGTH = "Protein length"
    EVENTS_PER_GENOME = "Mutation events per genome"
    EVENTS_PER_100_AA = "Mutation events per 100 AA"
    EVENTS_PER_100_AA_PER_GENOME = "Mutation events per 100 AA per genome"
    
    DISTINCT_PER_GENOME = "Distinct AA substitutions per genome"
    DISTINCT_PER_100_AA ="Distinct substitutions per 100 AA"


    GENOMES_MUTATED = "Genomes mutated"
    FREQUENCY = "Frequency (%)"
    POSITION_SUMMARY = "Protein position summary"
    MOST_COMMON_SUBSTITUTION = "Most common subsstitution"
    MOST_COMMON_COUNT = "Most common count"
    MOST_COMMON_FREQUENCY = "Most common frequency"
    DOMINANCE = "Dominance (%)"
    LABEL = "Label"

class MutationMetrics:

    SEQUENCES = "Sequences"
    WITH_MUTATIONS = "Sequences with mutations"
    WITHOUT_MUTATIONS = "Sequences without mutations"
    PERCENT_WITH_MUTATIONS = "Genomes with mutations (%)"
    MEAN = "Mean substitutions / genome"
    MEDIAN = "Median substitutions / genome"
    MAX = "Maximum substitutions"

class ClusterMetrics:

    SEQUENCES = "Sequences"
    TRANSMISSION_CLUSTERS = "Transmission clusters"
    CLUSTERED_SEQUENCES = "Clustered sequences"
    CLUSTERED_SEQUENCES_PERCENT = "Clustered sequences (%)"
    SINGLETON_SEQUENCES = "Singleton sequences"
    SINGLETON_SEQUENCES_PERCENT = "Singleton sequences (%)"
    LARGEST_CLUSTER_SIZE = "Largest cluster size"
    MEAN_CLUSTER_SIZE = "Mean cluster size"
    MEDIAN_CLUSTER_SIZE = "Median cluster size"

class ClusterColumns:

    PREFIX = "Cluster_"
    CLUSTER = "Cluster"
    SIZE = "Cluster size"
    YEAR = "Year"
    LINEAGE = "Lineage"
    COMMUNITY = "Autonomous community"

class Metadata:

    CORE = [

        Columns.ACCESSION,
        Columns.SUBTYPE,
        Columns.COLLECTION_DATE,
        Columns.LOCATION,
        Columns.SAMPLING_STRATEGY,
        Columns.COMMUNITY,
        Columns.CITY,
        Columns.INSTITUTION,
        Columns.ORIGINATING_LAB,
        Columns.SUBMITTER_LAB,
        Columns.TECHNOLOGY_NORMALIZED,
        Columns.ASSEMBLY_NORMALIZED,
        Columns.LINEAGE,
        Columns.GENOTYPE,
        Columns.AA_SUBSTITUTIONS,
    ]

    PATIENT = [

        Columns.HOST,
        Columns.GENDER_NORMALIZED,
        Columns.AGE_NORMALIZED,
        Columns.STATUS_NORMALIZED,
    ]

CATEGORY_NAMES = {
    Columns.LINEAGE: "virus_lineage",
    Columns.GENOTYPE: "virus_genotype",
}

class ReportColumns:

    YEAR = "Year"
    MONTH = "Month"
    PERCENTAGE = "Percentage"

    NUMERIC = {
        "A",
        "B",
        "Total",
        "Complete",
        "Missing",
        "Percentage",
        "Missing (%)",
    }

class Report:

    TABLE_EXTENSION = ".tsv"

    FIGURE_FORMATS = [
        "png",
        "svg",
    ]

    FIGURE_DPI = 300