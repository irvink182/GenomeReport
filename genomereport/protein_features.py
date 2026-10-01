###############################################################################
# RSV protein annotations
###############################################################################

G_FEATURES_A = [

    {
        "start": 1,
        "end": 37,
        "label": "Signal",
        "type": "signal",
        "level": "domain",
        "function": "Initial membrane anchoring.",
        "display": True,
    },

    {
        "start": 38,
        "end": 66,
        "label": "TM",
        "type": "transmembrane",
        "level": "domain",
        "function": "Hydrophobic membrane anchor.",
        "display": True,
    },

    {
        "start": 67,
        "end": 149,
        "label": "HVR1",
        "type": "hypervariable",
        "level": "domain",
        "function": "Immune evasion.",
        "display": True,
    },

    {
        "start": 151,
        "end": 190,
        "label": "CCD",
        "type": "domain",
        "level": "domain",
        "function": "Central conserved domain.",
        "display": True,
    },

    {
        "start": 182,
        "end": 186,
        "label": "CX3C",
        "type": "motif",
        "level": "feature",
        "function": "CX3CR1 receptor binding motif.",
        "display": True,
    },

    {
        "start": 191,
        "end": 321,
        "label": "HVR2",
        "type": "hypervariable",
        "level": "domain",
        "function": "Contains ON1 duplication.",
        "display": True,
    },

]

G_FEATURES_B = [

    {
        "start": 1,
        "end": 37,
        "label": "Signal",
        "type": "signal",
        "level": "domain",
        "function": "Initial membrane anchoring.",
        "display": True,
    },

    {
        "start": 38,
        "end": 66,
        "label": "TM",
        "type": "transmembrane",
        "level": "domain",
        "function": "Hydrophobic membrane anchor.",
        "display": True,
    },

    {
        "start": 67,
        "end": 149,
        "label": "HVR1",
        "type": "hypervariable",
        "level": "domain",
        "function": "Immune evasion.",
        "display": True,
    },

    {
        "start": 151,
        "end": 190,
        "label": "CCD",
        "type": "domain",
        "level": "domain",
        "function": "Central conserved domain.",
        "display": True,
    },

    {
        "start": 182,
        "end": 186,
        "label": "CX3C",
        "type": "motif",
        "level": "feature",
        "function": "CX3CR1 receptor binding motif.",
        "display": True,
    },

    {
        "start": 191,
        "end": 312,
        "label": "HVR2",
        "type": "hypervariable",
        "level": "domain",
        "function": "Contains ON1 duplication.",
        "display": True,
    },

]

F_FEATURES = [
    {
        "start": 1,
        "end": 25,
        "label": "Signal",
        "type": "signal",
        "level": "domain",
        "function": "Removed during maturation.",
        "display": True,
    },
    
    {
        "start": 26,
        "end": 109,
        "label": "F2",
        "type": "domain",
        "level": "domain",
        "function": "External globular head.",
        "display": True,
    },
            
    {
        "start": 110,
        "end": 136,
        "label": "Pep27",
        "type": "domain",
        "level": "domain",
        "function": "Furin cleavage peptide.",
        "display": True,
    },
    
    {
        "start": 137,
        "end": 574,
        "label": "F1",
        "type": "domain",
        "level": "domain",
        "function": "Fusion machinery.",
        "display": True,
    },
    
    {
        "start": 137,
        "end": 154,
        "label": "Fusion peptide",
        "type": "domain",
        "level": "feature",
        "function": "Target membrane insertion.",
        "display": True,
    },
    
    {
        "start": 62,
        "end": 69,
        "label": "Site Ø",
        "type": "epitope",
        "level": "feature",
        "target": "Nirsevimab",
        "function": "Nirsevimab epitope (part 1).",
        "display": True,
    },
    
    {
        "start": 196,
        "end": 212,
        "label": "Site Ø",
        "type": "epitope",
        "level": "feature",
        "target": "Nirsevimab",
        "function": "Nirsevimab epitope (part 2).",
        "display": True,
    },
    
    {
        "start": 250,
        "end": 277,
        "label": "Site II",
        "type": "epitope",
        "level": "feature",
        "target": "Palivizumab",
        "function": "Palivizumab epitope.",
        "display": True,
    },
    
    {
        "start": 422,
        "end": 468,
        "label": "Site IV",
        "type": "epitope",
        "level": "feature",
        "target": "Clesrovimab",
        "function": "Clesrovimab epitope.",
        "display": True,
    },
    
    {
        "start": 469,
        "end": 500,
        "label": "Trimer cavity",
        "type": "drug_target",
        "level": "feature",
        "target": "Ziresovir",
        "function": "Ziresovir binding site.",
        "display": True,
    },
    
    {
        "start": 525,
        "end": 548,
        "label": "TM",
        "type": "transmembrane",
        "level": "domain",
        "function": "Membrane anchor.",
        "display": True,
    },

]

L_FEATURES = [
    {
        "start": 1,
        "end": 550,
        "label": "N-term",
        "type": "domain",
        "level": "domain",
        "function": "RNA binding / initiation.",
        "display": True,
    },

    {
        "start": 551,
        "end": 1100,
        "label": "RdRp",
        "type": "domain",
        "level": "domain",
        "function": "RNA-dependent RNA polymerase.",
        "display": True,
    },

    {
        "start": 810,
        "end": 813,
        "label": "Motif C",
        "type": "motif",
        "level": "feature",
        "function": "Catalytic GDNQ motif.",
        "display": True,
    },
    
    {
        "start": 1101,
        "end": 1450,
        "label": "Capping",
        "type": "domain",
        "level": "domain",
        "function": "mRNA capping.",
        "display": True,
    },
    
    {
        "start": 1200,
        "end": 1250,
        "label": "Lumicitabine",
        "type": "drug_target",
        "level": "feature",
        "function": "Resistance-associated region.",
        "display": True,
    },
    
    {
        "start": 1600,
        "end": 1900,
        "label": "MTase",
        "type": "domain",
        "level": "domain",
        "function": "mRNA methyltransferase.",
        "display": True,
    },

]

PROTEIN_FEATURES = {

    "A": {

        "G": G_FEATURES_A,
        "F": F_FEATURES,
        "L": L_FEATURES,

    },

    "B": {

        "G": G_FEATURES_B,
        "F": F_FEATURES,
        "L": L_FEATURES,

    },
}