import pandas as pd

# Dataframe with GROMOS atom types and associated parameters
gromos_atp = pd.DataFrame(
    {
        "name": [
            "O",
            "OA",
            "N",
            "C",
            "CH1",
            "CH2",
            "CH3",
            "CH2r",
            "NT",
            "S",
            "NR",
            "OM",
            "NE",
            "NL",
            "NZ",
            "CH3p",
            "P",
            "OE",
            "CR1",
        ],
        "at.num": [8, 8, 7, 6, 6, 6, 6, 6, 7, 16, 7, 8, 7, 7, 7, 6, 15, 8, 6],
        "c12": [
            2.631580e-07,
            5.018430e-07,
            8.752940e-07,
            2.598570e-06,
            6.555574e-05,  # CH1
            1.543890e-05,
            8.595562e-06,
            1.193966e-05,
            2.596154e-06,
            2.724050e-06,
            1.506347e-06,
            1.724403e-07,
            8.752940e-07,
            8.752940e-07,
            8.752940e-07,
            8.736473e-06,
            3.893600e-06,
            3.558824e-07,
            6.29856e-06,
        ],
    }
)

# Dictionary mapping atom types from a force field to multiego representation
from_ff_to_multiego = {
    "OC1": "O1",
    "OC2": "O2",
    "OT1": "O1",
    "OT2": "O2",
}


def lj14_generator(df):
    """
    Generates types dictionary based on the provided DataFrame.

    Args:
    - df (pd.DataFrame): DataFrame containing atom types and parameters.

    Returns:
    - types_dict (dict): Dictionary containing different atom type combinations.
    """
    types_dict = {}
    types_dict["first_backbone_nitrogen"] = ((df["name"] == "N") & (df["type"] == "NL")).to_numpy()
    types_dict["backbone_nitrogen"] = ((df["name"] == "N") & (df["type"] != "NL")).to_numpy()
    types_dict["backbone_carbonyl"] = (df["name"] == "C").to_numpy()
    types_dict["backbone_oxygen"] = (df["name"] == "O").to_numpy()
    types_dict["ct_oxygen"] = ((df["name"] == "O1") | (df["name"] == "O2")).to_numpy()
    types_dict["sidechain_cb"] = (df["name"] == "CB").to_numpy()
    types_dict["sidechain_cgs"] = (
        (df["name"] == "CG")
        | (df["name"] == "CG1")
        | (df["name"] == "CG2")
        | (df["name"] == "SG")
        | (df["name"] == "OG")
        | (df["name"] == "OG1") & (df["resname"] != "PRO")
    ).to_numpy()

    return types_dict


# List of atom type combinations for LJ14 pairs
atom_type_combinations = [
    # Tuple of atom type combinations for LJ14 pairs
    ("backbone_carbonyl", "sidechain_cb", 0.275, None, 1),
    ("backbone_oxygen", "sidechain_cb", 0.1, None, 0),
    ("ct_oxygen", "sidechain_cb", 0.1, None, 0),
    ("backbone_nitrogen", "sidechain_cb", 0.65, None, -1),
    ("first_backbone_nitrogen", "backbone_nitrogen", None, 4.0e-6, 1),
    ("backbone_nitrogen", "backbone_nitrogen", 0.343, None, 1),
    ("backbone_carbonyl", "backbone_carbonyl", 0.5, None, -1),
    ("sidechain_cgs", "backbone_carbonyl", 0.078, None, 0),
    ("sidechain_cgs", "backbone_nitrogen", 0.087, None, 0),
    ("sidechain_cgs", "first_backbone_nitrogen", 0.087, None, 0),
]

# List of amino acids and nucleic acids
# TODO add capped termini
aminoacids_list = [
    "VAL",
    "ILE",
    "LEU",
    "GLU",
    "GLN",
    "ASP",
    "ASN",
    "HIS",
    "TRP",
    "PHE",
    "TYR",
    "ARG",
    "LYS",
    "SER",
    "THR",
    "MET",
    "ALA",
    "GLY",
    "PRO",
    "CYS",
    "ACE",
    "NME",
]
# TODO to check
nucleic_acid_list = ["A", "C", "G", "T"]
