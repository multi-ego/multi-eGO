import pandas as pd
import json
import sys

mg_OO_c12_rep = 1.5e-6
#mg_HH_c12_rep = 1.2e-8
mg_HH_c12_rep = 3e-9
mg_ON_c12_rep = 1.5e-6
#mg_NN_c12_rep = 2.5e-5
mg_NN_c12_rep = 5e-6
mg_HO_sigma = 0.169500
mg_eps_ch3 = 0.13
mg_eps_ch2 = 0.10
mg_eps_pol = 0.10
mg_eps_ch1 = 0.09

# Dataframe with atom types and associated parameters
gromos_atp = pd.DataFrame(
    {
        "name": [
            "O",
            "OM",
            "OA",
            "N",
            "NT",
            "NL",
            "NR",
            "NZ",
            "NE",
            "C",
            "CH",
            "CH1",
            "CAH",
            "CH2",
            "CAH2",
            "CH3",
            "CH2r",
            "S",
            "CH3p",
            "P",
            "OE",
            "CR1",
            "H",
            "C0",
        ],
        "at.num": [8, 8, 8, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 16, 6, 15, 8, 6, 1, 20],
        "rc_c12": [
            2.631580e-07,  # "O",
            1.724403e-07,  # "OM",
            5.018430e-07,  # "OA",
            8.752940e-07,  # "N",
            2.596154e-06,  # "NT",
            8.752940e-07,  # "NL",
            1.506347e-06,  # "NR",
            8.752940e-07,  # "NZ",
            8.752940e-07,  # "NE",
            2.598570e-06,  # "C",
            2.598570e-06,  # "CH"
            6.555574e-05,  # "CH1"
            6.555574e-05,  # "CAH"
            1.543890e-05,  # "CH2"
            1.543890e-05,  # "CAH2"
            8.595562e-06,  # "CH3"
            1.193966e-05,  # "CH2r"
            2.724050e-06,  # "S",
            8.736473e-06,  # "CH3p"
            3.893600e-06,  # "P",
            3.558824e-07,  # "OE",
            6.298560e-06,  # "CR1",
            9.148590e-10,  # "H",
            2.659360e-07,  # "C0",
        ],
        # 4*sig^12*eps 
        "mg_c12": [
            4.*0.276007**12 * mg_eps_pol,  # "O",
            4.*0.262585**12 * mg_eps_pol,  # "OM",
            4.*0.295484**12 * mg_eps_pol,  # "OA",
            4.*0.313647**12 * mg_eps_pol,  # "N",
            4.*0.357220**12 * mg_eps_pol,  # "NT",
            4.*0.313647**12 * mg_eps_pol,  # "NL",
            4.*0.334113**12 * mg_eps_pol,  # "NR",
            4.*0.313647**12 * mg_eps_pol,  # "NZ",
            4.*0.313647**12 * mg_eps_pol,  # "NE",
            4.*0.358118**12 * mg_eps_pol,  # "C",
            4.*0.358118**12 * mg_eps_ch3,  # "CH"
            4.*0.501918**12 * mg_eps_ch1,  # "CH1"
            4.*0.501918**12 * mg_eps_ch1,  # "CAH"
            4.*0.407038**12 * mg_eps_ch2,  # "CH2"
            4.*0.407038**12 * mg_eps_ch2,  # "CAH2"
            4.*0.374792**12 * mg_eps_ch3,  # "CH3"
            4.*0.395474**12 * mg_eps_ch2,  # "CH2r"
            4.*0.330769**12 * mg_eps_pol,  # "S",
            4.*0.374792**12 * mg_eps_ch3,  # "CH3p"
            4.*0.338557**12 * mg_eps_pol,  # "P",
            4.*0.284916**12 * mg_eps_pol,  # "OE",
            4.*0.374119**12 * mg_eps_pol,  # "CR1",
            0.0000000e-00 * mg_eps_pol,  # "H",
            0.0000000e-00 * mg_eps_pol,  # "C0",
        ],
        "mg_c6": [
            4.*0.276007**6 * mg_eps_pol,  # "O",
            4.*0.262585**6 * mg_eps_pol,  # "OM",
            4.*0.295484**6 * mg_eps_pol,  # "OA",
            4.*0.313647**6 * mg_eps_pol,  # "N",
            4.*0.357220**6 * mg_eps_pol,  # "NT",
            4.*0.313647**6 * mg_eps_pol,  # "NL",
            4.*0.334113**6 * mg_eps_pol,  # "NR",
            4.*0.313647**6 * mg_eps_pol,  # "NZ",
            4.*0.313647**6 * mg_eps_pol,  # "NE",
            4.*0.358118**6 * mg_eps_pol,  # "C",
            4.*0.358118**6 * mg_eps_ch3,  # "CH"
            4.*0.501918**6 * mg_eps_ch1,  # "CH1"
            4.*0.501918**6 * mg_eps_ch1,  # "CAH"
            4.*0.407038**6 * mg_eps_ch2,  # "CH2"
            4.*0.407038**6 * mg_eps_ch2,  # "CAH2"
            4.*0.374792**6 * mg_eps_ch3,  # "CH3"
            4.*0.395474**6 * mg_eps_ch2,  # "CH2r"
            4.*0.330769**6 * mg_eps_pol,  # "S",
            4.*0.374792**6 * mg_eps_ch3,  # "CH3p"
            4.*0.338557**6 * mg_eps_pol,  # "P",
            4.*0.284916**6 * mg_eps_pol,  # "OE",
            4.*0.374119**6 * mg_eps_pol,  # "CR1",
            0.0000000000 * mg_eps_pol,  # "H", # TODO
            0.0000000000 * mg_eps_pol,  # "C0",
        ],
    }
)

# Dictionary mapping atom types from a force field to multiego representation
from_ff_to_multiego = {
    "HN": "H",
    "OC1": "O1",
    "OC2": "O2",
    "OT1": "O1",
    "OT2": "O2",
    "Ca+2": "Cal",
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
    types_dict["backbone_calpha"] = (df["name"] == "CA").to_numpy()
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
    types_dict["sidechain_cds"] = (
        (df["name"] == "CD")
        | (df["name"] == "CD1")
        | (df["name"] == "CD2")
        | (df["name"] == "SD")
        | (df["name"] == "OD")
        | (df["name"] == "OD1")
        | (df["name"] == "OD2")
        | (df["name"] == "ND1")
        | (df["name"] == "ND2") & (df["resname"] != "PRO")
    ).to_numpy()

    return types_dict


# List of atom type combinations for LJ14 pairs
atom_type_combinations = [
    # Tuple of atom type combinations for LJ14 pairs
    ("backbone_carbonyl", "sidechain_cb", 0.275, 1.299682e-06, 1),
    ("backbone_oxygen", "sidechain_cb", 1, 1.5e-6, 0),
    ("ct_oxygen", "sidechain_cb", 1, 1.5e-6, 0),
    ("backbone_nitrogen", "sidechain_cb", 1, 2.7e-6, -1),
    ("first_backbone_nitrogen", "backbone_nitrogen", None, 4.0e-6, 1),
    ("backbone_nitrogen", "backbone_nitrogen", 0.343, None, 1),
    ("backbone_carbonyl", "backbone_carbonyl", 0.5, None, -1),
    ("sidechain_cgs", "backbone_carbonyl", 0.250, 1.2e-6, 0),
    ("sidechain_cgs", "backbone_nitrogen", 0.200, 5.5e-7, 0),
    ("sidechain_cgs", "first_backbone_nitrogen", 0.200, 5.5e-7, 0),
    ("sidechain_cds", "backbone_calpha", 0.100, 5e-7, 0),
]

# List of amino acids and nucleic acids
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


def parse_json(file_path):
    if file_path:
        try:
            with open(file_path, "r") as file:
                custom_dict = json.load(file)
                if not isinstance(custom_dict, dict):
                    raise ValueError("Error in reading the custom dictionary: Invalid dictionary format")
                return custom_dict
        except (json.JSONDecodeError, ValueError) as e:
            print(f"Error in reading the custom dictionary: {e}")
            sys.exit()
    else:
        return {}
