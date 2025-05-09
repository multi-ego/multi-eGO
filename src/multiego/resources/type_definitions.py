import pandas as pd
import json
import sys

mg_OO_c12_rep = 3e-6
mg_HH_c12_rep = 4.148590e-08
mg_ON_c12_rep = 9.799381e-06
mg_eps = 0.11
mg_eps_ch2 = 0.10
mg_eps_ch1 = 0.09

# Dataframe with GROMOS atom types and associated parameters
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
            "CH1a",
            "CH2",
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
            6.555574e-05,  # "CH1a"
            1.543890e-05,  # "CH2"
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
        "mg_c12": [
            1.0000000e-06 / 1.27911 * mg_eps,  # "O",
            7.4149321e-07 / 1.72504 * mg_eps,  # "OM",
            1.5055290e-06 / 0.84961 * mg_eps,  # "OA",
            2.3195290e-06 / 0.63980 * mg_eps,  # "N",
            5.0625000e-06 / 0.29314 * mg_eps,  # "NT",
            2.3195290e-06 / 0.63980 * mg_eps,  # "NL",
            3.3892810e-06 / 0.43786 * mg_eps,  # "NR",
            2.3195290e-06 / 0.63980 * mg_eps,  # "NZ",
            2.3195290e-06 / 0.63980 * mg_eps,  # "NE",
            4.9372840e-06 / 0.27741 * mg_eps,  # "C",
            4.9372840e-06 / 0.27741 * mg_eps,  # "CH"
            9.7022500e-05 / 0.09489 * mg_eps_ch1,  # "CH1"
            9.7022500e-05 / 0.09489 * mg_eps_ch1,  # "CAH"
            9.7022500e-05 / 0.09489 * mg_eps_ch1,  # "CH1a"
            3.3965584e-05 / 0.41050 * mg_eps_ch2,  # "CH2"
            2.6646244e-05 / 0.86710 * mg_eps,  # "CH3"
            2.8058209e-05 / 0.47920 * mg_eps_ch2,  # "CH2r"
            1.3075456e-05 / 1.90587 * mg_eps,  # "S",
            2.6646244e-05 / 0.86715 * mg_eps,  # "CH3p"
            2.2193521e-05 / 2.44674 * mg_eps,  # "P",
            1.2100000e-06 / 1.05711 * mg_eps,  # "OE",
            1.5116544e-05 / 0.50266 * mg_eps,  # "CR1",
            0.0000000e-00 / 1.00000 * mg_eps,  # "H",
            0.0000000e-00 / 1.00000 * mg_eps,  # "C0",
        ],
        "mg_c6": [
            0.0022619536 / 1.27911 * mg_eps,  # "O",
            0.0022619536 / 1.72504 * mg_eps,  # "OM",
            0.0022619536 / 0.84961 * mg_eps,  # "OA",
            0.0024364096 / 0.63980 * mg_eps,  # "N",
            0.0024364096 / 0.29314 * mg_eps,  # "NT",
            0.0024364096 / 0.63980 * mg_eps,  # "NL",
            0.0024364096 / 0.43786 * mg_eps,  # "NR",
            0.0024364096 / 0.63980 * mg_eps,  # "NZ",
            0.0024364096 / 0.63980 * mg_eps,  # "NE",
            0.0023406244 / 0.27741 * mg_eps,  # "C",
            0.0023406244 / 0.27741 * mg_eps,  # "CH"
            0.0060684100 / 0.09489 * mg_eps_ch1,  # "CH1"
            0.0060684100 / 0.09489 * mg_eps_ch1,  # "CAH"
            0.0060684100 / 0.09489 * mg_eps_ch1,  # "CH1a"
            0.0074684164 / 0.41054 * mg_eps_ch2,  # "CH2"
            0.0096138025 / 0.86715 * mg_eps,  # "CH3"
            0.0073342096 / 0.47928 * mg_eps_ch2,  # "CH2r"
            0.0099840064 / 1.90587 * mg_eps,  # "S",
            0.0096138025 / 0.86715 * mg_eps,  # "CH3p"
            0.0147379600 / 2.44674 * mg_eps,  # "P",
            0.0022619536 / 1.05711 * mg_eps,  # "OE",
            0.0055130625 / 0.50266 * mg_eps,  # "CR1",
            0.0000000000 / 1.00000 * mg_eps,  # "H", # TODO
            0.0000000000 / 1.00000 * mg_eps,  # "C0",
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
    ("backbone_oxygen", "sidechain_cb", 0.2, None, 0),
    ("ct_oxygen", "sidechain_cb", 0.2, None, 0),
    ("backbone_nitrogen", "sidechain_cb", 0.65, None, -1),
    ("first_backbone_nitrogen", "backbone_nitrogen", None, 4.0e-6, 1),
    ("backbone_nitrogen", "backbone_nitrogen", 0.343, None, 1),
    ("backbone_carbonyl", "backbone_carbonyl", 0.5, None, -1),
    ("sidechain_cgs", "backbone_carbonyl", 0.078, None, 0),
    ("sidechain_cgs", "backbone_nitrogen", 0.087, None, 0),
    ("sidechain_cgs", "first_backbone_nitrogen", 0.087, None, 0),
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
