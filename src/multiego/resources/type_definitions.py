import pandas as pd
import json
import sys

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
            "CH2",
            "CH3",
            "CH2r",
            "S",
            "CH3p",
            "P",
            "OE",
            "CR1",
        ],
        "at.num": [8, 8, 8, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 16, 6, 15, 8, 6],
        "bare_c12": [
            2.631580e-07,
            1.724403e-07,
            5.018430e-07,
            8.752940e-07,
            2.596154e-06,
            8.752940e-07,
            1.506347e-06,
            8.752940e-07,
            8.752940e-07,
            2.598570e-06,
            2.598570e-06,
            6.555574e-05,
            1.543890e-05,
            8.595562e-06,
            1.193966e-05,
            2.724050e-06,
            8.736473e-06,
            3.893600e-06,
            3.558824e-07,
            6.298560e-06,
        ],
        "c12": [
            1e-06 / 1.27911,
            7.4149321e-07 / 1.72504,
            1.505529e-06 / 0.84961,
            2.319529e-06 / 0.63980,
            5.0625e-06 / 0.29314,
            2.319529e-06 / 0.63980,
            3.389281e-06 / 0.43786,
            2.319529e-06 / 0.63980,
            2.319529e-06 / 0.63980,
            4.937284e-06 / 0.27741,
            4.937284e-06 / 0.27741,
            9.70225e-05 / 0.09489,
            3.3965584e-05 / 0.41054,
            2.6646244e-05 / 0.86715,
            2.8058209e-05 / 0.47928,
            1.3075456e-05 / 1.90587,
            2.6646244e-05 / 0.86715,
            2.2193521e-05 / 2.44674,
            1.21e-06 / 1.05711,
            1.5116544e-05 / 0.50266,
        ],
        "c6": [
            0.0022619536 / 1.27911,
            0.0022619536 / 1.72504,
            0.0022619536 / 0.84961,
            0.0024364096 / 0.63980,
            0.0024364096 / 0.29314,
            0.0024364096 / 0.63980,
            0.0024364096 / 0.43786,
            0.0024364096 / 0.63980,
            0.0024364096 / 0.63980,
            0.0023406244 / 0.27741,
            0.0023406244 / 0.27741,
            0.00606841 / 0.09489,
            0.0074684164 / 0.41054,
            0.0096138025 / 0.86715,
            0.0073342096 / 0.47928,
            0.0099840064 / 1.90587,
            0.0096138025 / 0.86715,
            0.01473796 / 2.44674,
            0.0022619536 / 1.05711,
            0.0055130625 / 0.50266,
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
