import pandas as pd
import json
import sys

mg_OO_c12_rep = 7.5e-7
mg_OMOM_c12_rep = 2.5e-5
mg_HH_c12_rep = 1.2e-8
mg_ON_c12_rep = 7.5e-7
mg_NN_c12_rep = 2.5e-5

mg_HO_sigma = 0.169500
mg_eps_HO = 0.12

eps_O = 0.085
eps_OM = 0.085
eps_OA = 0.13
eps_N = 0.085
eps_NT = 0.085
eps_NL = 0.085
eps_NR = 0.07
eps_NZ = 0.11
eps_NE = 0.085
eps_C = 0.085
eps_CH = 0.11
eps_CH1 = 0.085
eps_CAH = 0.085
eps_CH2 = 0.11
eps_CAH2 = 0.11
eps_CH3 = 0.12
eps_CH2r = 0.11
eps_S = 0.14
eps_CH3p = 0.00
eps_P = 0.00
eps_OE = 0.00
eps_CR1 = 0.00
eps_H = 0.00
eps_C0 = 0.00

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
            2.5 * 0.262134**12,  # "O",   2.631580e-07
            2.5 * 0.253061**12,  # "OM",  1.724403e-07
            2.5 * 0.276621**12,  # "OA",  5.018430e-07
            2.5 * 0.289746**12,  # "N",   8.752940e-07
            2.5 * 0.317224**12,  # "NT",  2.596154e-06
            2.5 * 0.289746**12,  # "NL",  8.752940e-07
            2.5 * 0.303155**12,  # "NR",  1.506347e-06
            2.5 * 0.289746**12,  # "NZ",  8.752940e-07
            2.5 * 0.289746**12,  # "NE",  8.752940e-07
            2.5 * 0.317248**12,  # "C",   2.598570e-06
            2.5 * 0.317248**12,  # "CH",  2.598570e-06
            2.5 * 0.415167**12,  # "CH1", 6.555574e-05
            2.5 * 0.415167**12,  # "CAH", 6.555574e-05
            2.5 * 0.368035**12,  # "CH2", 1.543890e-05
            2.5 * 0.368035**12,  # "CAH2",1.543890e-05
            2.5 * 0.350505**12,  # "CH3", 8.595562e-06
            2.5 * 0.360236**12,  # "CH2r",1.193966e-05
            2.5 * 0.318498**12,  # "S",   2.724050e-06
            2.5 * 0.350981**12,  # "CH3p",8.736473e-06
            2.5 * 0.328121**12,  # "P",   3.893600e-06
            2.5 * 0.268811**12,  # "OE",  3.558824e-07
            2.5 * 0.341540**12,  # "CR1", 6.298560e-06
            2.5 * 0.163538**12,  # "H",   9.148590e-10
            2.5 * 0.262363**12,  # "C0",  2.659360e-07
        ],
        # all to ch2 except ch e ch3
        "mg_c12": [
            4.0 * 0.27601**12 * eps_O,  # "O",  sig=0.27601
            4.0 * 0.26259**12 * eps_OM,  # "OM", sig=0.26259
            4.0 * 0.29548**12 * eps_OA,  # "OA", sig=0.29548
            4.0 * 0.31365**12 * eps_N,  # "N",  sig=0.31365
            4.0 * 0.35722**12 * eps_NT,  # "NT", sig=0.35722
            4.0 * 0.31365**12 * eps_NL,  # "NL", sig=0.31365
            4.0 * 0.33411**12 * eps_NR,  # "NR", sig=0.33411
            4.0 * 0.31365**12 * eps_NZ,  # "NZ", sig=0.31365
            4.0 * 0.31365**12 * eps_NE,  # "NE", sig=0.31365
            4.0 * 0.35812**12 * eps_C,  # "C",  sig=0.35812
            4.0 * 0.35812**12 * eps_CH,  # "CH", sig=0.35812
            4.0 * 0.44592**12 * eps_CH1,  # "CH1",  sig=0.50192
            4.0 * 0.44592**12 * eps_CAH,  # "CAH",  sig=0.50192
            4.0 * 0.40704**12 * eps_CH2,  # "CH2",  sig=0.40704
            4.0 * 0.40704**12 * eps_CAH2,  # "CAH2", sig=0.40704
            4.0 * 0.37479**12 * eps_CH3,  # "CH3",  sig=0.37479
            4.0 * 0.39547**12 * eps_CH2r,  # "CH2r", sig=0.39547
            4.0 * 0.33077**12 * eps_S,  # "S",    sig=0.33077
            4.0 * 0.37479**12 * eps_CH3p,  # "CH3p", sig=0.37479
            4.0 * 0.33856**12 * eps_P,  # "P",    sig=0.33856
            4.0 * 0.28492**12 * eps_OE,  # "OE",   sig=0.28492
            4.0 * 0.37412**12 * eps_CR1,  # "CR1",  sig=0.37412
            4.0 * 0.000000000 * eps_H,  # "H",
            4.0 * 0.000000000 * eps_C0,  # "C0",
        ],
        "mg_c6": [
            4.0 * 0.27601**6 * eps_O,  # "O",
            4.0 * 0.26259**6 * eps_OM,  # "OM",
            4.0 * 0.29548**6 * eps_OA,  # "OA",
            4.0 * 0.31365**6 * eps_N,  # "N",
            4.0 * 0.35722**6 * eps_NT,  # "NT",
            4.0 * 0.31365**6 * eps_NL,  # "NL",
            4.0 * 0.33411**6 * eps_NR,  # "NR",
            4.0 * 0.31365**6 * eps_NZ,  # "NZ",
            4.0 * 0.31365**6 * eps_NE,  # "NE",
            4.0 * 0.35812**6 * eps_C,  # "C",
            4.0 * 0.35812**6 * eps_CH,  # "CH"
            4.0 * 0.44592**6 * eps_CH1,  # "CH1"
            4.0 * 0.44592**6 * eps_CAH,  # "CAH"
            4.0 * 0.40704**6 * eps_CH2,  # "CH2"
            4.0 * 0.40704**6 * eps_CAH2,  # "CAH2"
            4.0 * 0.37479**6 * eps_CH3,  # "CH3"
            4.0 * 0.39547**6 * eps_CH2r,  # "CH2r"
            4.0 * 0.33077**6 * eps_S,  # "S",
            4.0 * 0.37479**6 * eps_CH3p,  # "CH3p"
            4.0 * 0.33856**6 * eps_P,  # "P",
            4.0 * 0.28492**6 * eps_OE,  # "OE",
            4.0 * 0.37412**6 * eps_CR1,  # "CR1",
            4.0 * 0.00000000 * eps_H,  # "H",
            4.0 * 0.00000000 * eps_C0,  # "C0",
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


# Special local repulsion which can be asymmetric
# List of atom type combinations for LJ14 pairs
atom_type_combinations = [
    # Tuple of atom type combinations for LJ14 pairs
    #   ("backbone_carbonyl", "sidechain_cb", 0.275, 1.299682e-06, 1),
    #   ("backbone_oxygen", "sidechain_cb", 1, 1.5e-6, 0),
    #   ("ct_oxygen", "sidechain_cb", None, 5.0e-7, 0),
    ("backbone_nitrogen", "sidechain_cb", None, 5.0e-7, -1),
    ("first_backbone_nitrogen", "backbone_nitrogen", None, 4.0e-6, 1),
    #    ("backbone_nitrogen", "backbone_nitrogen", 0.343, None, 1),
    #    ("backbone_carbonyl", "backbone_carbonyl", 0.5, None, -1),
    #    ("sidechain_cgs", "backbone_carbonyl", 0.250, 1.2e-6, 0),
    #    ("sidechain_cgs", "backbone_nitrogen", 0.200, 5.5e-7, 0),
    #    ("sidechain_cgs", "first_backbone_nitrogen", 0.200, 5.5e-7, 0),
    #    ("sidechain_cds", "backbone_calpha", 0.100, 5e-7, 0),
]

# Special non-local interactions different from basic mg combination rules
# PROTEIN
polar_sbtype = ["O", "OA", "OM", "N", "NL", "NT", "NR", "NZ", "NE", "C", "S"]
hyd_sbtype = ["CH3", "CH2", "CH2r", "CH1", "CAH", "CAH2"]
special_non_local = [
    {
        "atomtypes": (["O"], ["O"]),  # charged oxygen-oxygen repulsion
        "interaction": "rep",
        "sigma": None,  # not needed for repulsion
        "epsilon": mg_OO_c12_rep,
    },
    {
        "atomtypes": (["O"], ["OM"]),  # charged oxygen-oxygen repulsion
        "interaction": "rep",
        "sigma": None,  # not needed for repulsion
        "epsilon": mg_OO_c12_rep,
    },
    {
        "atomtypes": (["OM"], ["OM"]),  # charged oxygen-oxygen repulsion
        "interaction": "rep",
        "sigma": None,  # not needed for repulsion
        "epsilon": mg_OMOM_c12_rep,
    },
    {
        "atomtypes": (["NL"], ["NL"]),  # charged nitrogen-nitrogen repulsion
        "interaction": "rep",
        "sigma": None,  # not needed for repulsion
        "epsilon": mg_NN_c12_rep,
    },
    {
        "atomtypes": (["H"], ["H"]),  # hydrogen-hydrogen repulsion
        "interaction": "rep",
        "sigma": None,  # not needed for repulsion
        "epsilon": mg_HH_c12_rep,
    },
    {
        "atomtypes": (["O", "OM", "OA"], ["H"]),  # hydrogen bond attraction
        "interaction": "att",
        "sigma": mg_HO_sigma,
        "epsilon": mg_eps_HO,
    },
    {
        "atomtypes": (["NL"], ["NZ", "N", "NT", "NR", "NE", "C", "CH", "CH1", "CAH", "CH2", "CAH2", "CH3", "CH2r", "S"]),
        "interaction": "rep",
        "sigma": None,
        "epsilon": None,
    },
    {
        "atomtypes": (["OM"], ["C", "CH", "CH1", "CAH", "CH2", "CAH2", "CH3", "CH2r", "S"]),
        "interaction": "rep",
        "sigma": None,
        "epsilon": None,
    },
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
