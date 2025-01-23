args_dict = {
    "--system": {
        "type": str,
        "help": "Name of the system corresponding to system input folder.",
    },
    "--egos": {
        "type": str,
        "choices": ["mg", "production"],
        "help": "mg: creates a force-field for molten globule simulations."
        "production: creates a force-field combining random coil simulations and training simulations.",
    },
    "--reference": {
        "type": lambda x: x.split(","),
        "default": [],
        "help": "The folder including all the reference information needed to setup multi-eGO, "
        "corresponding to the subfolder to process.",
    },
    "--train": {
        "type": lambda x: x.split(","),
        "default": [],
        "help": "A list of the training simulations to be included in multi-eGO, "
        "corresponding to the subfolders to process and where the contacts are learned.",
    },
    "--input_refs": {
        "type": list,
        "default": [],
        "help": "A list of the training simulations to be included in multi-eGO, "
        "corresponding to the subfolders to process and where the contacts are learned.",
    },
    "--epsilon": {
        "type": float,
        "help": "Maximum interaction energy per contact.",
    },
    "--p_to_learn": {
        "type": float,
        "default": 0.9995,
        "help": "Fraction of training simulations to learn.",
    },
    "--epsilon_min": {
        "type": float,
        "default": 0.07,
        "help": "The minimum meaningful epsilon value.",
    },
    "--force_split": {
        "type": bool,
        "default": False,
        "action": "store_true",
        "help": "Split inter and intra-molecular interactions in the ffnonbonded and topology files.",
    },
    "--single_molecule": {
        "type": bool,
        "default": False,
        "action": "store_true",
        "help": "Enable optimisations valid if you are simulating a single molecule.",
    },
    "--custom_dict": {
        "type": str,
        "help": "Custom dictionary for special molecules",
    },
    "--custom_c12": {
        "type": str,
        "help": "Custom dictionary of c12 for special molecules",
    },
    "--no_header": {
        "type": bool,
        "default": False,
        "action": "store_true",
        "help": "Removes headers from the output files when set",
    },
    "--symmetry": {
        "default": [],
        "type": list,
        "help": "Symmetry list for the system",
    },
    "--symmetry_file": {
        "default": "",
        "type": str,
        "help": "Symmetry file for the system",
    },
    "--f": {
        "default": 1,
        "type": float,
        "help": "partition function normalization",
    },
    "--relative_c12d": {
        "default": 0.01,
        "type": float,
        "help": "Relative deviation from default to set new replulsive c12",
    },
    "--explicit_name": {
        "default": "",
        "type": str,
        "help": "Explicit name for the output directory stored in outputs/system",
    },
    "--config": {
        "default": "",
        "type": str,
        "help": "Configuration file for the system",
    },
}

# Global arguments
args_dict_global = {
    "--system": {
        "type": str,
        "help": "Name of the system corresponding to system input folder.",
    },
    "--egos": {
        "type": str,
        "choices": ["mg", "production"],
        "help": "mg: creates a force-field for molten globule simulations."
        "production: creates a force-field combining random coil simulations and training simulations.",
    },
    "--input_refs": {
        "type": list,
        "default": [],
        "help": "A list of the training simulations to be included in multi-eGO, "
        "corresponding to the subfolders to process and where the contacts are learned.",
    },
    "--force_split": {
        "type": bool,
        "default": False,
        "action": "store_true",
        "help": "Split inter and intra-molecular interactions in the ffnonbonded and topology files.",
    },
    "--single_molecule": {
        "type": bool,
        "default": False,
        "action": "store_true",
        "help": "Enable optimisations valid if you are simulating a single molecule.",
    },
    "--custom_dict": {
        "type": str,
        "help": "Custom dictionary for special molecules",
    },
    "--custom_c12": {
        "type": str,
        "help": "Custom dictionary of c12 for special molecules",
    },
    "--no_header": {
        "type": bool,
        "default": False,
        "action": "store_true",
        "help": "Removes headers from the output files when set",
    },
    "--symmetry": {
        "default": [],
        "type": list,
        "help": "Symmetry list for the system",
    },
    "--symmetry_file": {
        "default": "",
        "type": str,
        "help": "Symmetry file for the system",
    },
    "--relative_c12d": {
        "default": 0.001,
        "type": float,
        "help": "Relative deviation from default to set new replulsive c12",
    },
    "--explicit_name": {
        "default": "",
        "type": str,
        "help": "Explicit name for the output directory stored in outputs/system",
    },
    "--config": {
        "default": "",
        "type": str,
        "help": "Configuration file for the system",
    },
}

# Arguments needed in each reference group to define reference, trainings, spcify the matrix, and other variables
args_dict_single_reference = {
    "--reference": {
        "type": str,
        "default": "",
        "help": "The folder including all the reference information needed to setup multi-eGO, "
        "corresponding to the subfolder to process.",
        "required": True,
    },
    "--train": {
        "type": lambda x: x.split(","),
        "default": [],
        "help": "A list of the training simulations to be included in multi-eGO, "
        "corresponding to the subfolders to process and where the contacts are learned.",
        "required": True,
    },
    "--matrix": {
        "type": str,
        "help": "Matrix name",
        "required": True,
    },
    "--epsilon": {
        "type": float,
        "help": "Maximum interaction energy per contact.",
        "required": True,
    },
    "--p_to_learn": {
        "type": float,
        "default": 0.9995,
        "help": "Fraction of training simulations to learn.",
    },
    "--epsilon_min": {
        "type": float,
        "default": 0.07,
        "help": "The minimum meaningful epsilon value.",
    },
    "--f": {
        "default": 1,
        "type": float,
        "help": "partition function normalization",
    },
}
