# Arguments shared between the single-reference CLI mode and the config-file (global) mode.
# Add any new global argument here only — args_dict and args_dict_global are derived from this.
_args_dict_shared = {
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

# Single-reference CLI mode: extends the shared base with the per-reference arguments
# that can be provided directly on the command line.
args_dict = {
    **_args_dict_shared,
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
    "--epsilon": {
        "type": float,
        "help": "Maximum interaction energy per contact.",
    },
}

# Config-file (global) mode: the shared base only.
# Per-reference arguments (--reference, --train, --epsilon) are supplied inside
# each entry of args_dict_single_reference instead.
args_dict_global = _args_dict_shared

# Arguments needed in each reference group to define reference, trainings,
# specify the matrix, and other per-reference variables.
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
}
