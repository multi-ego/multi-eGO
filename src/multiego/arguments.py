args_dict = {
    "--system": {
        "type": str,
        "help": "Name of the system corresponding to system input folder.",
    },
    "--egos": {
        "type": str,
        "choices": ["rc", "production"],
        "help": "rc: creates a force-field for random coil simulations. "
        "production: creates a force-field combining random coil simulations and training simulations.",
    },
    "--epsilon": {
        "type": float,
        "help": "Maximum interaction energy per contact. The typical range is 0.2-0.4 kJ/mol",
    },
    "--reference": {
        "type": str,
        "default": "reference",
        "help": "The folder including all the reference information needed to setup multi-eGO, "
        "corresponding to the subfolder to process.",
    },
    "--train": {
        "type": lambda x: x.split(","),
        "default": [],
        "help": "A list of the training simulations to be included in multi-eGO, "
        "corresponding to the subfolders to process and where the contacts are learned.",
    },
    "--check": {
        "type": lambda x: [*x.split(",")],
        "default": [],
        "help": "A list of the simulations corresponding to the subfolders used to check "
        "whether the contacts learned are compatible with those provided in here.",
    },
    "--inter_epsilon": {
        "type": float,
        "help": "Maximum interaction energy per intermolecular contacts. The typical range is 0.2-0.4 kJ/mol",
    },
    "--inter_domain_epsilon": {
        "type": float,
        "help": "Maximum interaction energy per interdomain contacts. The typical range is 0.2-0.4 kJ/mol",
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
    "--multi_epsi_intra": {
        "type": str,
        "help": "Path to the input file specifying the intra epsilons",
    },
    "--multi_epsi_inter_domain": {
        "type": str,
        "help": "Path to the input file specifying the intra epsilons",
    },
    "--multi_epsi_inter": {
        "type": str,
        "help": "Path to the input file specifying the inter epsilons",
    },
    "--symmetry": {
        "default": "",
        "type": str,
        "help": "Symmetry file for the system",
    },
    "--f": {
        "default": 1,
        "type": float,
        "help": "partition function normalization",
    },
    "--inter_f": {
        "default": 1,
        "type": float,
        "help": "partition function normalization inter-molecular",
    },
    "--inter_domain_f": {
        "default": 1,
        "type": float,
        "help": "partition function normalization inter_domain",
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
