import os
import sys
import yaml
import argparse
from . import io

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
        "help": "mg: creates a force-field for molten globule simulations. "
        "production: creates a force-field combining reference simulations and training simulations.",
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
        "default": False,
        "action": "store_true",
        "help": "Split inter and intra-molecular interactions in the ffnonbonded and topology files.",
    },
    "--single_molecule": {
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
    "--learn_tolerance": {
        "default": 0.01,
        "type": float,
        "help": "Relative deviation from default to set new c6/c12",
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
    "--inputs_dir": {
        "default": None,
        "type": str,
        "help": "Directory containing the system input folders. "
        "Defaults to <root>/inputs when not specified. "
        "Not available in YAML config files.",
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


def build_parser():
    """
    Constructs and returns the argument parser for multi-eGO, registering
    all arguments from args_dict.

    Returns
    -------
    argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(
        prog="multiego.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""\
Generates a multi-eGO model based on one or more training simulations
and their corresponding reference simulations. In most cases one single
parameter is required, --epsilon, that sets the maximum interaction energy
for a contact pair.
""",
        epilog="""\
  example usage:

  1) generate a molten-globule prior model (reference data for intramolecular interactions)
     > python multiego.py --system GB1 --egos mg

  2) generate a production force-field using the reference data in the reference folder and
     the training data in the md_monomer folder; interaction energy is set to 0.3 kJ/mol
     > python multiego.py --system GB1 --egos production --train md_monomer --epsilon 0.3
""",
    )

    for arg, arg_spec in args_dict.items():
        # boolean flags must not carry a 'type' key when action is store_true/store_false
        spec = arg_spec.copy()
        if spec.get("action") in ("store_true", "store_false"):
            spec.pop("type", None)
        parser.add_argument(arg, **spec)

    return parser


def read_arguments(args, args_dict, args_dict_global, args_dict_single_reference):
    """
    Resolve final argument values from either a YAML config file or the
    command line, and normalise them into the ``input_refs`` list format
    expected by the rest of the pipeline.

    If ``args.config`` is set the YAML file takes precedence for any argument
    that was not explicitly overridden on the command line; otherwise the
    command-line values are converted into the same ``input_refs`` structure.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments as returned by ``build_parser().parse_args()``.
    args_dict : dict
        Full argument specification (single-reference CLI mode).
    args_dict_global : dict
        Argument specification for global (config-file) mode.
    args_dict_single_reference : dict
        Argument specification for per-reference blocks.

    Returns
    -------
    argparse.Namespace
        Resolved arguments with ``input_refs`` populated.
    """

    if args.config:

        config_yaml = read_config(args.config, args_dict)
        # check if yaml file is empty
        if not config_yaml:
            print("WARNING: Configuration file was parsed, but the dictionary is empty")
        else:
            args = combine_configurations(config_yaml, args, args_dict_global)
            # input_refs is a YAML-only block not handled by combine_configurations,
            # so we extract it directly from the raw YAML before passing to read_new_input
            args.input_refs = next(
                (e["input_refs"] for e in config_yaml if isinstance(e, dict) and "input_refs" in e),
                [],
            )
            args = read_new_input(args, args_dict_single_reference)

    # if config does not exist convert command line to input_ref dictionary format
    else:
        args.input_refs = []

        if args.egos == "production" and not args.reference:
            args.reference = ["reference"]

        args = convert_command_line_to_new_input(args, args_dict_single_reference)

    return args


def validate_args(args):
    """
    Validates parsed and fully resolved arguments, exiting with a clear
    message on any error. Should be called after read_arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments as returned by read_arguments.
    """
    if not args.system:
        print("ERROR: No system name found! Please provide a system name.")
        sys.exit()
    if not args.egos:
        print("ERROR: No egos mode found! Please provide an egos mode.")
        sys.exit()

    if args.p_to_learn < 0.9:
        print("WARNING: --p_to_learn should be large enough (suggested value is 0.9995)")

    if args.epsilon_min <= 0.0:
        print("ERROR: --epsilon_min must be greater than 0.")
        sys.exit()

    for ref in args.input_refs:
        required_keys = {"matrix", "epsilon", "train", "reference"}
        missing_keys = required_keys - set(ref)
        if missing_keys:
            raise ValueError(f"Missing required keys in {ref}:\n{missing_keys}")

        empty_required_keys = [key for key in required_keys if not ref[key]]
        if empty_required_keys:
            if "train" in empty_required_keys:
                print("ERROR: No training simulations found! Please provide a list of training simulations.")
                sys.exit()
            if "epsilon" in empty_required_keys:
                print("ERROR: No epsilon value found! Please provide an epsilon value.")
                sys.exit()
            raise ValueError(f"Empty values for required keys in {ref}.\nMissing {empty_required_keys}")

        if ref["epsilon"] < args.epsilon_min:
            print(f"ERROR: --epsilon ({ref['epsilon']}) must be greater-equal than --epsilon_min ({args.epsilon_min})")
            sys.exit()

    if args.symmetry_file and args.symmetry:
        print("ERROR: Both symmetry file and symmetry list provided. Please provide only one.")
        sys.exit()

    if args.custom_c12 is not None:
        custom_c12_dict = io.read_custom_c12_parameters(args.custom_c12)
        if custom_c12_dict is None or custom_c12_dict.empty:
            print("ERROR: Custom c12 parameter file was parsed, but the dictionary is empty")
            sys.exit()


def read_config(file, args_dict):
    """
    Reads a YAML file and returns its content as a dictionary.

    Parameters
    ----------
    file : str
        The path to the YAML file

    Returns
    -------
    list
        The parsed YAML content as a list of entries (scalars or single-key
        dicts), ready to be consumed by ``combine_configurations`` and
        ``read_new_input``.
    """
    with open(file, "r") as f:
        yml = yaml.safe_load(f)
    # check if the keys in the yaml file are valid
    # input_refs is a YAML-only key assembled by read_new_input, not an argparse argument
    yaml_only_keys = {"input_refs"}
    for element in yml:
        if not isinstance(element, dict):
            key = element
        else:
            key = list(element.keys())[0]
        if key in yaml_only_keys:
            continue
        if f"--{key}" not in args_dict:
            raise ValueError(f"ERROR: {key} in {file} is not a valid argument.")
    return yml


def read_new_input(args, args_dict_single_input):
    """
    Validate and normalise the ``input_refs`` list parsed from a YAML config.

    Each entry in ``args.input_refs`` is checked for unexpected keys and then
    merged with per-reference defaults from ``args_dict_single_input``.  Type
    coercion is applied to every value using the type callable from the spec.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments with ``args.input_refs`` populated from the YAML file.
        ``args.egos`` must be ``"production"``.
    args_dict_single_input : dict
        Per-reference argument specification (keys, types, and defaults).

    Returns
    -------
    argparse.Namespace
        The same ``args`` object with ``input_refs`` replaced by the
        validated and type-coerced list of per-reference dicts.
    """

    if args.egos != "production":
        raise ValueError("You should use 'input_refs' only with egos 'production'")

    # Check for invalid keys in input_dict
    input_refs = []
    valid_keys = {key.lstrip("--") for key in args_dict_single_input.keys()}
    for ref in args.input_refs:

        input_keys = set(ref.keys())

        # Raise an error if there are any unexpected keys
        unexpected_keys = input_keys - valid_keys
        if unexpected_keys:
            raise ValueError(f"Unexpected keys in {ref}: \n{unexpected_keys}")

        # Combine dictionaries with defaults
        combined_dict = {}
        for key, metadata in args_dict_single_input.items():
            stripped_key = key.lstrip("--")
            value = ref.get(stripped_key, metadata.get("default"))

            expected_type = metadata["type"]
            try:
                if value is not None:
                    value = expected_type(value)
                    ref[stripped_key] = value
            except (ValueError, TypeError):
                raise ValueError(
                    f"Invalid type for key '{stripped_key}'. Expected {expected_type}, got {type(value).__name__}."
                )

            combined_dict[stripped_key] = ref.get(stripped_key, metadata.get("default"))

        input_refs.append(combined_dict)

    args.input_refs = input_refs

    return args


def convert_command_line_to_new_input(args, args_dict_single_input):
    """
    Convert flat command-line arguments into the ``input_refs`` list format.

    Iterates over every reference folder in ``args.reference`` and every
    matrix file found inside it, building one per-reference dict per matrix.
    The result is stored in ``args.input_refs`` so the rest of the pipeline
    can treat CLI-mode and config-file-mode identically.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.  Must have ``system``, ``reference``,
        ``train``, ``epsilon``, and ``root_dir`` set.
    args_dict_single_input : dict
        Per-reference argument specification used to copy any additional
        per-reference keys from ``args`` into each entry.

    Returns
    -------
    argparse.Namespace
        The same ``args`` object with ``input_refs`` populated.
    """
    dict_input_ref = []
    appo = 0
    for reference in args.reference:

        matrices_intra = [m for m in os.listdir(f"{args.inputs_dir}/{args.system}/{reference}") if "intramat" in m]
        matrices_inter = [m for m in os.listdir(f"{args.inputs_dir}/{args.system}/{reference}") if "intermat" in m]
        matrices = matrices_intra + matrices_inter
        # check that in reference folder only 1 matrix of the same pair is present
        mat_type = [m.split(".ndx")[0] for m in matrices]
        if len(set(mat_type)) != len(mat_type):
            raise ValueError("In the reference folder, only one matrix of the same pair is allowed")

        for mat in mat_type:
            dict_input_ref.append({"reference": reference, "train": args.train, "matrix": mat, "epsilon": args.epsilon})
            for var in vars(args):
                if var in [key.lstrip("--") for key, _ in args_dict_single_input.items()]:
                    if var not in ["reference", "train", "epsilon"]:
                        dict_input_ref[appo].update({var: getattr(args, var)})
            appo += 1
    args.input_refs = dict_input_ref
    return args


def combine_configurations(yml, args, args_dict):
    """
    Combines the configuration from a YAML file with the command-line arguments. By overwriting
    values from the YAML configuration with the command-line arguments, the function ensures that
    the command-line arguments take precedence over the configuration file. Overwriting is done
    directly on the args namespace.

    Parameters
    ----------
    yml : dict
        The configuration from the YAML file
    args : dict
        The command-line arguments

    Returns
    -------
    dict
        The combined configuration
    """
    yaml_only_keys = {"input_refs"}
    for element in yml:
        if isinstance(element, dict):
            key, value = list(element.items())[0]
            if key in yaml_only_keys:
                continue
            value = args_dict[f"--{key}"]["type"](value)
            parse_key = f"--{key}"
            default_value = args_dict[parse_key]["default"] if "default" in args_dict[parse_key] else None

            if hasattr(args, key) and getattr(args, key) == default_value:
                setattr(args, key, value)
        else:
            if hasattr(args, element):
                setattr(args, element, True)

    return args
