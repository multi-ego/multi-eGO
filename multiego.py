import argparse
import sys
import os
import time
import gc

from src.multiego import bonded
from src.multiego import contacts
from src.multiego import io
from src.multiego import lj
from src.multiego import mg
from src.multiego import pairs
from src.multiego import generate_face
from src.multiego.arguments import args_dict
from src.multiego.arguments import args_dict_global
from src.multiego.arguments import args_dict_single_reference


def _build_parser():
    """
    Constructs the argument parser for multi-eGO.

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

  1) generate a random coil prior model to generate the reference data for a single domain intramolecular interactions
     > python multiego.py --system GB1 --egos rc

  2) generate a production simulation using the reference data in the reference folder and the training data in the md_monomer folder
     interaction energy is set to 0.3 kJ/mol
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


def _validate_args(args):
    """
    Validates parsed arguments and exits with a clear message on any error.

    Parameters
    ----------
    args : argparse.Namespace
        Arguments as returned by io.read_arguments.
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


def meGO_parsing():
    """
    Parses and validates command-line arguments, resolves symmetry and custom
    dictionaries, and initializes the meGO ensemble topology.

    Returns
    -------
    args : argparse.Namespace
        Fully resolved and validated arguments.
    mego_ensemble : MeGOEnsemble
        Initialized ensemble topology.
    custom_dict : dict
        Custom atom-name mapping dictionary (empty if not provided).
    """
    parser = _build_parser()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args, remaining = parser.parse_known_args()

    if remaining:
        print("Unknown arguments provided: " + str(remaining))
        parser.print_usage()
        sys.exit()

    args.root_dir = os.path.dirname(os.path.abspath(__file__))
    args = io.read_arguments(args, args_dict, args_dict_global, args_dict_single_reference)

    _validate_args(args)

    custom_dict = {}
    if args.custom_dict:
        custom_dict = io.parse_json(args.custom_dict)
        if custom_dict is None:
            print("ERROR: Custom dictionary was parsed, but the dictionary is empty")
            sys.exit()

    if args.symmetry_file:
        args.symmetry = io.read_symmetry_file(args.symmetry_file)
    elif args.symmetry:
        args.symmetry = io.parse_symmetry_list(args.symmetry)

    print(f"Running Multi-eGO: {args.egos}\n")
    print("- Processing Multi-eGO topology")
    mego_ensemble = contacts.init_meGO_ensemble(args, custom_dict)

    return args, mego_ensemble, custom_dict


def main():
    """
    Orchestrates the full multi-eGO model generation pipeline:
    argument parsing, bonded interactions, contact matrix processing,
    LJ parametrization, and output writing.
    """
    bt = time.time()
    generate_face.print_welcome()
    args, meGO_ensembles, custom_dict = meGO_parsing()

    st = time.time()
    print("- Done in:", st - bt, "seconds")
    print("- Checking for input files and folders")
    io.check_files_existence(args)
    if args.egos == "production":
        io.check_matrix_format(args)

    print("- Generating bonded interactions")
    meGO_ensembles = bonded.generate_bonded_interactions(meGO_ensembles)
    et = time.time()
    print("- Done in:", et - st, "seconds")
    st = et

    print("- Generating 1-4 data")
    pairs14, exclusion_bonds14 = bonded.generate_14_data(meGO_ensembles)
    et = time.time()
    print("- Done in:", et - st, "seconds")
    st = et

    if args.egos == "production":
        print("- Processing Multi-eGO contact matrices")
        meGO_ensembles, matrices = contacts.init_meGO_matrices(meGO_ensembles, args, custom_dict)
        et = time.time()
        print("- Done in:", et - st, "seconds")
        st = et

        print("- Initializing LJ dataset")
        train_dataset = lj.init_LJ_datasets(meGO_ensembles, matrices, pairs14, exclusion_bonds14, args)
        del matrices
        gc.collect()
        et = time.time()
        print("- Done in:", et - st, "seconds")
        st = et

        print("- Generating LJ dataset")
        meGO_LJ, meGO_LJ_14, stat_str = lj.generate_LJ(meGO_ensembles, train_dataset, args)
        del train_dataset
        gc.collect()
        et = time.time()
        print("- Done in:", et - st, "seconds")
        st = et

    elif args.egos == "mg":
        print("- Generating LJ dataset")
        meGO_LJ = mg.generate_MG_LJ(meGO_ensembles)
        stat_str = io.print_stats(meGO_LJ)
        meGO_LJ_14 = pairs14
        et = time.time()
        print("- Done in:", et - st, "seconds")
        st = et

    print("- Finalizing pairs and exclusions")
    meGO_LJ_14 = pairs.make_pairs_exclusion_topology(meGO_ensembles, meGO_LJ_14, args)
    et = time.time()
    print("- Done in:", et - st, "seconds")
    st = et

    print("- Writing Multi-eGO model")
    meGO_LJ = lj.sort_LJ(meGO_ensembles, meGO_LJ)
    io.write_model(meGO_ensembles, meGO_LJ, meGO_LJ_14, args, stat_str)
    et = time.time()
    print("- Done in:", et - st, "seconds")
    print("- Ran in:", et - bt, "seconds")

    generate_face.print_goodbye()


if __name__ == "__main__":
    main()
