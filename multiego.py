import argparse
import sys
import os
import pandas as pd
import time
import gc

from src.multiego import ensemble
from src.multiego import io
from tools.face_generator import generate_face
from src.multiego.resources.type_definitions import parse_json
from src.multiego.arguments import args_dict
from src.multiego.arguments import args_dict_global
from src.multiego.arguments import args_dict_single_reference


def meGO_parsing():
    """
    Parses command-line arguments for the multi-eGO model generation.

    Returns:
    argparse.Namespace: An object containing parsed arguments.

    This function sets up an argument parser using the argparse library to handle command-line arguments
    required for generating a multi-eGO model based on training simulations and reference simulations.
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

    for arg, arg_dict in args_dict.items():
        # necessary for the boolean flags
        if "action" in arg_dict.keys() and (arg_dict["action"] == "store_true" or arg_dict["action"] == "store_false"):
            arg_dict.pop("type")  # necessary for boolean flags
        parser.add_argument(arg, **arg_dict)

    args, remaining = parser.parse_known_args()
    args.root_dir = os.path.dirname(os.path.abspath(__file__))
    # multi_flag = False

    # Check if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # read arguments depending on how they were parsed command line/config file (with or without multi-reference )
    args = io.read_arguments(args, args_dict, args_dict_global, args_dict_single_reference)

    # TODO put all the checks in a function to have the code cleaner
    # check if the configuration file is provided or if system, and egos rc are provided or if system, egos production, train and epsilon are provided
    if not args.system:
        print("ERROR: No system name found! Please provide a system name.")
        sys.exit()
    if not args.egos:
        print("ERROR: No egos mode found! Please provide an egos mode.")
        sys.exit()

    # controls that all reference entries contains correct arguments
    for ref in args.input_refs:

        # Check for missing required keys (keys without a default)
        input_keys = set(ref)
        # required_keys = {key.lstrip("--") for key, value in args_dict_single_reference.items() if value.get("required", False)}
        required_keys = {"matrix", "epsilon", "train", "reference"}
        missing_keys = required_keys - input_keys
        if missing_keys:
            raise ValueError(f"Missing required keys in {ref}: \n{missing_keys}")

        # Check for empty required values
        empty_required_keys = [key for key in required_keys if not ref[key]]  # Checks for None, empty string, or empty list
        if empty_required_keys:
            if "train" in empty_required_keys:
                print("ERROR: No training simulations found! Please provide a list of training simulations.")
                sys.exit()
            if "epsilon" in empty_required_keys:
                print("ERROR: No epsilon value found! Please provide an epsilon value.")
                sys.exit()

            raise ValueError(f"Empty values for required keys in {ref}.\n Missing {empty_required_keys}")

        # Checks and warnings for each reference
        if ref["p_to_learn"] < 0.9:
            print("WARNING: --p_to_learn should be large enough (suggested value is 0.9995)")

        if ref["epsilon_min"] <= 0.0:
            print(f"ERROR: --epsilon_min ({ref['epsilon_min']}) must be greater than 0.")
            sys.exit()

        if ref["epsilon"] < ref["epsilon_min"]:
            raise ValueError(
                f"ERROR: in {ref}. \nEpsilon value epsi = {ref['epsilon']} is below the minimum meaningful value."
            )

    custom_dict = {}
    if args.custom_dict:
        custom_dict = parse_json(args.custom_dict)
        if custom_dict == None:
            print("ERROR: Custom dictionary was parsed, but the dictionary is empty")
            sys.exit()

    print(f"Running Multi-eGO: {args.egos}\n")
    print("- Processing Multi-eGO topology")
    mego_ensemble = ensemble.init_meGO_ensemble(args, custom_dict)

    if args.symmetry_file and args.symmetry:
        print("ERROR: Both symmetry file and symmetry list provided. Please provide only one.")
        sys.exit()
    if args.symmetry_file:
        args.symmetry = io.read_symmetry_file(args.symmetry_file)
    elif args.symmetry:
        args.symmetry = io.parse_symmetry_list(args.symmetry)

    custom_c12_dict = pd.DataFrame()
    if args.custom_c12 is not None:
        custom_c12_dict = io.read_custom_c12_parameters(args.custom_c12)
        if custom_c12_dict is None or custom_c12_dict.empty:
            print("ERROR: Custom c12 paramter file was parsed, but the dictionary is empty")
            sys.exit()

    if remaining:
        print("Unknown arguments provided: " + str(remaining))
        parser.print_usage()
        sys.exit()

    return args, mego_ensemble, custom_dict


def main():
    """
    Parses command-line arguments and generates a multi-eGO model by invoking various functions
    related to ensemble generation, LJ parameter computation, and writing the output.
    """

    bt = time.time()
    generate_face.print_welcome()
    args, meGO_ensembles, custom_dict = meGO_parsing()

    st = time.time()
    elapsed_time = st - bt
    print("- Done in:", elapsed_time, "seconds")
    print("- Checking for input files and folders")
    io.check_files_existence(args)
    if args.egos == "production":
        io.check_matrix_format(args)

    print("\t- Generating bonded interactions")
    meGO_ensembles = ensemble.generate_bonded_interactions(meGO_ensembles)
    print("\t- Generating 1-4 data")
    pairs14, exclusion_bonds14 = ensemble.generate_14_data(meGO_ensembles)
    et = time.time()
    elapsed_time = et - st
    st = et
    print("- Done in:", elapsed_time, "seconds")

    if args.egos == "production":
        print("- Processing Multi-eGO contact matrices")
        meGO_ensembles, matrices = ensemble.init_meGO_matrices(meGO_ensembles, args, custom_dict)
        et = time.time()
        elapsed_time = et - st
        st = et
        print("- Done in:", elapsed_time, "seconds")
        print("- Initializing LJ dataset")
        train_dataset = ensemble.init_LJ_datasets(meGO_ensembles, matrices, pairs14, exclusion_bonds14, args)
        # force memory cleaning to decrease footprint in case of large dataset
        del matrices
        gc.collect()
        et = time.time()
        elapsed_time = et - st
        st = et
        print("- Done in:", elapsed_time, "seconds")
        print("- Generate LJ dataset")
        meGO_LJ, meGO_LJ_14 = ensemble.generate_LJ(meGO_ensembles, train_dataset, args)
        # force memory cleaning to decrease footprint in case of large dataset
        del train_dataset
        gc.collect()
        et = time.time()
        elapsed_time = et - st
        st = et
        print("- Done in:", elapsed_time, "seconds")
    elif args.egos == "mg":
        print("- Generate the LJ dataset")
        meGO_LJ = ensemble.generate_OO_LJ(meGO_ensembles)
        io.print_stats(meGO_LJ)
        meGO_LJ_14 = pairs14
        et = time.time()
        elapsed_time = et - st
        st = et
        print("- Done in:", elapsed_time, "seconds")

    print("- Finalize pairs and exclusions")
    meGO_LJ_14 = ensemble.make_pairs_exclusion_topology(meGO_ensembles, meGO_LJ_14, args)
    et = time.time()
    elapsed_time = et - st
    st = et
    print("- Done in:", elapsed_time, "seconds")

    print("- Writing Multi-eGO model")
    meGO_LJ = ensemble.sort_LJ(meGO_ensembles, meGO_LJ)
    io.write_model(meGO_ensembles, meGO_LJ, meGO_LJ_14, args)
    et = time.time()
    elapsed_time = et - st
    print("- Done in:", elapsed_time, "seconds")
    print("- Ran in:", et - bt, "seconds")

    generate_face.print_goodbye()


if __name__ == "__main__":
    main()
