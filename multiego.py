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
    multi_flag = False

    # Check if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if args.config:
        config_yaml = io.read_config(args.config, args_dict)
        # check if yaml file is empty
        if not config_yaml:
            print("WARNING: Configuration file was parsed, but the dictionary is empty")
        else:
            args = io.combine_configurations(config_yaml, args, args_dict)

    # check if the configuration file is provided or if system, and egos rc are provided or if system, egos production, train and epsilon are provided
    if not args.system:
        print("ERROR: No system name found! Please provide a system name.")
        sys.exit()
    if not args.egos:
        print("ERROR: No egos mode found! Please provide an egos mode.")
        sys.exit()
    if args.egos == "production" and not args.train:
        print("ERROR: No training simulations found! Please provide a list of training simulations.")
        sys.exit()
    if args.egos == "production" and not (
        args.epsilon or args.multi_epsilon_intra or args.multi_epsilon_inter or args.inter_epsilon
    ):
        print("ERROR: No epsilon value found! Please provide an epsilon value.")
        sys.exit()
    if args.p_to_learn < 0.9:
        print("WARNING: --p_to_learn should be large enough (suggested value is 0.9995)")
    if args.epsilon_min <= 0.0:
        print("--epsilon_min (" + str(args.epsilon_min) + ") must be greater than 0.")
        sys.exit()

    if args.multi_epsilon_intra or args.multi_epsilon_inter_domain or args.multi_epsilon_inter:
        multi_flag = True

    custom_dict = {}
    if args.custom_dict:
        custom_dict = parse_json(args.custom_dict)
        if custom_dict == None:
            print("ERROR: Custom dictionary was parsed, but the dictionary is empty")
            sys.exit()

    print(f"Running Multi-eGO: {args.egos}\n")
    print("- Processing Multi-eGO topology")
    mego_ensemble = ensemble.init_meGO_ensemble(args, custom_dict)
    topol_names = [m for m in mego_ensemble["topology"].molecules]

    args.names = []
    for name in args.multi_epsilon_intra.keys():
        args.names.append(name)
    for name in args.multi_epsilon_inter_domain.keys():
        args.names.append(name)
    for name in args.multi_epsilon_inter.keys():
        args.names.append(name)
        for name in args.multi_epsilon_inter[name].keys():
            args.names.append(name)
    args.names = list(set(args.names))
    if sorted(args.names) != sorted(topol_names) and multi_flag:
        print("ERROR: The names of the molecules in the topology and the multi-epsilon files are different")
        sys.exit()
    elif not multi_flag:
        args.names = topol_names

    if args.egos == "production" and not args.reference:
        args.reference = ["reference"]

    if args.epsilon and not args.inter_epsilon:
        args.inter_epsilon = args.epsilon
    if args.epsilon and not args.inter_domain_epsilon:
        args.inter_domain_epsilon = args.epsilon
    if not args.multi_epsilon_intra:
        args.multi_epsilon_intra = {k: v for k, v in zip(args.names, [args.epsilon] * len(args.names))}
    if not args.multi_epsilon_inter_domain and args.inter_domain_epsilon:
        args.multi_epsilon_inter_domain = {k: v for k, v in zip(args.names, [args.inter_domain_epsilon] * len(args.names))}
    if not args.multi_epsilon_inter_domain and not args.inter_domain_epsilon:
        args.multi_epsilon_inter_domain = args.multi_epsilon_intra
    if not args.multi_epsilon_inter and args.inter_epsilon:
        args.multi_epsilon_inter = {k1: {k2: args.inter_epsilon for k2 in args.names} for k1 in args.names}

    # check all epsilons are set and greater than epsilon_min
    if args.egos == "production":
        for k, v in args.multi_epsilon_intra.items():
            if v < args.epsilon_min:
                print("ERROR: epsilon value for " + k + " is less than epsilon_min")
                sys.exit()
        for k, v in args.multi_epsilon_inter_domain.items():
            if v < args.epsilon_min:
                print("ERROR: epsilon value for " + k + " is less than epsilon_min")
                sys.exit()
        for k1, v1 in args.multi_epsilon_inter.items():
            for k2, v2 in v1.items():
                if v2 < args.epsilon_min:
                    print("ERROR: epsilon value for " + k1 + "-" + k2 + " is less than epsilon_min")
                    sys.exit()

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
        meGO_LJ_14 = pairs14
        et = time.time()
        elapsed_time = et - st
        st = et
        print("- Done in:", elapsed_time, "seconds")
    else:
        print("- Generate the LJ dataset")
        meGO_LJ = ensemble.generate_OO_LJ(meGO_ensembles)
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
    io.write_model(meGO_ensembles, meGO_LJ, meGO_LJ_14, args)
    et = time.time()
    elapsed_time = et - st
    print("- Done in:", elapsed_time, "seconds")
    print("- Ran in:", et - bt, "seconds")

    generate_face.print_goodbye()


if __name__ == "__main__":
    main()
