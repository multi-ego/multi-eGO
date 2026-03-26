import sys
import os
import time
import gc

from src.multiego import arguments
from src.multiego import bonded
from src.multiego import contacts
from src.multiego import generate_face
from src.multiego import io
from src.multiego import lj
from src.multiego import mg
from src.multiego import pairs


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
    parser = arguments.build_parser()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args, remaining = parser.parse_known_args()

    if remaining:
        print("Unknown arguments provided: " + str(remaining))
        parser.print_usage()
        sys.exit()

    args.root_dir = os.path.dirname(os.path.abspath(__file__))
    args = arguments.read_arguments(
        args, arguments.args_dict, arguments.args_dict_global, arguments.args_dict_single_reference
    )

    arguments.validate_args(args)

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
    pairs14, all_bd = bonded.generate_14_data(meGO_ensembles)
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
        train_dataset = lj.init_LJ_datasets(meGO_ensembles, matrices, pairs14, all_bd, args)
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
