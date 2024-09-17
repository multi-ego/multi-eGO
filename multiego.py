import argparse
import sys
import os
import parmed as pmd

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

    mego_topology = pmd.load_file(f"{args.root_dir}/inputs/{args.system}/topol.top")
    topol_names = [m for m in mego_topology.molecules]

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

    if args.egos != "rc" and not args.reference:
        args.reference = ["reference"]

    if args.epsilon and not args.inter_epsilon:
        args.inter_epsilon = args.epsilon
    if args.epsilon and not args.inter_domain_epsilon:
        args.inter_domain_epsilon = args.epsilon
    if not args.multi_epsilon_intra:
        args.multi_epsilon_intra = {k: v for k, v in zip(args.names, [args.epsilon] * len(args.names))}
    if not args.multi_epsilon_inter_domain:
        args.multi_epsilon_inter_domain = {k: v for k, v in zip(args.names, [args.inter_domain_epsilon] * len(args.names))}
    if not args.multi_epsilon_inter:
        args.multi_epsilon_inter = {k1: {k2: args.inter_epsilon for k2 in args.names} for k1 in args.names}

    # check all epsilons are set and greater than epsilon_min
    if args.egos != "rc":
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

    if args.custom_dict:
        custom_dict = parse_json(args.custom_dict)
        if custom_dict == None:
            print("WARNING: Custom dictionary was parsed, but the dictionary is empty")

    if remaining:
        print("Unknown arguments provided: " + str(remaining))
        parser.print_usage()
        sys.exit()

    return args


def get_meGO_LJ(meGO_ensemble, args):
    """
    This function generates Lennard-Jones (LJ) parameters for the multi-eGO ensemble based on the provided ensemble
    and command-line arguments.

    If the argument 'egos' is 'rc' (random coil), it generates basic LJ parameters for the ensemble.
    If 'egos' is 'production' or any other mode, it initializes LJ datasets, trains LJ parameters,
    and creates LJ pairs for 1-4 interactions.

    The resulting LJ parameters for 1-4 interactions are manipulated to get epsilon values,
    and a topology for exclusion pairs is created within the multi-eGO ensemble.
    The function returns two dataframes - meGO_LJ (LJ parameters) and meGO_LJ_14 (LJ parameters for 1-4 interactions).

    Parameters
    ----------
    meGO_ensemble : dict
        A dictionary containing the initialized multi-eGO ensemble.
    args : argparse.Namespace
        An object containing parsed arguments.

    Returns
    -------
    meGO_LJ : pandas.DataFrame
        A dataframe containing LJ parameters for the multi-eGO ensemble.
    meGO_LJ_14 : pandas.DataFrame
        A dataframe containing LJ parameters for 1-4 interactions in the multi-eGO ensemble.
    """
    pairs14, exclusion_bonds14 = ensemble.generate_14_data(meGO_ensemble)
    if args.egos == "rc":
        meGO_LJ = ensemble.generate_basic_LJ(meGO_ensemble, args)
        meGO_LJ_14 = pairs14
        meGO_LJ_14["epsilon"] = -meGO_LJ_14["c12"]
    else:
        train_dataset, check_dataset = ensemble.init_LJ_datasets(meGO_ensemble, pairs14, exclusion_bonds14, args)
        meGO_LJ, meGO_LJ_14 = ensemble.generate_LJ(meGO_ensemble, train_dataset, check_dataset, args)

    meGO_LJ_14 = ensemble.make_pairs_exclusion_topology(meGO_ensemble, meGO_LJ_14)

    return meGO_LJ, meGO_LJ_14


def main():
    """
    Parses command-line arguments and generates a multi-eGO model by invoking various functions
    related to ensemble generation, LJ parameter computation, and writing the output.
    """

    args = meGO_parsing()

    if not args.no_header:
        generate_face.print_welcome()

    print("- Checking for input files and folders")
    io.check_files_existence(args)

    print("- Initializing Multi-eGO model")
    meGO_ensembles = ensemble.init_meGO_ensemble(args)
    meGO_ensembles = ensemble.generate_bonded_interactions(meGO_ensembles)

    print("- Generating Multi-eGO model")
    meGO_LJ, meGO_LJ_14 = get_meGO_LJ(meGO_ensembles, args)

    print("- Writing Multi-eGO model")
    io.write_model(meGO_ensembles, meGO_LJ, meGO_LJ_14, args)

    generate_face.print_goodbye()


if __name__ == "__main__":
    main()
