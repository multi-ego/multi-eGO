import argparse
import sys
import os
import numpy as np

from src.multiego import ensemble
from src.multiego import io
from tools.face_generator import generate_face
from src.multiego.resources.type_definitions import parse_json


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
    # Required arguments
    required_args = parser.add_argument_group("Required arguments")
    required_args.add_argument(
        "--system",
        type=str,
        required=True,
        help="Name of the system corresponding to system input folder.",
    )
    required_args.add_argument(
        "--egos",
        choices=["rc", "production"],
        required=True,
        help="""\
            rc: creates a force-field for random coil simulations.
            production: creates a force-field combining random coil simulations and training simulations.
        """,
    )

    # Optional arguments
    optional_args = parser.add_argument_group("Optional arguments")
    optional_args.add_argument(
        "--epsilon",
        type=float,
        help="Maximum interaction energy per contact. The typical range is 0.2-0.4 kJ/mol",
    )
    optional_args.add_argument(
        "--reference",
        type=str,
        default="reference",
        help="""\
            The folder including all the reference information needed to setup multi-eGO,
            corresponding to the subfolder to process.
        """,
    )
    optional_args.add_argument(
        "--train",
        nargs="+",
        type=str,
        default=[],
        help="""\
            A list of the training simulations to be included in multi-eGO,
            corresponding to the subfolders to process and where the contacts are learned.
        """,
    )
    optional_args.add_argument(
        "--check",
        nargs="+",
        type=str,
        default=[],
        help="""\
            A list of the simulations corresponding to the subfolders used to check
            whether the contacts learned are compatible with those provided in here.
        """,
    )
    optional_args.add_argument("--out", type=str, default="", help="Suffix for the output directory name.")
    optional_args.add_argument(
        "--inter_epsilon",
        type=float,
        help="Maximum interaction energy per intermolecular contacts. The typical range is 0.2-0.4 kJ/mol",
    )
    optional_args.add_argument(
        "--inter_domain_epsilon",
        type=float,
        help="Maximum interaction energy per interdomain contacts. The typical range is 0.2-0.4 kJ/mol",
    )
    optional_args.add_argument(
        "--p_to_learn",
        type=float,
        default=0.9995,
        help="Fraction of training simulations to learn.",
    )
    optional_args.add_argument(
        "--epsilon_min",
        type=float,
        default=0.07,
        help="The minimum meaningful epsilon value.",
    )
    optional_args.add_argument(
        "--force_split",
        default=False,
        action="store_true",
        help="Split inter and intra-molecular interactions in the ffnonbonded and topology files.",
    )
    optional_args.add_argument(
        "--single_molecule",
        default=False,
        action="store_true",
        help="Enable optimisations valid if you are simulating a single molecule.",
    )
    optional_args.add_argument(
        "--custom_dict",
        type=str,
        help="Custom dictionary for special molecules",
    )
    optional_args.add_argument(
        "--custom_c12",
        type=str,
        help="Custom dictionary of c12 for special molecules",
    )
    optional_args.add_argument(
        "--no_header",
        action="store_true",
        help="Removes headers from the output files when set",
    )
    parser.add_argument("--multi_epsi_intra", type=str, help="Path to the input file specifying the intra epsilons")
    parser.add_argument("--multi_epsi_inter_domain", type=str, help="Path to the input file specifying the intra epsilons")
    parser.add_argument("--multi_epsi_inter", type=str, help="Path to the input file specifying the inter epsilons")
    optional_args.add_argument(
        "--symmetry",
        default="",
        type=str,
        help="Symmetry file for the system",
    )
    optional_args.add_argument(
        "--f",
        default=1,
        type=float,
        help="partition function normalization",
    )
    optional_args.add_argument(
        "--inter_f",
        default=1,
        type=float,
        help="partition function normalization inter-molecular",
    )
    optional_args.add_argument(
        "--inter_domain_f",
        default=1,
        type=float,
        help="partition function normalization inter_domain",
    )
    optional_args.add_argument(
        "--relative_c12d", default=0.001, type=float, help="Relative deviation from default to set new replulsive c12"
    )

    args, remaining = parser.parse_known_args()
    args.root_dir = os.path.dirname(os.path.abspath(__file__))

    # Set inter_epsilon default to epsilon if epsilon is provided
    if args.epsilon is not None and args.inter_epsilon is None:
        setattr(args, "inter_epsilon", args.epsilon)
    # Set inter_domain_epsilon default to epsilon if epsilon is provided
    if args.epsilon is not None and args.inter_domain_epsilon is None:
        setattr(args, "inter_domain_epsilon", args.epsilon)

    # checking the options provided in the commandline
    if args.egos != "rc" and not args.train:
        print("--egos=production requires the list of folders containing the training simulations using the --train flag")
        sys.exit()

    if (args.epsilon is None and args.multi_epsi_intra is None) and args.egos != "rc":
        print(
            "--epsilon or --multi_epsi_intra is required when using --egos=production. The typical range is between 0.2 and 0.4 kJ/mol"
        )
        sys.exit()

    if args.p_to_learn < 0.9:
        print("WARNING: --p_to_learn should be large enough (suggested value is 0.9995)")

    if args.egos != "rc" and args.epsilon_min <= 0.0:
        print("--epsilon_min (" + str(args.epsilon_min) + ") must be greater than 0.")
        sys.exit()

    # MULTI EPSILON CASES
    # TODO add the option to write in the multi epsi inter file Nones in order to remove interaction between systems
    # CHECK if either the single option or the multi option are provided. If both break
    if args.epsilon is not None and args.multi_epsi_intra is not None:
        print("""Choose either a single intra epsilon for the system or the multi-epsilon inter. Cannot choose both""")
        sys.exit()
    if args.inter_domain_epsilon is not None and args.multi_epsi_inter_domain is not None:
        print(
            "Choose either a single inter domain epsilon for the system or the multi-epsilon inter domain. Cannot choose both"
        )
        sys.exit()
    if args.inter_epsilon is not None and args.multi_epsi_inter is not None:
        print("Choose either a single inter epsilon for the system or the multi-epsilon inter. Cannot choose both")
        sys.exit()

    # CHECK if multi_epsi_inter_domain or multi_epsi_intra are parsed but not the multi_epsi_intra break
    if args.multi_epsi_intra is None and (args.multi_epsi_inter_domain is not None or args.multi_epsi_inter is not None):
        print(
            """--multi_epsi_inter_domain or --multi_epsi_inter where used, but --multi_epsi_intra was not parsed.
        In order to use the multi-epsilon option --multi_epsi_intra must be parsed. Please provide one or use the single epsslon case with:
            --epsilon
            --inter_domain_epsilon
            --inter_epsilon"""
        )
        exit()

    # if multi-epsi intra is parsed start overwrite other parameters
    if args.multi_epsi_intra is not None:
        setattr(args, "multi_mode", True)
        args.names, args.multi_epsilon = io.read_intra_file(args.multi_epsi_intra)

        # INTER-DOMAIN
        # multi_epsi inter domain
        if args.multi_epsi_inter_domain is not None:
            args.names_inter_domain, args.multi_epsilon_inter_domain = io.read_intra_file(args.multi_epsi_inter_domain)

        # multi_inter domain None but inter domain parsed --> ERROR
        if args.multi_epsi_inter_domain is None and args.inter_domain_epsilon is not None:
            print(
                """Inter domain option should be parsed with --multi_epsi_inter_domain if --multi_epsi_intra is used and not with --inter_domain_epsilon
            Choose either multiple epsilon options:
                   --multi_epsi_intra PATH_TO_FILE
                   --multi_epsi_inter_domain PATH_TO_FILE
            Or the single epsilon options:
                   --epsilon VALUE
                   --inter_domain_epsilon VALUE
                   """
            )
            exit()

        # CASE: multi intra but no multi inter domain --> set multi_inter_domain as multi_intra
        if args.multi_epsi_inter_domain is None and args.inter_domain_epsilon is None and args.multi_epsi_intra is not None:
            setattr(args, "names_inter_domain", args.names)
            setattr(args, "multi_epsilon_inter_domain", args.multi_epsilon)

        # INTER
        if args.multi_epsi_inter:
            args.names_inter, args.multi_epsilon_inter = io.read_inter_file(args.multi_epsi_inter)

        # No multi_epsilon_inter, no inter_epsilon --> set multi_epsilon_inter as one of the multi_epsi_intra (should not be needed if it's not defined explicetily)
        if args.multi_epsi_inter is None and args.inter_epsilon is None and args.multi_epsi_intra is not None:
            print(
                """--multi intra mode activated, but no information for inter epsilon was set.
Please set also the inter molecular interaction using one of the following options:
                  -inter_epsilon VALUE
                  -multi_epsi_inter PATH_TO_FILE """
            )
            exit()

        # No multi_epsilon_inter, inter_epsilon --> set multi_epsilon_inter as inter_epsilon
        if args.multi_epsi_inter is None and args.inter_epsilon is not None:
            setattr(args, "names_inter", np.array(args.names))
            setattr(
                args, "multi_epsilon_inter", np.zeros((len(args.multi_epsilon), len(args.multi_epsilon))) + args.inter_epsilon
            )

        # Multi-case checks:
        if args.multi_epsi_inter is not None and args.multi_epsi_intra is not None:
            if np.any(np.array(args.names) != np.array(args.names_inter)):
                print(
                    f"""ERROR: the names of the molecules in the files {args.multi_epsi_intra} and {args.multi_epsi_inter} are different.
                    The names of the molecules must be consistent with each other and with those in the topology"""
                )
                exit()

        # if multi_inter and no multi intra break
        if args.multi_epsi_inter is not None and args.multi_epsi_intra is None:
            print(
                """if multi_epsi_inter is used, also multi_epsi must be used. define also the set of epsilons via --multi_epsi_intra """
            )

        # if multi_inter_domain and no multi intra break
        if args.multi_epsi_inter_domain is not None and args.multi_epsi_intra is None:
            print(
                """--if multi_epsi_inter_domain is used, also multi_epsi must be used. define also the set of epsilons via --multi_epsi_intra """
            )
    else:
        setattr(args, "multi_mode", False)

    # CHECK all epsilons are greater than epsilon_min
    if args.epsilon is not None:
        if args.egos != "rc" and args.epsilon <= args.epsilon_min:
            print("--epsilon (" + str(args.epsilon) + ") must be greater than --epsilon_min (" + str(args.epsilon_min) + ")")
            sys.exit()

        if args.egos != "rc" and args.inter_domain_epsilon <= args.epsilon_min:
            print(
                "--inter_domain_epsilon ("
                + str(args.inter_domain_epsilon)
                + ") must be greater than --epsilon_min ("
                + str(args.epsilon_min)
                + ")"
            )
            sys.exit()

        if args.egos != "rc" and args.inter_epsilon <= args.epsilon_min:
            print(
                "--inter_epsilon ("
                + str(args.inter_epsilon)
                + ") must be greater than --epsilon_min ("
                + str(args.epsilon_min)
                + ")"
            )
            sys.exit()

    elif args.multi_mode is not None:
        if args.egos != "rc" and np.min(args.multi_epsilon) <= args.epsilon_min:
            print(
                f"all epsilons in {args.multi_epsi_intra} must be greater than --epsilon_min (" + str(args.epsilon_min) + ")"
            )
            sys.exit()

        if args.egos != "rc" and np.min(args.multi_epsilon_inter_domain) <= args.epsilon_min:
            print(
                f"all epsilons in {args.multi_epsi_inter_domain} must be greater than --epsilon_min ("
                + str(args.epsilon_min)
                + ")"
            )
            sys.exit()

        if args.egos != "rc" and np.min(args.multi_epsilon_inter) <= args.epsilon_min:
            print(
                f"all epsilons in {args.multi_epsi_inter} must be greater than --epsilon_min (" + str(args.epsilon_min) + ")"
            )
            sys.exit()

    if args.custom_dict:
        custom_dict = parse_json(args.custom_dict)
        if custom_dict == None:
            print("WARNING: Custom dictionary was parsed, but the dictionary is empty")

    if remaining:
        print("Unknown arguments provided: " + str(remaining))
        parser.print_usage()
        sys.exit()

    return args


def init_meGO_ensembles(args):
    """
    Initializes a multi-eGO ensemble based on the provided arguments.

    Args:
    args (argparse.Namespace): Parsed command-line arguments.

    Returns:
    meGO_ensemble: Initialized multi-eGO ensemble.

    This function initializes a multi-eGO ensemble by utilizing the provided arguments.
    It uses these arguments to create and configure the initial ensemble,
    generating bonded interactions within the ensemble.
    The resulting meGO_ensemble is returned for further processing.
    """
    meGO_ensemble = ensemble.init_meGO_ensemble(args)
    meGO_ensemble = ensemble.generate_bonded_interactions(meGO_ensemble)

    return meGO_ensemble


def get_meGO_LJ(meGO_ensemble, args):
    """
    Generates Lennard-Jones (LJ) parameters for the multi-eGO ensemble based on the provided ensemble and arguments.

    Args:
    meGO_ensemble: Initialized multi-eGO ensemble.
    args (argparse.Namespace): Parsed command-line arguments.

    Returns:
    tuple: A tuple containing two dataframes - meGO_LJ and meGO_LJ_14.

    This function generates Lennard-Jones (LJ) parameters for the multi-eGO ensemble based on the provided ensemble
    and command-line arguments.

    If the argument 'egos' is 'rc' (random coil), it generates basic LJ parameters for the ensemble.
    If 'egos' is 'production' or any other mode, it initializes LJ datasets, trains LJ parameters,
    and creates LJ pairs for 1-4 interactions.

    The resulting LJ parameters for 1-4 interactions are manipulated to get epsilon values,
    and a topology for exclusion pairs is created within the multi-eGO ensemble.
    The function returns two dataframes - meGO_LJ (LJ parameters) and meGO_LJ_14 (LJ parameters for 1-4 interactions).
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
    Main function that processes command-line arguments and generates a multi-eGO model.

    Parses command-line arguments and generates a multi-eGO model by invoking various functions
    related to ensemble generation, LJ parameter computation, and writing the output.

    Command-line Arguments:
    --system: Name of the system corresponding to the system input folder.
    --egos: Type of EGO. 'rc' for creating a force-field for random coil simulations,
            'production' for creating a force-field combining random coil simulations and training simulations.
    --epsilon: Maximum interaction energy per contact.
    --reference: The folder including all the reference information needed to setup multi-eGO, corresponding to the subfolder to process.
    --train: A list of the simulations to be included in multi-eGO, corresponding to the subfolders to process and where the contacts are learned.
    --check: Contacts from a simulation or a structure used to check whether the contacts learned are compatible with the structures provided.
    --out: Suffix for the output directory name.
    --inter_epsilon: Maximum interaction energy per intermolecular contacts.
    --inter_domain_epsilon: Maximum interaction energy per interdomain contacts.
    --p_to_learn: Amount of the simulation to learn.
    --epsilon_min: The minimum meaningful epsilon value.
    --no_header: Removes headers from output when set.
    """

    args = meGO_parsing()

    if not args.no_header:
        generate_face.print_welcome()

    print("- Checking for input files and folders")
    io.check_files_existence(args)

    print("- Initializing Multi-eGO model")
    meGO_ensembles = init_meGO_ensembles(args)

    print("- Generating Multi-eGO model")
    meGO_LJ, meGO_LJ_14 = get_meGO_LJ(meGO_ensembles, args)

    print("- Writing Multi-eGO model")
    io.write_model(meGO_ensembles, meGO_LJ, meGO_LJ_14, args)

    generate_face.print_goodbye()


if __name__ == "__main__":
    main()
