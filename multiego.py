import argparse
import sys
import os

from src.multiego import ensemble
from src.multiego import io
from tools.face_generator import generate_face


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
     in this case multiego expect the following input structure:
     └── input
         └── GB1
             └── reference
                 ├── topol.top
                 └── multi-eGO_basic.ff

  2) generate a production simulation using the reference data in the reference folder and the training data in the md_monomer folder
     interaction energy is set to 0.3 kJ/mol
     > python multiego.py --system GB1 --egos production --train md_monomer --epsilon 0.3
     in this case multiego expect the following input structure:
     └── input
         └── GB1
             ├── reference
             │   ├── topol.top
             │   ├── intramat_1_1.ndx
             │   └── multi-eGO_basic.ff
             └── md_monomer
                 ├── topol.top
                 ├── intramat_1_1.ndx
                 └── all-atom.ff

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
        "--no_header",
        action="store_true",
        help="Removes headers from the output files when set",
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
    if args.egos != "rc" and args.train is None:
        print("--egos=production requires the list of folders containing the training simulations using the --train flag")
        sys.exit()

    if args.epsilon is None and args.egos != "rc":
        print("--epsilon is required when using --egos=production. The typical range is between 0.2 and 0.4 kJ/mol")
        sys.exit()

    if args.p_to_learn < 0.9:
        print("WARNING: --p_to_learn should be large enough (suggested value is 0.9995)")

    if args.egos != "rc" and args.epsilon_min <= 0.0:
        print("--epsilon_min (" + str(args.epsilon_min) + ") must be greater than 0.")
        sys.exit()

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

    if remaining:
        print("Unknown arguments provided: " + str(remaining))
        parser.print_usage()
        sys.exit()

    if not args.no_header:
        generate_face.print_welcome()

    print("- Checking for input files and folders")
    md_ensembles_list = [args.reference] + args.train + args.check
    io.check_files_existence(args.egos, args.system, args.root_dir, md_ensembles_list)

    # Initializing Multi-eGO ensemble, which will gather all the multiego.ensemble contact etc.
    print("- Initializing Multi-eGO ensemble")
    meGO_ensemble = ensemble.init_meGO_ensemble(args)
    meGO_ensemble = ensemble.generate_bonded_interactions(meGO_ensemble)

    print("- Generating the model")
    pairs14, exclusion_bonds14 = ensemble.generate_14_data(meGO_ensemble)
    if args.egos == "rc":
        meGO_LJ = ensemble.generate_basic_LJ(meGO_ensemble)
        meGO_LJ_14 = pairs14
        meGO_LJ_14["epsilon"] = -meGO_LJ_14["c12"]
    else:
        train_dataset, check_dataset = ensemble.init_LJ_datasets(meGO_ensemble, pairs14, exclusion_bonds14)
        meGO_LJ, meGO_LJ_14 = ensemble.generate_LJ(meGO_ensemble, train_dataset, check_dataset, args)

    meGO_LJ_14 = ensemble.make_pairs_exclusion_topology(meGO_ensemble, meGO_LJ_14)

    print("- Writing Multi-eGO model")
    io.write_model(meGO_ensemble, meGO_LJ, meGO_LJ_14, args)

    generate_face.print_goodbye()


if __name__ == "__main__":
    main()
