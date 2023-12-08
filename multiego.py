import argparse
import sys
import os

from src.multiego import ensemble
from src.multiego import io
from src.multiego.util import float_range
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
    --reference_from: The folder including all the reference information needed to setup multi-eGO, corresponding to the subfolder to process.
    --train_from: A list of the simulations to be included in multi-eGO, corresponding to the subfolders to process and where the contacts are learned.
    --check_with: Contacts from a simulation or a structure used to check whether the contacts learned are compatible with the structures provided.
    --out: Suffix for the output directory name.
    --inter_epsilon: Maximum interaction energy per intermolecular contacts.
    --inter_domain_epsilon: Maximum interaction energy per interdomain contacts.
    --p_to_learn: Amount of the simulation to learn.
    --epsilon_min: The minimum meaningful epsilon value.
    --no_header: Removes headers from output when set.
    """
    parser = argparse.ArgumentParser(description="Generate a multi-eGO model based on provided parameters.")
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
        choices=[float_range.FloatRange(0.0, 1000.0)],
        help="Maximum interaction energy per contact.",
    )
    optional_args.add_argument(
        "--reference_from",
        type=str,
        default="reference",
        help="""\
            The folder including all the reference information needed to setup multi-eGO,
            corresponding to the subfolder to process.
        """,
    )
    optional_args.add_argument(
        "--train_from",
        nargs="+",
        type=str,
        default=[],
        help="""\
            A list of the training simulations to be included in multi-eGO,
            corresponding to the subfolders to process and where the contacts are learned.
        """,
    )
    optional_args.add_argument(
        "--check_with",
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
        help="Maximum interaction energy per intermolecular contacts.",
    )
    optional_args.add_argument(
        "--inter_domain_epsilon",
        type=float,
        help="Maximum interaction energy per interdomain contacts.",
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
    if args.egos != "rc" and args.train_from is None:
        print(
            "--egos=production require the definition of simulation folders containing the simulations to learn contacts from using --train_from flag"
        )
        sys.exit()

    if args.egos == "production" and not args.train_from:
        print(
            "--egos=production requires the definition of the intramolecular and intermolecular ensembles by using --train_from"
        )
        sys.exit()

    if args.epsilon is None and args.egos != "rc":
        print("--epsilon is required when using --egos=production")
        sys.exit()

    if args.p_to_learn < 0.9:
        print("WARNING: --p_to_learn should be high enough")

    if args.egos != "rc" and args.epsilon <= args.epsilon_min:
        print("--epsilon must be greater than --epsilon_min")
        sys.exit()

    if args.egos != "rc" and args.inter_domain_epsilon <= args.epsilon_min:
        print("--inter_domain_epsilon must be greater than --epsilon_min")
        sys.exit()

    if args.egos != "rc" and args.inter_epsilon <= args.epsilon_min:
        print("--inter_epsilon must be greater than --epsilon_min")
        sys.exit()

    if not os.path.exists(f"{args.root_dir}/outputs"):
        os.mkdir(f"{args.root_dir}/outputs")

    if not args.no_header:
        generate_face.print_wellcome()

    output_dir = io.create_output_directories(args)

    print("- Checking for input files and folders")
    md_ensembles_list = ["reference"] + args.train_from + args.check_with
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

    io.write_model(meGO_ensemble, meGO_LJ, meGO_LJ_14, args, output_dir, args.out)

    generate_face.print_goodbye()


if __name__ == "__main__":
    main()
