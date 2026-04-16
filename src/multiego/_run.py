"""
Core pipeline logic shared by the root-level ``multiego.py`` launcher script
and the ``multiego`` console command installed by pip.

Neither caller should be imported directly; both invoke ``main(root_dir)``.
"""

import os
import sys
import gc
import time

from multiego import arguments
from multiego import contacts
from multiego import generate_face
from multiego import fileio as io
from multiego import lj
from multiego import mg
from multiego import pairs
from multiego import _term
from multiego.ensemble_data import MeGOEnsemble


def meGO_parsing(root_dir):
    """
    Parses and validates command-line arguments, resolves symmetry and custom
    dictionaries, and initializes the meGO ensemble topology.

    Parameters
    ----------
    root_dir : str
        Directory that contains the ``inputs/`` and ``outputs/`` folders.
        Pass ``os.path.dirname(os.path.abspath(__file__))`` when calling from
        the repo-root launcher script, or ``os.getcwd()`` when calling from
        the pip-installed console command.

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

    args.root_dir = root_dir
    if not args.inputs_dir:
        if args.config:
            # Config lives at {inputs_dir}/{system}/config.yml, so inputs_dir
            # is two levels above the config file.
            args.inputs_dir = os.path.dirname(os.path.dirname(os.path.abspath(args.config)))
        else:
            args.inputs_dir = os.path.join(root_dir, "inputs")
    if not args.outputs_dir:
        args.outputs_dir = os.path.join(root_dir, "outputs")
    args = arguments.read_arguments(
        args, arguments.args_dict, arguments.args_dict_global, arguments.args_dict_single_reference
    )

    arguments.validate_args(args)

    custom_dict = io.parse_json(args.custom_dict) if args.custom_dict else {}

    if args.symmetry_file:
        args.symmetry = io.read_symmetry_file(args.symmetry_file)
    elif args.symmetry:
        args.symmetry = io.parse_symmetry_list(args.symmetry)

    print(f"{_term._c(_term._BOLD)}Running Multi-eGO: {args.egos}{_term._c(_term._RESET)}\n")
    _term.header("Initialising Multi-eGO from")
    mego_ensemble = MeGOEnsemble.from_topology(args, custom_dict)

    return args, mego_ensemble, custom_dict


def main(root_dir):
    """
    Orchestrates the full multi-eGO model generation pipeline:
    argument parsing, bonded interactions, contact matrix processing,
    LJ parametrization, and output writing.

    Parameters
    ----------
    root_dir : str
        Directory that contains the ``inputs/`` and ``outputs/`` folders.
    """
    bt = time.time()
    generate_face.print_welcome()
    try:
        _main(root_dir, bt)
    except (ValueError, RuntimeError, FileNotFoundError) as e:
        _term.error(str(e))
        sys.exit(1)


def _main(root_dir, bt):
    args, meGO_ensembles, custom_dict = meGO_parsing(root_dir)

    st = time.time()
    _term.section_timing(st - bt)
    _term.header("Checking for input files and folders")
    io.check_files_existence(args)
    meGO_LJ_14 = None
    if args.egos == "production":
        io.check_matrix_format(args)
        _term.header("Processing Multi-eGO contact matrices")
        meGO_ensembles, matrices = contacts.init_meGO_matrices(meGO_ensembles, args, custom_dict)
        et = time.time()
        _term.section_timing(et - st)
        st = et

        _term.header("Initializing LJ dataset")
        with _term.spinner("building LJ dataset"):
            train_dataset = lj.init_LJ_datasets(meGO_ensembles, matrices, args)
        del matrices
        gc.collect()
        et = time.time()
        _term.section_timing(et - st)
        st = et

        _term.header("Generating LJ dataset")
        meGO_LJ, meGO_LJ_14, stat_str = lj.generate_LJ(meGO_ensembles, train_dataset, args)
        del train_dataset
        gc.collect()
        et = time.time()
        _term.section_timing(et - st)
        st = et

    elif args.egos == "mg":
        _term.header("Generating LJ dataset")
        with _term.spinner("building MG LJ"):
            meGO_LJ = mg.generate_MG_LJ(meGO_ensembles)
        stat_str = io.print_stats(meGO_LJ)
        et = time.time()
        _term.section_timing(et - st)
        st = et

    _term.header("Pairs and exclusions")
    with _term.spinner("computing pairs and exclusions"):
        meGO_LJ_14 = pairs.make_pairs_exclusion_topology(meGO_ensembles, args, meGO_LJ_14=meGO_LJ_14)
    et = time.time()
    _term.section_timing(et - st)
    st = et

    _term.header("Writing Multi-eGO model")
    with _term.spinner("writing files"):
        output_dir = io.write_model(meGO_ensembles, meGO_LJ, meGO_LJ_14, args, stat_str)
    _term.success(f"  Output written to {output_dir}")
    et = time.time()
    _term.section_timing(et - st)

    print()
    _term.rule(f"done in {et - bt:.2f} s")

    generate_face.print_goodbye()
