import multiego.ensemble
import multiego.io
import multiego.util.float_range

import argparse
import pandas as pd
import sys
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='TODO!')
    parser.add_argument('--system', type=str, required=True, help='Name of the system corresponding to system input folder.')
    parser.add_argument('--egos', choices=['rc', 'production'], required=True, help='Type of EGO.\n rc: creates a force-field for random coil simulations.\n production: creates a force-field combining random coil simulations and training simulations.')
    parser.add_argument('--epsilon', type=float, choices=[multiego.util.float_range.FloatRange(0.0, 1.0)], help='Maximum interaction energy per contact.')
    # This is to use epsilon as default for inter molecular epsilon and ligand epsilon
    args, remaining = parser.parse_known_args()

    # Default arguments
    parser.add_argument('--train_from', nargs='+', type=str, default=[], help='A list of the simulations to be included in multi-eGO, corresponding to the subfolders to process and where the contacts are learned')
    parser.add_argument('--check_with', nargs='+', type=str, default=[], help='Those are contacts from a simulation or a structure used to check whether the contacts learned are compatible with the structures provided in here')
    parser.add_argument('--md_threshold', type=float, default=0.001, help='Minimum population for attractive interactions in training simulations.')
    parser.add_argument('--rc_threshold', type=float, default=0.000020, help='Minimum population for repulsive interactions in reference simulations.')
    parser.add_argument('--inter_epsilon', type=float, default=args.epsilon, help='Maximum interaction energy per intermolecular contacts.')
    args = parser.parse_args()

    # checking the options provided in the commandline
    if args.egos != 'rc' and args.train_from is None:
        print('--egos=production require the definition of simulation folders containing the simulations to learn contacts from using --train_from flag')
        sys.exit()

    if args.egos == 'production' and not args.train_from:
        print('--egos=production requires the definition of the intramolecular and intermolecular ensembles by using --train_from')
        sys.exit()

    if args.epsilon is None and args.egos != 'rc':
        print('--epsilon is required when using --egos=production')
        sys.exit()

    if not os.path.exists('outputs'): os.mkdir('outputs')
    output_dir = multiego.io.create_output_directories(args)

    print('- Checking for input files and folders')
    md_ensembles_list = ['reference']+args.train_from+args.check_with
    multiego.io.check_files_existence(args.egos, args.system, md_ensembles_list)

    # Initializing Multi-eGO ensemble, which will gather all the multiego.ensemble contact etc.
    print('- Initializing Multi-eGO ensemble')
    meGO_ensemble = multiego.ensemble.init_meGO_ensemble(args.system, args.egos, args.train_from, args.check_with)
    meGO_ensemble = multiego.ensemble.generate_bonded_interactions(meGO_ensemble)

    print('- Generating the model')
    pairs14, exclusion_bonds14 = multiego.ensemble.generate_14_data(meGO_ensemble)
    if args.egos == 'rc':
        meGO_LJ = pd.DataFrame()
        meGO_LJ_14 = pairs14
        meGO_LJ_14['epsilon'] = - meGO_LJ_14['c12']
    else:
       train_dataset, check_dataset = multiego.ensemble.init_LJ_datasets(meGO_ensemble, pairs14, exclusion_bonds14) 
       meGO_LJ, meGO_LJ_14  = multiego.ensemble.generate_LJ(meGO_ensemble, train_dataset, check_dataset, args)

    meGO_LJ_14 = multiego.ensemble.make_pairs_exclusion_topology(meGO_ensemble, meGO_LJ_14)
    
    multiego.io.write_model(meGO_ensemble, meGO_LJ, meGO_LJ_14, args, output_dir)
