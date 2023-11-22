import argparse
import sys
import os

from src.multiego import ensemble
from src.multiego import io
from src.multiego.util import float_range

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='TODO!')
    parser.add_argument('--system', type=str, required=True, help='Name of the system corresponding to system input folder.')
    parser.add_argument('--egos', choices=['rc', 'production'], required=True, help='Type of EGO.\n rc: creates a force-field for random coil simulations.\n production: creates a force-field combining random coil simulations and training simulations.')
    parser.add_argument('--epsilon', type=float, choices=[float_range.FloatRange(0.0, 1000.0)], help='Maximum interaction energy per contact.')
    parser.add_argument('--no_header', action='store_true', help='Removes headers from output_files when set')
    # This is to use epsilon as default for inter molecular epsilon and ligand epsilon
    args, remaining = parser.parse_known_args()

    # Default arguments
    parser.add_argument('--train_from', nargs='+', type=str, default=[], help='A list of the simulations to be included in multi-eGO, corresponding to the subfolders to process and where the contacts are learned')
    parser.add_argument('--check_with', nargs='+', type=str, default=[], help='Those are contacts from a simulation or a structure used to check whether the contacts learned are compatible with the structures provided in here')
    
    parser.add_argument('--p_to_learn', type=float, default=0.9995, help='Amount of the simulation to learn.')
    #parser.add_argument('--fraction', type=float, default=0.2, help='Minimum fraction of the maximum interaction energy per contact.')
    parser.add_argument('--epsilon_min', type=float, default=0.07, help='The minimum meaningfull epsilon value.')

    parser.add_argument('--inter_epsilon', type=float, default=args.epsilon, help='Maximum interaction energy per intermolecular contacts.')
    parser.add_argument('--inter_domain_epsilon', type=float, default=args.epsilon, help='Maximum interaction energy per interdomain contacts.')
    parser.add_argument('--out', type=str, default='', help='Suffix for the output directory name.')
    args = parser.parse_args()

    args.root_dir = os.path.dirname(os.path.abspath(__file__))

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

    if args.p_to_learn<0.9:
        print('WARNING: --p_to_learn should be high enough')

    #if args.fraction<0.1 or args.fraction>0.3:
    #    print('--fraction should be between 0.1 and 0.3')
    #    sys.exit()

    if args.egos != 'rc' and args.epsilon <= args.epsilon_min:
        print('--epsilon must be greater than --epsilon_min')
        sys.exit()
    
    if args.egos != 'rc' and args.inter_domain_epsilon <= args.epsilon_min:
        print('--inter_domain_epsilon must be greater than --epsilon_min')
        sys.exit()

    if args.egos != 'rc' and args.inter_epsilon <= args.epsilon_min:
        print('--inter_epsilon must be greater than --epsilon_min')
        sys.exit()

    if not os.path.exists(f'{args.root_dir}/outputs'): os.mkdir(f'{args.root_dir}/outputs')
    output_dir = io.create_output_directories(args)

    print('- Checking for input files and folders')
    md_ensembles_list = ['reference']+args.train_from+args.check_with
    io.check_files_existence(args.egos, args.system, args.root_dir, md_ensembles_list)

    # Initializing Multi-eGO ensemble, which will gather all the multiego.ensemble contact etc.
    print('- Initializing Multi-eGO ensemble')
    meGO_ensemble = ensemble.init_meGO_ensemble(args)
    meGO_ensemble = ensemble.generate_bonded_interactions(meGO_ensemble)

    print('- Generating the model')
    pairs14, exclusion_bonds14 = ensemble.generate_14_data(meGO_ensemble)
    if args.egos == 'rc':
        meGO_LJ = ensemble.generate_basic_LJ(meGO_ensemble)
        meGO_LJ_14 = pairs14
        meGO_LJ_14['epsilon'] = -meGO_LJ_14['c12']
    else:
       train_dataset, check_dataset = ensemble.init_LJ_datasets(meGO_ensemble, pairs14, exclusion_bonds14) 
       meGO_LJ, meGO_LJ_14  = ensemble.generate_LJ(meGO_ensemble, train_dataset, check_dataset, args)

    meGO_LJ_14 = ensemble.make_pairs_exclusion_topology(meGO_ensemble, meGO_LJ_14)
    
    io.write_model(meGO_ensemble, meGO_LJ, meGO_LJ_14, args, output_dir, args.out)
