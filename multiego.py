from vanessa import float_range, check_files_existance
from Ensemble import read_simulations
from Multi_eGO_Ensemble import Multi_eGO_Ensemble
import argparse
import concurrent.futures
import sys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Metti una descrizione caruccina, tipo sul come nominare i file.')
    parser.add_argument('--system', type=str, required=True, help='Name of the system corresponding to system input folder.')
    parser.add_argument('--egos', choices=['rc', 'production'], required=True, help='Type of EGO.\n rc: creates a force-field for random coil simulations.\n production: creates a force-field combining random coil simulations and training simulations.')
    parser.add_argument('--epsilon', type=float_range(0.0, 1.0), help='Maximum interaction energy per contact.')
    # This is to use epsilon as default for inter molecular epsilon and ligand epsilon
    args, remaining = parser.parse_known_args()

    # Default arguments
    parser.add_argument('--train_from', nargs='+', type=str, default=[], help='A list of the simulations to be included in multi-eGO, corresponding to the subfolders to process and where the contacts are learned')
    parser.add_argument('--check_with', nargs='+', type=str, default=[], help='Those are contacts from a simulation or a structure used to check whether the contacts learned are compatible with the structures provided in here')
    parser.add_argument('--md_threshold', type=float, default=0.001, help='Minimum population for attractive interactions in training simulations.')
    parser.add_argument('--rc_threshold', type=float, default=0.000010, help='Minimum population for repulsive interactions in reference simulations.')
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

    print('Checking the presence of directories, .top, and .ndx files')
    md_ensembles_list = ['reference']+args.train_from+args.check_with
    check_files_existance(args.egos, args.system, md_ensembles_list)
    # TODO qui potrei aggiungere un print che mi dice tutte le cartelle che sta leggendo prima di farlo

    # Reading and preparing all the simulations defined in --train_from
    with concurrent.futures.ThreadPoolExecutor() as executor:
        # submit the read_subfolders function to be executed in parallel
        results = [executor.submit(read_simulations, args, simulation) for simulation in md_ensembles_list]
        ensembles = [result.result() for result in concurrent.futures.as_completed(results)]

    # Initializing Multi-eGO ensemble, which will gather all the Ensemble contact etc.
    print('- Initializing Multi-eGO ensemble')
    multiego_ensemble = Multi_eGO_Ensemble(args)

    for ensemble in ensembles:
        multiego_ensemble.add_ensemble_from(ensemble)

    multiego_ensemble.check_topology_conversion()
    multiego_ensemble.generate_bonded_interactions()
    multiego_ensemble.generate_LJ_potential()
    multiego_ensemble.write_model()
