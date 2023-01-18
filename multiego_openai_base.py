from vanessa import float_range, check_files_existance
from Ensemble import read_simulations
from Multi_eGO_Ensemble import Multi_eGO_Ensemble
import argparse
import concurrent.futures
import sys
import psutil


# Start tracking the CPU usage
psutil.cpu_percent(percpu=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Metti una descrizione caruccina, tipo sul come nominare i file.')
    parser.add_argument('--protein', type=str, required=True, help='Name of the proteina corresponding to the master folder containing subfolders.')
    parser.add_argument('--md_ensembles', nargs='+', type=str, help='A list of the simulations to be included in multi-eGO, corresponding to the subfolders to process.')
    parser.add_argument('--egos', choices=['all', 'split', 'rc'], required=True, help='Type of EGOs.')
    parser.add_argument('--epsilon', type=float_range(0.0, 1.0), help='Define a custom Epsilon value for the LJ parametrization from 0 to 1.')
    parser.add_argument('--intra', type=str, help='From the list of simulation defined in md_ensembles, choose from which one the intramolecular contacts are defined.')
    parser.add_argument('--inter', type=str, help='From the list of simulation defined in md_ensembles, choose from which one the intermolecular contacts are defined.')
    # This is to use epsilon as default for inter molecular epsilon and ligand epsilon
    args, remaining = parser.parse_known_args()

    # Default arguments
    parser.add_argument('--md_threshold', type=float, default=0.001, help='Contacts in intramat.ndx or intermat.ndx below this trehsold are dropped.')
    parser.add_argument('--rc_threshold', type=float, default=0.0001, help='Contacts in intramat.ndx or intermat.ndx below this trehsold are dropped.')
    parser.add_argument('--inter_epsilon', type=float, default=args.epsilon, help='Contacts in intramat.ndx or intermat.ndx below this trehsold are dropped.')
    args = parser.parse_args()

    # checking the options provided in the commandline
    if args.egos == 'split' and not all(arg is not None for arg in (args.intra, args.inter)):
        print('--egos=split requires the definition of the intramolecular and intermolecular ensembles')
        sys.exit()

    if args.egos != 'rc' and args.md_ensembles is None:
        print('--egos=split or egos=all requires the definition of simulation folder using --md_ensembles flag')
        sys.exit()

    if args.epsilon is None and args.egos != 'rc':
        print('--epsilon is required unless --egos=rc')
        sys.exit()

    print('Checking the presence of directories and .top, .ndx files')
    md_ensembles_list = ['reference']+args.md_ensembles
    check_files_existance(args.protein, md_ensembles_list)

    # Reading and preparing all the simulations defined in --md_ensembles
    with concurrent.futures.ThreadPoolExecutor() as executor:
        # submit the read_subfolders function to be executed in parallel
        results = [executor.submit(read_simulations, args, simulation) for simulation in md_ensembles_list]
        ensembles = [result.result() for result in concurrent.futures.as_completed(results)]

    # Initializing Multi-eGO ensemble, which will gather all the Ensemble contact etc.
    print('- Initializing Multi-eGO ensemble')
    multiego_ensemble = Multi_eGO_Ensemble(args)

    for ensemble in ensembles:
        multiego_ensemble.add_ensemble_from(ensemble)
    
    multiego_ensemble.generate_LJ_potential()




    # Get the CPU usage after the code has finished executing
    cpu_percent = psutil.cpu_percent(percpu=True)
    
    # Print the CPU usage
    print(f"CPU usage: {cpu_percent}%")
