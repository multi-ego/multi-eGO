import os
import pandas as pd
import sys, getopt
#from read_input_openai import search_and_read_csvs
from write_output import write_LJ, write_atomtypes_atp, write_topology, write_ligand_topology
from topology_parser import read_topology
from greta import ensemble, multiego_ensemble, convert_topology, from_ff_to_multiego
pd.options.mode.chained_assignment = None  # default='warn'



import argparse



# Non si legge più il pdb
# Non serve più la chain e il residue number per la lettura delle matrici


def main(argv):

    parameters = {
        # native pair distance cut-off, used only when learning from structures
        'distance_cutoff':5.5,
        # this is the minimum probability for a pair to be considered
        'md_threshold':0.001, 
        # this is the minimum probability for the random-coil matrix
        'rc_threshold':0.0001,
        # this is the interaction energy, by default this value is propagated to all following epsilons
        'epsilon_md':0.300,
        # This is the interaction energy of the amyloid cross beta
        'epsilon_amyl':0.300,
        # This is the protein-ligand interaction energy 
        'epsilon_ligand':0.300,
        # Does the model include the interaction with a ligand
        'ligand':False,
        # This is to reduce the kds when taking the ligand from another FF
        'ligand_reduction':6.75, # 2*1.5*1.5*1.5
        # This is the project name used to read inputs and write outputs
        'protein':None,
        # This is the force-field mode: random coil, single reference, multi reference
        'egos':None,
        # Input files folder
        'input_folders':None,
        # Output files folder
        'output_folder':None
    }

    print('\n\nMulti-eGO (codename: GRETA)\n')

    #readall=0
    md_ensembles_list = []
    pdb_ensembles_list = []

# TODO Vanessa: qui si fa un elenco degli ensemble in base alle flag.
# Ad esempio se abbiamo due MD sarebbe --MD=[MD1, MD2] oppure --MD=MD1 --MD=MD2
# Ciclare per ogni input in modo da fare un ensemble multiego in questo modo abbiamo due ensemble soliti ma possono esserne usati di piu',
# come diverse fibrille o informazioni parziali sui complessi.
# Dunque niente più egos.
# https://stackoverflow.com/questions/32761999/how-to-pass-an-entire-list-as-command-line-argument-in-python
# Questo vuol dire che c'è da rinominare le cartelle di input in qualche modo standard.


    parser = argparse.ArgumentParser(description='Multiego')

    # define the arguments
    parser.add_argument('--protein', type=str, required=True, help='Name of the directory corresponding to the protein')
    parser.add_argument('--md_ensembles', nargs='+', type=str, required=True, help='List of MD ensembles')
    parser.add_argument('--egos', choices=['all', 'split', 'rc'], required=True, help='Type of EGOs')
    parser.add_argument('--epsilon', type=float, help='Epsilon value')
    #parser.add_argument('--epsilon_amyloid', type=float, help='Epsilon value for amyloid')
    #parser.add_argument('--ligand', action='store_true', help='Include ligand')
    #parser.add_argument('--epsilon_ligand', type=float, help='Epsilon value for ligand')
    parser.add_argument('--intra', type=str, help='Intramolecular ensemble for split EGOs')
    parser.add_argument('--inter', type=str, help='Intermolecular ensemble for split EGOs')

    args = parser.parse_args()

    # checking the options provided in the commandline
    if args.egos == 'split' and not all(arg is not None for arg in (args.intra, args.inter)):
        print('--egos=split requires the definition of the intramolecular and intermolecular ensembles')
        sys.exit()

    if args.epsilon is None and args.egos != 'rc':
        print('--epsilon is required unless --egos=rc')
        sys.exit()




    print("Multi-eGO initializing")
    print(f"Reading folder {args.protein}")

    # The reference is always included. The starting conf.gro and topol.top are created with gmx pdb2gmx using multi-ego-basic.ff included.
    ensembles_to_read = ['reference']

    # Defining which ensemble will be added in multi-eGO
    if args.intra:
        ensembles_to_read.append(args.intra)
        print(f'Including intramolecular contacts from {args.intra}')
    if args.inter:
        ensembles_to_read.append(args.inter)
        print(f'Including intermolecular contacts from {args.inter}')

    print(f"The following ensembles will be included in multi-eGO force field: {ensembles_to_read}")
    
    # Here we insert the path for all the files needed. Luigi.py should follow the same nomenclature.

    # Adding topol.top path
    file_list = {}
    for ensemble in ensembles_to_read:
        # find file
        
        file_list[f'{ensemble}'] = {'topology':f"{args.protein}/{ensemble}/topol.top"}
    

    print(file_list)

    
    
    
    
    exit()



















    # We store all the arguments in a dictionary, but maybe this step can be removed
    parameters['protein'] = args.protein
    parameters['egos'] = args.egos
    parameters['epsilon_md'] = args.epsilon
    parameters['epsilon_amyl'] = args.epsilon_amyloid if args.epsilon_amyloid else args.epsilon
    #parameters['epsilon_ligand'] = args.epsilon_ligand if args.epsilon_ligand else args.epsilon
    #parameters['ligand'] = args.ligand
    parameters['intra'] = args.intra
    parameters['inter'] = args.inter
    parameters['input_folders'] = f"inputs/{args.protein}"









    # Definition of the matrices files to read obtained by gmx clustsize
    atomic_contact_matrices_to_read = [f"{parameters['input_folders']}/reference/"]

    # Folders to save the output files generated by the script
    if parameters['egos'] == 'rc':
        parameters['output_folder'] = f"outputs/{args.protein}_{parameters['egos']}"
    else:
        epsilon_string = f"e{parameters['epsilon_md']}_{parameters['epsilon_amyl']}"
        #if parameters['ligand'] == True:
        #    parameters['output_folder'] = f"outputs/{args.protein}_{parameters['egos']}_{epsilon_string}_ligand{parameters['epsilon_ligand']}"
        #else:
        parameters['output_folder'] = f"outputs/{args.protein}_{parameters['egos']}_{epsilon_string}"

    print('- Creating a multi-eGO force-field using the following parameters:')
    for k,v in parameters.items():
        if v == None: continue
        print('\t{:<20}: {:<20}'.format(k,v))
    
    try:
        os.mkdir(parameters['output_folder'])
    except OSError as error:
        pass



    # TODO ego_native diventa una inizializzazione di multi-eGO:
    # Bisogna partire da una topologia pdb2gmx con il nostro force field e poi verranno aggiunte cose

    # initialize multiego
    print('- Initializing multi-eGO')
    multi_ego = multiego_ensemble(parameters)

    


    print('- Reading all files')


    
    
    
    


    # Read di partenza per top e pdb di multiego
    














    file_paths = find_files(ensemble='reference', parameters=parameters)
    reference = ensemble(parameters=parameters, ensemble_parameters=file_paths, name='reference')
    reference.prepare_ensemble()
    reference.get_parsed_topology()
    
    if parameters['egos'] != 'rc':
        print('- Adding Random Coil probability matrix to multi-eGO ensemble')
        # Multi-eGO always require the random coil probability
        reference.add_random_coil()
        multi_ego.add_structure_based_contacts('random_coil', reference.atomic_mat_random_coil)
            
    multi_ego.add_ensemble_top(reference)
    multi_ego.add_parsed_topology(reference)
    
    # TODO delete reference saving the information of convert topology
    reference_atoms_size = reference.atoms_size
    #del reference

    if parameters['egos'] == 'all':
        print('\t- egos = all: intra and inter molecular contacts will be learned from all the ensembles')
    elif parameters['egos'] == 'split':
        print('\t- egos = split: intra and inter molecular contacts will be learned from the corresponding ensembles')
    elif parameters['egos'] == 'rc':
        print('\t- egos = rc: no molecular contacts will be parametrized')
    else:
        print("egos is not either 'all', 'split', or 'rc'")
        exit()

    print('- Creating md ensembles')
    for ensemble_name in md_ensembles_list: # TODO qui perche' per forza MD? gli posso dire di non leggere il nat-all e dovrebbe stare a posto e fare il PDB LJ     
        file_paths = find_files(ensemble_name, parameters=parameters)
        md_ensemble = ensemble(parameters=parameters, ensemble_parameters=file_paths, name=ensemble_name)
        md_ensemble.prepare_ensemble()
        convert_topology(from_ff_to_multiego, reference.ensemble_top, md_ensemble.ensemble_top)
        md_ensemble.add_MD_contacts()
        print('- Adding MD probability matrix to multi-eGO ensemble')
        multi_ego.add_structure_based_contacts(ensemble_name, md_ensemble.atomic_mat_MD)



    print('- Generating multi-eGO LJ')
    multi_ego.generate_multiego_LJ()

    print('- Generating pairs and exclusions for multi-eGO topology')
    multi_ego.generate_pairs_exclusions()
    print('- Generating writable')
    multi_ego.generate_outputs_toWrite()  

    print('- Writing multi-eGO Force-Field')
    write_atomtypes_atp(multi_ego)        
    write_LJ(multi_ego)
    write_topology(multi_ego)
    if parameters['ligand'] == True:
        write_ligand_topology(multi_ego)

    print('- Force-Field files saved in ' + parameters['output_folder'])
    print('\nGRETA completed! Carlo is happy!\t\^o^/\n')


if __name__ == "__main__":
   main(sys.argv[1:])
 
