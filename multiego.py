import os
import pandas as pd
import sys, getopt
from read_input import find_files
from write_output import write_LJ, write_atomtypes_atp, write_topology, write_ligand_topology
from topology_parser import read_topology
from greta import ensemble, multiego_ensemble
pd.options.mode.chained_assignment = None  # default='warn'

def main(argv):

    parameters = {
        # native pair distance cut-off, used only when learning from structures
        'distance_cutoff':5.5,
        # neighbor aminoacid to exclude < x, used only when learning from structures
        'distance_residue':2,
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
        # Acid FFnonbondend it only works on the native pairs
        'acid_ff':False,
        # Default behavior is to train from a simulation
        'ensemble':True,
        # Does the model include the interaction with a ligand
        'ligand':False,
        # This is to reduce the kds when taking the ligand from another FF
        'ligand_reduction':6.75, # 2*1.5*1.5*1.5
        # This is the project name used to read inputs and write outputs
        'protein':None,
        # This is the force-field mode: random coil, single reference, multi reference
        'egos':None,
        # Input files folder
        'input_folder':None,
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

    try:
        opts, args = getopt.getopt(argv,"",["protein=", "md_ensembles=", "PDB1=", "PDB2=", "egos=", "epsilon1=", "epsilon2=", "ligand", "epsilon_ligand=", "noensemble", "help"])
        #opts, args = getopt.getopt(argv,"",["protein=", "egos=", "epsilon=", "epsilon_amyloid=", "ligand", "epsilon_ligand=", "noensemble", "help"])
    except getopt.GetoptError:
        print('multiego.py --MD1=path/to/folder --MD2=path/to/folder --egos=<all|split|rc> --epsilon=0.x (not used with --egos=rc) --ligand (optional) --epsilon_amyloid=0.x (optional) --epsilon_ligand=0.x (optional) --noensemble (optional)')
        sys.exit(2)
    if(len(opts)==0):
            print('multiego.py --protein=<protein> --egos=<single|merge|rc> --epsilon=0.x (not used with --egos=rc) --ligand (optional) --epsilon_amyloid=0.x (optional) --epsilon_ligand=0.x (optional) --noensemble (optional)')
            sys.exit()

    for opt, arg in opts:
        if opt == '--help':
            print('multiego.py --protein=<protein> --egos=<single|merge|rc> --epsilon=0.x (not used with --egos=rc) --ligand (optional) --epsilon_amyloid=0.x (optional) --epsilon_ligand=0.x (optional) --noensemble (optional)')
            sys.exit()

        elif opt in ("--protein"):
            if not arg:
                print('Provide a protein name')
                sys.exit()
            else:
                parameters['protein'] = arg

        elif opt in ("--egos"):
            if arg in ('all', 'split', 'rc'):
                parameters['egos'] = arg
            else:
                print('--egos accepts <all|split|rc> options')
                # TODO mettere una guida
                sys.exit()

        elif opt in ("--md_ensembles"):
            if not arg:
                print('Provide a path for the first MD simulation')
                sys.exit()
            md_ensembles_list = arg.split(',')

            for n, e in enumerate(md_ensembles_list, start=1):
                parameters[f'md_ensemble{n}'] = e

        elif opt in ("--epsilon"):
            arg = float(arg)
            if arg > 1 or arg < 0:
                # TODO mettere una guida
                print('Epsilon values must be chosen between 0 and 1')
                sys.exit()
            else:
                epsilon_list = arg.split(',')
                parameters['epsilon1'] = float(arg)
                parameters['epsilon2'] = float(arg)
                parameters['epsilon_ligand'] = float(arg)
                #readall +=1
        elif opt in ("--epsilon2"):
            # if set this overwrite the epsilon_md value
            arg = float(arg)
            if arg > 1 or arg < 0:
                print('Epsilon values must be chosen between 0 and 1')
                sys.exit()
            else:
                parameters['epsilon2'] = float(arg)
        elif opt in ("--ligand"):
            parameters['ligand'] = True
        elif opt in ("--epsilon_ligand"):
            # if set this overwrite the epsilon_md value
            arg = float(arg)
            if arg > 1 or arg < 0:
                print('Epsilon values must be chosen between 0 and 1')
            else:
                parameters['epsilon_ligand'] = float(arg)
        
        elif opt in ("--noensemble"):
            parameters['ensemble'] = False 
  
    # TODO figure out valid parameter combinations
    # check if input parameter combination is valid

    parameters['input_folder'] = f"inputs/{parameters['protein']}"

    # Folders to save the output files generated by the script
    if parameters['egos'] == 'rc':
        parameters['output_folder'] = f"outputs/{parameters['protein']}_{parameters['egos']}"
    else:
        epsilon_string = f"e{parameters['epsilon_md']}_{parameters['epsilon_amyl']}"
        if parameters['ligand'] == True:
            parameters['output_folder'] = f"outputs/{parameters['protein']}_{parameters['egos']}_{epsilon_string}_ligand{parameters['epsilon_ligand']}"
        else:
            parameters['output_folder'] = f"outputs/{parameters['protein']}_{parameters['egos']}_{epsilon_string}"

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
    

    print('- Adding the reference structure')
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


    print('- Creating md ensembles')
    for ensemble_name in md_ensembles_list:        
        file_paths = find_files(ensemble_name, parameters=parameters)
        md_ensemble = ensemble(parameters=parameters, ensemble_parameters=file_paths, name=ensemble_name)
        md_ensemble.prepare_ensemble()
        md_ensemble.assign_chains(reference_atoms_size)
        md_ensemble.add_MD_contacts()
        md_ensemble.convert_topology(reference)
        print(f'- The following contacts were converted: {md_ensemble.conversion_dict}')
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

    exit()

    #if parameters['egos'] != 'rc':
    #    print('- Adding Random Coil probability matrix to multi-eGO ensemble')
    #    # Multi-eGO always require the random coil probability
    #    ego_native.add_random_coil()
    #    multi_ego.add_structure_based_contacts(random_coil = ego_native.atomic_mat_random_coil)
    if parameters['egos'] != 'rc':

        if parameters['ensemble'] == True:
            print('- Ensemble = True: creating MD ensemble')         
            ego_md_parameters = {
                'topology_file':f"{parameters['input_folder']}/topol_md.top",
                'structure_file': f"{parameters['input_folder']}/native_md.pdb",
                'mdmat_contacts_file': f"{parameters['input_folder']}/plainMD_contacts.ndx",
            }

            ego_md = ensemble(parameters = parameters, ensemble_parameters=ego_md_parameters)
            ego_md.prepare_ensemble()
            ego_md.add_MD_contacts()
            ego_md.convert_topology(ego_native)
            print(f'- The following contacts were converted: {ego_md.conversion_dict}')
            print('- Adding MD probability matrix to multi-eGO ensemble')
            multi_ego.add_structure_based_contacts(atomic_mat_plainMD = ego_md.atomic_mat_MD)
        
        else:
            print('- Adding Structure-Based contact pairs to multi-eGO ensemble')
            # TODO ci potrebbe piacere il fatto di avere SB e MD insieme della nativa?
            ego_native.get_structure_pairs()
            multi_ego.add_structure_based_contacts(native_pairs=ego_native.structure_pairs)

        if parameters['egos'] == 'merge':
            print('- Merge = True: creating fibril ensemble')
            #TODO fix the fibril temp folder -> in general fix all the inputs with an automated script
            ego_fibril_parameters = {
                #'topology_file':f"{parameters['input_folder']}/fibril_2m5m/topol.top",
                'topology_file':f"{parameters['input_folder']}/fibril_temp/topol.top",
                'structure_file': f"{parameters['input_folder']}/fibril.pdb",
            }
            ego_fibril = ensemble(parameters = parameters, ensemble_parameters=ego_fibril_parameters)
            ego_fibril.prepare_ensemble()
            
            print('- Matching fibril topology to native topology')
            ego_fibril.match_native_topology(ego_native.sbtype_idx_dict)
            
            print('- Making fibril Structure-Based contact pairs')
            ego_fibril.get_structure_pairs(ego_native)
            
            print('- Adding fibril Structure-Based contact pairs to multi-eGO ensemble')
            multi_ego.add_structure_based_contacts(fibril_pairs = ego_fibril.structure_pairs)

        if parameters['ligand']:
            print('- Ligand = True: creating ligand ensemble')
            
            top = read_topology(f'{parameters["input_folder"]}/topol.top')
            ego_ligand_parameters = {
                'topology_file':f"{parameters['input_folder']}/topol_ligand.top",
                'structure_file': f"{parameters['input_folder']}/topol_native_ligand.pdb",
                'mdmat_contacts_file': f"{parameters['input_folder']}/ligandMD_contacts.ndx",
                'itp_file': f'{parameters["input_folder"]}/topol_ligand.itp',
                'prm_file': f'{parameters["input_folder"]}/topol_ligand.prm'
            }
            ego_ligand = ensemble(parameters=parameters, ensemble_parameters=ego_ligand_parameters)
            ego_ligand.prepare_ensemble()
            ego_ligand.add_MD_contacts()
            ego_ligand.convert_topology(ego_native)
            print(f'- The following contacts were converted: {ego_ligand.conversion_dict}')
            
            print('- Extracting ligand ensemble')
            ego_ligand.get_ligand_ensemble()
            ego_ligand.ligand_MD_LJ_pairs()
            
            #multi_ego.add_ensemble_top(ego_ligand)
            multi_ego.add_parsed_ligand_topology(ego_ligand)
            print('- Adding MD probability matrix to multi-eGO ensemble')
            multi_ego.add_structure_based_contacts(ligand_MD_pairs = ego_ligand.ligand_atomic_mat_MD)
        
    elif parameters['egos'] == 'rc':
        pass

    else: # one should never get here
        print("I dont' understand --egos=",parameters['egos'])
        exit()
    
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
