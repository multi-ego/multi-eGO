import multiego.resources.type_definitions
import multiego.io
import multiego.topology

import glob
import pandas as pd
import parmed
import os
import numpy as np
import warnings

def assign_molecule_type(molecule_type_dict, molecule_name, molecule_topology):
    '''
    Decides if the molecule type of the system is 'protein', 'nucleic_acid' or 'other'
    and writes said information into molecule_type_dict before returning it

    Parameters
    ----------
    molecule_type_dict : dict
        Contains the molecule type information per system
    molecule_name : str
        The name of system
    molecule_topology : parmed.Topology
        The topology of the molecule, which will be used to figure out the molecule_type

    Returns
    -------
    molecule_type_dict : dict
        Updated molecule_type_dict with the added new system name
    '''
    first_aminoacid = molecule_topology.residues[0].name
    if first_aminoacid in multiego.resources.type_definitions.aminoacids_list:
        molecule_type_dict[molecule_name] = 'protein'
    elif first_aminoacid in multiego.resources.type_definitions.nucleic_acid_list:
        molecule_type_dict[molecule_name] = 'nucleic_acid'
    else:
        molecule_type_dict[molecule_name] = 'other'

    return molecule_type_dict


def initialize_topology(topology):
    '''
    '''
    # In a single topology different type of molecules can be present (e.g. protein, ligand).
    # For each molecule the different atomtypes are saved.
    columns_to_drop = ['nb_idx', 'solvent_radius', 'screen', 'occupancy', 'bfactor',
                       'altloc', 'join', 'irotat', 'rmin', 'rmin_14', 'epsilon_14', 'tree']
    ensemble_topology_dataframe, new_number, col_molecule, new_resnum, ensemble_molecules_idx_sbtype_dictionary, temp_number_c12_dict, temp_number_c6_dict = pd.DataFrame(), [], [], [], {}, {}, {}

    molecule_type_dict = {}
    # I needed to add this for loop as by creating the topology dataframe by looping over molecules, the c12 information is lost
    for atom in topology.atoms:
        temp_number_c12_dict[str(atom.idx+1)] = atom.epsilon*4.184
        temp_number_c6_dict[str(atom.idx+1)] = atom.sigma*0.1

    for molecule_number, (molecule_name, molecule_topology) in enumerate(topology.molecules.items(), 1):
        molecule_type_dict = assign_molecule_type(
            molecule_type_dict, molecule_name, molecule_topology[0])
        ensemble_molecules_idx_sbtype_dictionary[f'{str(molecule_number)}_{molecule_name}'] = {}
        ensemble_topology_dataframe = pd.concat(
            [ensemble_topology_dataframe, molecule_topology[0].to_dataframe()], axis=0)
        for atom in molecule_topology[0].atoms:
            new_number.append(str(atom.idx+1))
            col_molecule.append(f'{molecule_number}_{molecule_name}')
            new_resnum.append(str(atom.residue.number))

    ensemble_topology_dataframe['number'] = new_number
    ensemble_topology_dataframe['molecule'] = col_molecule
    ensemble_topology_dataframe['molecule_number'] = col_molecule
    ensemble_topology_dataframe[['molecule_number', 'molecule_name']] = ensemble_topology_dataframe.molecule.str.split('_', expand=True)
    ensemble_topology_dataframe['resnum'] = new_resnum
    ensemble_topology_dataframe['cgnr'] = ensemble_topology_dataframe['resnum']
    ensemble_topology_dataframe['ptype'] = 'A'
    ensemble_topology_dataframe = ensemble_topology_dataframe.replace({'name': multiego.resources.type_definitions.from_ff_to_multiego})
    ensemble_topology_dataframe['sb_type'] = ensemble_topology_dataframe['name'] + '_' + ensemble_topology_dataframe['molecule_name'] + '_' + ensemble_topology_dataframe['resnum'].astype(str)
    ensemble_topology_dataframe.rename(columns={'epsilon': 'c12'}, inplace=True)

    ensemble_topology_dataframe['charge'] = 0.
    ensemble_topology_dataframe['c6'] = ensemble_topology_dataframe['number'].map(temp_number_c6_dict)
    ensemble_topology_dataframe['c12'] = ensemble_topology_dataframe['number'].map(temp_number_c12_dict)
    ensemble_topology_dataframe['molecule_type'] = ensemble_topology_dataframe['molecule_name'].map(molecule_type_dict)

    for molecule in ensemble_molecules_idx_sbtype_dictionary.keys():
        temp_topology_dataframe = ensemble_topology_dataframe.loc[ensemble_topology_dataframe['molecule'] == molecule]
        number_sbtype_dict = temp_topology_dataframe[['number', 'sb_type']].set_index('number')['sb_type'].to_dict()
        ensemble_molecules_idx_sbtype_dictionary[molecule] = number_sbtype_dict
    sbtype_c12_dict = ensemble_topology_dataframe[['sb_type', 'c12']].set_index('sb_type')['c12'].to_dict()
    sbtype_name_dict = ensemble_topology_dataframe[['sb_type', 'name']].set_index('sb_type')['name'].to_dict()
    sbtype_moltype_dict = ensemble_topology_dataframe[['sb_type', 'molecule_type']].set_index('sb_type')['molecule_type'].to_dict()

    return ensemble_topology_dataframe, ensemble_molecules_idx_sbtype_dictionary, sbtype_c12_dict, sbtype_name_dict, sbtype_moltype_dict, molecule_type_dict


def initialize_molecular_contacts(contact_matrix, path, ensemble_molecules_idx_sbtype_dictionary, simulation):
    '''
    This function is called "initialize_molecular_contacts" and it takes three arguments:
     1) contact_matrices: a dictionary of contact matrices, where the keys are the file names (intramat_1_1.ndx) and the values are the contents of the files in the form of a pandas dataframe.
     2) ensemble_molecules_idx_sbtype_dictionary: a dictionary that associates the atom number with the structure-based type (sbtype) for each molecule in the ensemble.
     3) simulation: a string that represents the source of the simulation (e.g. "reference" or "native_MD").

    The function does the following:
     - Initializes an empty pandas dataframe called "ensemble_contact_matrix" that will be used to store the processed contact matrices.
     - Initializes an empty dictionary called "molecule_names_dictionary" that will be used to associate a molecule number with its name.
     - Loops through the keys of the ensemble_molecules_idx_sbtype_dictionary, and split the key by '_' and store the first element of the split as the type of contact matrix (e.g. "intra" or "inter"), second and third as the number of molecule.
     - Loops through the contact_matrices dictionary, and for each matrix:
        - Rename the column 'molecule_name_ai' and 'molecule_name_aj' by adding the name of the molecule to the 'molecule_number_ai' and 'molecule_number_aj' columns respectively.
        - Map the 'ai' and 'aj' columns, containing the atom number, to the corresponding sbtype from the ensemble_molecules_idx_sbtype_dictionary by using the name of the molecule.
        - If the file name starts with 'intramat' set the 'same_chain' column as True, if it starts with 'intermat' set it as False, otherwise print an error message and exit the script
        - Remove all the lines containing H atoms as the final model only contains heavy-atoms.
        - Concatenate all the dataframes contained in the simulation folder

    Parameters
    ----------
    contact_matrix : pd.DataFrame
        Contains the contact informations read from intra-/intermat
    path:

    ensemble_molecules_inx_sbtype_dictionary : dict
        Associates atom indices to atoms named according to multi-eGO conventions
    simulation : str
        The simulation classified equivalent to the input folder

    Returns
    -------
    contact_matrix : pd.DataFrame
        A contact matrix containing contact data for each of the different simulations
    '''
    print('\t\t-', f'Initializing {simulation} contact matrix')
    molecule_names_dictionary = {}
    for molecule_name in ensemble_molecules_idx_sbtype_dictionary.keys():
        name = molecule_name.split('_')
        molecule_names_dictionary[str(name[0])] = name[1]

    counter = 1
    name = path.split('/')[-1].split('_')
    # Renaming stuff
    contact_matrix['molecule_name_ai'] = contact_matrix['molecule_number_ai'].astype(str) + '_' + contact_matrix['molecule_number_ai'].map(molecule_names_dictionary)
    contact_matrix['molecule_name_aj'] = contact_matrix['molecule_number_aj'].astype(str) + '_' + contact_matrix['molecule_number_aj'].map(molecule_names_dictionary)
    contact_matrix['ai'] = contact_matrix['ai'].map(ensemble_molecules_idx_sbtype_dictionary[contact_matrix['molecule_name_ai'][0]])
    contact_matrix['aj'] = contact_matrix['aj'].map(ensemble_molecules_idx_sbtype_dictionary[contact_matrix['molecule_name_aj'][0]])

    contact_matrix = contact_matrix[~contact_matrix['ai'].astype(str).str.startswith('H')]
    contact_matrix = contact_matrix[~contact_matrix['aj'].astype(str).str.startswith('H')]

    contact_matrix = contact_matrix[['molecule_name_ai', 'ai', 'molecule_name_aj', 'aj', 'distance_m', 'distance', 'distance_14', 'probability', 'cutoff']]
    if name[0] == 'intramat': contact_matrix['same_chain'] = True
    elif name[0] == 'intermat': contact_matrix['same_chain'] = False
    else: raise Exception('There might be an error in the contact matrix naming. It must be intermat_X_X or intramat_X_X')

    contact_matrix['source'] = simulation
    contact_matrix['file'] = '_'.join(name)
    contact_matrix[['idx_ai', 'idx_aj']] = contact_matrix[['ai', 'aj']]
    contact_matrix.set_index(['idx_ai', 'idx_aj'], inplace=True)

    return contact_matrix


def init_meGO_ensemble(args):
    '''
    TODO
    '''
    # we initialize the reference topology
    reference_path = f'inputs/{args.system}/reference'
    ensemble_type = reference_path.split('/')[-1]
    print('\t-', f'Initializing {ensemble_type} ensemble topology')
    topology_path = f'{reference_path}/topol.top'
    if not os.path.isfile(topology_path): raise FileNotFoundError(f"{topology_path} not found.")
    print('\t-', f'Reading {topology_path}')
    # ignore the dihedral type overriding in parmed
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        reference_topology = parmed.load_file(topology_path)

    topology_dataframe, molecules_idx_sbtype_dictionary, sbtype_c12_dict, sbtype_name_dict, sbtype_moltype_dict, molecule_type_dict = initialize_topology(reference_topology)

    reference_contact_matrices = {}
    if args.egos != 'rc':
        matrix_paths = glob.glob(f'{reference_path}/int??mat_?_?.ndx')
        if matrix_paths == []: 
            raise FileNotFoundError('.ndx files must be named as intramat_X_X.ndx or intermat_1_1.ndx')
        for path in matrix_paths:
            name = path.replace('inputs/', '')
            name = name.replace('/', '_')
            name = name.replace('.ndx', '')
            reference_contact_matrices[name] = initialize_molecular_contacts(multiego.io.read_molecular_contacts(path), path, molecules_idx_sbtype_dictionary, 'reference')
            reference_contact_matrices[name] = reference_contact_matrices[name].add_prefix('rc_')

    ensemble = {}
    ensemble['topology'] = reference_topology
    ensemble['topology_dataframe'] = topology_dataframe
    ensemble['molecules_idx_sbtype_dictionary'] = molecules_idx_sbtype_dictionary
    ensemble['sbtype_c12_dict'] = sbtype_c12_dict
    ensemble['sbtype_name_dict'] = sbtype_name_dict
    ensemble['sbtype_moltype_dict'] = sbtype_moltype_dict
    ensemble['sbtype_number_dict'] = ensemble['topology_dataframe'][['sb_type', 'number']].set_index('sb_type')['number'].to_dict()
    ensemble['molecule_type_dict'] = molecule_type_dict
    ensemble['reference_matrices'] = reference_contact_matrices
    ensemble['train_matrix_tuples'] = []
    ensemble['check_matrix_tuples'] = []
  
    if args.egos == 'rc':
        return ensemble

    reference_set = set(ensemble['topology_dataframe']['name'].to_list())

    # now we process the train contact matrices
    train_contact_matrices = {}
    train_topology_dataframe = pd.DataFrame()
    for simulation in args.train_from: 
        print('\t-', f'Initializing {simulation} ensemble topology')
        simulation_path = f'inputs/{args.system}/{simulation}'
        topology_path = f'{simulation_path}/topol.top'
        if not os.path.isfile(topology_path): raise FileNotFoundError(f"{topology_path} not found.")
        print('\t-', f'Reading {topology_path}')
        # ignore the dihedral type overriding in parmed
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            topology = parmed.load_file(topology_path)

        print('\t-', f'{simulation} topology contains: {topology.molecules}')
        temp_topology_dataframe, molecules_idx_sbtype_dictionary, _, _, _, _ = initialize_topology(topology)
        train_topology_dataframe = pd.concat([train_topology_dataframe, temp_topology_dataframe], axis=0, ignore_index=True) 
        matrix_paths = glob.glob(f'{simulation_path}/int??mat_?_?.ndx')
        if matrix_paths == []: 
            raise FileNotFoundError('.ndx files must be named as intramat_X_X.ndx or intermat_1_1.ndx')
        for path in matrix_paths:
            name = path.replace('inputs/', '')
            name = name.replace('/', '_')
            name = name.replace('.ndx', '')
            train_contact_matrices[name] = initialize_molecular_contacts(multiego.io.read_molecular_contacts(path), path, molecules_idx_sbtype_dictionary, simulation)
            ref_name = reference_path+'_'+path.split('/')[-1]
            ref_name = ref_name.replace('inputs/', '')
            ref_name = ref_name.replace('/', '_')
            ref_name = ref_name.replace('.ndx', '')
            ensemble['train_matrix_tuples'].append((name, ref_name))

    ensemble['train_matrices'] = train_contact_matrices

    for number, molecule in enumerate(ensemble['topology'].molecules, 1):
        comparison_dataframe = train_topology_dataframe.loc[train_topology_dataframe['molecule'] == f'{number}_{molecule}']
        if not comparison_dataframe.empty:
            comparison_set = set(comparison_dataframe[~comparison_dataframe['name'].astype(str).str.startswith('H')]['name'].to_list())

    difference_set = comparison_set.difference(reference_set)
    if difference_set:
        print(f'The following atomtypes are not converted:\n{difference_set} \nYou MUST add them in "from_ff_to_multiego" dictionary to properly merge all the contacts.')
        exit()
 
    # now we process the check contact matrices
    check_contact_matrices = {}
    check_topology_dataframe = pd.DataFrame()
    for simulation in args.check_with: 
        print('\t-', f'Initializing {simulation} ensemble topology')
        simulation_path = f'inputs/{args.system}/{simulation}'
        topology_path = f'{simulation_path}/topol.top'
        if not os.path.isfile(topology_path): raise FileNotFoundError(f"{topology_path} not found.")
        print('\t-', f'Reading {topology_path}')
        # ignore the dihedral type overriding in parmed
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            topology = parmed.load_file(topology_path)

        print('\t-', f'{simulation} topology contains: {topology.molecules}')
 
        temp_topology_dataframe, molecules_idx_sbtype_dictionary, _, _, _, _ = initialize_topology(topology)
        check_topology_dataframe = pd.concat([check_topology_dataframe, temp_topology_dataframe], axis=0, ignore_index=True) 

        matrix_paths = glob.glob(f'{simulation_path}/int??mat_?_?.ndx')
        if matrix_paths == []: 
            raise FileNotFoundError('.ndx files must be named as intramat_X_X.ndx or intermat_1_1.ndx')
        for path in matrix_paths:
            name = path.replace('inputs/', '')
            name = name.replace('/', '_')
            name = name.replace('.ndx', '')
            check_contact_matrices[name] = initialize_molecular_contacts(multiego.io.read_molecular_contacts(path), path, molecules_idx_sbtype_dictionary, simulation)
            ref_name = reference_path+'_'+path.split('/')[-1]
            ref_name = ref_name.replace('inputs/', '')
            ref_name = ref_name.replace('/', '_')
            ref_name = ref_name.replace('.ndx', '')
            ensemble['check_matrix_tuples'].append((name, ref_name))

    ensemble['check_matrices'] = check_contact_matrices

    if not check_topology_dataframe.empty:
        for number, molecule in enumerate(ensemble['topology'].molecules, 1):
            comparison_dataframe = check_topology_dataframe.loc[check_topology_dataframe['molecule'] == f'{number}_{molecule}']
            if not comparison_dataframe.empty:
                comparison_set = set(comparison_dataframe[~comparison_dataframe['name'].astype(str).str.startswith('H')]['name'].to_list())
 
        difference_set = comparison_set.difference(reference_set)
        if difference_set:
            print(f'The following atomtypes are not converted:\n{difference_set} \nYou MUST add them in "from_ff_to_multiego" dictionary to properly merge all the contacts.')
            exit()

    return ensemble 
    

def generate_bonded_interactions(meGO_ensemble):
    '''
    Generates the bonded interactions and stores them in meGO_ensemble

    Parameters
    ----------
    meGO_ensemble : dict
        The meGO_ensemble object containing all the relevant system information
    
    Returns
    -------
    meGO_ensemble : dict
        The updated meGO_ensemble object with updated/added bonded parameters
    '''
    if 'meGO_bonded_interactions' not in meGO_ensemble.keys(): meGO_ensemble['meGO_bonded_interactions'] = {}
    if 'bond_pairs' not in meGO_ensemble.keys(): meGO_ensemble['bond_pairs'] = {}
    if 'user_pairs' not in meGO_ensemble.keys(): meGO_ensemble['user_pairs'] = {}

    for molecule, topol in meGO_ensemble['topology'].molecules.items():      
        meGO_ensemble['meGO_bonded_interactions'][molecule] = {
            'bonds' : multiego.topology.get_bonds(topol[0].bonds),
            'angles' : multiego.topology.get_angles(topol[0].angles),
            'dihedrals' : multiego.topology.get_dihedrals(topol[0].dihedrals),
            'impropers' : multiego.topology.get_impropers(topol[0].impropers),
            'pairs' : multiego.topology.get_pairs(topol[0].adjusts)
        }
        # The following bonds are used in the parametrization of LJ 1-4
        meGO_ensemble['bond_pairs'][molecule] = multiego.topology.get_bond_pairs(topol[0].bonds)
        meGO_ensemble['user_pairs'][molecule] = multiego.topology.get_pairs(topol[0].adjusts)

    return meGO_ensemble


def generate_14_data(meGO_ensemble):
    # First of all we generate the random-coil 1-4 interactions:
    pairs14 = pd.DataFrame()
    exclusion_bonds14 = pd.DataFrame()
    for molecule, bond_pair in meGO_ensemble['bond_pairs'].items():
        reduced_topology = meGO_ensemble['topology_dataframe'].loc[meGO_ensemble['topology_dataframe']['molecule_name'] == molecule][['number', 'sb_type', 'resnum', 'name', 'type', 'resname', 'molecule_type']].copy()

        reduced_topology['number'] = reduced_topology['number'].astype(str)
        reduced_topology['resnum'] = reduced_topology['resnum'].astype(int)
        # Dictionaries definitions to map values
        atnum_type_dict = reduced_topology.set_index('sb_type')['number'].to_dict()
        type_atnum_dict = reduced_topology.set_index('number')['sb_type'].to_dict()

        # Building the exclusion bonded list
        # exclusion_bonds are all the interactions within 3 bonds
        # p14 are specifically the interactions at exactly 3 bonds
        exclusion_bonds, tmp_p14 = multiego.topology.get_14_interaction_list(reduced_topology, bond_pair)
        # split->convert->remerge:
        tmp_ex = pd.DataFrame(columns = ['ai', 'aj', 'exclusion_bonds'])
        tmp_ex['exclusion_bonds'] = exclusion_bonds
        tmp_ex[['ai','aj']] = tmp_ex['exclusion_bonds'].str.split('_', expand=True)
        tmp_ex['ai'] = tmp_ex['ai'].map(type_atnum_dict)
        tmp_ex['aj'] = tmp_ex['aj'].map(type_atnum_dict)
        tmp_ex['1-4'] = '1_2_3' 
        tmp_ex['same_chain'] = True
        tmp_ex.loc[(tmp_ex['exclusion_bonds'].isin(tmp_p14)), '1-4'] = '1_4'
        exclusion_bonds14 = pd.concat([exclusion_bonds14, tmp_ex], axis=0, sort=False, ignore_index=True)

        # Adding the c12 for 1-4 interactions
        reduced_topology['c12'] = reduced_topology['sb_type'].map(meGO_ensemble['sbtype_c12_dict'])

        pairs = pd.DataFrame()
        if meGO_ensemble['molecule_type_dict'][molecule] == 'protein':
            pairs = multiego.topology.protein_LJ14(reduced_topology)
            pairs['ai'] = pairs['ai'].map(type_atnum_dict)
            pairs['aj'] = pairs['aj'].map(type_atnum_dict)
            pairs['rep'] = pairs['c12']
            pairs['same_chain'] = True
        else:
            pairs['ai'] = meGO_ensemble['user_pairs'][molecule].ai.astype(str)
            pairs['aj'] = meGO_ensemble['user_pairs'][molecule].aj.astype(str)
            pairs['ai'] = pairs['ai'].map(type_atnum_dict)
            pairs['aj'] = pairs['aj'].map(type_atnum_dict)
            nonprotein_c12 = []
            for test in meGO_ensemble['user_pairs'][molecule].type:
                nonprotein_c12.append(float(test.epsilon)*4.184)
            pairs['c12'] = nonprotein_c12
            pairs['rep'] = pairs['c12']
            pairs['same_chain'] = True

        pairs14 = pd.concat([pairs14, pairs], axis=0, sort=False, ignore_index=True)

    return pairs14, exclusion_bonds14


def init_LJ_datasets(meGO_ensemble, pairs14, exclusion_bonds14):

    # we cycle over train matrices to pair them with reference matrices and then we add 1-4 assignments and defaults c12s and concatenate everything
    train_dataset = pd.DataFrame()
    for (name, ref_name) in meGO_ensemble['train_matrix_tuples']:
        #sysname_train_from_intramat_1_1 <-> sysname_reference_intramat_1_1
        if not ref_name in meGO_ensemble['reference_matrices'].keys():
            raise RuntimeError('say something')

        temp_merged = pd.merge(meGO_ensemble['train_matrices'][name], meGO_ensemble['reference_matrices'][ref_name], left_index=True, right_index=True, how='outer')
        train_dataset = pd.concat([train_dataset, temp_merged], axis=0, sort = False, ignore_index = True)

    # This is to FLAG 1-2, 1-3, 1-4 cases:
    train_dataset = pd.merge(train_dataset, exclusion_bonds14[["ai", "aj", "same_chain", "1-4"]], how="left", on=["ai", "aj", "same_chain"])
    train_dataset.loc[(train_dataset['ai']==train_dataset['aj'])&(train_dataset['same_chain']==True), '1-4'] = '0'
    train_dataset['1-4'] = train_dataset['1-4'].fillna('1>4')
    # This is to set the correct default C12 values taking into account specialised 1-4 values (including the special 1-5 O-O)
    train_dataset = pd.merge(train_dataset, pairs14[["ai", "aj", "same_chain", "rep"]], how="left", on=["ai", "aj", "same_chain"])
    train_dataset.loc[(train_dataset['rep'].notna())&(train_dataset['rep']!=0.), '1-4'] = '1_4'
    train_dataset.loc[(train_dataset['1-4']=="0"), 'rep'] = 0.
    train_dataset.loc[(train_dataset['1-4']=="1_2_3"), 'rep'] = 0.
    train_dataset.loc[(train_dataset['1-4']=="1_4")&(train_dataset['rep'].isnull()), 'rep'] = 0.
    train_dataset['rep'] = train_dataset['rep'].fillna(np.sqrt(train_dataset['ai'].map(meGO_ensemble['sbtype_c12_dict'])*train_dataset['aj'].map(meGO_ensemble['sbtype_c12_dict'])))

    # we cycle over check matrices to pair them with reference matrices and then we add 1-4 assignments and defaults c12s and concatenate everything
    check_dataset = pd.DataFrame()
    for (name, ref_name) in meGO_ensemble['check_matrix_tuples']:
        #sysname_check_from_intramat_1_1 <-> sysname_reference_intramat_1_1
        if not ref_name in meGO_ensemble['reference_matrices'].keys():
            raise RuntimeError('say something')

        temp_merged = pd.merge(meGO_ensemble['check_matrices'][name], meGO_ensemble['reference_matrices'][ref_name], left_index=True, right_index=True, how='outer')
        check_dataset = pd.concat([check_dataset, temp_merged], axis=0, sort = False, ignore_index = True)

    if not check_dataset.empty:
        # This is to FLAG 1-2, 1-3, 1-4 cases:
        check_dataset = pd.merge(check_dataset, exclusion_bonds14[["ai", "aj", "same_chain", "1-4"]], how="left", on=["ai", "aj", "same_chain"])
        check_dataset.loc[(check_dataset['ai']==check_dataset['aj'])&(check_dataset['same_chain']==True), '1-4'] = '0'
        check_dataset['1-4'] = check_dataset['1-4'].fillna('1>4')
        # This is to set the correct default C12 values taking into account specialised 1-4 values (including the special 1-5 O-O)
        check_dataset = pd.merge(check_dataset, pairs14[["ai", "aj", "same_chain", "rep"]], how="left", on=["ai", "aj", "same_chain"])
        check_dataset.loc[(check_dataset['rep'].notna())&(check_dataset['rep']!=0.), '1-4'] = '1_4'
        check_dataset.loc[(check_dataset['1-4']=="0"), 'rep'] = 0.
        check_dataset.loc[(check_dataset['1-4']=="1_2_3"), 'rep'] = 0.
        check_dataset.loc[(check_dataset['1-4']=="1_4")&(check_dataset['rep'].isnull()), 'rep'] = 0.
        check_dataset['rep'] = check_dataset['rep'].fillna(np.sqrt(check_dataset['ai'].map(meGO_ensemble['sbtype_c12_dict'])*check_dataset['aj'].map(meGO_ensemble['sbtype_c12_dict'])))

    return train_dataset, check_dataset


def generate_LJ(meGO_ensemble, train_dataset, check_dataset, parameters):
    '''
    This function reads the probabilities obtained using gmx_clustsize from the ensembles defined in the command line.
    The random coil probabilities are used to reweight the explicit water ones.
    Intra and inter molecular contacts are splitted as different rules are applied during the reweighting.
    For each atom contact the sigma and epsilon are obtained.

    Parameters
    ----------
    meGO_ensemble : dict
        Contains the relevant meGO data such as interactions and statistics
    parameters : dict
        Contains the command-line parsed parameters

    Returns
    -------
    meGO_atomic_contacts_merged : pd.DataFrame
        Contains the non-bonded atomic contacts associated to LJ parameters and statistics
    meGO_LJ_14 : pd.DataFrame
        Contains the 1-4 (paris-exclusions) atomic contacts associated to LJ parameters and statistics
    '''

    # this is the minimum probability ratio for attractive contacts and is set so that epsilon should not be smaller than 0.1*epsilon
    limit_rc = 1./parameters.rc_threshold**0.1
    # This keep only significat attractive/repulsive interactions
    meGO_LJ = train_dataset.loc[(train_dataset['probability']>parameters.md_threshold)|((train_dataset['probability']<=parameters.md_threshold)&(train_dataset['rc_probability']>limit_rc*np.maximum(train_dataset['probability'],parameters.rc_threshold)))].copy()

    # Add sigma, add epsilon reweighted, add c6 and c12
    meGO_LJ['sigma'] = (meGO_LJ['distance']) / (2.**(1./6.))
    meGO_LJ['epsilon'] = np.nan 

    # The index has been reset as here I have issues with multiple index duplicates. The same contact is kept twice: one for intra and one for inter.
    # The following pandas functions cannot handle multiple rows with the same index although it has been defined the "same_chain" filter.
    meGO_LJ.reset_index(inplace=True)

    # Epsilon reweight based on probability
    # Paissoni Equation 2.1
    # Attractive intramolecular
    meGO_LJ.loc[(meGO_LJ['probability']>limit_rc*np.maximum(meGO_LJ['rc_probability'],parameters.rc_threshold))&(meGO_LJ['same_chain']==True), 'epsilon'] = -(parameters.epsilon/np.log(parameters.rc_threshold))*(np.log(meGO_LJ['probability']/np.maximum(meGO_LJ['rc_probability'],parameters.rc_threshold)))
    # Attractive intermolecular
    meGO_LJ.loc[(meGO_LJ['probability']>limit_rc*np.maximum(meGO_LJ['rc_probability'],parameters.rc_threshold))&(meGO_LJ['same_chain']==False), 'epsilon'] = -(parameters.inter_epsilon/np.log(parameters.rc_threshold))*(np.log(meGO_LJ['probability']/np.maximum(meGO_LJ['rc_probability'],parameters.rc_threshold)))

    # Probability only Repulsive intramolecular
    meGO_LJ.loc[(np.maximum(meGO_LJ['probability'],parameters.rc_threshold)<1./limit_rc*meGO_LJ['rc_probability'])&(meGO_LJ['same_chain']==True)&(meGO_LJ['probability']<=parameters.md_threshold), 'epsilon'] = -(parameters.epsilon/np.log(parameters.rc_threshold))*meGO_LJ['cutoff']**12*np.log(np.maximum(meGO_LJ['probability'],parameters.rc_threshold)/meGO_LJ['rc_probability'])-meGO_LJ['rep']
    # Probability only Repulsive intermolecular
    meGO_LJ.loc[(np.maximum(meGO_LJ['probability'],parameters.rc_threshold)<1./limit_rc*meGO_LJ['rc_probability'])&(meGO_LJ['same_chain']==False)&(meGO_LJ['probability']<=parameters.md_threshold), 'epsilon'] = -(parameters.inter_epsilon/np.log(parameters.rc_threshold))*meGO_LJ['cutoff']**12*np.log(np.maximum(meGO_LJ['probability'],parameters.rc_threshold)/meGO_LJ['rc_probability'])-meGO_LJ['rep']

    # Full Repulsive intramolecular
    meGO_LJ.loc[(meGO_LJ['probability']<1./limit_rc*np.maximum(meGO_LJ['rc_probability'],parameters.rc_threshold))&(meGO_LJ['probability']>parameters.md_threshold)&(meGO_LJ['rc_probability']>parameters.md_threshold)&(meGO_LJ['same_chain']==True)&((meGO_LJ['rc_distance_m']-meGO_LJ['distance_m'])<0.), 'epsilon'] = -(parameters.epsilon/np.log(parameters.rc_threshold))*meGO_LJ['distance_m']**12*np.log(meGO_LJ['probability']/meGO_LJ['rc_probability'])-(meGO_LJ['rep']*(meGO_LJ['distance_m']**12/meGO_LJ['rc_distance_m']**12))
    # Full Repulsive intermolecular
    meGO_LJ.loc[(meGO_LJ['probability']<1./limit_rc*np.maximum(meGO_LJ['rc_probability'],parameters.rc_threshold))&(meGO_LJ['probability']>parameters.md_threshold)&(meGO_LJ['rc_probability']>parameters.md_threshold)&(meGO_LJ['same_chain']==False)&((meGO_LJ['rc_distance_m']-meGO_LJ['distance_m'])<0.), 'epsilon'] = -(parameters.inter_epsilon/np.log(parameters.rc_threshold))*meGO_LJ['distance_m']**12*np.log(meGO_LJ['probability']/meGO_LJ['rc_probability'])-(meGO_LJ['rep']*(meGO_LJ['distance_m']**12/meGO_LJ['rc_distance_m']**12))

    # mild c12s (for small probability ratios): smaller
    meGO_LJ.loc[(meGO_LJ['probability']<=limit_rc*np.maximum(meGO_LJ['rc_probability'],parameters.rc_threshold))&(meGO_LJ['probability']>=meGO_LJ['rc_probability'])&((meGO_LJ['rc_distance_m']-meGO_LJ['distance_m'])>0.)&(meGO_LJ['probability']>parameters.md_threshold)&(meGO_LJ['rc_probability']>parameters.md_threshold), 'epsilon'] = -meGO_LJ['rep']*(meGO_LJ['distance_m']/meGO_LJ['rc_distance_m'])**12 
    # mild c12s (for small probability ratios): larger 
    meGO_LJ.loc[(meGO_LJ['probability']>=1./limit_rc*np.maximum(meGO_LJ['rc_probability'],parameters.rc_threshold))&(meGO_LJ['probability']<=meGO_LJ['rc_probability'])&((meGO_LJ['rc_distance_m']-meGO_LJ['distance_m'])<0.)&(meGO_LJ['rc_probability']>parameters.md_threshold)&(meGO_LJ['probability']>parameters.md_threshold), 'epsilon'] = -meGO_LJ['rep']*(meGO_LJ['distance_m']/meGO_LJ['rc_distance_m'])**12

    # This set the default (random coil) value for selected 1-4 interactions
    meGO_LJ.loc[(meGO_LJ['1-4'] == '1_4'), 'epsilon'] = -meGO_LJ['rep']
    # Rescale c12 1-4 interactions 
    meGO_LJ.loc[(np.abs(meGO_LJ['rc_distance_14']-meGO_LJ['distance_14'])>0.)&(meGO_LJ['rc_probability']>parameters.md_threshold)&(meGO_LJ['1-4']=="1_4")&(meGO_LJ['same_chain']==True), 'epsilon'] = -meGO_LJ['rep']*(meGO_LJ['distance_14']/meGO_LJ['rc_distance_14'])**12

    # where the probability is less than md_threshold we do not trust the sigma estimate and so we set it to cutoff, so that then this low probability contact are
    # better accounted for when merging contact from multiple simulations
    meGO_LJ.loc[(meGO_LJ['probability']<parameters.md_threshold), 'sigma'] = meGO_LJ['cutoff'] / (2.**(1./6.))

    # Here we are reindexing like before
    meGO_LJ[['idx_ai', 'idx_aj']] = meGO_LJ[['ai', 'aj']]
    meGO_LJ.set_index(['idx_ai', 'idx_aj'], inplace=True)

    # clean NaN and zeros
    meGO_LJ.dropna(subset=['epsilon'], inplace=True)
    meGO_LJ = meGO_LJ[meGO_LJ.epsilon != 0]

    # remove unnecessary fields
    meGO_LJ = meGO_LJ[['molecule_name_ai', 'ai', 'molecule_name_aj', 'aj', 
    'probability', 'same_chain', 'source', 'file', 'rc_probability', 'rc_file', 'sigma', 'epsilon', '1-4', 'distance', 'distance_m', 'cutoff', 'rep']]
    # Inverse pairs calvario
    # this must list ALL COLUMNS!
    inverse_meGO_LJ = meGO_LJ[['molecule_name_aj', 'aj', 'molecule_name_ai', 'ai',
    'probability', 'same_chain', 'source', 'file', 'rc_probability', 'rc_file', 'sigma', 'epsilon', '1-4', 'distance', 'distance_m', 'cutoff', 'rep']].copy()
    inverse_meGO_LJ.columns = ['molecule_name_ai', 'ai', 'molecule_name_aj', 'aj',
    'probability', 'same_chain', 'source', 'file', 'rc_probability', 'rc_file', 'sigma', 'epsilon', '1-4', 'distance', 'distance_m', 'cutoff', 'rep']
    # The contacts are duplicated before cleaning due to the inverse pairs and the sigma calculation requires a simmetric dataframe
    meGO_LJ = pd.concat([meGO_LJ, inverse_meGO_LJ], axis=0, sort=False, ignore_index=True)

    # add a flag to identify learned contacts vs check ones
    meGO_LJ['learned'] = 1

    # process check_dataset
    # the idea here is that if a contact would be attractive/mild c12 smaller in the check dataset
    # and the same contact is instead repulsive in the training, then the repulsive term can be rescaled
    #
    if not check_dataset.empty:
        # Remove low probability ones
        meGO_check_contacts = check_dataset.loc[((check_dataset['probability']>parameters.md_threshold)&(check_dataset['probability']>=check_dataset['rc_probability']))|(check_dataset['1-4']=="1_4")].copy()
        meGO_check_contacts = meGO_check_contacts.loc[~((meGO_check_contacts['same_chain']==True)&(meGO_check_contacts['rep']==0))]
        meGO_check_contacts['sigma'] = (meGO_check_contacts['distance']) / (2.**(1./6.))
        meGO_check_contacts['learned'] = 0
        # set the epsilon of these contacts using the default c12 repulsive term
        meGO_check_contacts['epsilon'] = -meGO_check_contacts['rep']
        meGO_LJ = pd.concat([meGO_LJ, meGO_check_contacts], axis=0, sort=False, ignore_index=True)
        meGO_LJ.drop_duplicates(inplace=True, ignore_index = True)
        # this calculates the increase in energy (if any) to form the "check" contact 
        energy_at_check_dist = meGO_LJ.groupby(by=['ai', 'aj', 'same_chain'])[['distance', 'distance_14', 'epsilon', 'source', 'same_chain', '1-4']].apply(check_LJ, parameters)
        meGO_LJ = pd.merge(meGO_LJ, energy_at_check_dist.rename('energy_at_check_dist'), how="inner", on=["ai", "aj", "same_chain"])
        # now we should keep only those check_with contacts that are unique and whose energy_at_check_dist is large 
        meGO_LJ.sort_values(by = ['ai', 'aj', 'same_chain', 'learned'], ascending = [True, True, True, False], inplace = True)
        # Cleaning the duplicates
        meGO_LJ = meGO_LJ.drop_duplicates(subset = ['ai', 'aj', 'same_chain'], keep = 'first')
        # rescale problematic contacts
        meGO_LJ['epsilon'] *= meGO_LJ['energy_at_check_dist']
        meGO_LJ.drop('energy_at_check_dist', axis=1, inplace=True)

    # Here we create a copy of contacts to be added in pairs-exclusion section in topol.top.
    # All contacts should be applied intermolecularly, but intermolecular specific contacts are not used intramolecularly.
    # meGO_LJ_14 will be handled differently to overcome this issue.
    meGO_LJ = meGO_LJ[meGO_LJ.epsilon != 0]
    meGO_LJ_14 = meGO_LJ.copy()

    # Here we sort all the atom pairs based on the distance and we keep the closer ones.
    # Sorting the pairs prioritising intermolecular interactions and shorter length ones
    # When distances are the same than we choose the most likely 
    meGO_LJ.sort_values(by = ['ai', 'aj', 'same_chain', 'sigma', 'probability'], ascending = [True, True, True, True, False], inplace = True)
    # Cleaning the duplicates
    meGO_LJ = meGO_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    # Removing the reverse duplicates
    cols = ['ai', 'aj']
    meGO_LJ[cols] = np.sort(meGO_LJ[cols].values, axis=1)
    meGO_LJ = meGO_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')

    # Pairs prioritise intramolecular interactions
    meGO_LJ_14.sort_values(by = ['ai', 'aj', 'same_chain', 'sigma', 'probability'], ascending = [True, True, False, True, False], inplace = True)
    meGO_LJ_14 = meGO_LJ_14.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    meGO_LJ_14[cols] = np.sort(meGO_LJ_14[cols].values, axis=1)
    meGO_LJ_14 = meGO_LJ_14.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')

    # where meGO_LJ_14 is the same of meGO_LJ and same_chain is yes that the line can be dropped
    # that is I want to keep lines with same_chain no or lines with same chain yes that have same_chain no in meGO_LJ
    test = pd.merge(meGO_LJ_14, meGO_LJ, how="right", on=["ai", "aj"])
    meGO_LJ_14 = test.loc[(test['same_chain_x']==False)|((test['same_chain_x']==True)&(test['same_chain_y']==False))]
    meGO_LJ_14 = meGO_LJ_14.drop(columns = ['sigma_y', 'epsilon_y', 'same_chain_y', 'probability_y', 'rc_probability_y', 'source_y', '1-4_y', 'cutoff_y', 'rep_y'])
    meGO_LJ_14.rename(columns = {'sigma_x': 'sigma', 'probability_x': 'probability', 'rc_probability_x': 'rc_probability', 'epsilon_x': 'epsilon', 'same_chain_x': 'same_chain', 'source_x': 'source', '1-4_x': '1-4', 'cutoff_x': 'cutoff', 'rep_x': 'rep'}, inplace = True)

    # copy 1-4 interactions into meGO_LJ_14
    copy14 = meGO_LJ.loc[(meGO_LJ['1-4']=='1_4')]
    meGO_LJ_14 = pd.concat([meGO_LJ_14,copy14], axis=0, sort = False, ignore_index = True)
    # remove them from the default force-field
    meGO_LJ = meGO_LJ.loc[(meGO_LJ['1-4']!='1_4')]

    meGO_LJ['c6'] = 4 * meGO_LJ['epsilon'] * (meGO_LJ['sigma'] ** 6)
    meGO_LJ['c12'] = abs(4 * meGO_LJ['epsilon'] * (meGO_LJ['sigma'] ** 12))
    meGO_LJ.loc[(meGO_LJ['epsilon']<0.), 'c6'] = 0.
    meGO_LJ.loc[(meGO_LJ['epsilon']<0.), 'c12'] = -meGO_LJ['epsilon']
    # remove attractive self interactions for proteins' side-chains and scale their c12s
    meGO_LJ.loc[(meGO_LJ['ai']==meGO_LJ['aj'])&(meGO_LJ['ai'].map(meGO_ensemble['sbtype_moltype_dict'])=="protein")&(meGO_LJ['epsilon']>0.)&~((meGO_LJ['ai'].map(meGO_ensemble['sbtype_name_dict'])=="CA")|(meGO_LJ['ai'].map(meGO_ensemble['sbtype_name_dict'])=="N")|(meGO_LJ['ai'].map(meGO_ensemble['sbtype_name_dict'])=="C")|(meGO_LJ['ai'].map(meGO_ensemble['sbtype_name_dict'])=="O")), 'c12'] = np.minimum(meGO_LJ['distance']**12, meGO_LJ['rep'])
    meGO_LJ.loc[(meGO_LJ['ai']==meGO_LJ['aj'])&(meGO_LJ['ai'].map(meGO_ensemble['sbtype_moltype_dict'])=="protein")&(meGO_LJ['epsilon']>0.)&~((meGO_LJ['ai'].map(meGO_ensemble['sbtype_name_dict'])=="CA")|(meGO_LJ['ai'].map(meGO_ensemble['sbtype_name_dict'])=="N")|(meGO_LJ['ai'].map(meGO_ensemble['sbtype_name_dict'])=="C")|(meGO_LJ['ai'].map(meGO_ensemble['sbtype_name_dict'])=="O")), 'c6'] = 0 


    meGO_LJ_14['c6'] = 4 * meGO_LJ_14['epsilon'] * (meGO_LJ_14['sigma'] ** 6)
    meGO_LJ_14['c12'] = abs(4 * meGO_LJ_14['epsilon'] * (meGO_LJ_14['sigma'] ** 12))
    # repulsive interactions have just a very large C12
    meGO_LJ_14.loc[(meGO_LJ_14['epsilon']<0.), 'c6'] = 0.
    meGO_LJ_14.loc[(meGO_LJ_14['epsilon']<0.), 'c12'] = -meGO_LJ_14['epsilon']  

    meGO_LJ['type'] = 1
    meGO_LJ['number_ai'] = meGO_LJ['ai'].map(meGO_ensemble['sbtype_number_dict'])
    meGO_LJ['number_aj'] = meGO_LJ['aj'].map(meGO_ensemble['sbtype_number_dict'])
    meGO_LJ['number_ai'] = meGO_LJ['number_ai'].astype(int)
    meGO_LJ['number_aj'] = meGO_LJ['number_aj'].astype(int)

    meGO_LJ = meGO_LJ[['ai', 'aj', 'type', 'c6', 'c12', 'sigma', 'epsilon', 'probability', 'rc_probability', 'molecule_name_ai',  'molecule_name_aj', 'same_chain', 'source', 'file', 'rc_file', 'number_ai', 'number_aj', 'cutoff']]
    # Here we want to sort so that ai is smaller than aj
    inv_meGO = meGO_LJ[['aj', 'ai', 'type', 'c6', 'c12', 'sigma', 'epsilon', 'probability', 'rc_probability', 'molecule_name_aj',  'molecule_name_ai', 'same_chain', 'source', 'file', 'rc_file', 'number_aj', 'number_ai', 'cutoff']].copy()
    inv_meGO.columns = ['ai', 'aj', 'type', 'c6', 'c12', 'sigma', 'epsilon', 'probability', 'rc_probability', 'molecule_name_ai',  'molecule_name_aj', 'same_chain', 'source', 'file', 'rc_file', 'number_ai', 'number_aj', 'cutoff'] 
    meGO_LJ = pd.concat([meGO_LJ,inv_meGO], axis=0, sort = False, ignore_index = True)
    meGO_LJ = meGO_LJ[meGO_LJ['number_ai']<=meGO_LJ['number_aj']]
    meGO_LJ.sort_values(by = ['number_ai', 'number_aj'], inplace = True)
    meGO_LJ = meGO_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')

    return meGO_LJ, meGO_LJ_14


def check_LJ(test, parameters):
    energy = 1.
    if len(test) == 1:
        # this is the case where we have a contact from check and default c12s from train
        if (len(test.loc[test.source.isin(parameters.check_with)])) == 1:
            # default c12
            eps = -test.loc[(test.source.isin(parameters.check_with))].iloc[0]['epsilon']
            #distance from check
            dist_check = test.loc[(test.source.isin(parameters.check_with))].iloc[0]['distance']
            if dist_check < eps**(1./12.):
                energy = ((dist_check)**12)/eps
    else:
        # this is the special case for 1-4 interactions
        if((test.loc[test.source.isin(parameters.check_with)]).iloc[0]['1-4'] == "1_4"):
            #distance from check
            dist_check = test.loc[(test.source.isin(parameters.check_with))].iloc[0]['distance_14']
            #distance from train 
            dist_train = test.loc[~(test.source.isin(parameters.check_with))].iloc[0]['distance_14']
            if dist_check < dist_train:
                energy = (dist_check/dist_train)**12

        # this is the case where we can a contact defined in both check and train
        else:
            #distance from check
            dist_check = test.loc[(test.source.isin(parameters.check_with))].iloc[0]['distance']
            #distance from train 
            dist_train = test.loc[~(test.source.isin(parameters.check_with))].iloc[0]['distance']
            #epsilon from train
            eps = test.loc[~(test.source.isin(parameters.check_with))].iloc[0]['epsilon']
            if dist_check < dist_train and eps < 0:
                #energy = -eps/(dist_check)**12+eps/(dist_train)**12
                energy = (dist_check/dist_train)**12
 
    return energy


def make_pairs_exclusion_topology(meGO_ensemble, meGO_LJ_14):
    '''
    This function prepares the [ exclusion ] and [ pairs ] section to output to topology.top
    
    Parameters
    ----------
    meGO_ensemlbe : dict
        The meGO_ensemble object contains all the relevant system information
    meGO_LJ_14 : pd.DataFrame
        Contains the contact information for the 1-4 interactions

    Returns
    -------
    pairs_molecule_dict : dict
        Contains the "write out"-ready pairs-exclusions interactions for each molecule  
    '''
    pairs_molecule_dict = {}
    for molecule, bond_pair in meGO_ensemble['bond_pairs'].items():
        reduced_topology = meGO_ensemble['topology_dataframe'].loc[meGO_ensemble['topology_dataframe']['molecule_name'] == molecule][['number', 'sb_type', 'resnum', 'name', 'type', 'resname', 'molecule_type']].copy()

        reduced_topology['number'] = reduced_topology['number'].astype(str)
        reduced_topology['resnum'] = reduced_topology['resnum'].astype(int)

        atnum_type_dict = reduced_topology.set_index('sb_type')['number'].to_dict()
        # Building the exclusion bonded list
        # exclusion_bonds are all the interactions within 3 bonds
        # p14 are specifically the interactions at exactly 3 bonds
        exclusion_bonds, p14 = multiego.topology.get_14_interaction_list(reduced_topology, bond_pair) 
        pairs = pd.DataFrame()
        if not meGO_LJ_14.empty:
            # pairs from greta does not have duplicates because these have been cleaned before
            pairs = meGO_LJ_14[['ai', 'aj', 'c6', 'c12', 'epsilon', 'same_chain', 'probability', 'rc_probability', 'source']].copy()
            pairs['c12_ai'] = pairs['ai']
            pairs['c12_aj'] = pairs['aj']

            # The exclusion list was made based on the atom number
            pairs['ai'] = pairs['ai'].map(atnum_type_dict)
            pairs['aj'] = pairs['aj'].map(atnum_type_dict)
            pairs['check'] = pairs['ai'] + '_' + pairs['aj']
            # Here the drop the contacts which are already defined by GROMACS, including the eventual 1-4 exclusion defined in the LJ_pairs
            pairs['remove'] = ''
            pairs.loc[(pairs['check'].isin(exclusion_bonds)), 'remove'] = 'Yes'
            pairs.loc[(pairs['check'].isin(p14)&(pairs['same_chain']==True)), 'remove'] = 'No'
            mask = pairs.remove == 'Yes'
            pairs = pairs[~mask]
            pairs['c12_ai'] = pairs['c12_ai'].map(meGO_ensemble['sbtype_c12_dict'])
            pairs['c12_aj'] = pairs['c12_aj'].map(meGO_ensemble['sbtype_c12_dict'])
            pairs['func'] = 1
            # Intermolecular interactions are excluded 
            pairs.loc[(pairs['same_chain'] == False), 'c6'] = 0.
            pairs.loc[(pairs['same_chain'] == False), 'c12'] = np.sqrt(pairs['c12_ai'] * pairs['c12_aj'])  
            # this is a safety check 
            pairs = pairs[pairs['c12']>0.]
            pairs.drop(columns = ['same_chain', 'c12_ai', 'c12_aj', 'check', 'remove', 'epsilon'], inplace = True)
            pairs = pairs[['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source']]

            # Drop NaNs. This is an issue when adding the ligand ensemble.
            pairs.dropna(inplace=True)
            pairs['ai'] = pairs['ai'].astype(int)
            pairs['aj'] = pairs['aj'].astype(int)  
    
            # Here we want to sort so that ai is smaller than aj
            inv_pairs = pairs[['aj', 'ai', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source']].copy()
            inv_pairs.columns = ['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source']
            pairs = pd.concat([pairs,inv_pairs], axis=0, sort = False, ignore_index = True)
            pairs = pairs[pairs['ai']<pairs['aj']]
            pairs.drop_duplicates(inplace=True, ignore_index = True)
            pairs.sort_values(by = ['ai', 'aj'], inplace = True)

        pairs_molecule_dict[molecule] = pairs

    return pairs_molecule_dict
