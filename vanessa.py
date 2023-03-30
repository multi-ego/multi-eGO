from atomtypes_definitions import gromos_atp, from_ff_to_multiego, aminoacids_list, nucleic_acid_list
import argparse
import pandas as pd
import os
import glob
import numpy as np
import sys

def read_molecular_contacts(path):
    '''
    This function add column headers and add a column based whether the contacts are intramolecular or intermolecular
    '''

    print('\t-', f"Reading {path}")
    contact_matrix = pd.read_csv(path, header=None, sep='\s+')
    contact_matrix.columns = ['molecule_number_ai', 'ai', 'molecule_number_aj', 'aj', 'distance_m', 'distance_NMR', 'probability', 'flag', 'sigma']
    contact_matrix.drop(columns=['sigma'], inplace=True)
    contact_matrix.columns = ['molecule_number_ai', 'ai', 'molecule_number_aj', 'aj', 'distance_m', 'distance', 'probability', 'flag']
    contact_matrix['molecule_number_ai'] = contact_matrix['molecule_number_ai'].astype(str)
    contact_matrix['ai'] = contact_matrix['ai'].astype(str)
    contact_matrix['molecule_number_aj'] = contact_matrix['molecule_number_aj'].astype(str)
    contact_matrix['aj'] = contact_matrix['aj'].astype(str)

    return contact_matrix


def initialize_ensemble_topology(topology, simulation):
    '''
    '''
    # In a single topology different type of molecules can be present (e.g. protein, ligand).
    # For each molecule the different atomtypes are saved.
    print(
        '\t-', f'Reading {simulation} topology containing: {topology.molecules}')
    columns_to_drop = ['nb_idx', 'solvent_radius', 'screen', 'occupancy', 'bfactor',
                       'altloc', 'join', 'irotat', 'rmin', 'rmin_14', 'epsilon_14', 'tree']
    ensemble_topology_dataframe, new_number, col_molecule, new_resnum, ensemble_molecules_idx_sbtype_dictionary, temp_number_c12_dict = pd.DataFrame(), [], [], [], {}, {}

    # I needed to add this for loop as by creating the topology dataframe by looping over molecules, the c12 information is lost
    for atom in topology.atoms:    
        temp_number_c12_dict[str(atom.idx+1)] = atom.epsilon*4.184
        
    for molecule_number, (molecule_name, molecule_topology) in enumerate(topology.molecules.items(), 1):
        molecule_type_dict = assign_molecule_type(molecule_name, molecule_topology[0])
        ensemble_molecules_idx_sbtype_dictionary[f'{str(molecule_number)}_{molecule_name}'] = {}
        ensemble_topology_dataframe = pd.concat([ensemble_topology_dataframe, molecule_topology[0].to_dataframe()], axis=0)
        for atom in molecule_topology[0].atoms:
            new_number.append(str(atom.idx+1))
            col_molecule.append(f'{molecule_number}_{molecule_name}')
            new_resnum.append(str(atom.residue.number))
    del molecule_name

    ensemble_topology_dataframe['number'] = new_number
    ensemble_topology_dataframe['molecule'] = col_molecule
    ensemble_topology_dataframe['molecule_number'] = col_molecule
    ensemble_topology_dataframe[['molecule_number', 'molecule_name']] = ensemble_topology_dataframe.molecule.str.split('_', expand=True)
    ensemble_topology_dataframe['resnum'] = new_resnum
    ensemble_topology_dataframe['cgnr'] = ensemble_topology_dataframe['resnum']
    ensemble_topology_dataframe['ptype'] = 'A'
    ensemble_topology_dataframe = ensemble_topology_dataframe.replace(
        {'name': from_ff_to_multiego})
    ensemble_topology_dataframe['sb_type'] = ensemble_topology_dataframe['name'] + '_' +ensemble_topology_dataframe['molecule_name']+ '_' + ensemble_topology_dataframe['resnum'].astype(str)
    ensemble_topology_dataframe.rename(columns = {'epsilon':'c12'}, inplace=True)
    
    ensemble_topology_dataframe['charge'] = 0.
    ensemble_topology_dataframe['c6'] = 0.
    ensemble_topology_dataframe['c12'] = ensemble_topology_dataframe['number'].map(temp_number_c12_dict)
    ensemble_topology_dataframe['molecule_type'] = ensemble_topology_dataframe['molecule_name'].map(molecule_type_dict)
    
    for molecule in ensemble_molecules_idx_sbtype_dictionary.keys():
        temp_topology_dataframe = ensemble_topology_dataframe.loc[ensemble_topology_dataframe['molecule'] == molecule]
        number_sbtype_dict = temp_topology_dataframe[['number', 'sb_type']].set_index('number')['sb_type'].to_dict()
        ensemble_molecules_idx_sbtype_dictionary[molecule] = number_sbtype_dict
    sbtype_c12_dict = ensemble_topology_dataframe[['sb_type', 'c12']].set_index('sb_type')['c12'].to_dict()

    return ensemble_topology_dataframe, ensemble_molecules_idx_sbtype_dictionary, sbtype_c12_dict, molecule_type_dict


def get_bonds(topology):
    bonds_dataframe = pd.DataFrame({
        'ai': [bonds.atom1.idx + 1 for bonds in topology],
        'aj': [bonds.atom2.idx + 1 for bonds in topology],
        'funct': [bonds.funct for bonds in topology],
        'req': [bonds.type.req for bonds in topology],
        'k': [bonds.type.k for bonds in topology]
    })
    # Conversion from KCal/mol/A^2 to KJ/mol/nm^2 and from Amber to Gromos
    bonds_dataframe['req'] = bonds_dataframe['req']/10.
    bonds_dataframe['k'] = bonds_dataframe['k']*4.184*100*2
    bonds_dataframe['k'] = bonds_dataframe['k'].map(lambda x:'{:.6e}'.format(x))
    return bonds_dataframe


def get_bond_pairs(topology):
    ai, aj = [], []
    for bonds in topology:
        ai.append(bonds.atom1.idx + 1)
        aj.append(bonds.atom2.idx + 1)
    bond_tuple = list([(str(ai), str(aj)) for ai, aj in zip(ai, aj)])
    return bond_tuple


def get_angles(topology):
    angles_dataframe = pd.DataFrame({
        'ai' : [angle.atom1.idx + 1 for angle in topology],
        'aj' : [angle.atom2.idx + 1 for angle in topology],
        'ak' : [angle.atom3.idx + 1 for angle in topology],
        'funct' : [angle.funct for angle in topology],
        'theteq' : [angle.type.theteq for angle in topology],
        'k' : [angle.type.k for angle in topology]
    })
    angles_dataframe['k'] = angles_dataframe['k']*4.184*2
    angles_dataframe['k'] = angles_dataframe['k'].map(lambda x:'{:.6e}'.format(x))
    return angles_dataframe


def get_dihedrals(topology):
    dihedrals_dataframe = pd.DataFrame({
        'ai' : [dihedral.atom1.idx + 1 for dihedral in topology],
        'aj' : [dihedral.atom2.idx + 1 for dihedral in topology],
        'ak' : [dihedral.atom3.idx + 1 for dihedral in topology],
        'al' : [dihedral.atom4.idx + 1 for dihedral in topology],
        'funct' : [dihedral.funct for dihedral in topology],
        'phase' : [dihedral.type.phase for dihedral in topology],
        'phi_k' : [dihedral.type.phi_k for dihedral in topology],
        'per' : [dihedral.type.per for dihedral in topology]
    })
    dihedrals_dataframe['phi_k'] = dihedrals_dataframe['phi_k']*4.184
    return dihedrals_dataframe


def get_impropers(topology):
    impropers_dataframe = pd.DataFrame({
        'ai' : [improper.atom1.idx + 1 for improper in topology],
        'aj' : [improper.atom2.idx + 1 for improper in topology],
        'ak' : [improper.atom3.idx + 1 for improper in topology],
        'al' : [improper.atom4.idx + 1 for improper in topology],
        'funct' : [improper.funct for improper in topology],
        'psi_eq' : [improper.type.psi_eq for improper in topology],
        'psi_k' : [improper.type.psi_k for improper in topology]
    })
    impropers_dataframe['psi_k'] = impropers_dataframe['psi_k']*4.184*2
    return impropers_dataframe


def get_pairs(topology):
    pairs_dataframe = pd.DataFrame({
        'ai' : [pair.atom1.idx + 1 for pair in topology],
        'aj' : [pair.atom2.idx + 1 for pair in topology],
        'funct' : [pair.funct for pair in topology],
        'type' : [pair.type for pair in topology],
    })

    # TODO change unit measures
    return pairs_dataframe
  

def assign_molecule_type(molecule_name, molecule_topology):
    molecule_type_dict = {}
    first_aminoacid = molecule_topology.residues[0].name

    if first_aminoacid in aminoacids_list:
        molecule_type_dict[molecule_name] = 'protein'
    elif first_aminoacid in nucleic_acid_list:
        molecule_type_dict[molecule_name] = 'nucleic_acid'
    else:
        molecule_type_dict[molecule_name] = 'other'
    return molecule_type_dict


def initialize_molecular_contacts(contact_matrices, ensemble_molecules_idx_sbtype_dictionary, simulation):
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
    '''
    print('\t\t-', f'Initializing {simulation} contact matrix')
    ensemble_contact_matrix = pd.DataFrame()
    molecule_names_dictionary = {}
    for molecule_name in ensemble_molecules_idx_sbtype_dictionary.keys():
        name = molecule_name.split('_')
        molecule_names_dictionary[str(name[0])] = name[1]
    
    counter = 1
    original_contact_matrices = contact_matrices.copy()
    for file_name, contact_matrix in original_contact_matrices.items():

        #for file_name, contact_matrix in contact_matrices.items():
        name = file_name.split('_')
        name = name[:-1] + name[-1].split('.')[:-1]
        
        # Renaming stuff
        contact_matrix['molecule_name_ai'] = contact_matrix['molecule_number_ai'].astype(
            str) + '_' + contact_matrix['molecule_number_ai'].map(molecule_names_dictionary)
        contact_matrix['molecule_name_aj'] = contact_matrix['molecule_number_aj'].astype(
            str) + '_' + contact_matrix['molecule_number_aj'].map(molecule_names_dictionary)
        contact_matrix['ai'] = contact_matrix['ai'].map(
            ensemble_molecules_idx_sbtype_dictionary[contact_matrix['molecule_name_ai'][0]])
        contact_matrix['aj'] = contact_matrix['aj'].map(
            ensemble_molecules_idx_sbtype_dictionary[contact_matrix['molecule_name_aj'][0]])

        contact_matrix = contact_matrix[~contact_matrix['ai'].astype(
            str).str.startswith('H')]
        contact_matrix = contact_matrix[~contact_matrix['aj'].astype(
            str).str.startswith('H')]

        contact_matrix = contact_matrix[[
            'molecule_name_ai', 'ai', 'molecule_name_aj', 'aj', 'distance_m', 'distance', 'probability', 'flag']]
        if name[0] == 'intramat':
            contact_matrix['same_chain'] = True
        elif name[0] == 'intermat':
            contact_matrix['same_chain'] = False
        else:
            raise Exception(
                'There might be an error in the contact matrix naming. It must be intermat_X_X or intramat_X_X')

        contact_matrix['source'] = simulation
        contact_matrix['file'] = '_'.join(name)
        ensemble_contact_matrix = pd.concat([ensemble_contact_matrix, contact_matrix], axis=0)
        ensemble_contact_matrix[['idx_ai', 'idx_aj']] = ensemble_contact_matrix[['ai', 'aj']]
        ensemble_contact_matrix.set_index(['idx_ai', 'idx_aj'], inplace=True)
    return ensemble_contact_matrix
    

def check_LJ(test, parameters):
    if len(test) == 1: 
        return 0. 
    else:
        #distance comes from check
        dist_check = test.loc[(test.source.isin(parameters.check_with))].iloc[0]['distance']
        dist_train = test.loc[~(test.source.isin(parameters.check_with))].iloc[0]['distance']
        #epsilon from train
        if dist_check < dist_train:
            eps = test.loc[~(test.source.isin(parameters.check_with))].iloc[0]['epsilon']
        else:
            eps = 0.
 
        if eps < 0. :
            return -eps/(dist_check)**12
        else:
            return 0.


def parametrize_LJ(topology_dataframe, molecule_type_dict, bond_tuple, pairs_tuple, type_c12_dict, meGO_atomic_contacts, reference_atomic_contacts, check_atomic_contacts, sbtype_number_dict, parameters):
    '''
    This function reads the probabilities obtained using gmx_clustsize from the ensembles defined in the command line.
    The random coil probabilities are used to reweight the explicit water ones.
    Intra and inter molecular contacts are splitted as different rules are applied during the reweighting.
    For each atom contact the sigma and epsilon are obtained.
    '''

    # First of all we generate the random-coil 1-4 interactions:
    pairs14 = pd.DataFrame()
    exclusion_bonds14 = pd.DataFrame()

    for molecule, bond_pair in bond_tuple.items():
        reduced_topology = topology_dataframe.loc[topology_dataframe['molecule_name'] == molecule][['number', 'sb_type', 'resnum', 'name', 'type', 'resname']].copy()

        reduced_topology['number'] = reduced_topology['number'].astype(str)
        reduced_topology['resnum'] = reduced_topology['resnum'].astype(int)
        # Dictionaries definitions to map values
        atnum_type_dict = reduced_topology.set_index('sb_type')['number'].to_dict()
        type_atnum_dict = reduced_topology.set_index('number')['sb_type'].to_dict()

        # Building the exclusion bonded list
        # exclusion_bonds are all the interactions within 3 bonds
        # p14 are specifically the interactions at exactly 3 bonds
        exclusion_bonds, tmp_p14 = list14(reduced_topology, bond_pair)
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
        reduced_topology['c12'] = reduced_topology['sb_type'].map(type_c12_dict)
        # TODO: this should be done with a check on the molecule type (PROTEIN/NUCLEIC/LIPID/OTHER)
        pairs = pd.DataFrame()
        if molecule_type_dict[molecule] == 'protein':
            pairs = protein_LJ14(reduced_topology)
            pairs['ai'] = pairs['ai'].map(type_atnum_dict)
            pairs['aj'] = pairs['aj'].map(type_atnum_dict)
            pairs['rep'] = pairs['c12']
            pairs['same_chain'] = True
        else:
            pairs['ai'] = pairs_tuple[molecule].ai.astype(str)
            pairs['aj'] = pairs_tuple[molecule].aj.astype(str)
            pairs['ai'] = pairs['ai'].map(type_atnum_dict)
            pairs['aj'] = pairs['aj'].map(type_atnum_dict)
            nonprotein_c12 = []
            for test in pairs_tuple[molecule].type:
                nonprotein_c12.append(float(test.epsilon)*4.184)
            pairs['c12'] = nonprotein_c12
            pairs['rep'] = pairs['c12']
            pairs['same_chain'] = True
        pairs14 = pd.concat([pairs14, pairs], axis=0, sort=False, ignore_index=True)

    if parameters.egos != 'rc':
        meGO_atomic_contacts_merged = pd.merge(meGO_atomic_contacts, reference_atomic_contacts, left_index=True, right_index=True, how='outer')
        # TODO arrabbiati quando nel fibril non c'e' un RC intermol che senno' qui da' errore
        meGO_atomic_contacts_merged = meGO_atomic_contacts_merged.loc[meGO_atomic_contacts_merged['same_chain'] == meGO_atomic_contacts_merged['rc_same_chain']]

        # This is to FLAG 1-2, 1-3, 1-4 cases:
        meGO_atomic_contacts_merged = pd.merge(meGO_atomic_contacts_merged, exclusion_bonds14[["ai", "aj", "same_chain", "1-4"]], how="left", on=["ai", "aj", "same_chain"])
        meGO_atomic_contacts_merged['1-4'] = meGO_atomic_contacts_merged['1-4'].fillna('1>4')
        # This is to set the correct default C12 values taking into account specialised 1-4 values (including the special 1-5 O-O)
        meGO_atomic_contacts_merged = pd.merge(meGO_atomic_contacts_merged, pairs14[["ai", "aj", "same_chain", "rep"]], how="left", on=["ai", "aj", "same_chain"])
        meGO_atomic_contacts_merged['1-4'].loc[(meGO_atomic_contacts_merged['rep'].notna())&(meGO_atomic_contacts_merged['rep']!=0.)] = '1_4'
        meGO_atomic_contacts_merged['rep'].loc[(meGO_atomic_contacts_merged['1-4']=="1_2_3")] = 0.
        meGO_atomic_contacts_merged['rep'].loc[(meGO_atomic_contacts_merged['1-4']=="1_4")&(meGO_atomic_contacts_merged['rep'].isnull())] = 0.
        meGO_atomic_contacts_merged['rep'] = meGO_atomic_contacts_merged['rep'].fillna(np.sqrt(meGO_atomic_contacts_merged['ai'].map(type_c12_dict)*meGO_atomic_contacts_merged['aj'].map(type_c12_dict)))

        # This keep only significat attractive/repulsive interactions
        meGO_atomic_contacts_merged = meGO_atomic_contacts_merged.loc[(meGO_atomic_contacts_merged['probability']>parameters.md_threshold)]

        # Add sigma, add epsilon reweighted, add c6 and c12
        meGO_atomic_contacts_merged['sigma'] = (meGO_atomic_contacts_merged['distance']) / (2.**(1/6))
        meGO_atomic_contacts_merged['epsilon'] = np.nan 

        # The index has been reset as here I have issues with multiple index duplicates. The same contact is kept twice: one for intra and one for inter.
        # The following pandas functions cannot handle multiple rows with the same index although it has been defined the "same_chain" filter.
        meGO_atomic_contacts_merged.reset_index(inplace=True)

        # this is the minimum probability ratio for attractive contacts and is set so that epsilon should not be smaller than 0.1*epsilon
        limit_rc = 1./parameters.rc_threshold**0.1

        # Epsilon reweight based on probability
        # Paissoni Equation 2.1
        # Attractive intramolecular
        meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']>limit_rc*np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold))&(meGO_atomic_contacts_merged['same_chain']==True)] = -(parameters.epsilon/np.log(parameters.rc_threshold))*(np.log(meGO_atomic_contacts_merged['probability']/np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold)))

        # Attractive intermolecular
        meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']>limit_rc*np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold))&(meGO_atomic_contacts_merged['same_chain']==False)] = -(parameters.inter_epsilon/np.log(parameters.rc_threshold))*(np.log(meGO_atomic_contacts_merged['probability']/np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold)))
        
        # Repulsive intramolecular
        meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']<1./limit_rc*np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold))&(meGO_atomic_contacts_merged['rc_probability']>parameters.md_threshold)&(meGO_atomic_contacts_merged['same_chain']==True)] = parameters.epsilon*meGO_atomic_contacts_merged['distance_m']**12*np.log(meGO_atomic_contacts_merged['probability']/meGO_atomic_contacts_merged['rc_probability'])-(meGO_atomic_contacts_merged['rep']*(meGO_atomic_contacts_merged['distance_m']**12/meGO_atomic_contacts_merged['rc_distance_m']**12))

        # Repulsive intermolecular
        meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']<1./limit_rc*np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold))&(meGO_atomic_contacts_merged['rc_probability']>parameters.md_threshold)&(meGO_atomic_contacts_merged['same_chain']==False)] = parameters.inter_epsilon*meGO_atomic_contacts_merged['distance_m']**12*np.log(meGO_atomic_contacts_merged['probability']/meGO_atomic_contacts_merged['rc_probability'])-(meGO_atomic_contacts_merged['rep']*(meGO_atomic_contacts_merged['distance_m']**12/meGO_atomic_contacts_merged['rc_distance_m']**12))

        # remove repulsive with values close to the default ones
        meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['epsilon']<0.)&(np.abs(meGO_atomic_contacts_merged['epsilon'])<meGO_atomic_contacts_merged['rep'])] = np.nan
 
        # This set the default (random coil) value for selected 1-4 interactions
        meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['1-4'] == '1_4')] = -meGO_atomic_contacts_merged['rep']

        # mild c12s (for small probability ratios): smaller
        meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']<=limit_rc*np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold))&(meGO_atomic_contacts_merged['probability']>meGO_atomic_contacts_merged['rc_probability'])&((meGO_atomic_contacts_merged['rc_distance_m']-meGO_atomic_contacts_merged['distance_m'])>0.)] = -meGO_atomic_contacts_merged['rep']*(meGO_atomic_contacts_merged['distance_m']/meGO_atomic_contacts_merged['rc_distance_m'])**12 
        # mild c12s (for small probability ratios): larger 
        meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']>=1./limit_rc*np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold))&(meGO_atomic_contacts_merged['probability']<=meGO_atomic_contacts_merged['rc_probability'])&((meGO_atomic_contacts_merged['rc_distance_m']-meGO_atomic_contacts_merged['distance_m'])<0.)&(meGO_atomic_contacts_merged['rc_probability']>parameters.md_threshold)] = -meGO_atomic_contacts_merged['rep']*(meGO_atomic_contacts_merged['distance_m']/meGO_atomic_contacts_merged['rc_distance_m'])**12

        # Rescale c12 1-4 interactions 
        meGO_atomic_contacts_merged['epsilon'].loc[(np.abs(meGO_atomic_contacts_merged['rc_distance_m']-meGO_atomic_contacts_merged['distance_m'])>0.)&(meGO_atomic_contacts_merged['rc_probability']>parameters.md_threshold)&(meGO_atomic_contacts_merged['1-4']=="1_4")&(meGO_atomic_contacts_merged['same_chain']==True)] = -meGO_atomic_contacts_merged['rep']*(meGO_atomic_contacts_merged['distance_m']/meGO_atomic_contacts_merged['rc_distance_m'])**12

        # Here we are reindexing like before
        meGO_atomic_contacts_merged[['idx_ai', 'idx_aj']] = meGO_atomic_contacts_merged[['ai', 'aj']]
        meGO_atomic_contacts_merged.set_index(['idx_ai', 'idx_aj'], inplace=True)
        
        # clean NaN and zeros
        meGO_atomic_contacts_merged.dropna(subset=['epsilon'], inplace=True)
        meGO_atomic_contacts_merged = meGO_atomic_contacts_merged[meGO_atomic_contacts_merged.epsilon != 0]

        # remove unnecessary fields
        meGO_atomic_contacts_merged = meGO_atomic_contacts_merged[['molecule_name_ai', 'ai', 'molecule_name_aj', 'aj', 
        'probability', 'same_chain', 'source', 'flag', 'file', 
        'rc_probability', 'rc_file', 'sigma', 'epsilon', '1-4', 'distance']]
        # Inverse pairs calvario
        # this must list ALL COLUMNS!
        inverse_meGO_atomic_contacts_merged = meGO_atomic_contacts_merged[['molecule_name_aj', 'aj', 'molecule_name_ai', 'ai',
        'probability', 'same_chain', 'source', 'flag', 'file',
        'rc_probability', 'rc_file', 'sigma', 'epsilon', '1-4', 'distance']].copy()
        inverse_meGO_atomic_contacts_merged.columns = ['molecule_name_ai', 'ai', 'molecule_name_aj', 'aj',
        'probability', 'same_chain', 'source', 'flag', 'file',
        'rc_probability', 'rc_file', 'sigma', 'epsilon', '1-4', 'distance']
        # The contacts are duplicated before cleaning due to the inverse pairs and the sigma calculation requires a simmetric dataframe
        meGO_atomic_contacts_merged = pd.concat([meGO_atomic_contacts_merged, inverse_meGO_atomic_contacts_merged], axis=0, sort=False, ignore_index=True)

        # process check_atomic_contacts
        if not check_atomic_contacts is None:
            if not check_atomic_contacts.empty:
                check_atomic_contacts.drop('distance_m', axis=1, inplace=True)
                check_atomic_contacts = check_atomic_contacts.loc[(check_atomic_contacts['probability']>limit*parameters.rc_threshold)] 
                meGO_atomic_contacts_merged = pd.concat([meGO_atomic_contacts_merged, check_atomic_contacts], axis=0, sort=False, ignore_index=True)
                meGO_atomic_contacts_merged.drop_duplicates(inplace=True, ignore_index = True)
                energy_at_check_dist = meGO_atomic_contacts_merged.groupby(by=['ai', 'aj', 'same_chain'])[['distance', 'epsilon', 'source', 'same_chain']].apply(check_LJ, parameters)
                meGO_atomic_contacts_merged = pd.merge(meGO_atomic_contacts_merged, energy_at_check_dist.rename('energy_at_check_dist'), how="inner", on=["ai", "aj", "same_chain"])
                ## remove check_with contacts 
                meGO_atomic_contacts_merged=meGO_atomic_contacts_merged[~meGO_atomic_contacts_merged.source.isin(parameters.check_with)]
                ## rescale problematic contacts
                meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['same_chain']==True)&(meGO_atomic_contacts_merged['energy_at_check_dist']>parameters.epsilon)] *= parameters.epsilon/meGO_atomic_contacts_merged['energy_at_check_dist']
                meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['same_chain']==False)&(meGO_atomic_contacts_merged['energy_at_check_dist']>parameters.inter_epsilon)] *= parameters.inter_epsilon/meGO_atomic_contacts_merged['energy_at_check_dist']
                meGO_atomic_contacts_merged.drop('energy_at_check_dist', axis=1, inplace=True)

        # Here we create a copy of contacts to be added in pairs-exclusion section in topol.top.
        # All contacts should be applied intermolecularly, but intermolecular specific contacts are not used intramolecularly.
        # meGO_LJ_14 will be handled differently to overcome this issue.
        meGO_LJ_14 = meGO_atomic_contacts_merged.copy()

        # Here we sort all the atom pairs based on the distance and we keep the closer ones.
        # Sorting the pairs prioritising intermolecular interactions and shorter length ones
        meGO_atomic_contacts_merged.sort_values(by = ['ai', 'aj', 'same_chain', 'sigma'], ascending = [True, True, True, True], inplace = True)
        # Cleaning the duplicates
        meGO_atomic_contacts_merged = meGO_atomic_contacts_merged.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
        # Removing the reverse duplicates
        cols = ['ai', 'aj']
        meGO_atomic_contacts_merged[cols] = np.sort(meGO_atomic_contacts_merged[cols].values, axis=1)
        meGO_atomic_contacts_merged = meGO_atomic_contacts_merged.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
        
        # Pairs prioritise intramolecular interactions
        meGO_LJ_14.sort_values(by = ['ai', 'aj', 'same_chain', 'sigma'], ascending = [True, True, False, True], inplace = True)
        meGO_LJ_14 = meGO_LJ_14.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
        meGO_LJ_14[cols] = np.sort(meGO_LJ_14[cols].values, axis=1)
        meGO_LJ_14 = meGO_LJ_14.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
        
        # where meGO_LJ_14 is the same of meGO_atomic_contacts_merged and same_chain is yes that the line can be dropped
        # that is I want to keep lines with same_chain no or lines with same chain yes that have same_chain no in meGO_atomic_contacts_merged
        test = pd.merge(meGO_LJ_14, meGO_atomic_contacts_merged, how="right", on=["ai", "aj"])
        meGO_LJ_14 = test.loc[(test['same_chain_x']==False)|((test['same_chain_x']==True)&(test['same_chain_y']==False))]
        meGO_LJ_14.drop(columns = ['sigma_y', 'epsilon_y', 'same_chain_y', 'probability_y', 'rc_probability_y', 'source_y', '1-4_y'], inplace = True)
        meGO_LJ_14.rename(columns = {'sigma_x': 'sigma', 'probability_x': 'probability', 'rc_probability_x': 'rc_probability', 'epsilon_x': 'epsilon', 'same_chain_x': 'same_chain', 'source_x': 'source', '1-4_x': '1-4'}, inplace = True)

        # copy 1-4 interactions into meGO_LJ_14
        copy14 = meGO_atomic_contacts_merged.loc[(meGO_atomic_contacts_merged['1-4']=='1_4')]
        meGO_LJ_14 = pd.concat([meGO_LJ_14,copy14], axis=0, sort = False, ignore_index = True)
        
        meGO_atomic_contacts_merged['c6'] = 4 * meGO_atomic_contacts_merged['epsilon'] * (meGO_atomic_contacts_merged['sigma'] ** 6)
        meGO_atomic_contacts_merged['c12'] = abs(4 * meGO_atomic_contacts_merged['epsilon'] * (meGO_atomic_contacts_merged['sigma'] ** 12))
        meGO_atomic_contacts_merged['c6'].loc[(meGO_atomic_contacts_merged['epsilon']<0.)] = 0.
        meGO_atomic_contacts_merged['c12'].loc[(meGO_atomic_contacts_merged['epsilon']<0.)] = -meGO_atomic_contacts_merged['epsilon']

        meGO_LJ_14['c6'] = 4 * meGO_LJ_14['epsilon'] * (meGO_LJ_14['sigma'] ** 6)
        meGO_LJ_14['c12'] = abs(4 * meGO_LJ_14['epsilon'] * (meGO_LJ_14['sigma'] ** 12))
        # repulsive interactions have just a very large C12
        meGO_LJ_14['c6'].loc[(meGO_LJ_14['epsilon']<0.)] = 0.
        meGO_LJ_14['c12'].loc[(meGO_LJ_14['epsilon']<0.)] =-meGO_LJ_14['epsilon']  

        # TODO aggiungere nel print anche il contributo di ogni ensemble
        print(f'''
        \t- LJ parameterization completed with a total of {len(meGO_atomic_contacts_merged)} contacts.
        \t- The average epsilon is {meGO_atomic_contacts_merged['epsilon'].loc[meGO_atomic_contacts_merged['epsilon']>0.].mean()}
        \t- The maximum epsilon is {meGO_atomic_contacts_merged['epsilon'].max()}
        ''')

        meGO_atomic_contacts_merged['type'] = 1
        meGO_atomic_contacts_merged['number_ai'] = meGO_atomic_contacts_merged['ai'].map(sbtype_number_dict)
        meGO_atomic_contacts_merged['number_aj'] = meGO_atomic_contacts_merged['aj'].map(sbtype_number_dict)
        meGO_atomic_contacts_merged['number_ai'] = meGO_atomic_contacts_merged['number_ai'].astype(int)
        meGO_atomic_contacts_merged['number_aj'] = meGO_atomic_contacts_merged['number_aj'].astype(int)

        meGO_atomic_contacts_merged = meGO_atomic_contacts_merged[['ai', 'aj', 'type', 'c6', 'c12', 'sigma', 'epsilon', 'probability', 'rc_probability', 'flag', 'molecule_name_ai',  'molecule_name_aj', 'same_chain', 'source', 'file', 'rc_file', 'number_ai', 'number_aj']]
        # Here we want to sort so that ai is smaller than aj
        inv_meGO = meGO_atomic_contacts_merged[['aj', 'ai', 'type', 'c6', 'c12', 'sigma', 'epsilon', 'probability', 'rc_probability', 'flag', 'molecule_name_aj',  'molecule_name_ai', 'same_chain', 'source', 'file', 'rc_file', 'number_aj', 'number_ai']].copy()
        inv_meGO.columns = ['ai', 'aj', 'type', 'c6', 'c12', 'sigma', 'epsilon', 'probability', 'rc_probability', 'flag', 'molecule_name_ai',  'molecule_name_aj', 'same_chain', 'source', 'file', 'rc_file', 'number_ai', 'number_aj'] 
        meGO_atomic_contacts_merged = pd.concat([meGO_atomic_contacts_merged,inv_meGO], axis=0, sort = False, ignore_index = True)
        meGO_atomic_contacts_merged = meGO_atomic_contacts_merged[meGO_atomic_contacts_merged['number_ai']<=meGO_atomic_contacts_merged['number_aj']]
        meGO_atomic_contacts_merged.sort_values(by = ['number_ai', 'number_aj'], inplace = True)
        meGO_atomic_contacts_merged = meGO_atomic_contacts_merged.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')

    else:
        meGO_atomic_contacts_merged = pd.DataFrame()
        meGO_LJ_14 = pairs
        meGO_LJ_14['epsilon'] = -meGO_LJ_14['c12'] 
        
    return meGO_atomic_contacts_merged, meGO_LJ_14


def make_pairs_exclusion_topology(topology_dataframe, bond_tuple, type_c12_dict, parameters, meGO_LJ_14):
    '''
    This function prepares the [ exclusion ] and [ pairs ] section to paste in topology.top
    '''
    pairs_molecule_dict = {}
    for molecule, bond_pair in bond_tuple.items():
        reduced_topology = topology_dataframe.loc[topology_dataframe['molecule_name'] == molecule][['number', 'sb_type', 'resnum', 'name', 'type', 'resname', 'molecule_type']].copy()

        reduced_topology['number'] = reduced_topology['number'].astype(str)
        reduced_topology['resnum'] = reduced_topology['resnum'].astype(int)

        atnum_type_dict = reduced_topology.set_index('sb_type')['number'].to_dict()
        # Building the exclusion bonded list
        # exclusion_bonds are all the interactions within 3 bonds
        # p14 are specifically the interactions at exactly 3 bonds
        exclusion_bonds, p14 = list14(reduced_topology, bond_pair) 
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
            pairs['c12_ai'] = pairs['c12_ai'].map(type_c12_dict)
            pairs['c12_aj'] = pairs['c12_aj'].map(type_c12_dict)
            pairs['func'] = 1
            # Intermolecular interactions are excluded 
            pairs['c6'].loc[(pairs['same_chain'] == False)] = 0.
            pairs['c12'].loc[(pairs['same_chain'] == False)] = np.sqrt(pairs['c12_ai'] * pairs['c12_aj'])  
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


def list14(reduced_topology, bond_pair):
    # Building the exclusion bonded list
    # exclusion_bonds are all the interactions within 3 bonds
    # p14 are specifically the interactions at exactly 3 bonds
    ex, ex14, p14, exclusion_bonds = [], [], [], []
    for atom in reduced_topology['number'].to_list():
        for t in bond_pair:
            if t[0] == atom:
                first = t[1]
                ex.append(t[1])
            elif t[1] == atom:
                first = t[0]
                ex.append(t[0])
            else: continue
            for tt in bond_pair:
                if (tt[0] == first) & (tt[1] != atom):
                    second = tt[1]
                    ex.append(tt[1])
                elif (tt[1] == first) & (tt[0] != atom):
                    second = tt[0]
                    ex.append(tt[0])
                else: continue
                for ttt in bond_pair:
                    if (ttt[0] == second) & (ttt[1] != first):
                        ex.append(ttt[1])
                        ex14.append(ttt[1])

                    elif (ttt[1] == second) & (ttt[0] != first):
                        ex.append(ttt[0])
                        ex14.append(ttt[0])
        for e in ex:
            exclusion_bonds.append((str(str(atom) + '_' + str(e))))
            exclusion_bonds.append((str(str(e) + '_' + str(atom))))
        ex = []
        for e in ex14:
            p14.append((str(str(atom) + '_' + str(e))))
            p14.append((str(str(e) + '_' + str(atom))))
        ex14 = []

    return exclusion_bonds, p14


def protein_LJ14(reduced_topology):
    # Here we make a dictionary of the atoms used for local geometry 
    backbone_nitrogen = reduced_topology.loc[reduced_topology['name'] == 'N']
    backbone_carbonyl = reduced_topology.loc[reduced_topology['name'] == 'C']
    backbone_oxygen = reduced_topology.loc[reduced_topology['name']=='O']
    ct_oxygen = reduced_topology.loc[(reduced_topology['name']=='O1')|(reduced_topology['name']=='O2')]
    sidechain_cb = reduced_topology.loc[reduced_topology['name'] == 'CB']
    pro_cd = reduced_topology.loc[(reduced_topology['name'] == 'CD')&(reduced_topology['resname'] == 'PRO')]
    sidechain_cgs = reduced_topology.loc[(reduced_topology['name'] == 'CG')|(reduced_topology['name'] == 'CG1')|(reduced_topology['name'] == 'CG2')|(reduced_topology['name'] == 'SG')|(reduced_topology['name'] == 'OG')|(reduced_topology['name'] == 'OG1')&(reduced_topology['resname'] != 'PRO')]
    
    pairs = pd.DataFrame()  
   
    # For proline CD take the CB, N of the previous residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=pro_cd, atomtype2=sidechain_cb, constant=2.715402e-06, shift=-1)], axis=0, sort=False, ignore_index=True)
    
    # For backbone carbonyl take the CB of the next residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=backbone_carbonyl, atomtype2=sidechain_cb, prefactor=0.275, shift=+1)], axis=0, sort=False, ignore_index=True)
    
    # For backbone oxygen take the CB of the same residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=backbone_oxygen, atomtype2=sidechain_cb, prefactor=0.1)], axis=0, sort=False, ignore_index=True)
    
    # now we add the pair between the last CB and the two OCT ones
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=ct_oxygen, atomtype2=sidechain_cb, prefactor=0.1)], axis=0, sort=False, ignore_index=True)
    
    # For each backbone nitrogen take the CB of the previuos residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=backbone_nitrogen, atomtype2=sidechain_cb, prefactor=0.65, shift=-1)], axis=0, sort=False, ignore_index=True)
    
    # For each backbone nitrogen take the N of the next residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=backbone_nitrogen, atomtype2=backbone_nitrogen, prefactor=0.343, shift=+1)], axis=0, sort=False, ignore_index=True)
    
    # For each backbone oxygen take the O of the next residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=backbone_oxygen, atomtype2=backbone_oxygen, prefactor=11.4, shift=+1)], axis=0, sort=False, ignore_index=True)
    
    # now we add the pair between the penultimate oxygen and the two CT ones
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=ct_oxygen, atomtype2=backbone_oxygen, prefactor=11.4, shift=-1)], axis=0, sort=False, ignore_index=True)
    
    # For each backbone carbonyl take the carbonyl of the next residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=backbone_carbonyl, atomtype2=backbone_carbonyl, prefactor=0.5, shift=-1)], axis=0, sort=False, ignore_index=True)
    
    # For each backbone carbonyl take the CGs of the same residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=sidechain_cgs, atomtype2=backbone_carbonyl, prefactor=0.078)], axis=0, sort=False, ignore_index=True)
    
    # For each backbone nitrogen take the CGs of the same residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=sidechain_cgs, atomtype2=backbone_nitrogen, prefactor=0.087)], axis=0, sort=False, ignore_index=True)
    # make it symmetric 
    inv_LJ = pairs[['aj', 'ai', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source']].copy()
    inv_LJ.columns = ['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source']
    pairs = pd.concat([pairs, inv_LJ], axis=0, sort = False, ignore_index = True)
    pairs['ai'] = pairs['ai'].astype(str)
    pairs['aj'] = pairs['aj'].astype(str)

    return pairs


def create_pairs_14_dataframe(atomtype1, atomtype2, c6 = 0.0, shift = 0, prefactor = None, constant = None):
    pairs_14_ai, pairs_14_aj, pairs_14_c6, pairs_14_c12 = [], [], [], []
    for index, line_atomtype1 in atomtype1.iterrows():
        line_atomtype2 = atomtype2.loc[(atomtype2['resnum'] == line_atomtype1['resnum']+shift)].squeeze(axis=None)
        if not line_atomtype2.empty:
            pairs_14_ai.append(line_atomtype1['number'])
            pairs_14_aj.append(line_atomtype2['number'])
            pairs_14_c6.append(c6)
            if constant is not None:
                pairs_14_c12.append(constant)
            elif prefactor is not None:
                pairs_14_c12.append(prefactor*np.sqrt(line_atomtype1['c12']*line_atomtype2['c12']))
            else:
                print('You need to chose prefactor or constant')
                exit()

    pairs_14 = pd.DataFrame(columns=['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source'])
    pairs_14['ai'] = pairs_14_ai
    pairs_14['aj'] = pairs_14_aj
    pairs_14['func'] = 1
    pairs_14['c6'] = pairs_14_c6
    pairs_14['c12'] = pairs_14_c12
    pairs_14['source'] = '1-4'
    pairs_14['probability'] = 1.0
    pairs_14['rc_probability'] = 1.0

    return pairs_14


# Errors
def float_range(min_value, max_value):
    def check_float_range(value):
        try:
            value = float(value)
        except ValueError:
            raise argparse.ArgumentTypeError("Invalid float value")
        if not min_value <= value <= max_value:
            raise argparse.ArgumentTypeError(f"Float value must be between {min_value} and {max_value}")
        return value
    return check_float_range


def check_files_existance(ego, protein, md_ensembles):
    for ensemble in md_ensembles:
        ensemble = f'inputs/{protein}/{ensemble}'
        if not os.path.exists(ensemble):
            raise FileNotFoundError(f"Folder {ensemble}/ does not exist.")
        else:
            top_files = glob.glob(f'{ensemble}/*.top')
            if not top_files:
                raise FileNotFoundError(f"No .top files found in {ensemble}/")
            ndx_files = glob.glob(f'{ensemble}/*.ndx')
            if not ndx_files and not ego=="rc":
                raise FileNotFoundError(f"No .ndx files found in {ensemble}/")
