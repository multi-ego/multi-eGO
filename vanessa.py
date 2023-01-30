from atomtypes_definitions import gromos_atp, from_ff_to_multiego
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
    contact_matrix.columns = ['molecule_number_ai', 'ai', 'molecule_number_aj', 'aj', 'distance', 'distance_NMR', 'probability']
    contact_matrix.drop(columns=['distance'], inplace=True)
    contact_matrix.columns = ['molecule_number_ai', 'ai', 'molecule_number_aj', 'aj', 'distance', 'probability']
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
    ensemble_topology_dataframe = topology.to_dataframe()
    columns_to_drop = ['nb_idx', 'solvent_radius', 'screen', 'occupancy', 'bfactor',
                       'altloc', 'join', 'irotat', 'rmin', 'rmin_14', 'epsilon_14', 'tree']
    ensemble_topology_dataframe.drop(columns=columns_to_drop, inplace=True)
    new_number, col_molecule, new_resnum, ensemble_molecules_idx_sbtype_dictionary = [], [], [], {}
    for molecule_number, (molecule_name, molecule_topology) in enumerate(topology.molecules.items(), 1):
        ensemble_molecules_idx_sbtype_dictionary[f'{str(molecule_number)}_{molecule_name}'] = {}
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
    ensemble_topology_dataframe['sb_type'] = ensemble_topology_dataframe['name'] + '_'+ensemble_topology_dataframe['molecule_name']+'_' + ensemble_topology_dataframe['resnum'].astype(str)
    ensemble_topology_dataframe.rename(columns = {'epsilon':'c12'}, inplace=True)#, 'sb_type':'type', 'residue_number':'resnr'}, inplace=True)
    
    ensemble_topology_dataframe['charge'] = 0.
    ensemble_topology_dataframe['c6'] = 0.
    ensemble_topology_dataframe['c6'] = ensemble_topology_dataframe['c6'].map(lambda x:'{:.6e}'.format(x))
    # TODO c'è ancora da sistemare il c12 che ha numeri a caso
    ensemble_topology_dataframe['c12'] = ensemble_topology_dataframe['c12']*4.184
    ensemble_topology_dataframe['c12'] = ensemble_topology_dataframe['c12'].map(lambda x:'{:.6e}'.format(x))

    for molecule in ensemble_molecules_idx_sbtype_dictionary.keys():
        temp_topology_dataframe = ensemble_topology_dataframe.loc[ensemble_topology_dataframe['molecule'] == molecule]
        number_sbtype_dict = temp_topology_dataframe[['number', 'sb_type']].set_index('number')['sb_type'].to_dict()
        ensemble_molecules_idx_sbtype_dictionary[molecule] = number_sbtype_dict

    return ensemble_topology_dataframe, ensemble_molecules_idx_sbtype_dictionary


def get_bonds(topology):
    ai, aj, funct, k, req = [], [], [], [], []
    bonds_dataframe = pd.DataFrame()
    for bonds in topology:
        ai.append(bonds.atom1.idx + 1)
        aj.append(bonds.atom2.idx + 1)
        funct.append(bonds.funct)
        k.append(bonds.type.k)
        req.append(bonds.type.req)
    bonds_dataframe['ai'] = ai
    bonds_dataframe['aj'] = aj
    bonds_dataframe['funct'] = funct
    bonds_dataframe['k'] = k
    bonds_dataframe['req'] = req

    # TODO convert cose
    #bonds_dataframe
    return bonds_dataframe


def get_angles(topology):
    ai, aj, ak, funct, k, theteq = [], [], [], [], [], []
    angles_dataframe = pd.DataFrame()
    for angle in topology:
        ai.append(angle.atom1.idx + 1)
        aj.append(angle.atom2.idx + 1)
        ak.append(angle.atom3.idx + 1)
        funct.append(angle.funct)
        k.append(angle.type.k)
        theteq.append(angle.type.theteq)
    angles_dataframe['ai'] = ai
    angles_dataframe['aj'] = aj
    angles_dataframe['ak'] = ak
    angles_dataframe['funct'] = funct
    angles_dataframe['k'] = k
    angles_dataframe['theteq'] = theteq

    # TODO convert cose
    #angles_dataframe
    return angles_dataframe


def get_dihedrals(topology):
    ai, aj, ak, al, funct, phi_k, per, phase = [], [], [], [], [], [], [], []
    dihedrals_dataframe = pd.DataFrame()
    for dihedral in topology:
        ai.append(dihedral.atom1.idx + 1)
        aj.append(dihedral.atom2.idx + 1)
        ak.append(dihedral.atom3.idx + 1)
        al.append(dihedral.atom4.idx + 1)
        funct.append(dihedral.funct)
        phi_k.append(dihedral.type.phi_k)
        per.append(dihedral.type.per)
        phase.append(dihedral.type.phase)
    dihedrals_dataframe['ai'] = ai
    dihedrals_dataframe['aj'] = aj
    dihedrals_dataframe['ak'] = ak
    dihedrals_dataframe['al'] = al
    dihedrals_dataframe['funct'] = funct
    dihedrals_dataframe['phi_k'] = phi_k
    dihedrals_dataframe['per'] = per
    dihedrals_dataframe['phase'] = phase

    # TODO convert cose
    #dihedrals_dataframe
    return dihedrals_dataframe

def get_impropers(topology):
    ai, aj, ak, al, funct, psi_k, psi_eq = [], [], [], [], [], [], []
    impropers_dataframe = pd.DataFrame()
    for improper in topology:
        ai.append(improper.atom1.idx + 1)
        aj.append(improper.atom2.idx + 1)
        ak.append(improper.atom3.idx + 1)
        al.append(improper.atom4.idx + 1)
        funct.append(improper.funct)
        psi_k.append(improper.type.psi_k)
        psi_eq.append(improper.type.psi_eq)
    impropers_dataframe['ai'] = ai
    impropers_dataframe['aj'] = aj
    impropers_dataframe['ak'] = ak
    impropers_dataframe['al'] = al
    impropers_dataframe['funct'] = funct
    impropers_dataframe['psi_k'] = psi_k
    impropers_dataframe['psi_eq'] = psi_eq

    # TODO convert cose
    #impropers_dataframe
    return impropers_dataframe


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
        - Set the 'probability' column less than 0.000001 to 0.000001
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
        print(contact_matrix, file_name, simulation)
        contact_matrix['molecule_name_ai'] = contact_matrix['molecule_number_ai'].astype(
            str) + '_' + contact_matrix['molecule_number_ai'].map(molecule_names_dictionary)
        contact_matrix['molecule_name_aj'] = contact_matrix['molecule_number_aj'].astype(
            str) + '_' + contact_matrix['molecule_number_aj'].map(molecule_names_dictionary)
        contact_matrix['ai'] = contact_matrix['ai'].map(
            ensemble_molecules_idx_sbtype_dictionary[contact_matrix['molecule_name_ai'][0]])
        print(contact_matrix, file_name, simulation)
        contact_matrix['aj'] = contact_matrix['aj'].map(
            ensemble_molecules_idx_sbtype_dictionary[contact_matrix['molecule_name_aj'][0]])


        # TODO drop nan per far andare avanti lo script MA E' DA TOGLIERE ASSOLUTAMENTE
        # TODO
        # TODO
        # TODO
        # TODO
        contact_matrix.dropna(inplace=True)
        # TODO
        # TODO
        # TODO
        # TODO
        # TODO

        contact_matrix = contact_matrix[~contact_matrix['ai'].astype(
            str).str.startswith('H')]
        contact_matrix = contact_matrix[~contact_matrix['aj'].astype(
            str).str.startswith('H')]
        contact_matrix['probability'].loc[contact_matrix['probability'] < (
            0.000001)] = 0.000001

        contact_matrix = contact_matrix[[
            'molecule_name_ai', 'ai', 'molecule_name_aj', 'aj', 'distance', 'probability']]
        if name[0] == 'intramat':
            contact_matrix['same_chain'] = True
        elif name[0] == 'intermat':
            contact_matrix['same_chain'] = False
        else:  # TODO che questa magari si puo' togliere se c'e' un controllo sui file
            print(
                'There might be an error in the contact matrix naming. It must be inter_X_X or intra_X_X')
            sys.exit

        contact_matrix['source'] = simulation
        contact_matrix['file'] = '_'.join(name)
        ensemble_contact_matrix = pd.concat([ensemble_contact_matrix, contact_matrix], axis=0)
        ensemble_contact_matrix[['idx_ai', 'idx_aj']
                                ] = ensemble_contact_matrix[['ai', 'aj']]
        ensemble_contact_matrix.set_index(['idx_ai', 'idx_aj'], inplace=True)
    return ensemble_contact_matrix
    

def parametrize_LJ(meGO_atomic_contacts, reference_atomic_contacts, check_atomic_contacts, sbtype_c12_dict, parameters): # TODO probabilmente qui simulation e' inutile
    '''
    This function reads the probabilities obtained using gmx_clustsize from the ensembles defined in the command line.
    The random coil probabilities are used to reweight the explicit water ones.
    Intra and inter molecular contacts are splitted as different rules are applied during the reweighting.
    For each atom contact the sigma and epsilon are obtained.
    '''
    #print(meGO_atomic_contacts)
    # TODO aggiungere il dizionario per convertire l'atomtype
    meGO_atomic_contacts_merged = pd.merge(meGO_atomic_contacts, reference_atomic_contacts, left_index=True, right_index=True, how='outer')
    # TODO controlla che abbia senso: qui voglio togliere quando same_chain e rc_same_chain sono diversi. Perché o sono intra o sono inter.

    #print(meGO_atomic_contacts)
    #print(reference_atomic_contacts)
    print(meGO_atomic_contacts_merged)

    meGO_atomic_contacts_merged = meGO_atomic_contacts_merged.loc[meGO_atomic_contacts_merged['same_chain'] == meGO_atomic_contacts_merged['rc_same_chain']]
    #meGO_atomic_contacts_merged.to_csv(f'analysis/meGO_atomic_contacts_merged')
    #meGO_atomic_contacts_merged.drop(columns = ['rc_ai', 'rc_aj'], inplace=True)

    meGO_atomic_contacts_merged = meGO_atomic_contacts_merged.loc[(meGO_atomic_contacts_merged['probability']>parameters.md_threshold)]
    # Add sigma, add epsilon reweighted, add c6 and c12
    meGO_atomic_contacts_merged['sigma'] = (meGO_atomic_contacts_merged['distance']) / (2**(1/6))
    meGO_atomic_contacts_merged['epsilon'] = np.nan 

    # Epsilon reweight based on probability
    # Paissoni Equation 2.1
    # Attractive
    #print(meGO_atomic_contacts_merged)
    meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']>1.2*meGO_atomic_contacts_merged['rc_probability'])&(meGO_atomic_contacts_merged['same_chain']=='Yes')] = -(parameters.epsilon/np.log(parameters.rc_threshold))*(np.log(meGO_atomic_contacts_merged['probability']/np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold)))
    meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']>1.2*meGO_atomic_contacts_merged['rc_probability'])&(meGO_atomic_contacts_merged['same_chain']== 'No')] = -(parameters.epsilon/np.log(parameters.md_threshold))*(np.log(meGO_atomic_contacts_merged['probability']/np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.md_threshold)))
    # Repulsive
    meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']<(1./1.2)*meGO_atomic_contacts_merged['rc_probability'])] = np.log(meGO_atomic_contacts_merged['probability']/meGO_atomic_contacts_merged['rc_probability'])*(np.minimum(meGO_atomic_contacts_merged['distance'],meGO_atomic_contacts_merged['rc_distance'])**12)

    # clean NaN and zeros 
    meGO_atomic_contacts_merged.dropna(inplace=True)
    meGO_atomic_contacts_merged = meGO_atomic_contacts_merged[meGO_atomic_contacts_merged.epsilon != 0]

    # Inverse pairs calvario
    # this must list ALL COLUMNS!
    inverse_meGO_atomic_contacts_merged = meGO_atomic_contacts_merged[['molecule_name_ai', 'aj', 'molecule_name_aj', 'ai', 'distance',
     'probability', 'same_chain', 'source', 'file', 'rc_molecule_name_ai',
     'rc_ai', 'rc_molecule_name_aj', 'rc_aj', 'rc_distance',
     'rc_probability', 'rc_same_chain', 'rc_source', 'rc_file', 'sigma',
     'epsilon']].copy()
    inverse_meGO_atomic_contacts_merged.columns = ['molecule_name_ai', 'aj', 'molecule_name_aj', 'ai', 'distance',
     'probability', 'same_chain', 'source', 'file', 'rc_molecule_name_ai',
     'rc_ai', 'rc_molecule_name_aj', 'rc_aj', 'rc_distance',
     'rc_probability', 'rc_same_chain', 'rc_source', 'rc_file', 'sigma',
     'epsilon']
    
    # The contacts are duplicated before cleaning due to the inverse pairs and the sigma calculation requires a simmetric dataframe
    meGO_atomic_contacts_merged = pd.concat([meGO_atomic_contacts_merged, inverse_meGO_atomic_contacts_merged, check_atomic_contacts], axis=0, sort=False, ignore_index=True)

    # Here we change all the sigma for each combination of ai and aj and based on intra or inter molecular contacts. For all those pairs we keep the sigma min.
    meGO_atomic_contacts_merged['new_sigma'] = meGO_atomic_contacts_merged.groupby(by=['ai', 'aj', 'same_chain'])['sigma'].transform('min')
    # If the sigma minimum is shorter than the original sigma in the case of a repulsive contact it is removed
    meGO_atomic_contacts_merged = meGO_atomic_contacts_merged.loc[~((meGO_atomic_contacts_merged['new_sigma']<meGO_atomic_contacts_merged['sigma'])&(meGO_atomic_contacts_merged['epsilon']<0))]
    # And update the sigma, which includes also the sigma check. Sigma check is a list of contacts from other ensembles to prevent the parametrization of long sigmas.
    meGO_atomic_contacts_merged['sigma'] = meGO_atomic_contacts_merged['new_sigma']
    meGO_atomic_contacts_merged.drop(columns=['new_sigma'], inplace=True)

    # Here we create a copy of contacts to be added in pairs-exclusion section in topol.top. All contacts should be applied intermolecularly, but some are not compatible intramolecularly.
    # meGO_pairs will be handled differently to overcome this issue.
    meGO_pairs = meGO_atomic_contacts_merged.copy()

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
    meGO_pairs.sort_values(by = ['ai', 'aj', 'same_chain', 'sigma'], ascending = [True, True, False, True], inplace = True)
    meGO_pairs = meGO_pairs.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    meGO_pairs[cols] = np.sort(meGO_pairs[cols].values, axis=1)
    meGO_pairs = meGO_pairs.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    
    # where meGO_pairs is the same of meGO_atomic_contacts_merged and same_chain is yes that the line can be dropped
    # that is I want to keep lines with same_chain no or lines with same chain yes that have same_chain no in meGO_atomic_contacts_merged
    test = pd.merge(meGO_pairs, meGO_atomic_contacts_merged, how="right", on=["ai", "aj"])
    #print(meGO_pairs)
    meGO_pairs = test.loc[(test['same_chain_x']=='No')|((test['same_chain_x']=='Yes')&(test['same_chain_y']=='No'))]
    #print(meGO_pairs)
    meGO_pairs.drop(columns = ['sigma_y', 'epsilon_y', 'same_chain_y', 'probability_y', 'rc_probability_y', 'source_y'], inplace = True)
    meGO_pairs.rename(columns = {'sigma_x': 'sigma', 'probability_x': 'probability', 'rc_probability_x': 'rc_probability', 'epsilon_x': 'epsilon', 'same_chain_x': 'same_chain', 'source_x': 'source'}, inplace = True)
    # now we copy the lines with negative epsilon from greta to pairs because we want repulsive interactions only intramolecularly
    # and we use same-chain as a flag to keep track of them
    meGO_atomic_contacts_merged['same_chain'].loc[(meGO_atomic_contacts_merged['epsilon']<0.0)&(meGO_atomic_contacts_merged['same_chain']=='Yes')] = 'Move'
    meGO_pairs = pd.concat([meGO_pairs, meGO_atomic_contacts_merged.loc[(meGO_atomic_contacts_merged['same_chain']=='Move')]], axis=0, sort=False, ignore_index = True)
    # and we remove the same lines from meGO_atomic_contacts_merged
    meGO_atomic_contacts_merged = meGO_atomic_contacts_merged[(meGO_atomic_contacts_merged['same_chain']!='Move')]

    
    
    print(meGO_atomic_contacts_merged)
    exit()
    
    
    meGO_atomic_contacts_merged.insert(2, 'type', 1)
    meGO_atomic_contacts_merged.insert(3, 'c6', '')
    meGO_atomic_contacts_merged['c6'] = 4 * meGO_atomic_contacts_merged['epsilon'] * (meGO_atomic_contacts_merged['sigma'] ** 6)
    meGO_atomic_contacts_merged.insert(4, 'c12', '')
    meGO_atomic_contacts_merged['c12'] = abs(4 * meGO_atomic_contacts_merged['epsilon'] * (meGO_atomic_contacts_merged['sigma'] ** 12))
    # repulsive interactions have just a very large C12
    meGO_atomic_contacts_merged['c12ij'] = np.sqrt(meGO_atomic_contacts_merged['ai'].map(sbtype_c12_dict)*meGO_atomic_contacts_merged['aj'].map(sbtype_c12_dict))
    meGO_atomic_contacts_merged['c6'].loc[(meGO_atomic_contacts_merged['epsilon']<0.)] = 0.
    meGO_atomic_contacts_merged['c12'].loc[(meGO_atomic_contacts_merged['epsilon']<0.)] = np.maximum(-meGO_atomic_contacts_merged['epsilon'],meGO_atomic_contacts_merged['c12ij'])

    meGO_pairs.insert(2, 'type', 1)
    meGO_pairs.insert(3, 'c6', '')
    meGO_pairs['c6'] = 4 * meGO_pairs['epsilon'] * (meGO_pairs['sigma'] ** 6)
    meGO_pairs.insert(4, 'c12', '')
    meGO_pairs['c12'] = abs(4 * meGO_pairs['epsilon'] * (meGO_pairs['sigma'] ** 12))
    # repulsive interactions have just a very large C12
    meGO_pairs['c12ij'] = np.sqrt(meGO_pairs['ai'].map(sbtype_c12_dict)*meGO_pairs['aj'].map(sbtype_c12_dict))
    meGO_pairs['c6'].loc[(meGO_pairs['epsilon']<0.)] = 0.
    meGO_pairs['c12'].loc[(meGO_pairs['epsilon']<0.)] = np.maximum(-meGO_pairs['epsilon'],meGO_pairs['c12ij'])  

    print('\t- Cleaning and Merging Complete, pairs count:', len(meGO_atomic_contacts_merged))
    print('\t- LJ Merging completed: ', len(meGO_atomic_contacts_merged))
    print("\t\t- Average epsilon is ", meGO_atomic_contacts_merged['epsilon'].mean())
    print("\t\t- Maximum epsilon is ", meGO_atomic_contacts_merged['epsilon'].max())

    return meGO_atomic_contacts_merged, meGO_pairs











    return meGO_atomic_contacts_merged


def merge_and_clean_LJ(meGO_LJ_potential):#, meGO_LJ_potential, type_c12_dict, type_q_dict, parameters):
    '''
    This function...
    '''






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


def check_files_existance(protein, md_ensembles):
    for ensemble in md_ensembles:
        ensemble = f'inputs/{protein}/{ensemble}'
        if not os.path.exists(ensemble):
            raise FileNotFoundError(f"Folder {ensemble} does not exist.")
        else:
            top_files = glob.glob(f'{ensemble}/*.top')
            if not top_files:
                raise FileNotFoundError(f"No .top files found in {ensemble}.")
            ndx_files = glob.glob(f'{ensemble}/*.ndx')
            if not ndx_files:
                raise FileNotFoundError(f"No .ndx files found in {ensemble}.")
