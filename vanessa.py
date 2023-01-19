from atomtypes_definitions import gromos_atp, from_ff_to_multiego
import argparse
import pandas as pd
import os
import glob
import numpy as np


def read_molecular_contacts(path, is_same_chain):
    '''
    This function add column headers and add a column based whether the contacts are intramolecular or intermolecular
    '''

    print('\t-', f"Reading {path}")
    contact_matrix = pd.read_csv(path, header=None, sep='\s+')
    contact_matrix.columns = ['ai', 'aj', 'distance', 'distance_NMR', 'probability']
    contact_matrix.drop(columns=['distance'], inplace=True)
    contact_matrix.columns = ['ai', 'aj', 'distance', 'probability']
    if is_same_chain:
        contact_matrix['same_chain'] = 'Yes'
    else:
        contact_matrix['same_chain'] = 'No'
    return contact_matrix


def initialize_ensemble_topology(topology, simulation):
    '''
    This function prepare the ensemble topology starting from .top and .gro/.pdb.
    First, it removes the water, solvent and sodium (other atomtypes will be added when required).
    Renumbers the atom ids starting from 1
    '''


    # In a single topology different type of molecules can be present (e.g. protein, ligand).
    # For each molecule the different atomtypes are saved.

    print('\t-', f'Reading {simulation} topology containing:')
    for molecule in topology.molecules:
        print('\t\t', f'{molecule}')

    #print(topology)
    #print(topology.molecules)
    coso = topology.molecules
    topology_dataframe = topology.to_dataframe()
    print(topology_dataframe)
    #print(type(coso))
    

    # TODO questo funziona, bisogna segnarsi il nome ed il numero di molecole all'interno della topologia. Controllare la matrice con più molecole diverse.
    for c, o in coso.items():
        print(c, o[0].to_dataframe())
    #topology = topology["!((:TIP3)|(:SOL)|(:WAT)|(:NA))"]

    # To read the matrix you do not really need all this stuff as all the molecules are defined in the reference topology

    # Dropping some unused columns
    topology_dataframe.drop(columns=['nb_idx', 'solvent_radius', 'screen', 'bfactor', 'occupancy', 'altloc', 'join', 'irotat', 'rmin', 'rmin_14', 'epsilon', 'epsilon_14'], inplace=True)

    # All atomnumbers are -1, here we renumber
    topology_dataframe['number'] = list(range(1, len(topology.atoms)+1))
    # The count starts from 1
    topology_dataframe['resid'] = topology_dataframe['resid'] + 1
    topology_dataframe['resnum'] = topology_dataframe['resid']
    topology_dataframe['cgnr'] = topology_dataframe['resid']
    topology_dataframe['ptype'] = 'A'
    # Excluded volume definition
    topology_dataframe['c6'] = '0.00000e+00'
    topology_dataframe['c12'] = topology_dataframe['type'].map(gromos_atp['c12'])

    topology_dataframe.rename(columns={
        'number':'atom_number',
        'type':'atom_type',
        'resnum':'residue_number',
        'resname':'residue',
        'name':'atom',
    }, inplace=True)

    print('\t\t-', f'Applying topology fixes to {simulation}')
    # Removing an extra H from PRO
    pd.options.mode.chained_assignment = None 
    mask = ((topology_dataframe['residue'] == "PRO") & (topology_dataframe['atom_type'] == 'N'))
    topology_dataframe['mass'][mask] = topology_dataframe['mass'][mask].astype(float).sub(1)
    # Adding an extra H to the N terminal
    mask = ((topology_dataframe['residue_number'] == topology_dataframe['residue_number'].min()) & (topology_dataframe['atom_type'] == 'N'))
    topology_dataframe['mass'][mask] = topology_dataframe['mass'][mask].astype(float).add(2)

    # Aromatic carbons dictionary
    aromatic_carbons_dict = {
        'PHE': ['CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TYR': ['CD1', 'CD2', 'CE1', 'CE2'],
        'HIS': ['CE1', 'CD2'],
        'TRP': ['CD1', 'CE3', 'CZ2', 'CZ3', 'CH2']
    }

    for resname, atomnames in aromatic_carbons_dict.items():
        for atom in atomnames:
            mask = ((topology_dataframe['residue'] == resname) & (topology_dataframe['atom'] == atom))
            topology_dataframe['mass'][mask] = topology_dataframe['mass'][mask].astype(float).add(1)

    print('\t\t-',f'Converting {simulation} atomtypes from different forcefields')
    topology_dataframe = topology_dataframe.replace({'atom':from_ff_to_multiego})

    print('\t\t- Defining multi-eGO atomtypes')
    topology_dataframe['sb_type'] = topology_dataframe['atom'] + '_' + topology_dataframe['residue_number'].astype(str)

    return topology, topology_dataframe


def create_ensemble_dictionaries(topology_dataframe, simulation):
    print('\t-', f'Building {simulation} ensemble dictionaries')
    number_of_atoms = len(topology_dataframe)
    print('\t-', f'Number of atoms in {simulation}: {number_of_atoms}')
    sbtype_idx_dict = topology_dataframe[['atom_number', 'sb_type']].set_index('sb_type')['atom_number'].to_dict()
    print('\t-',f'Created {simulation} dictionary of structure-based atomtypes and atom number')
    idx_sbtype_dict = topology_dataframe[['atom_number', 'sb_type']].set_index('atom_number')['sb_type'].to_dict()
    print('\t-',f'Created {simulation} dictionary of atom number and structure-based atomtypes')
    type_c12_dict = topology_dataframe[['sb_type', 'c12']].rename(columns={'sb_type':'; type'}).set_index('; type')['c12'].to_dict()
    print('\t-',f'Created {simulation} dictionary of structure-based atomtypes and c12')
    type_q_dict = topology_dataframe[['sb_type', 'charge']].rename(columns={'sb_type':'; type'}).set_index('; type')['charge'].to_dict()
    print('\t-',f'Created {simulation} dictionary of structure-based atomtypes and charge')

    return number_of_atoms, sbtype_idx_dict, idx_sbtype_dict, type_c12_dict, type_q_dict


def initialize_molecular_contacts(contact_matrix, idx_sbtype_dict, simulation):

        contact_matrix['idx_ai'] = contact_matrix['ai'].map(idx_sbtype_dict)
        contact_matrix['idx_aj'] = contact_matrix['aj'].map(idx_sbtype_dict)
        contact_matrix.drop(['ai', 'aj'], axis=1, inplace=True)
        contact_matrix.rename(columns={'idx_ai':'ai', 'idx_aj':'aj'}, inplace=True)
        contact_matrix = contact_matrix[~contact_matrix['ai'].astype(str).str.startswith('H')]
        contact_matrix = contact_matrix[~contact_matrix['aj'].astype(str).str.startswith('H')]
        contact_matrix['probability'].loc[contact_matrix['probability'] < (0.000001)] = 0.000001 
        contact_matrix['source'] = simulation

        #contact_matrix.replace({'ai':idx_sbtype_dict}, inplace=True)
        #contact_matrix.replace({'aj':idx_sbtype_dict}, inplace=True)
        #contact_matrix = contact_matrix[~contact_matrix['ai'].astype(str).str.startswith('H')]
        #contact_matrix = contact_matrix[~contact_matrix['aj'].astype(str).str.startswith('H')]
        #contact_matrix['probability'].loc[contact_matrix['probability'] < (0.000001)] = 0.000001 
        #contact_matrix['source'] = simulation

        # Reference simulation only contains random coil contacts
        if simulation == 'reference':
            contact_matrix[['type_ai', 'residue_ai']] = contact_matrix.ai.str.split("_", expand = True)
            contact_matrix[['type_aj', 'residue_aj']] = contact_matrix.aj.str.split("_", expand = True)
            contact_matrix['residue_ai'] = contact_matrix['residue_ai'].astype(int)
            contact_matrix['residue_aj'] = contact_matrix['residue_aj'].astype(int)
            contact_matrix.drop(columns=['type_ai', 'type_aj'], inplace=True)
            #contact_matrix['probability'].loc[contact_matrix['probability'] < (0.000001)] = 0.000001
            new_colnames = []
            for colname in contact_matrix.columns:
                new_colnames.append(f'rc_{colname}')
            contact_matrix.columns = new_colnames
        return contact_matrix


def parametrize_ensemble_LJ(ensemble_contact_matrix, reference_contact_matrix, parameters, simulation): # TODO probabilmente qui simulation e' inutile
    '''
    This function reads the probabilities obtained using gmx_clustsize from the ensembles defined in the command line.
    The random coil probabilities are used to reweight the explicit water ones.
    Intra and inter molecular contacts are splitted as different rules are applied during the reweighting.
    For each atom contact the sigma and epsilon are obtained.
    '''

    
    # TODO qui probabilmente si possono tenere due file separati
    # Definining the structure-based atomtypes as the DataFrame index
    ensemble_contact_matrix[['idx_ai', 'idx_aj']] = ensemble_contact_matrix[['ai', 'aj']]
    ensemble_contact_matrix.set_index(['idx_ai', 'idx_aj'], inplace = True)

    reference_contact_matrix[['idx_ai', 'idx_aj']] = reference_contact_matrix[['rc_ai', 'rc_aj']]
    reference_contact_matrix.set_index(['idx_ai', 'idx_aj'], inplace = True)

    intra_ensemble_contact_matrix = ensemble_contact_matrix.loc[ensemble_contact_matrix['same_chain'] == 'Yes']
    intra_reference_contact_matrix = reference_contact_matrix.loc[reference_contact_matrix['rc_same_chain'] == 'Yes']
    inter_ensemble_contact_matrix = ensemble_contact_matrix.loc[ensemble_contact_matrix['same_chain'] == 'No']
    inter_reference_contact_matrix = reference_contact_matrix.loc[reference_contact_matrix['rc_same_chain'] == 'No']

    intra_ensemble_contacts_reweighted = pd.DataFrame()
    inter_ensemble_contacts_reweighted = pd.DataFrame()

    if not intra_ensemble_contact_matrix.empty:
        intra_ensemble_contacts_reweighted = reweight_contacts(intra_ensemble_contact_matrix, intra_reference_contact_matrix, parameters, simulation)
        print(f'\t\t- Intramolecular contact pairs added after removing duplicates: ', len(intra_ensemble_contacts_reweighted))
        print("\t\t\t- Epsilon input:", parameters.inter_epsilon) 
        print("\t\t\t- Average epsilon is", intra_ensemble_contacts_reweighted['epsilon'].loc[(intra_ensemble_contacts_reweighted['epsilon']>0.)].mean())
        print("\t\t\t- Maximum epsilon is", intra_ensemble_contacts_reweighted['epsilon'].max())

    if not inter_ensemble_contact_matrix.empty:
        inter_ensemble_contacts_reweighted = reweight_contacts(inter_ensemble_contact_matrix, inter_reference_contact_matrix, parameters, simulation)
        print(f'\t\t- Intermolecular contact pairs added after removing duplicates: ', len(inter_ensemble_contacts_reweighted))
        print("\t\t\t- Epsilon input:", parameters.inter_epsilon) 
        print("\t\t\t- Average epsilon is", inter_ensemble_contacts_reweighted['epsilon'].loc[(inter_ensemble_contacts_reweighted['epsilon']>0.)].mean())
        print("\t\t\t- Maximum epsilon is", inter_ensemble_contacts_reweighted['epsilon'].max())

    reweighted_contact_matrix = pd.concat([intra_ensemble_contacts_reweighted, inter_ensemble_contacts_reweighted], axis=0, sort = False, ignore_index = True)

    return reweighted_contact_matrix


def reweight_contacts(ensemble_contact_matrix, reference_contact_matrix, parameters, simulation):
    '''
    This function merge the contact matrix from any ensemble read and concat the random coil information as a comparison to reweight the LJ pairs accordingly.
    
    In the case of intramolecular contacts the rc_threshold is applied
    In the case of intermolecular contacts the md_threshold is applied TODO ma non dovrebbe corrispondere anche qui la rc con la inter_rc?

    '''

    ensemble_contact_matrix.to_csv(f'analysis/{simulation}_before_reweight')
    
    merged_ensemble_matrix = pd.concat([ensemble_contact_matrix, reference_contact_matrix], axis=1)
    merged_ensemble_matrix.drop(columns=['rc_ai', 'rc_aj'], inplace=True)
    merged_ensemble_matrix = merged_ensemble_matrix.loc[(merged_ensemble_matrix['probability'] > parameters.md_threshold)]
    # TODO questa non ricordo cosa faccia
    #merged_ensemble_matrix = merged_ensemble_matrix.loc[(merged_ensemble_matrix['probability']>parameters['md_threshold'])|((merged_ensemble_matrix['probability']<merged_ensemble_matrix['rc_probability'])&(merged_ensemble_matrix['probability']>parameters['rc_threshold']))]

    merged_ensemble_matrix['sigma'] = merged_ensemble_matrix['distance']/(2**(1/6))
    merged_ensemble_matrix['epsilon'] = np.nan 

    # Reweighting the epsilon based on probability using the Paissoni Equation 2.1
    # Attractive pairs. Those are the general rule where the probability from the MD is reduced based on the RC. High RC probabilities define a contact which does not require a strong LJ parametrization as it can be sampled without any contact but it is due to local geometry.
    merged_ensemble_matrix['epsilon'].loc[(merged_ensemble_matrix['probability']>1.2*merged_ensemble_matrix['rc_probability'])&(merged_ensemble_matrix['same_chain']=='Yes')] = -(parameters.inter_epsilon/np.log(parameters.rc_threshold))*(np.log(merged_ensemble_matrix['probability']/np.maximum(merged_ensemble_matrix['rc_probability'],parameters.rc_threshold)))
    merged_ensemble_matrix['epsilon'].loc[(merged_ensemble_matrix['probability']>1.2*merged_ensemble_matrix['rc_probability'])&(merged_ensemble_matrix['same_chain']=='No')] = -(parameters.inter_epsilon/np.log(parameters.md_threshold))*(np.log(merged_ensemble_matrix['probability']/np.maximum(merged_ensemble_matrix['rc_probability'],parameters.md_threshold)))
    # Repulsive pairs define contacts with higher RC probability compared to MD.
    merged_ensemble_matrix['epsilon'].loc[(merged_ensemble_matrix['probability']<(1./1.2)*merged_ensemble_matrix['rc_probability'])] = np.log(merged_ensemble_matrix['probability']/merged_ensemble_matrix['rc_probability'])*(np.minimum(merged_ensemble_matrix['distance'],merged_ensemble_matrix['rc_distance'])**12)

    # clean NaN and zeros 
    merged_ensemble_matrix.dropna(inplace=True)
    merged_ensemble_matrix = merged_ensemble_matrix[merged_ensemble_matrix.epsilon != 0]

    # Retrieving the ai and aj information from the index and reindexing back again
    merged_ensemble_matrix = merged_ensemble_matrix.reset_index()
    merged_ensemble_matrix[['ai', 'aj']] = merged_ensemble_matrix[['idx_ai', 'idx_aj']]
    merged_ensemble_matrix.set_index(['idx_ai', 'idx_aj'], inplace = True)

    # Here we just change the column order to be more readable in case of debugging
    merged_ensemble_matrix = merged_ensemble_matrix[['ai', 'aj', 'sigma', 'epsilon', 'distance', 'probability', 'same_chain', 'source', 'rc_distance', 'rc_probability', 'rc_same_chain', 'rc_source', 'rc_residue_ai', 'rc_residue_aj']]

    # Inverse pairs calvario
    # THE FOLLOWING FUNCTION MUST LIST ALL COLUMNS!
    # Here we make a copy by changing order between ai with aj
    inverse_ensemble_matrix = merged_ensemble_matrix[['aj', 'ai', 'sigma', 'epsilon', 'distance', 'probability', 'same_chain', 'source', 'rc_distance', 'rc_probability', 'rc_same_chain', 'rc_source', 'rc_residue_ai', 'rc_residue_aj']].copy()
    # Here we rename aj as ai and viceversa
    inverse_ensemble_matrix.columns = ['ai', 'aj', 'sigma', 'epsilon', 'distance', 'probability', 'same_chain', 'source', 'rc_distance', 'rc_probability', 'rc_same_chain', 'rc_source', 'rc_residue_ai', 'rc_residue_aj']
    # Here we concat the two dataframes. All the contacts are doubled but with the copy as ai and aj inverted. This is due to the struggle of inverse pairs removal
    merged_ensemble_matrix = pd.concat([merged_ensemble_matrix, inverse_ensemble_matrix], axis=0, sort = False, ignore_index = True)

    # TODO queste due suggerite da chatGPT
    #inverse_ensemble_matrix = merged_ensemble_matrix.reindex(columns=['aj', 'ai', 'sigma', 'epsilon', 'distance', 'probability', 'same_chain', 'source', 'rc_distance', 'rc_probability', 'rc_same_chain', 'rc_source', 'rc_residue_ai', 'rc_residue_aj'])
    #merged_ensemble_matrix = merged_ensemble_matrix.append(inverse_ensemble_matrix, ignore_index=True)

    # Here we sort all pairs based on sigma
    merged_ensemble_matrix.sort_values(by = ['ai', 'aj', 'sigma'], inplace = True)
    # Cleaning the duplicates
    merged_ensemble_matrix = merged_ensemble_matrix.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    # Removing the reverse duplicates
    cols = ['ai', 'aj']
    merged_ensemble_matrix[cols] = np.sort(merged_ensemble_matrix[cols].values, axis=1)
    merged_ensemble_matrix = merged_ensemble_matrix.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    # And finally reindex, just to be sure ^_____^
    merged_ensemble_matrix[['idx_ai', 'idx_aj']] = merged_ensemble_matrix[['ai', 'aj']]
    merged_ensemble_matrix.set_index(['idx_ai', 'idx_aj'], inplace=True)   

    return merged_ensemble_matrix


def merge_and_clean_LJ(topology, meGO_LJ_potential, type_c12_dict, type_q_dict, parameters):
    '''
    This function...
    '''
    print(topology)


    pass

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
