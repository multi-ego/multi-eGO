import src.atomtypes_definitions as atomtypes_definitions
from src.io import IO 
from src.topology import Topology

import glob
import pandas as pd
import parmed
import os
import numpy as np

class EnsembleClass:
    def __init__(self):
        pass

    def assign_molecule_type(self, molecule_type_dict, molecule_name, molecule_topology):
        first_aminoacid = molecule_topology.residues[0].name

        if first_aminoacid in atomtypes_definitions.aminoacids_list:
            molecule_type_dict[molecule_name] = 'protein'
        elif first_aminoacid in atomtypes_definitions.nucleic_acid_list:
            molecule_type_dict[molecule_name] = 'nucleic_acid'
        else:
            molecule_type_dict[molecule_name] = 'other'
        return molecule_type_dict


    def initialize_ensemble_topology(self, topology, simulation):
        '''
        '''
        # In a single topology different type of molecules can be present (e.g. protein, ligand).
        # For each molecule the different atomtypes are saved.
        print(
            '\t-', f'Reading {simulation} topology containing: {topology.molecules}')
        columns_to_drop = ['nb_idx', 'solvent_radius', 'screen', 'occupancy', 'bfactor',
                           'altloc', 'join', 'irotat', 'rmin', 'rmin_14', 'epsilon_14', 'tree']
        ensemble_topology_dataframe, new_number, col_molecule, new_resnum, ensemble_molecules_idx_sbtype_dictionary, temp_number_c12_dict = pd.DataFrame(), [], [], [], {}, {}

        molecule_type_dict = {}
        # I needed to add this for loop as by creating the topology dataframe by looping over molecules, the c12 information is lost
        for atom in topology.atoms:
            temp_number_c12_dict[str(atom.idx+1)] = atom.epsilon*4.184

        for molecule_number, (molecule_name, molecule_topology) in enumerate(topology.molecules.items(), 1):
            molecule_type_dict = self.assign_molecule_type(
                molecule_type_dict, molecule_name, molecule_topology[0])
            ensemble_molecules_idx_sbtype_dictionary[f'{str(molecule_number)}_{molecule_name}'] = {
            }
            ensemble_topology_dataframe = pd.concat(
                [ensemble_topology_dataframe, molecule_topology[0].to_dataframe()], axis=0)
            for atom in molecule_topology[0].atoms:
                new_number.append(str(atom.idx+1))
                col_molecule.append(f'{molecule_number}_{molecule_name}')
                new_resnum.append(str(atom.residue.number))
        del molecule_name

        ensemble_topology_dataframe['number'] = new_number
        ensemble_topology_dataframe['molecule'] = col_molecule
        ensemble_topology_dataframe['molecule_number'] = col_molecule
        ensemble_topology_dataframe[['molecule_number', 'molecule_name']
                                    ] = ensemble_topology_dataframe.molecule.str.split('_', expand=True)
        ensemble_topology_dataframe['resnum'] = new_resnum
        ensemble_topology_dataframe['cgnr'] = ensemble_topology_dataframe['resnum']
        ensemble_topology_dataframe['ptype'] = 'A'
        ensemble_topology_dataframe = ensemble_topology_dataframe.replace(
            {'name': atomtypes_definitions.from_ff_to_multiego})
        ensemble_topology_dataframe['sb_type'] = ensemble_topology_dataframe['name'] + '_' + \
            ensemble_topology_dataframe['molecule_name'] + '_' + \
            ensemble_topology_dataframe['resnum'].astype(str)
        ensemble_topology_dataframe.rename(
            columns={'epsilon': 'c12'}, inplace=True)

        ensemble_topology_dataframe['charge'] = 0.
        ensemble_topology_dataframe['c6'] = 0.
        ensemble_topology_dataframe['c12'] = ensemble_topology_dataframe['number'].map(
            temp_number_c12_dict)
        ensemble_topology_dataframe['molecule_type'] = ensemble_topology_dataframe['molecule_name'].map(
            molecule_type_dict)

        for molecule in ensemble_molecules_idx_sbtype_dictionary.keys():
            temp_topology_dataframe = ensemble_topology_dataframe.loc[
                ensemble_topology_dataframe['molecule'] == molecule]
            number_sbtype_dict = temp_topology_dataframe[[
                'number', 'sb_type']].set_index('number')['sb_type'].to_dict()
            ensemble_molecules_idx_sbtype_dictionary[molecule] = number_sbtype_dict
        sbtype_c12_dict = ensemble_topology_dataframe[[
            'sb_type', 'c12']].set_index('sb_type')['c12'].to_dict()

        return ensemble_topology_dataframe, ensemble_molecules_idx_sbtype_dictionary, sbtype_c12_dict, molecule_type_dict


    # def read_files(self):
    #     '''
    #     This function search the .top files and the contact matrices .ndx for each folder defined in the command line.
    #     It is important to call properly the names.
    #     '''
    #     top_path = glob.glob(f'inputs/{self.simulation_path}/*.top')
    #     print('\t-', f'Reading {top_path[0]}')
    #     self.topology = parmed.load_file(top_path[0])
    #     # Reading contact matrix created using gmx_clustsize
    #     # Reference requires both intra and inter molecular matrices
    #     if self.args.egos == 'rc':
    #         self.ensemble_contact_matrices = pd.DataFrame()
    #     else:
    #         matrix = None
    #         for matrix in glob.glob(f'inputs/{self.simulation_path}/int??mat_?_?.ndx'):
    #             name = matrix.replace(f'inputs/{self.simulation_path}/', '')
    #             self.ensemble_contact_matrices[name] = io.read_molecular_contacts(
    #                 matrix)
    #         if not matrix:
    #             raise FileNotFoundError('.ndx files must be named as intramat_X_X.ndx or intermat_1_1.ndx')


    def initialize_molecular_contacts(self, contact_matrices, ensemble_molecules_idx_sbtype_dictionary, simulation):
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
                'molecule_name_ai', 'ai', 'molecule_name_aj', 'aj', 'distance_m', 'distance', 'distance_14', 'probability', 'cutoff']]
            if name[0] == 'intramat':
                contact_matrix['same_chain'] = True
            elif name[0] == 'intermat':
                contact_matrix['same_chain'] = False
            else:
                raise Exception(
                    'There might be an error in the contact matrix naming. It must be intermat_X_X or intramat_X_X')

            contact_matrix['source'] = simulation
            contact_matrix['file'] = '_'.join(name)
            ensemble_contact_matrix = pd.concat(
                [ensemble_contact_matrix, contact_matrix], axis=0)
            ensemble_contact_matrix[['idx_ai', 'idx_aj']
                                    ] = ensemble_contact_matrix[['ai', 'aj']]
            ensemble_contact_matrix.set_index(['idx_ai', 'idx_aj'], inplace=True)
        return ensemble_contact_matrix


    def initialize_ensemble(self, simulation_path, egos):
        '''
        This function creates an ensemble dictionary containing various data.
        Contents of the dictionary include:
            - the simulation name from which the data was calculated
            - the path to said simunlation
            - the parmed.topology of the used for the simulation
            - e cazzi vari

        Parameters
        ----------
        simulation_path : string
            The path to the directory of the simulation
        egos : string
            The type of simulation. Either 'rc' for random coil or 'production' for a 
            full-fledged multi-eGO simulation

        Returns
        -------
        ensemble : dict
            The dictionary containing the simulation data
        '''
        ensemble_type = simulation_path.split('/')[-1]
        print('\t-', f'Initializing {ensemble_type} ensemble topology')
        # top_path = glob.glob(f'inputs/{self.simulation_path}/*.top')
        topology_path = f'inputs/{simulation_path}/topol.top'
        if not os.path.isfile(topology_path):
            raise FileNotFoundError(f"{topology_path} not found.")
        print('\t-', f'Reading {topology_path}')
        topology = parmed.load_file(topology_path)
        ensemble_contact_matrices = {}
        if egos == 'rc':
            pass
        else:
            matrix_paths = glob.glob(f'inputs/{simulation_path}/int??mat_?_?.ndx')
            if matrix_paths == []: 
                raise FileNotFoundError('.ndx files must be named as intramat_X_X.ndx or intermat_1_1.ndx')
            for path in matrix_paths:
                name = path.replace(f'inputs/{simulation_path}/', '')
                ensemble_contact_matrices[name] = IO.read_molecular_contacts(path)

        ensemble_topology_dataframe, ensemble_molecules_idx_sbtype_dictionary, sbtype_c12_dict, molecule_type_dict = self.initialize_ensemble_topology(topology, ensemble_type)
        atomic_contacts = self.initialize_molecular_contacts(ensemble_contact_matrices, ensemble_molecules_idx_sbtype_dictionary, ensemble_type)

        ensemble = {}
        ensemble['simulation'] = ensemble_type
        ensemble['simulation_path'] = simulation_path
        ensemble['topology'] = topology
        ensemble['ensemble_topology_dataframe'] = ensemble_topology_dataframe
        ensemble['ensemble_molecules_idx_sbtype_dictionary'] = ensemble_molecules_idx_sbtype_dictionary
        ensemble['ensemble_contact_matrices'] = ensemble_contact_matrices
        ensemble['sbtype_c12_dict'] = sbtype_c12_dict
        ensemble['molecule_type_dict'] = molecule_type_dict
        ensemble['atomic_contacts'] = atomic_contacts

        return ensemble
    
    def add_ensemble_from(self, meGO_ensemble, ensemble, check_with):
        # TODO da aggiungere il dizionario di conversione delle topologie!!!
        print('\t-', f'Adding topology from {ensemble["simulation"]}')

        if ensemble['simulation'] == 'reference':
            # This defines as the reference structure and eventual molecules will be added
            if not 'reference_topology_dataframe' in ensemble.keys(): meGO_ensemble['reference_topology_dataframe'] = pd.DataFrame()
            meGO_ensemble['reference_topology'] = ensemble['topology']
            meGO_ensemble['reference_topology_dataframe'] = pd.concat([meGO_ensemble['reference_topology_dataframe'], ensemble['ensemble_topology_dataframe']], axis=0, ignore_index=True)
            meGO_ensemble['sbtype_c12_dict'] = ensemble['sbtype_c12_dict'] # WARNING redundant?
            meGO_ensemble['sbtype_number_dict'] = meGO_ensemble['reference_topology_dataframe'][['sb_type', 'number']].set_index('sb_type')['number'].to_dict()
            meGO_ensemble['reference_atomic_contacts'] = ensemble['atomic_contacts'].add_prefix('rc_')
            meGO_ensemble['molecule_type_dict'] = ensemble['molecule_type_dict'] # WARNING redundant?
        
        elif ensemble['simulation'] in check_with:
            if not 'check_atomic_contacts' in meGO_ensemble.keys(): meGO_ensemble['check_atomic_contacts'] = pd.DataFrame()
            meGO_ensemble['check_atomic_contacts'] = pd.concat([meGO_ensemble['check_atomic_contacts'], ensemble['atomic_contacts']], axis=0, ignore_index=True)
        
        else:
            if not 'meGO_topology_dataframe' in meGO_ensemble.keys(): meGO_ensemble['meGO_topology_dataframe'] = pd.DataFrame()
            if not 'meGO_atomic_contacts' in meGO_ensemble.keys(): meGO_ensemble['meGO_atomic_contacts'] = pd.DataFrame()
            meGO_ensemble['meGO_topology_dataframe'] = pd.concat([meGO_ensemble['meGO_topology_dataframe'], ensemble['ensemble_topology_dataframe']], axis=0, ignore_index=True)
            meGO_ensemble['meGO_atomic_contacts'] = pd.concat([meGO_ensemble['meGO_atomic_contacts'], ensemble['atomic_contacts']], axis=0)

        return meGO_ensemble

    def check_topology_conversion(self, meGO_ensemble, egos):
        '''
        This function is required to check the different atomtypes between different force fields.
        The atom types MUST match otherwise a proper ffnobonded cannot be created.
        


        This function is called "check_topology_conversion" and it is a method of a class.
        It does the following:
            Initializes an empty set called "reference_set" and fills it with the unique 'name' values of the "reference_topology_dataframe" attribute of the class.
            Loops through the "molecules" attribute of the "reference_topology" attribute of the class.
            For each molecule, it creates a new DataFrame called "comparison_dataframe" by filtering the "meGO_topology_dataframe" attribute of the class by the molecule number and name (e.g. "1_protein").
            If the "comparison_dataframe" is not empty, it creates a new set called "comparison_set" and fills it with the unique 'name' values of the "comparison_dataframe", after removing all rows with 'name' value starting with H
            Lastly, it finds the difference between the two sets (comparison_set and reference_set) and store it in difference_set.
            If difference_set is not empty, it prints a message indicating that the atomtypes in difference_set are not converted and that they must be added to the "from_ff_to_multiego" dictionary to properly merge all the contacts and exits the program.
        This function is checking if there are any atom types present in the reference topology that are not present in the meGO topology and if there are any it exits the program.
        '''
        if egos != 'rc':
            reference_set = set(meGO_ensemble['reference_topology_dataframe']['name'].to_list())
            for number, molecule in enumerate(meGO_ensemble['reference_topology'].molecules, 1):
                comparison_dataframe = meGO_ensemble['meGO_topology_dataframe'].loc[meGO_ensemble['meGO_topology_dataframe']['molecule'] == f'{number}_{molecule}']
                if not comparison_dataframe.empty:
                    comparison_set = set(comparison_dataframe[~comparison_dataframe['name'].astype(str).str.startswith('H')]['name'].to_list())
            difference_set = comparison_set.difference(reference_set)
            if difference_set:
                print(f'The following atomtypes are not converted:\n{difference_set} \nYou MUST add them in "from_ff_to_multiego" dictionary to properly merge all the contacts.')
                exit()

    def generate_bonded_interactions(self, meGO_ensemble):#, reference_topology):
        # meGO_bonded_interactions = {}
        # bond_pairs = {}
        # user_pairs = {}
        
        if 'meGO_bonded_interactions' not in meGO_ensemble.keys(): meGO_ensemble['meGO_bonded_interactions'] = {}
        if 'bond_pairs' not in meGO_ensemble.keys(): meGO_ensemble['bond_pairs'] = {}
        if 'user_pairs' not in meGO_ensemble.keys(): meGO_ensemble['user_pairs'] = {}

        for molecule, topol in meGO_ensemble['reference_topology'].molecules.items():      
            meGO_ensemble['meGO_bonded_interactions'][molecule] = {
                'bonds' : Topology.get_bonds(topol[0].bonds),
                'angles' : Topology.get_angles(topol[0].angles),
                'dihedrals' : Topology.get_dihedrals(topol[0].dihedrals),
                'impropers' : Topology.get_impropers(topol[0].impropers),
                'pairs' : Topology.get_pairs(topol[0].adjusts)
            }
            # The following bonds are used in the parametrization of LJ 1-4
            meGO_ensemble['bond_pairs'][molecule] = Topology.get_bond_pairs(topol[0].bonds)
            meGO_ensemble['user_pairs'][molecule] = Topology.get_pairs(topol[0].adjusts)

        return meGO_ensemble


    # def parametrize_LJ(self, topology_dataframe, molecule_type_dict, bond_tuple, pairs_tuple, type_c12_dict, meGO_atomic_contacts, reference_atomic_contacts, check_atomic_contacts, sbtype_number_dict, parameters):
    def parametrize_LJ(self, meGO_ensemble, parameters):
        '''
        This function reads the probabilities obtained using gmx_clustsize from the ensembles defined in the command line.
        The random coil probabilities are used to reweight the explicit water ones.
        Intra and inter molecular contacts are splitted as different rules are applied during the reweighting.
        For each atom contact the sigma and epsilon are obtained.
        '''

        # First of all we generate the random-coil 1-4 interactions:
        pairs14 = pd.DataFrame()
        exclusion_bonds14 = pd.DataFrame()

        for molecule, bond_pair in meGO_ensemble['bond_pairs'].items():
            reduced_topology = meGO_ensemble['reference_topology_dataframe'].loc[meGO_ensemble['reference_topology_dataframe']['molecule_name'] == molecule][['number', 'sb_type', 'resnum', 'name', 'type', 'resname', 'molecule_type']].copy()

            reduced_topology['number'] = reduced_topology['number'].astype(str)
            reduced_topology['resnum'] = reduced_topology['resnum'].astype(int)
            # Dictionaries definitions to map values
            atnum_type_dict = reduced_topology.set_index('sb_type')['number'].to_dict()
            type_atnum_dict = reduced_topology.set_index('number')['sb_type'].to_dict()

            # Building the exclusion bonded list
            # exclusion_bonds are all the interactions within 3 bonds
            # p14 are specifically the interactions at exactly 3 bonds
            exclusion_bonds, tmp_p14 = Topology.get_14_interaction_list(reduced_topology, bond_pair)
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
                pairs = Topology.protein_LJ14(reduced_topology)
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

        if parameters.egos != 'rc':
            meGO_atomic_contacts_merged = pd.merge(meGO_ensemble['meGO_atomic_contacts'], meGO_ensemble['reference_atomic_contacts'], left_index=True, right_index=True, how='outer')
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
            meGO_atomic_contacts_merged['rep'] = meGO_atomic_contacts_merged['rep'].fillna(np.sqrt(meGO_atomic_contacts_merged['ai'].map(meGO_ensemble['sbtype_c12_dict'])*meGO_atomic_contacts_merged['aj'].map(meGO_ensemble['sbtype_c12_dict'])))

            # this is the minimum probability ratio for attractive contacts and is set so that epsilon should not be smaller than 0.1*epsilon
            limit_rc = 1./parameters.rc_threshold**0.1
            # This keep only significat attractive/repulsive interactions
            meGO_atomic_contacts_merged = meGO_atomic_contacts_merged.loc[(meGO_atomic_contacts_merged['probability']>parameters.md_threshold)|((meGO_atomic_contacts_merged['probability']<=parameters.md_threshold)&(meGO_atomic_contacts_merged['rc_probability']>limit_rc*np.maximum(meGO_atomic_contacts_merged['probability'],parameters.rc_threshold)))]
            # This cleans contacts with just a small tailed histogram
            meGO_atomic_contacts_merged = meGO_atomic_contacts_merged.loc[(meGO_atomic_contacts_merged['distance']<meGO_atomic_contacts_merged['cutoff']-0.02)|(meGO_atomic_contacts_merged['probability']>0.01)] 

            # Add sigma, add epsilon reweighted, add c6 and c12
            meGO_atomic_contacts_merged['sigma'] = (meGO_atomic_contacts_merged['distance']) / (2.**(1./6.))
            meGO_atomic_contacts_merged['epsilon'] = np.nan 

            # The index has been reset as here I have issues with multiple index duplicates. The same contact is kept twice: one for intra and one for inter.
            # The following pandas functions cannot handle multiple rows with the same index although it has been defined the "same_chain" filter.
            meGO_atomic_contacts_merged.reset_index(inplace=True)


            # Epsilon reweight based on probability
            # Paissoni Equation 2.1

            # Attractive intramolecular
            meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']>limit_rc*np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold))&(meGO_atomic_contacts_merged['same_chain']==True)] = -(parameters.epsilon/np.log(parameters.rc_threshold))*(np.log(meGO_atomic_contacts_merged['probability']/np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold)))
            # Attractive intermolecular
            meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']>limit_rc*np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold))&(meGO_atomic_contacts_merged['same_chain']==False)] = -(parameters.inter_epsilon/np.log(parameters.rc_threshold))*(np.log(meGO_atomic_contacts_merged['probability']/np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold)))

            # Probability only Repulsive intramolecular
            meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']<1./limit_rc*np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold))&(meGO_atomic_contacts_merged['same_chain']==True)&(meGO_atomic_contacts_merged['probability']<=parameters.md_threshold)] = -(parameters.epsilon/np.log(parameters.rc_threshold))*meGO_atomic_contacts_merged['cutoff']**12*np.log(np.maximum(meGO_atomic_contacts_merged['probability'],parameters.rc_threshold)/meGO_atomic_contacts_merged['rc_probability'])-meGO_atomic_contacts_merged['rep']
            # Probability only Repulsive intermolecular
            meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']<1./limit_rc*np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold))&(meGO_atomic_contacts_merged['same_chain']==False)&(meGO_atomic_contacts_merged['probability']<=parameters.md_threshold)] = -(parameters.inter_epsilon/np.log(parameters.rc_threshold))*meGO_atomic_contacts_merged['cutoff']**12*np.log(np.maximum(meGO_atomic_contacts_merged['probability'],parameters.rc_threshold)/meGO_atomic_contacts_merged['rc_probability'])-meGO_atomic_contacts_merged['rep']

            # Full Repulsive intramolecular
            meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']<1./limit_rc*np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold))&(meGO_atomic_contacts_merged['probability']>parameters.md_threshold)&(meGO_atomic_contacts_merged['rc_probability']>parameters.md_threshold)&(meGO_atomic_contacts_merged['same_chain']==True)&((meGO_atomic_contacts_merged['rc_distance_m']-meGO_atomic_contacts_merged['distance_m'])<0.)] = -(parameters.epsilon/np.log(parameters.rc_threshold))*meGO_atomic_contacts_merged['distance_m']**12*np.log(meGO_atomic_contacts_merged['probability']/meGO_atomic_contacts_merged['rc_probability'])-(meGO_atomic_contacts_merged['rep']*(meGO_atomic_contacts_merged['distance_m']**12/meGO_atomic_contacts_merged['rc_distance_m']**12))
            # Full Repulsive intermolecular
            meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']<1./limit_rc*np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold))&(meGO_atomic_contacts_merged['probability']>parameters.md_threshold)&(meGO_atomic_contacts_merged['rc_probability']>parameters.md_threshold)&(meGO_atomic_contacts_merged['same_chain']==False)&((meGO_atomic_contacts_merged['rc_distance_m']-meGO_atomic_contacts_merged['distance_m'])<0.)] = -(parameters.inter_epsilon/np.log(parameters.rc_threshold))*meGO_atomic_contacts_merged['distance_m']**12*np.log(meGO_atomic_contacts_merged['probability']/meGO_atomic_contacts_merged['rc_probability'])-(meGO_atomic_contacts_merged['rep']*(meGO_atomic_contacts_merged['distance_m']**12/meGO_atomic_contacts_merged['rc_distance_m']**12))

            # mild c12s (for small probability ratios): smaller
            meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']<=limit_rc*np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold))&(meGO_atomic_contacts_merged['probability']>meGO_atomic_contacts_merged['rc_probability'])&((meGO_atomic_contacts_merged['rc_distance_m']-meGO_atomic_contacts_merged['distance_m'])>0.)&(meGO_atomic_contacts_merged['probability']>parameters.md_threshold)&(meGO_atomic_contacts_merged['rc_probability']>parameters.md_threshold)] = -meGO_atomic_contacts_merged['rep']*(meGO_atomic_contacts_merged['distance_m']/meGO_atomic_contacts_merged['rc_distance_m'])**12 
            # mild c12s (for small probability ratios): larger 
            meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['probability']>=1./limit_rc*np.maximum(meGO_atomic_contacts_merged['rc_probability'],parameters.rc_threshold))&(meGO_atomic_contacts_merged['probability']<=meGO_atomic_contacts_merged['rc_probability'])&((meGO_atomic_contacts_merged['rc_distance_m']-meGO_atomic_contacts_merged['distance_m'])<0.)&(meGO_atomic_contacts_merged['rc_probability']>parameters.md_threshold)&(meGO_atomic_contacts_merged['probability']>parameters.md_threshold)] = -meGO_atomic_contacts_merged['rep']*(meGO_atomic_contacts_merged['distance_m']/meGO_atomic_contacts_merged['rc_distance_m'])**12

            # This set the default (random coil) value for selected 1-4 interactions
            meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['1-4'] == '1_4')] = -meGO_atomic_contacts_merged['rep']
            # Rescale c12 1-4 interactions 
            meGO_atomic_contacts_merged['epsilon'].loc[(np.abs(meGO_atomic_contacts_merged['rc_distance_14']-meGO_atomic_contacts_merged['distance_14'])>0.)&(meGO_atomic_contacts_merged['rc_probability']>parameters.md_threshold)&(meGO_atomic_contacts_merged['1-4']=="1_4")&(meGO_atomic_contacts_merged['same_chain']==True)] = -meGO_atomic_contacts_merged['rep']*(meGO_atomic_contacts_merged['distance_14']/meGO_atomic_contacts_merged['rc_distance_14'])**12

            # Here we are reindexing like before
            meGO_atomic_contacts_merged[['idx_ai', 'idx_aj']] = meGO_atomic_contacts_merged[['ai', 'aj']]
            meGO_atomic_contacts_merged.set_index(['idx_ai', 'idx_aj'], inplace=True)

            # clean NaN and zeros
            meGO_atomic_contacts_merged.dropna(subset=['epsilon'], inplace=True)
            meGO_atomic_contacts_merged = meGO_atomic_contacts_merged[meGO_atomic_contacts_merged.epsilon != 0]

            # remove unnecessary fields
            meGO_atomic_contacts_merged = meGO_atomic_contacts_merged[['molecule_name_ai', 'ai', 'molecule_name_aj', 'aj', 
            'probability', 'same_chain', 'source', 'file', 'rc_probability', 'rc_file', 'sigma', 'epsilon', '1-4', 'distance', 'distance_m']]
            # Inverse pairs calvario
            # this must list ALL COLUMNS!
            inverse_meGO_atomic_contacts_merged = meGO_atomic_contacts_merged[['molecule_name_aj', 'aj', 'molecule_name_ai', 'ai',
            'probability', 'same_chain', 'source', 'file', 'rc_probability', 'rc_file', 'sigma', 'epsilon', '1-4', 'distance', 'distance_m']].copy()
            inverse_meGO_atomic_contacts_merged.columns = ['molecule_name_ai', 'ai', 'molecule_name_aj', 'aj',
            'probability', 'same_chain', 'source', 'file', 'rc_probability', 'rc_file', 'sigma', 'epsilon', '1-4', 'distance', 'distance_m']
            # The contacts are duplicated before cleaning due to the inverse pairs and the sigma calculation requires a simmetric dataframe
            meGO_atomic_contacts_merged = pd.concat([meGO_atomic_contacts_merged, inverse_meGO_atomic_contacts_merged], axis=0, sort=False, ignore_index=True)

            # process check_atomic_contacts
            if 'check_atomic_contacts' is meGO_ensemble.keys():
                if not meGO_ensemble['check_atomic_contacts'] == {}:
                    # Remove low probability ones
                    meGO_ensemble['check_atomic_contacts'] = meGO_ensemble['check_atomic_contacts'].loc[(meGO_ensemble['check_atomic_contacts']['probability']>limit_rc*parameters.md_threshold)] 
                    meGO_atomic_contacts_merged = pd.concat([meGO_atomic_contacts_merged, meGO_ensemble['check_atomic_contacts']], axis=0, sort=False, ignore_index=True)
                    meGO_atomic_contacts_merged.drop_duplicates(inplace=True, ignore_index = True)
                    # this calculates the increase in energy (if any) to reach the "check" configuration
                    energy_at_check_dist = meGO_atomic_contacts_merged.groupby(by=['ai', 'aj', 'same_chain'])[['distance', 'epsilon', 'source', 'same_chain']].apply(self.check_LJ, parameters)
                    meGO_atomic_contacts_merged = pd.merge(meGO_atomic_contacts_merged, energy_at_check_dist.rename('energy_at_check_dist'), how="inner", on=["ai", "aj", "same_chain"])
                    ## remove check_with contacts 
                    meGO_atomic_contacts_merged=meGO_atomic_contacts_merged[~meGO_atomic_contacts_merged.source.isin(parameters.check_with)]
                    ## rescale problematic contacts
                    meGO_atomic_contacts_merged['epsilon'].loc[(meGO_atomic_contacts_merged['same_chain']==True)&(meGO_atomic_contacts_merged['energy_at_check_dist']>parameters.epsilon)&(meGO_atomic_contacts_merged['1-4']!="1_4")] *= parameters.epsilon/meGO_atomic_contacts_merged['energy_at_check_dist']
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
            # remove them from the default force-field
            meGO_atomic_contacts_merged = meGO_atomic_contacts_merged.loc[(meGO_atomic_contacts_merged['1-4']!='1_4')]

            meGO_atomic_contacts_merged['c6'] = 4 * meGO_atomic_contacts_merged['epsilon'] * (meGO_atomic_contacts_merged['sigma'] ** 6)
            meGO_atomic_contacts_merged['c12'] = abs(4 * meGO_atomic_contacts_merged['epsilon'] * (meGO_atomic_contacts_merged['sigma'] ** 12))
            meGO_atomic_contacts_merged['c6'].loc[(meGO_atomic_contacts_merged['epsilon']<0.)] = 0.
            meGO_atomic_contacts_merged['c12'].loc[(meGO_atomic_contacts_merged['epsilon']<0.)] = -meGO_atomic_contacts_merged['epsilon']

            meGO_LJ_14['c6'] = 4 * meGO_LJ_14['epsilon'] * (meGO_LJ_14['sigma'] ** 6)
            meGO_LJ_14['c12'] = abs(4 * meGO_LJ_14['epsilon'] * (meGO_LJ_14['sigma'] ** 12))
            # repulsive interactions have just a very large C12
            meGO_LJ_14['c6'].loc[(meGO_LJ_14['epsilon']<0.)] = 0.
            meGO_LJ_14['c12'].loc[(meGO_LJ_14['epsilon']<0.)] =-meGO_LJ_14['epsilon']  

            meGO_atomic_contacts_merged['type'] = 1
            meGO_atomic_contacts_merged['number_ai'] = meGO_atomic_contacts_merged['ai'].map(meGO_ensemble['sbtype_number_dict'])
            meGO_atomic_contacts_merged['number_aj'] = meGO_atomic_contacts_merged['aj'].map(meGO_ensemble['sbtype_number_dict'])
            meGO_atomic_contacts_merged['number_ai'] = meGO_atomic_contacts_merged['number_ai'].astype(int)
            meGO_atomic_contacts_merged['number_aj'] = meGO_atomic_contacts_merged['number_aj'].astype(int)

            meGO_atomic_contacts_merged = meGO_atomic_contacts_merged[['ai', 'aj', 'type', 'c6', 'c12', 'sigma', 'epsilon', 'probability', 'rc_probability', 'molecule_name_ai',  'molecule_name_aj', 'same_chain', 'source', 'file', 'rc_file', 'number_ai', 'number_aj']]
            # Here we want to sort so that ai is smaller than aj
            inv_meGO = meGO_atomic_contacts_merged[['aj', 'ai', 'type', 'c6', 'c12', 'sigma', 'epsilon', 'probability', 'rc_probability', 'molecule_name_aj',  'molecule_name_ai', 'same_chain', 'source', 'file', 'rc_file', 'number_aj', 'number_ai']].copy()
            inv_meGO.columns = ['ai', 'aj', 'type', 'c6', 'c12', 'sigma', 'epsilon', 'probability', 'rc_probability', 'molecule_name_ai',  'molecule_name_aj', 'same_chain', 'source', 'file', 'rc_file', 'number_ai', 'number_aj'] 
            meGO_atomic_contacts_merged = pd.concat([meGO_atomic_contacts_merged,inv_meGO], axis=0, sort = False, ignore_index = True)
            meGO_atomic_contacts_merged = meGO_atomic_contacts_merged[meGO_atomic_contacts_merged['number_ai']<=meGO_atomic_contacts_merged['number_aj']]
            meGO_atomic_contacts_merged.sort_values(by = ['number_ai', 'number_aj'], inplace = True)
            meGO_atomic_contacts_merged = meGO_atomic_contacts_merged.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')

        else:
            meGO_atomic_contacts_merged = pd.DataFrame()
            meGO_LJ_14 = pairs14
            meGO_LJ_14['epsilon'] = -meGO_LJ_14['c12'] 

        return meGO_atomic_contacts_merged, meGO_LJ_14


    def make_pairs_exclusion_topology(self, meGO_ensemble, meGO_LJ_14):#topology_dataframe, bond_tuple, type_c12_dict, meGO_LJ_14, parameters):
        '''
        This function prepares the [ exclusion ] and [ pairs ] section to paste in topology.top
        '''
        pairs_molecule_dict = {}
        for molecule, bond_pair in meGO_ensemble['bond_pairs'].items():
            reduced_topology = meGO_ensemble['reference_topology_dataframe'].loc[meGO_ensemble['reference_topology_dataframe']['molecule_name'] == molecule][['number', 'sb_type', 'resnum', 'name', 'type', 'resname', 'molecule_type']].copy()

            reduced_topology['number'] = reduced_topology['number'].astype(str)
            reduced_topology['resnum'] = reduced_topology['resnum'].astype(int)

            atnum_type_dict = reduced_topology.set_index('sb_type')['number'].to_dict()
            # Building the exclusion bonded list
            # exclusion_bonds are all the interactions within 3 bonds
            # p14 are specifically the interactions at exactly 3 bonds
            exclusion_bonds, p14 = Topology.get_14_interaction_list(reduced_topology, bond_pair) 
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

    def generate_LJ_potential(self, meGO_ensemble, parameters):
        '''
        This function merges all the LJ contacts from ensembles into a single one.
        All contacts are reweighted based on the RC probability.
        Duplicates are removed.
        '''
        # meGO_LJ_potential, meGO_LJ_14 = self.parametrize_LJ(self.reference_topology_dataframe, self.molecule_type_dict, self.bond_pairs, self.user_pairs, self.sbtype_c12_dict, self.meGO_atomic_contacts, self.reference_atomic_contacts, self.check_atomic_contacts, self.sbtype_number_dict, self.parameters)
        meGO_LJ_potential, meGO_LJ_14 = self.parametrize_LJ(meGO_ensemble, parameters)
        # meGO_LJ_14 = self.make_pairs_exclusion_topology(meGO_ensemble['reference_topology_dataframe'], meGO_ensemble['bond_pairs'], meGO_ensemble['sbtype_c12_dict'], meGO_LJ_14, parameters)
        meGO_LJ_14 = self.make_pairs_exclusion_topology(meGO_ensemble, meGO_LJ_14)

        return meGO_LJ_potential, meGO_LJ_14

Ensemble = EnsembleClass()
