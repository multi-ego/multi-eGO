import os
import pandas as pd
from vanessa import get_bonds, get_bond_pairs, get_angles, get_dihedrals, get_impropers, parametrize_LJ, make_pairs_exclusion_topology
from write_output import get_name, write_topology, write_nonbonded

class Multi_eGO_Ensemble:

    def __init__(self, args):
        self.parameters = args
        self.output_folder = f'outputs/{get_name(args)}'
        try:
            os.mkdir('outputs')
        except OSError as error:
            pass
        try:
            os.mkdir(self.output_folder)
        except OSError as error:
            pass

        # TODO inserire la possibilità di leggere un RC fatto da più di una molecola
        self.reference_topology = None
        self.reference_topology_dataframe = pd.DataFrame()
        self.meGO_bonded_interactions = {}
        self.bond_pairs = {}
        self.reference_atomic_contacts = pd.DataFrame()

        self.meGO_topology_dataframe = pd.DataFrame()
        self.meGO_atomic_contacts = pd.DataFrame()
        self.check_atomic_contacts = None
        self.meGO_LJ_potential = pd.DataFrame()
        self.meGO_LJ_14 = pd.DataFrame()
        self.sbtype_c12_dict = None
        self.sbtype_number_dict = None

        #reference_n_residues = None
        #topology_dataframe = pd.DataFrame()
        
        # Contact matrices definition
        #contact_matrices_dictionary = {}

    def add_ensemble_from(self, ensemble):
        '''
        This method...
        This function inserts in a dictionary all the contact matrices learned from different simulations
        '''
        # TODO da aggiungere il dizionario di conversione delle topologie!!!
        print('\t-', f'Adding topology from {ensemble.simulation}')

        if ensemble.simulation == 'reference':
            # This defines as the reference structure and eventual molecules will be added
            self.reference_topology = ensemble.topology
            self.reference_topology_dataframe = pd.concat([self.reference_topology_dataframe, ensemble.ensemble_topology_dataframe], axis=0, ignore_index=True)
            self.sbtype_c12_dict = ensemble.sbtype_c12_dict
            self.sbtype_number_dict = self.reference_topology_dataframe[['sb_type', 'number']].set_index('sb_type')['number'].to_dict()
            self.reference_atomic_contacts = ensemble.atomic_contacts.add_prefix('rc_')
        
        elif ensemble.simulation in self.parameters.check_with:
            self.check_atomic_contacts = pd.concat([self.check_atomic_contacts, ensemble.atomic_contacts], axis=0, ignore_index=True)
        
        else:
            self.meGO_topology_dataframe = pd.concat([self.meGO_topology_dataframe, ensemble.ensemble_topology_dataframe], axis=0, ignore_index=True)
            self.meGO_atomic_contacts = pd.concat([self.meGO_atomic_contacts, ensemble.atomic_contacts], axis=0)


    def check_topology_conversion(self):
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
        if self.parameters.egos != 'rc':
            reference_set = set(self.reference_topology_dataframe['name'].to_list())
            for number, molecule in enumerate(self.reference_topology.molecules, 1):
                comparison_dataframe = self.meGO_topology_dataframe.loc[self.meGO_topology_dataframe['molecule'] == f'{number}_{molecule}']
                if not comparison_dataframe.empty:
                    comparison_set = set(comparison_dataframe[~comparison_dataframe['name'].astype(str).str.startswith('H')]['name'].to_list())
            difference_set = comparison_set.difference(reference_set)
            if difference_set:
                print(f'The following atomtypes are not converted:\n{difference_set} \nYou MUST add them in "from_ff_to_multiego" dictionary to properly merge all the contacts.')
                exit()


    def generate_bonded_interactions(self):
        '''
        '''
        for molecule, topol in self.reference_topology.molecules.items():
            self.meGO_bonded_interactions[molecule] = {
                'bonds' : get_bonds(topol[0].bonds),
                'angles' : get_angles(topol[0].angles),
                'dihedrals' : get_dihedrals(topol[0].dihedrals),
                'impropers' : get_impropers(topol[0].impropers)
            }
            # The following bonds are used in the parametrization of LJ 1-4
            self.bond_pairs[molecule] = get_bond_pairs(topol[0].bonds)


    def generate_LJ_potential(self):
        '''
        This function merges all the LJ contacts from ensembles into a single one.
        All contacts are reweighted based on the RC probability.
        Duplicates are removed.
        '''
        # TODO qui si mette parametrize_LJ che equivale a MD_LJ_pairs
        # All contacts are reweighted by the random coil probability both as intra and intermolecular and added to the LJ pairs of multi-eGO.
        self.meGO_LJ_potential, self.meGO_LJ_14 = parametrize_LJ(self.reference_topology_dataframe, self.meGO_atomic_contacts, self.reference_atomic_contacts, self.check_atomic_contacts, self.sbtype_c12_dict, self.sbtype_number_dict, self.parameters)
        self.meGO_LJ_14 = make_pairs_exclusion_topology(self.reference_topology_dataframe, self.bond_pairs, self.sbtype_c12_dict, self.parameters, self.meGO_LJ_14)
        # TODO controllare bene i numeri che sono troppi pairs in uscita!!! (i c12 in input potrebbero essere il problema)
        #print(self.meGO_LJ_14)

    def write_model(self):
        '''        
        '''
        print('- Writing Multi-eGO model')
        write_topology(self.reference_topology_dataframe, self.meGO_bonded_interactions, self.meGO_LJ_14, self.parameters, self.output_folder)
        write_nonbonded(self.reference_topology_dataframe, self.meGO_LJ_potential, self.parameters, self.output_folder)
        
        print('\n- The model is baked with the following parameters:\n')
        for argument, value in vars(self.parameters).items():
            if type(value) is list:
                print('\t- {:<20} = {:<20}'.format(argument, ", ".join(value)))
            elif type(value) is not str:
                print('\t- {:<20} = {:<20}'.format(argument, str(value)))
            else:
                print('\t- {:<20} = {:<20}'.format(argument, value))
        print(f'\nAnd it can be found in the following folder:\n{self.output_folder}')
        print('\nNessuno è più basito, nessuno è più sorpreso. Ognuno di voi ha capito tutto.\nCarlo is happy!\t\^o^/\n')
