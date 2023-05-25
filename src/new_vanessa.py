
class VanessaClass:
    def __init__(self):
        pass

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
                'impropers' : get_impropers(topol[0].impropers),
                'pairs' : get_pairs(topol[0].adjusts)
            }
            # The following bonds are used in the parametrization of LJ 1-4
            self.bond_pairs[molecule] = get_bond_pairs(topol[0].bonds)
            self.user_pairs[molecule] = get_pairs(topol[0].adjusts)

    def generate_LJ_potential(self):
        '''
        This function merges all the LJ contacts from ensembles into a single one.
        All contacts are reweighted based on the RC probability.
        Duplicates are removed.
        '''
        self.meGO_LJ_potential, self.meGO_LJ_14 = parametrize_LJ(self.reference_topology_dataframe, self.molecule_type_dict, self.bond_pairs, self.user_pairs, self.sbtype_c12_dict, self.meGO_atomic_contacts, self.reference_atomic_contacts, self.check_atomic_contacts, self.sbtype_number_dict, self.parameters)
        self.meGO_LJ_14 = make_pairs_exclusion_topology(self.reference_topology_dataframe, self.bond_pairs, self.sbtype_c12_dict, self.parameters, self.meGO_LJ_14)
