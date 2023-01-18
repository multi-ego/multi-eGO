import pandas as pd
from vanessa import parametrize_ensemble_LJ, merge_and_clean_LJ
import parmed as pmd
from parmed import Structure
from parmed.gromacs import GromacsTopologyFile


class Multi_eGO_Ensemble:
    # TODO inserire la possibilità di leggere un RC fatto da più di una molecola

    meGO_topology = None
    reference_n_residues = None
    topology_dataframe = pd.DataFrame()
    
    # Contact matrices definition
    contact_matrices_dictionary = {}
    reference_contact_matrix = pd.DataFrame()
    meGO_LJ_potential = pd.DataFrame()
    meGO_LJ_14 = pd.DataFrame()

    def __init__(self, args):
        self.parameters = args

    def add_ensemble_from(self, ensemble):
        '''
        This method...
        This function inserts in a dictionary all the contact matrices learned from different simulations
        '''

        print('\t-', f'Adding topology and contacts from {ensemble.simulation}')

        if ensemble.simulation == 'reference':
            # This defines as the reference structure and eventual molecules will be added
            self.meGO_topology = ensemble.topology.copy(Structure)
            self.reference_n_residues = len(self.meGO_topology.residues)

        
        # Here we add only the non prortein without H atoms. To be properly tested with ligands
        # I  do not know how to select "protein" with ParmED so I am using the residues
        # TODO what if I have multiple molecules in different reference files?

        ensemble_to_add = ensemble.topology.copy(Structure)
        not_protein_not_H = ensemble_to_add[f'!((:1-{self.reference_n_residues})|(@%H=))']        
        self.meGO_topology = self.meGO_topology + not_protein_not_H

        # This part is for contact matrix
        if ensemble.simulation == 'reference':
            print('\t- Adding random coil probabilities from reference')
            self.reference_contact_matrix = ensemble.atomic_contacts
        else:
            print('\t-', f'Adding {ensemble.simulation} contacts in multi-eGO ensemble')
            self.contact_matrices_dictionary[f'{ensemble.simulation}'] = ensemble.atomic_contacts
        #print(self.contact_matrices_dictionary)
    

    def create_topology_dictionaries(self):
        
        # TODO from topology create another dataframe, probably the Ensemble.py is not needed anymore?
        self.meGO_topology 
        type_c12_dict = topology_dataframe[['sb_type', 'c12']].rename(columns={'sb_type':'; type'}).set_index('; type')['c12'].to_dict()



    def generate_LJ_potential(self):
        '''
        This function merges all the LJ contacts from ensembles into a single one.
        All contacts are reweighted based on the RC probability.
        Duplicates are removed.
        '''

        # TODO qui si mette parametrize_ensemble_LJ che equivale a MD_LJ_pairs
        # All contacts are reweighted by the random coil probability both as intra and intermolecular and added to the LJ pairs of multi-eGO.
        for simulation, contact_matrix in self.contact_matrices_dictionary.items():
            print('\t-', f'Parametrizing {simulation} probabilities')
            self.meGO_LJ_potential = pd.concat([self.meGO_LJ_potential, parametrize_ensemble_LJ(contact_matrix, self.reference_contact_matrix, self.parameters, simulation)], axis=0, sort=False, ignore_index=True)

        # TODO if LJ empty then RC, to be inserted here

        self.meGO_LJ_potential, self.meGO_LJ_14 = merge_and_clean_LJ(self.meGO_topology, self.meGO_LJ_potential, self.type_c12_dict, self.type_q_dict, self.parameters)

        print(self.meGO_LJ_potential)




    