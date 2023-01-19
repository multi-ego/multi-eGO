from parmed import load_file
from parmed.topologyobjects import ParameterWarning
from vanessa import read_molecular_contacts, initialize_ensemble_topology, initialize_molecular_contacts, create_ensemble_dictionaries
import os
import pandas as pd
import warnings
import glob
import sys

warnings.filterwarnings("ignore", category=ParameterWarning)


def read_simulations(args, simulation):
    '''
    This function creates an Ensemble for each simulation defined by the folder name in --md_ensembles.
    '''
    
    simulation_path = os.path.join(args.protein, simulation)
    ensemble = Ensemble(simulation, simulation_path, args)
    ensemble.read_files()
    ensemble.initialize_ensemble()
    return ensemble


class Ensemble:
    def __init__(self, simulation, simulation_path, args):
        self.args = args
        self.simulation = simulation
        self.simulation_path = simulation_path
        self.topology = None
        self.ensemble_molecules = None
        self.topology_dataframe = pd.DataFrame()
        self.intramolecular_contacts = pd.DataFrame()
        self.intermolecular_contacts = pd.DataFrame()
        self.atomic_contacts = pd.DataFrame()
        self.structure = None

        print(f"- Creating the {simulation} ensemble")


    def read_files(self):
        '''
        This function search the .top files and the contact matrices .ndx for each folder defined in the command line.
        It is important to call properly the names.
        '''

        top_path = glob.glob(f'inputs/{self.simulation_path}/*.top')
        print('\t-', f'Reading {top_path[0]}')
        self.topology = load_file(top_path[0])
        # Reading contact matrix created using gmx_clustsize
        # Reference requires both intra and inter molecular matrices

        #if self.simulation == 'reference':
        intramolecular_path = os.path.join('inputs/', self.simulation_path + '/intramat.ndx')
        self.intramolecular_contacts = read_molecular_contacts(intramolecular_path, is_same_chain=True)
        intermolecular_path = os.path.join('inputs/', self.simulation_path + '/intermat.ndx')
        self.intermolecular_contacts = read_molecular_contacts(intermolecular_path, is_same_chain=False)

        #if self.simulation == self.args.intra:
        #    intramolecular_path = os.path.join('inputs/', self.simulation_path + '/intramat.ndx')
        #    self.intramolecular_contacts = read_molecular_contacts(intramolecular_path, is_same_chain=True)

        #if self.simulation == self.args.inter:
        #    intermolecular_path = os.path.join('inputs/', self.simulation_path + '/intermat.ndx')
        #    self.intermolecular_contacts = read_molecular_contacts(intermolecular_path, is_same_chain=False)
    
    
    def initialize_ensemble(self):
        '''
        After reading the input files, this function builds the ensemble topology and contact matrices
        '''

        # TODO here I need to move all the ensemble part from greta.py by reading the topology and structure dataframes
        print('\t-', f'Initializing {self.simulation} ensemble topology')
        self.topology, self.topology_dataframe = initialize_ensemble_topology(self.topology, self.simulation)

        # TODO questi dizionari probabilmente non servono più dentro la classe ensemble, controlla poi con il topology parser
        self.number_of_atoms, self.sbtype_idx_dict, self.idx_sbtype_dict, self.type_c12_dict, self.type_q_dict = create_ensemble_dictionaries(self.topology_dataframe, self.simulation)
    
        print('\t-', f'{self.simulation} structure-based topology created')
        print('\t-', f'Initializing {self.simulation} ensemble molecular contacts')

        self.atomic_contacts = pd.concat([self.atomic_contacts, initialize_molecular_contacts(self.intramolecular_contacts, self.idx_sbtype_dict, self.simulation), initialize_molecular_contacts(self.intermolecular_contacts, self.idx_sbtype_dict, self.simulation)], axis=0)
