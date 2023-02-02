from parmed import load_file
from parmed.topologyobjects import ParameterWarning
from vanessa import read_molecular_contacts, initialize_ensemble_topology, initialize_molecular_contacts
import os
import pandas as pd
import warnings
import glob
import sys
import copy

warnings.simplefilter("ignore")
warnings.filterwarnings("ignore", category=ParameterWarning)

def read_simulations(args, simulation):
    '''
    This function creates an Ensemble for each simulation defined by the folder name in --md_ensembles.
    '''
    simulation_path = os.path.join(args.protein, simulation)
    ensemble_tuple = (args, simulation, simulation_path)
    ensemble = Ensemble(ensemble_tuple)
    return copy.deepcopy(ensemble)


class Ensemble:
    args = None
    simulation = None
    simulation_path = None
    topology = None
    ensemble_topology_dataframe = None
    ensemble_molecules_idx_sbtype_dictionary = None
    ensemble_contact_matrices = {}
    sbtype_c12_dict = None
    atomic_contacts = None

    def __init__(self, ensemble_tuple):
        self.args = ensemble_tuple[0]
        self.simulation = ensemble_tuple[1]
        self.simulation_path = ensemble_tuple[2]
        print(f"- Creating the {ensemble_tuple[1]} ensemble")


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
        for matrix in glob.glob(f'inputs/{self.simulation_path}/*.ndx'):
            name = matrix.replace(f'inputs/{self.simulation_path}/', '')
            print(name)
            self.ensemble_contact_matrices[name] = read_molecular_contacts(matrix)
        print(self.ensemble_contact_matrices)

    def initialize_ensemble(self):
        '''
        After reading the input files, this function builds the ensemble topology and contact matrices
        '''
        print('\t-', f'Initializing {self.simulation} ensemble topology')
        self.ensemble_topology_dataframe, self.ensemble_molecules_idx_sbtype_dictionary, self.sbtype_c12_dict = initialize_ensemble_topology(self.topology, self.simulation)
        # TODO qua ci si ritorna dopo con una matrice scritta giusta
        self.atomic_contacts = initialize_molecular_contacts(self.ensemble_contact_matrices, self.ensemble_molecules_idx_sbtype_dictionary, self.simulation)
        self.ensemble_contact_matrices = {}