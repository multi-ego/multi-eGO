import pandas as pd
import MDAnalysis as mda
import warnings
from topology_parser import topology_parser
import parmed as pmd

warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')
warnings.filterwarnings("ignore", category=DeprecationWarning) 

pd.options.mode.chained_assignment = None  # default='warn'

def read_pdbs(parameters, flag):
    if not flag:
        directory = f"{parameters['input_folder']}/native.pdb"

    else:
        directory = f"{parameters['input_folder']}/fibril.pdb"

    print('\tReading ', directory)        
    pdb = mda.Universe(directory, guess_bonds = True)

    return pdb


def plainMD_mdmat(parameters, ensemble_parameters, idx_sbtype_dict):
    # Reading PlainMD contacts
#    if ensemble_parameters['is_ligand']:
#        contact_map_file = f'{parameters["input_folder"]}/ligandMD_contacts.ndx'
#    else:
#        contact_map_file = f'{parameters["input_folder"]}/plainMD_contacts.ndx'
    contact_map_file = ensemble_parameters['mdmat_contacts_file']
    print('\tReading ', contact_map_file)        
    atomic_mat_plainMD = pd.read_csv(contact_map_file, header=None, sep = '\s+')
    atomic_mat_plainMD.columns = ['residue_ai', 'ai', 'residue_aj', 'aj', 'distance', 'distance_NMR', 'probability']
    atomic_mat_plainMD.drop(columns=['distance'], inplace=True)
    atomic_mat_plainMD.columns = ['residue_ai', 'ai', 'residue_aj', 'aj', 'distance', 'probability']
    atomic_mat_plainMD = atomic_mat_plainMD.replace({'ai':idx_sbtype_dict})
    atomic_mat_plainMD = atomic_mat_plainMD.replace({'aj':idx_sbtype_dict})
    atomic_mat_plainMD[['type_ai', 'residue_ai']] = atomic_mat_plainMD.ai.str.split("_", expand = True)
    atomic_mat_plainMD[['type_aj', 'residue_aj']] = atomic_mat_plainMD.aj.str.split("_", expand = True)
    atomic_mat_plainMD['residue_ai'] = atomic_mat_plainMD['residue_ai'].astype(int)
    atomic_mat_plainMD['residue_aj'] = atomic_mat_plainMD['residue_aj'].astype(int)
    atomic_mat_plainMD.drop(atomic_mat_plainMD[abs(atomic_mat_plainMD['residue_aj'] - atomic_mat_plainMD['residue_ai']) < parameters['distance_residue']].index, inplace=True)
    atomic_mat_plainMD.drop(columns=['type_ai', 'type_aj'], inplace=True)

    return atomic_mat_plainMD


def random_coil_mdmat(parameters, idx_sbtype_dict):
    # Reading Random Coil contacts
    contact_map_file = f'{parameters["input_folder"]}/random_coil_contacts.ndx'
    print('\tReading ', contact_map_file)        
    atomic_mat_random_coil = pd.read_csv(contact_map_file, header=None, sep = '\s+')
    atomic_mat_random_coil.columns = ['residue_ai', 'ai', 'residue_aj', 'aj', 'distance', 'distance_NMR', 'probability']
    atomic_mat_random_coil.drop(columns=['distance_NMR'], inplace=True)

    atomic_mat_random_coil = atomic_mat_random_coil.replace({'ai':idx_sbtype_dict})
    atomic_mat_random_coil = atomic_mat_random_coil.replace({'aj':idx_sbtype_dict})
    atomic_mat_random_coil[['type_ai', 'residue_ai']] = atomic_mat_random_coil.ai.str.split("_", expand = True)
    atomic_mat_random_coil[['type_aj', 'residue_aj']] = atomic_mat_random_coil.aj.str.split("_", expand = True)
    atomic_mat_random_coil['residue_ai'] = atomic_mat_random_coil['residue_ai'].astype(int)
    atomic_mat_random_coil['residue_aj'] = atomic_mat_random_coil['residue_aj'].astype(int)
    atomic_mat_random_coil.drop(atomic_mat_random_coil[abs(atomic_mat_random_coil['residue_aj'] - atomic_mat_random_coil['residue_ai']) < parameters["distance_residue"]].index, inplace=True)
    atomic_mat_random_coil.drop(columns=['type_ai', 'type_aj'], inplace=True)
    atomic_mat_random_coil['probability'].loc[atomic_mat_random_coil['probability'] < (parameters['rc_threshold'])] = parameters['rc_threshold']

    new_colnames = []
    for colname in atomic_mat_random_coil.columns:
        new_colnames.append(f'rc_{colname}')
    atomic_mat_random_coil.columns = new_colnames

    return atomic_mat_random_coil
