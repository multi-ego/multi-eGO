import pandas as pd
import MDAnalysis as mda
import warnings
from topology_parser import topology_atoms, topology_bonds

warnings.filterwarnings('ignore', category=UserWarning, module='MDAnalysis')
warnings.filterwarnings("ignore", category=DeprecationWarning) 

pd.options.mode.chained_assignment = None  # default='warn'

def read_pdbs(parameters, flag):
    if not flag:
        directory = f"{parameters['input_folder']}/native.pdb"

    else:
        directory = f"{parameters['input_folder']}/fibril.pdb"
        
    pdb = mda.Universe(directory, guess_bonds = True)

    return pdb

def read_topology(parameters):
    # Read the topology created from pbd2gmx with gromos-primefull
    topol_atoms = topology_atoms(f'{parameters["input_folder"]}/topol.top')
    topol_bonds = topology_bonds(f'{parameters["input_folder"]}/topol.top')
	
    return topol_atoms, topol_bonds

def plainMD_mdmat(parameters):
    # Reading PlainMD contacts
    atomic_mat_plainMD = pd.read_csv(f'{parameters["input_folder"]}/plainMD_contacts.ndx', header=None, sep = '\s+')
    atomic_mat_plainMD.columns = ['residue_ai', 'ai', 'residue_aj', 'aj', 'distance', 'distance_NMR', 'probability']
    atomic_mat_plainMD.drop(columns=['distance'], inplace=True)
    atomic_mat_plainMD.columns = ['residue_ai', 'ai', 'residue_aj', 'aj', 'distance', 'probability']
    reference_plainMD_structure = f'{parameters["input_folder"]}/plainMD-noh.gro'

    plainMD = mda.Universe(reference_plainMD_structure)
    peptides = plainMD.select_atoms('all')
    plain_atomtypes_dict = {}
    for atom in peptides:
        plain_atomtypes_dict[atom.id] = str(atom.name) + '_' + str(atom.resnum)

    atomic_mat_plainMD = atomic_mat_plainMD.replace({'ai':plain_atomtypes_dict})
    atomic_mat_plainMD = atomic_mat_plainMD.replace({'aj':plain_atomtypes_dict})
    atomic_mat_plainMD[['type_ai', 'residue_ai']] = atomic_mat_plainMD.ai.str.split("_", expand = True)
    atomic_mat_plainMD[['type_aj', 'residue_aj']] = atomic_mat_plainMD.aj.str.split("_", expand = True)
    atomic_mat_plainMD['residue_ai'] = atomic_mat_plainMD['residue_ai'].astype(int)
    atomic_mat_plainMD['residue_aj'] = atomic_mat_plainMD['residue_aj'].astype(int)
    atomic_mat_plainMD.drop(atomic_mat_plainMD[abs(atomic_mat_plainMD['residue_aj'] - atomic_mat_plainMD['residue_ai']) < parameters['distance_residue']].index, inplace=True)
    atomic_mat_plainMD.drop(columns=['type_ai', 'type_aj'], inplace=True)

    return atomic_mat_plainMD


def random_coil_mdmat(parameters):
	# Reading Random Coil contacts
	atomic_mat_random_coil = pd.read_csv(f'{parameters["input_folder"]}/random_coil_contacts.ndx', header=None, sep = '\s+')
	atomic_mat_random_coil.columns = ['residue_ai', 'ai', 'residue_aj', 'aj', 'distance', 'distance_NMR', 'probability']
	atomic_mat_random_coil.drop(columns=['distance_NMR'], inplace=True)

	reference_random_coil_structure = f'{parameters["input_folder"]}/native.pdb'

	random_coil = mda.Universe(reference_random_coil_structure)
	peptides = random_coil.select_atoms('all')
	random_coil_atomtypes_dict = {}
	for atom in peptides:
		random_coil_atomtypes_dict[atom.id] = str(atom.name) + '_' + str(atom.resnum)

	atomic_mat_random_coil = atomic_mat_random_coil.replace({'ai':random_coil_atomtypes_dict})
	atomic_mat_random_coil = atomic_mat_random_coil.replace({'aj':random_coil_atomtypes_dict})
	atomic_mat_random_coil[['type_ai', 'residue_ai']] = atomic_mat_random_coil.ai.str.split("_", expand = True)
	atomic_mat_random_coil[['type_aj', 'residue_aj']] = atomic_mat_random_coil.aj.str.split("_", expand = True)
	atomic_mat_random_coil['residue_ai'] = atomic_mat_random_coil['residue_ai'].astype(int)
	atomic_mat_random_coil['residue_aj'] = atomic_mat_random_coil['residue_aj'].astype(int)
	atomic_mat_random_coil.drop(atomic_mat_random_coil[abs(atomic_mat_random_coil['residue_aj'] - atomic_mat_random_coil['residue_ai']) < parameters["distance_residue"]].index, inplace=True)
	atomic_mat_random_coil.drop(columns=['type_ai', 'type_aj'], inplace=True)
	pd.options.mode.chained_assignment = None
	atomic_mat_random_coil['probability'].loc[atomic_mat_random_coil['probability'] < (parameters['ratio_threshold']/10)] = parameters['ratio_threshold']/10
	pd.options.mode.chained_assignment = 'warn' 

	new_colnames = []
	for colname in atomic_mat_random_coil.columns:
		new_colnames.append(f'rc_{colname}')
	atomic_mat_random_coil.columns = new_colnames

	return atomic_mat_random_coil
