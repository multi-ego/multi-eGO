import pandas as pd
import MDAnalysis
from protein_configuration import distance_residue, ratio_treshold


#TODO integration in GRETA

# Reading PlainMD contacts
atomic_mat_plainMD = pd.read_csv('inputs/native_ABeta/plainMD_contacts.ndx', header=None, sep = '\s+')
atomic_mat_plainMD.columns = ['residue_ai', 'ai', 'residue_aj', 'aj', 'distance', 'distance_NMR', 'probability']
atomic_mat_plainMD.drop(columns=['distance'], inplace=True)
atomic_mat_plainMD.columns = ['residue_ai', 'ai', 'residue_aj', 'aj', 'distance', 'probability']
atomic_mat_plainMD['distance'] = atomic_mat_plainMD['distance']

plainMD_directory = '/home/emanuele/ABeta/markov'
#plainMD_directory = 'inputs/native_%s/native.pdb' %(protein)
reference_plainMD_structure = f'{plainMD_directory}/reduced-noh.gro'

plainMD = MDAnalysis.Universe(reference_plainMD_structure)
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
atomic_mat_plainMD.drop(atomic_mat_plainMD[abs(atomic_mat_plainMD['residue_aj'] - atomic_mat_plainMD['residue_ai']) < distance_residue].index, inplace=True)
#atomic_mat_plainMD.drop(columns=['residue_ai', 'residue_aj', 'type_ai', 'type_aj'], inplace=True)
atomic_mat_plainMD.drop(columns=['type_ai', 'type_aj'], inplace=True)

# Reading Random Coil contacts
atomic_mat_random_coil = pd.read_csv('inputs/native_ABeta/random_coil_contacts.ndx', header=None, sep = '\s+')
atomic_mat_random_coil.columns = ['residue_ai', 'ai', 'residue_aj', 'aj', 'distance', 'probability']
atomic_mat_random_coil['distance'] = atomic_mat_random_coil['distance']

random_coil_directory = '/home/emanuele/ABeta/random_coil/monomer_test/native_278K'
#plainMD_directory = 'inputs/native_%s/native.pdb' %(protein)

reference_random_coil_structure = f'{random_coil_directory}/box_abeta_greta.gro'
random_coil = MDAnalysis.Universe(reference_random_coil_structure)
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
atomic_mat_random_coil.drop(atomic_mat_random_coil[abs(atomic_mat_random_coil['residue_aj'] - atomic_mat_random_coil['residue_ai']) < distance_residue].index, inplace=True)
#atomic_mat_random_coil.drop(columns=['residue_ai', 'residue_aj', 'type_ai', 'type_aj'], inplace=True)
atomic_mat_random_coil.drop(columns=['type_ai', 'type_aj'], inplace=True)


# QUESTA CI PIACE
atomic_mat_random_coil['probability'].loc[atomic_mat_random_coil['probability'] < (ratio_treshold/10)] = ratio_treshold/10
