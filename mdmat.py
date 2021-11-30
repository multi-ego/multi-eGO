import pandas as pd
import MDAnalysis
from protein_configuration import distance_residue


#TODO integration in GRETA

columns=['residue_ai', 'ai', 'residue_aj', 'aj', 'distance', 'probability']
mdmat_plainMD = pd.read_csv('inputs/native_ABeta/plainMD_contacts.ndx', header=None, sep = '\s+')
mdmat_plainMD.columns = columns
mdmat_plainMD['ai'] = mdmat_plainMD['ai']+1
mdmat_plainMD['aj'] = mdmat_plainMD['aj']+1
mdmat_plainMD['distance'] = mdmat_plainMD['distance']*10

plainMD_directory = '/home/emanuele/ABeta/markov'
#plainMD_directory = 'inputs/native_%s/native.pdb' %(protein)
reference_plainMD_structure = f'{plainMD_directory}/reduced-noh.gro'
reference_plainMD_trajectory = f'{plainMD_directory}/reduced-noh.xtc'

plainMD = MDAnalysis.Universe(reference_plainMD_structure, reference_plainMD_trajectory)
peptides = plainMD.select_atoms('all')
atomtypes_dict = {}
for atom in peptides:
    atomtypes_dict[atom.id] = str(atom.name) + '_' + str(atom.resnum)

mdmat_plainMD = mdmat_plainMD.replace({'ai':atomtypes_dict})
mdmat_plainMD = mdmat_plainMD.replace({'aj':atomtypes_dict})
mdmat_plainMD[['type_ai', 'residue_ai']] = mdmat_plainMD.ai.str.split("_", expand = True)
mdmat_plainMD[['type_aj', 'residue_aj']] = mdmat_plainMD.aj.str.split("_", expand = True)
mdmat_plainMD['residue_ai'] = mdmat_plainMD['residue_ai'].astype(int)
mdmat_plainMD['residue_aj'] = mdmat_plainMD['residue_aj'].astype(int)
mdmat_plainMD.drop(mdmat_plainMD[abs(mdmat_plainMD['residue_aj'] - mdmat_plainMD['residue_ai']) < distance_residue].index, inplace=True)
mdmat_plainMD.drop(columns=['residue_ai', 'residue_aj', 'type_ai', 'type_aj'], inplace=True)
mdmat_random_coil = pd.read_csv('inputs/native_ABeta/random_coil_contacts.ndx', header=None, sep = '\s+')
mdmat_random_coil.columns = columns
mdmat_random_coil['ai'] = mdmat_random_coil['ai']+1
mdmat_random_coil['aj'] = mdmat_random_coil['aj']+1
mdmat_random_coil['distance'] = mdmat_random_coil['distance']*10

plainMD_directory = '/home/emanuele/ABeta/markov'
reference_plainMD_structure = f'{plainMD_directory}/reduced-noh.gro'
reference_plainMD_trajectory = f'{plainMD_directory}/reduced-noh.xtc'

plainMD = MDAnalysis.Universe(reference_plainMD_structure, reference_plainMD_trajectory)
peptides = plainMD.select_atoms('all')

atomtypes_dict = {}
for atom in peptides:
    atomtypes_dict[atom.id] = str(atom.name) + '_' + str(atom.resnum)

mdmat_random_coil = mdmat_random_coil.replace({'ai':atomtypes_dict})
mdmat_random_coil = mdmat_random_coil.replace({'aj':atomtypes_dict})
mdmat_random_coil[['type_ai', 'residue_ai']] = mdmat_random_coil.ai.str.split("_", expand = True)
mdmat_random_coil[['type_aj', 'residue_aj']] = mdmat_random_coil.aj.str.split("_", expand = True)
mdmat_random_coil['residue_ai'] = mdmat_random_coil['residue_ai'].astype(int)
mdmat_random_coil['residue_aj'] = mdmat_random_coil['residue_aj'].astype(int)
mdmat_random_coil.drop(mdmat_random_coil[abs(mdmat_random_coil['residue_aj'] - mdmat_random_coil['residue_ai']) < distance_residue].index, inplace=True)
mdmat_random_coil.drop(columns=['residue_ai', 'residue_aj', 'type_ai', 'type_aj'], inplace=True)
