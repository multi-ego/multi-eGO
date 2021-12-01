import pandas as pd
import MDAnalysis
from protein_configuration import distance_residue


#TODO integration in GRETA

# Reading PlainMD contacts
columns=['residue_ai', 'ai', 'residue_aj', 'aj', 'distance', 'probability']
mdmat_plainMD = pd.read_csv('inputs/native_ABeta/plainMD_contacts.ndx', header=None, sep = '\s+')
mdmat_plainMD.columns = columns
#mdmat_plainMD['ai'] = mdmat_plainMD['ai']+1
#mdmat_plainMD['aj'] = mdmat_plainMD['aj']+1
mdmat_plainMD['distance'] = mdmat_plainMD['distance']*10

plainMD_directory = '/home/emanuele/ABeta/markov'
#plainMD_directory = 'inputs/native_%s/native.pdb' %(protein)
reference_plainMD_structure = f'{plainMD_directory}/reduced-noh.gro'

plainMD = MDAnalysis.Universe(reference_plainMD_structure)
peptides = plainMD.select_atoms('all')
plain_atomtypes_dict = {}
for atom in peptides:
    plain_atomtypes_dict[atom.id] = str(atom.name) + '_' + str(atom.resnum)

mdmat_plainMD = mdmat_plainMD.replace({'ai':plain_atomtypes_dict})
mdmat_plainMD = mdmat_plainMD.replace({'aj':plain_atomtypes_dict})
mdmat_plainMD[['type_ai', 'residue_ai']] = mdmat_plainMD.ai.str.split("_", expand = True)
mdmat_plainMD[['type_aj', 'residue_aj']] = mdmat_plainMD.aj.str.split("_", expand = True)
mdmat_plainMD['residue_ai'] = mdmat_plainMD['residue_ai'].astype(int)
mdmat_plainMD['residue_aj'] = mdmat_plainMD['residue_aj'].astype(int)
mdmat_plainMD.drop(mdmat_plainMD[abs(mdmat_plainMD['residue_aj'] - mdmat_plainMD['residue_ai']) < distance_residue].index, inplace=True)
mdmat_plainMD.drop(columns=['residue_ai', 'residue_aj', 'type_ai', 'type_aj'], inplace=True)


# Reading Random Coil contacts
mdmat_random_coil = pd.read_csv('inputs/native_ABeta/random_coil_contacts.ndx', header=None, sep = '\s+')
mdmat_random_coil.columns = columns
#mdmat_random_coil['ai'] = mdmat_random_coil['ai']+1
#mdmat_random_coil['aj'] = mdmat_random_coil['aj']+1
mdmat_random_coil['distance'] = mdmat_random_coil['distance']*10

# kkkdkd
probability_min_rc = mdmat_random_coil[mdmat_random_coil.probability != 0].min()['probability']
mdmat_random_coil['probability'].loc[mdmat_random_coil['probability'] == 0] = probability_min_rc/2


random_coil_directory = '/home/emanuele/ABeta/random_coil/monomer_test/native_278K'
#plainMD_directory = 'inputs/native_%s/native.pdb' %(protein)

reference_random_coil_structure = f'{random_coil_directory}/box_abeta_greta.gro'
random_coil = MDAnalysis.Universe(reference_random_coil_structure)
peptides = random_coil.select_atoms('all')
random_coil_atomtypes_dict = {}
for atom in peptides:
    random_coil_atomtypes_dict[atom.id] = str(atom.name) + '_' + str(atom.resnum)


mdmat_random_coil = mdmat_random_coil.replace({'ai':random_coil_atomtypes_dict})
mdmat_random_coil = mdmat_random_coil.replace({'aj':random_coil_atomtypes_dict})
mdmat_random_coil[['type_ai', 'residue_ai']] = mdmat_random_coil.ai.str.split("_", expand = True)
mdmat_random_coil[['type_aj', 'residue_aj']] = mdmat_random_coil.aj.str.split("_", expand = True)
mdmat_random_coil['residue_ai'] = mdmat_random_coil['residue_ai'].astype(int)
mdmat_random_coil['residue_aj'] = mdmat_random_coil['residue_aj'].astype(int)
mdmat_random_coil.drop(mdmat_random_coil[abs(mdmat_random_coil['residue_aj'] - mdmat_random_coil['residue_ai']) < distance_residue].index, inplace=True)
mdmat_random_coil.drop(columns=['residue_ai', 'residue_aj', 'type_ai', 'type_aj'], inplace=True)

