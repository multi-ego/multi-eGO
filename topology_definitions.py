import itertools
import pandas as pd
from read_input import read_top
import numpy as np

# Import the topology informations
topology = read_top()
protein = topology.molecules[0]
exclusion_list_gromologist = []


# Unfortunately those are strings which must be separated
# BONDS
atom_types, atom_resids = protein.list_bonds(by_resid=True)
ai_type, aj_type, ai_resid, aj_resid = [], [], [], []

for atyp in atom_types:
    atyp_split = atyp.split(' ')
    ai_type.append(atyp_split[0])
    aj_type.append(atyp_split[1])

for ares in atom_resids:
    ares_split = ares.split(' ')
    ai_resid.append(ares_split[0])
    aj_resid.append(ares_split[1])

topology_bonds = pd.DataFrame(np.column_stack([ai_type, ai_resid, aj_type, aj_resid]), columns=['ai_type', 'ai_resid','aj_type', 'aj_resid'])
topology_bonds['ai'] = topology_bonds['ai_type'] + '_' + topology_bonds['ai_resid'].astype(str)
topology_bonds['aj'] = topology_bonds['aj_type'] + '_' + topology_bonds['aj_resid'].astype(str)
topology_bonds.drop(['ai_type', 'ai_resid','aj_type', 'aj_resid'], axis=1, inplace=True)

# ANGLES
atom_types, atom_resids = protein.list_angles(by_resid=True)
ai_type, aj_type, ak_type, ai_resid, aj_resid, ak_resid = [], [], [], [], [], []

for atyp in atom_types:
    atyp_split = atyp.split(' ')
    ai_type.append(atyp_split[0])
    aj_type.append(atyp_split[1])
    ak_type.append(atyp_split[2])

for ares in atom_resids:
    ares_split = ares.split(' ')
    ai_resid.append(ares_split[0])
    aj_resid.append(ares_split[1])
    ak_resid.append(ares_split[2])

topology_angles = pd.DataFrame(np.column_stack([ai_type, ai_resid, aj_type, aj_resid, ak_type, ak_resid]), columns=['ai_type', 'ai_resid', 'aj_type', 'aj_resid', 'ak_type', 'ak_resid'])
topology_angles['ai'] = topology_angles['ai_type'] + '_' + topology_angles['ai_resid'].astype(str)
topology_angles['aj'] = topology_angles['aj_type'] + '_' + topology_angles['aj_resid'].astype(str)
topology_angles['ak'] = topology_angles['ak_type'] + '_' + topology_angles['ak_resid'].astype(str)
topology_angles.drop(['ai_type', 'ai_resid','aj_type', 'aj_resid', 'ak_type', 'ak_resid'], axis=1, inplace=True)

# DIHEDRALS
atom_types, atom_resids = protein.list_dihedrals(by_resid=True)
ai_type, aj_type, ak_type, al_type, ai_resid, aj_resid, ak_resid, al_resid = [], [], [], [], [], [], [], []

for atyp in atom_types:
    atyp_split = atyp.split(' ')
    ai_type.append(atyp_split[0])
    aj_type.append(atyp_split[1])
    ak_type.append(atyp_split[2])
    al_type.append(atyp_split[3])

for ares in atom_resids:
    ares_split = ares.split(' ')
    ai_resid.append(ares_split[0])
    aj_resid.append(ares_split[1])
    ak_resid.append(ares_split[2])
    al_resid.append(ares_split[3])

topology_dihedrals = pd.DataFrame(np.column_stack([ai_type, ai_resid, aj_type, aj_resid, ak_type, ak_resid, al_type, al_resid]), columns=['ai_type', 'ai_resid', 'aj_type', 'aj_resid', 'ak_type', 'ak_resid', 'al_type', 'al_resid'])
topology_dihedrals['ai'] = topology_dihedrals['ai_type'] + '_' + topology_dihedrals['ai_resid'].astype(str)
topology_dihedrals['aj'] = topology_dihedrals['aj_type'] + '_' + topology_dihedrals['aj_resid'].astype(str)
topology_dihedrals['ak'] = topology_dihedrals['ak_type'] + '_' + topology_dihedrals['ak_resid'].astype(str)
topology_dihedrals['al'] = topology_dihedrals['al_type'] + '_' + topology_dihedrals['al_resid'].astype(str)
topology_dihedrals.drop(['ai_type', 'ai_resid','aj_type', 'aj_resid', 'ak_type', 'ak_resid', 'al_type', 'al_resid'], axis=1, inplace=True)

# IMPROPERS
atom_types, atom_resids = protein.list_impropers(by_resid=True)
ai_type, aj_type, ak_type, al_type, ai_resid, aj_resid, ak_resid, al_resid = [], [], [], [], [], [], [], []

for atyp in atom_types:
    atyp_split = atyp.split(' ')
    ai_type.append(atyp_split[0])
    aj_type.append(atyp_split[1])
    ak_type.append(atyp_split[2])
    al_type.append(atyp_split[3])

for ares in atom_resids:
    ares_split = ares.split(' ')
    ai_resid.append(ares_split[0])
    aj_resid.append(ares_split[1])
    ak_resid.append(ares_split[2])
    al_resid.append(ares_split[3])

topology_impropers = pd.DataFrame(np.column_stack([ai_type, ai_resid, aj_type, aj_resid, ak_type, ak_resid, al_type, al_resid]), columns=['ai_type', 'ai_resid', 'aj_type', 'aj_resid', 'ak_type', 'ak_resid', 'al_type', 'al_resid'])
topology_impropers['ai'] = topology_impropers['ai_type'] + '_' + topology_impropers['ai_resid'].astype(str)
topology_impropers['aj'] = topology_impropers['aj_type'] + '_' + topology_impropers['aj_resid'].astype(str)
topology_impropers['ak'] = topology_impropers['ak_type'] + '_' + topology_impropers['ak_resid'].astype(str)
topology_impropers['al'] = topology_impropers['al_type'] + '_' + topology_impropers['al_resid'].astype(str)
topology_impropers.drop(['ai_type', 'ai_resid','aj_type', 'aj_resid', 'ak_type', 'ak_resid', 'al_type', 'al_resid'], axis=1, inplace=True)


# Exclusion list creation

# It takes the column names of each dataframe (ai, aj, ak, al) and make combinations pairwise.
# From the dataframes takes the ai aj columns and marges making the atom combination and assigning to a list

for c in list(itertools.combinations((topology_bonds.columns.values.tolist()), 2)):
    exclusion_list_gromologist.append((topology_bonds[c[0]] + '_' + topology_bonds[c[1]]).to_list())
    exclusion_list_gromologist.append((topology_bonds[c[1]] + '_' + topology_bonds[c[0]]).to_list())

for c in list(itertools.combinations((topology_angles.columns.values.tolist()), 2)):
    exclusion_list_gromologist.append((topology_angles[c[0]] + '_' + topology_angles[c[1]]).to_list())
    exclusion_list_gromologist.append((topology_angles[c[1]] + '_' + topology_angles[c[0]]).to_list())

for c in list(itertools.combinations((topology_dihedrals.columns.values.tolist()), 2)):
    exclusion_list_gromologist.append((topology_dihedrals[c[0]] + '_' + topology_dihedrals[c[1]]).to_list())
    exclusion_list_gromologist.append((topology_dihedrals[c[1]] + '_' + topology_dihedrals[c[0]]).to_list())

for c in list(itertools.combinations((topology_impropers.columns.values.tolist()), 2)):
    exclusion_list_gromologist.append((topology_impropers[c[0]] + '_' + topology_impropers[c[1]]).to_list())
    exclusion_list_gromologist.append((topology_impropers[c[1]] + '_' + topology_impropers[c[0]]).to_list())

# The resulting list is a list of lists and this makes a unique flat one
exclusion_list_gromologist = [item for sublist in exclusion_list_gromologist for item in sublist]
#print(len(flat_exclusion_list_gromologist)) # 2238

exclusion_list_gromologist = list(set(exclusion_list_gromologist))
#print(len(exclusion_list_gromologist)) # 514