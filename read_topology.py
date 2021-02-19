import gromologist
import pandas as pd
from read_input import read_top

# Import the topology informations
topology = read_top()
protein = topology.molecules[0]


# Unfortunately those are strings which must be separated
# BONDS
atom_types, atom_resids = protein.list_bonds(by_resid=True)
topology_bonds = pd.DataFrame(columns=['ai', 'aj', 'ai_type', 'ai_resid','aj_type', 'aj_resid'])
ai_type = []
aj_type = []
ai_resid = []
aj_resid = []

for atyp in atom_types:
    atyp_split = atyp.split(' ')
    ai_type.append(atyp_split[0])
    aj_type.append(atyp_split[1])

for ares in atom_resids:
    ares_split = ares.split(' ')
    ai_resid.append(ares_split[0])
    aj_resid.append(ares_split[1])

topology_bonds['ai_type'] = ai_type
topology_bonds['aj_type'] = aj_type
topology_bonds['ai_resid'] = ai_resid
topology_bonds['aj_resid'] = aj_resid
topology_bonds['ai'] = topology_bonds['ai_type'] + '_' + topology_bonds['ai_resid'].astype(str)
topology_bonds['aj'] = topology_bonds['aj_type'] + '_' + topology_bonds['aj_resid'].astype(str)
topology_bonds.drop(['ai_type', 'ai_resid','aj_type', 'aj_resid'], axis=1, inplace=True)

# ANGLES
atom_types, atom_resids = protein.list_angles(by_resid=True)
topology_angles = pd.DataFrame(columns=['ai', 'aj', 'ak', 'ai_type', 'ai_resid', 'aj_type', 'aj_resid', 'ak_type', 'ak_resid'])
ai_type = []
aj_type = []
ak_type = []
ai_resid = []
aj_resid = []
ak_resid = []

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

topology_angles['ai_type'] = ai_type
topology_angles['aj_type'] = aj_type
topology_angles['ak_type'] = ak_type
topology_angles['ai_resid'] = ai_resid
topology_angles['aj_resid'] = aj_resid
topology_angles['ak_resid'] = ak_resid
topology_angles['ai'] = topology_angles['ai_type'] + '_' + topology_angles['ai_resid'].astype(str)
topology_angles['aj'] = topology_angles['aj_type'] + '_' + topology_angles['aj_resid'].astype(str)
topology_angles['ak'] = topology_angles['ak_type'] + '_' + topology_angles['ak_resid'].astype(str)
topology_angles.drop(['ai_type', 'ai_resid','aj_type', 'aj_resid', 'ak_type', 'ak_resid'], axis=1, inplace=True)

# DIHEDRALS
atom_types, atom_resids = protein.list_dihedrals(by_resid=True)
topology_dihedrals = pd.DataFrame(columns=['ai', 'aj', 'ak', 'al', 'ai_type', 'ai_resid', 'aj_type', 'aj_resid', 'ak_type', 'ak_resid', 'al_type', 'al_resid'])
ai_type = []
aj_type = []
ak_type = []
al_type = []
ai_resid = []
aj_resid = []
ak_resid = []
al_resid = []

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

topology_dihedrals['ai_type'] = ai_type
topology_dihedrals['aj_type'] = aj_type
topology_dihedrals['ak_type'] = ak_type
topology_dihedrals['al_type'] = al_type
topology_dihedrals['ai_resid'] = ai_resid
topology_dihedrals['aj_resid'] = aj_resid
topology_dihedrals['ak_resid'] = ak_resid
topology_dihedrals['al_resid'] = al_resid
topology_dihedrals['ai'] = topology_dihedrals['ai_type'] + '_' + topology_dihedrals['ai_resid'].astype(str)
topology_dihedrals['aj'] = topology_dihedrals['aj_type'] + '_' + topology_dihedrals['aj_resid'].astype(str)
topology_dihedrals['ak'] = topology_dihedrals['ak_type'] + '_' + topology_dihedrals['ak_resid'].astype(str)
topology_dihedrals['al'] = topology_dihedrals['al_type'] + '_' + topology_dihedrals['al_resid'].astype(str)
topology_dihedrals.drop(['ai_type', 'ai_resid','aj_type', 'aj_resid', 'ak_type', 'ak_resid', 'al_type', 'al_resid'], axis=1, inplace=True)

# IMPROPERS
atom_types, atom_resids = protein.list_impropers(by_resid=True)
topology_impropers = pd.DataFrame(columns=['ai', 'aj', 'ak', 'al', 'ai_type', 'ai_resid', 'aj_type', 'aj_resid', 'ak_type', 'ak_resid', 'al_type', 'al_resid'])
ai_type = []
aj_type = []
ak_type = []
al_type = []
ai_resid = []
aj_resid = []
ak_resid = []
al_resid = []

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

topology_impropers['ai_type'] = ai_type
topology_impropers['aj_type'] = aj_type
topology_impropers['ak_type'] = ak_type
topology_impropers['al_type'] = al_type
topology_impropers['ai_resid'] = ai_resid
topology_impropers['aj_resid'] = aj_resid
topology_impropers['ak_resid'] = ak_resid
topology_impropers['al_resid'] = al_resid
topology_impropers['ai'] = topology_impropers['ai_type'] + '_' + topology_impropers['ai_resid'].astype(str)
topology_impropers['aj'] = topology_impropers['aj_type'] + '_' + topology_impropers['aj_resid'].astype(str)
topology_impropers['ak'] = topology_impropers['ak_type'] + '_' + topology_impropers['ak_resid'].astype(str)
topology_impropers['al'] = topology_impropers['al_type'] + '_' + topology_impropers['al_resid'].astype(str)
topology_impropers.drop(['ai_type', 'ai_resid','aj_type', 'aj_resid', 'ak_type', 'ak_resid', 'al_type', 'al_resid'], axis=1, inplace=True)



print(topology_bonds) # 87 rows
print(topology_angles) # 120 rows
print(topology_dihedrals) # 69 rows
print(topology_impropers) # 43 rows