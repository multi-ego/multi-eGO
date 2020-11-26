import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PDB_small
from MDAnalysis.analysis import dihedrals, distances, contacts
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis.dihedrals import Dihedral
from read_input import read_gro_dihedrals
import pandas as pd
import itertools

pdb_pep = mda.Universe('GRETA/native/pep.pdb') # Da spostare su read pdb
native_dihedrals = read_gro_dihedrals()


########################## DIHEDRALS

phi_dihedrals = []
for index, row in native_dihedrals.iterrows():
    
    atom_selection = pdb_pep.atoms[row['ai'] - 1] + pdb_pep.atoms[row['aj'] - 1] + pdb_pep.atoms[row['ak'] - 1] + pdb_pep.atoms[row['al'] - 1]
    phi = atom_selection.dihedral.value()
    phi_dihedrals.append(phi)

native_dihedrals['func'] = 9
native_dihedrals['phi'] = phi_dihedrals
native_dihedrals['kd'] = ''
native_dihedrals['mult'] = ''

#print(native_dihedrals)



################################ PAIRS

u = mda.Universe(PDB_small)

ca = u.select_atoms('name CA')
n_ca = len(ca)
print(n_ca)

self_distances = distances.self_distance_array(ca.positions)
print(self_distances)



### HA DEL POTENZIALE

atom_sel = pdb_pep.select_atoms('all')
print(atom_sel)
print(len(atom_sel))

self_distances = distances.self_distance_array(atom_sel.positions)
print(self_distances)
print(len(self_distances))

atomtype = []
for atom in atom_sel:
    atp = (str(atom.name) + '_' + str(atom.resnum))
    atomtype.append(atp)
    #print(atp)

#print(atomtype)

contatti = list(itertools.combinations(atomtype, 2))
ai = []
aj = []
for n in range(0, len(contatti)):
    i = contatti[n][0]
    ai.append(i)
    j = contatti[n][1]
    aj.append(j)

#print(ai)
#print(len(ai))
#print(aj)
#print(len(aj))

native_LJ = pd.DataFrame(columns = [';ai', 'aj', 'distance', 'sigma', 'epsilon'])

native_LJ[';ai'] = ai
native_LJ['aj'] = aj
native_LJ['distance'] = self_distances

native_LJ = native_LJ[native_LJ.distance < 6]

native_LJ['sigma'] = (native_LJ['distance']/10) / (2**(1/6))
native_LJ['epsilon'] = 1

print(native_LJ)
#print(native_LJ.to_string())
