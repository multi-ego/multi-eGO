from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.vectors import Vector
from Bio.PDB.vectors import calc_dihedral
import pandas as pd


##############
# QUESTA PARTE DA SPOSTARE IN READ INPUT
##############



# PDB file reading for LJ parametrization

parser = PDBParser(PERMISSIVE=1)
structure_id = 'pep'
filename = 'GRETA/native/pep.pdb'
structure = parser.get_structure(structure_id, filename)

gro_dihedrals = pd.read_csv('GRETA/native/gro_dihedrals', sep = "\\s+", header = None)
gro_dihedrals.columns = ["ai", "aj", "ak", "al", "func", "def"]

print(gro_dihedrals.to_string())
# Reading PEP dihedrals



##################
# LJ
##################

LJ_pep = pd.DataFrame(columns = [';ai', 'aj', 'sigma', 'epsilon'])

pdb_ai = []
pdb_aj = []
pdb_sigma = []
pdb_epsilon = []


#prova = list(structure.get_atoms())
#print(prova)
#a0 = I.get(0)
#print(a0)

for atom1 in structure.get_atoms():
    for atom2 in structure.get_atoms():
        if atom2.get_serial_number() < atom1.get_serial_number(): 
            continue
        dist = atom2-atom1
        if dist < 6 and dist > 0:
            
            # Writing the ai and aj columns as types for pairs
            atomtype1 = atom1.name+'_'+str(atom1.get_parent().id[1])
            atomtype2 = atom2.name+'_'+str(atom2.get_parent().id[1])

            pdb_ai.append(atomtype1)
            pdb_aj.append(atomtype2)
          
            sigma = (dist/10) / (2**(1/6))
            epsilon = 1

            pdb_sigma.append(sigma)
            pdb_epsilon.append(epsilon)

            #values = [atomtype1, atomtype2, sigma, epsilon]
            #print(values)

LJ_pep[';ai'] = pdb_ai
LJ_pep['aj'] = pdb_aj
LJ_pep['sigma'] = pdb_sigma
LJ_pep['epsilon'] = pdb_epsilon

#print(LJ_pep) #sembra ragionevole

###########################
# DIHEDRALS
###########################

#import pandas as pd
#import numpy as np

#df = pd.DataFrame({'c1': [10, 11, 12], 'c2': [100, 110, 120]})

#for index, row in df.iterrows():
#    print(row['c1'], row['c2'])

atoms = structure.get_atoms()
#print(atoms)


#v1 = atom1.get_vector()

v1 = []

# E qui mi sono creato v1

for ato in structure.get_atoms():
    v = ato.get_vector()
    v1.append(v)

for index, row in gro_dihedrals.iterrows():

    dihedral_coso = calc_dihedral(v1[row['ai'] - 1], v1[row['aj'] - 1], v1[row['ak'] - 1], v1[row['al'] - 1])


#print(v1)
#print(len(v1))

#print(v1[84])


#print(prova_dihe)

#pdb_dihedrals = calc_dihedral(atom1.get_vector(), atom2.get_vector(), atom3.get_vector(), atom4.get_vector())

#avrà sei colonne, quartetti func angolo e kd = 1


# mi calcolo i vettori per ogni atomo e poi faccio un append su qualcosa. Poi semplicemente accedo ai valori che corrispondono a cose

# tipo l'atomo 1 sarà il primo elemento dunque 0. 


