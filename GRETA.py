from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.vectors import Vector
from Bio.PDB.vectors import calc_dihedral
import pandas as pd
import numpy as np
import math
from read_input import read_pdb, read_gro_dihedrals


##############
# QUESTA PARTE DA SPOSTARE IN READ INPUT
##############

# PDB file reading for LJ parametrization

#parser = PDBParser(PERMISSIVE=1)
#structure_id = 'pep'
#filename = 'GRETA/native/pep.pdb'
#structure = parser.get_structure(structure_id, filename)


#gro_dihedrals = pd.read_csv('GRETA/native/gro_dihedrals', sep = "\\s+", header = None)
#gro_dihedrals.columns = ["ai", "aj", "ak", "al", "func", "def"]



native_structure, fibril_structure = read_pdb()
native_dihedrals = read_gro_dihedrals()


# Reading PEP dihedrals

##################
# LJ
##################
native_LJ = pd.DataFrame(columns = [';ai', 'aj', 'sigma', 'epsilon'])

pdb_ai = []
pdb_aj = []
pdb_sigma = []
pdb_epsilon = []


#prova = list(structure.get_atoms())
#print(prova)
#a0 = I.get(0)
#print(a0)

for atom1 in native_structure.get_atoms():
    for atom2 in native_structure.get_atoms():
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

native_LJ[';ai'] = pdb_ai
native_LJ['aj'] = pdb_aj
native_LJ['sigma'] = pdb_sigma
native_LJ['epsilon'] = pdb_epsilon

print(native_LJ)

#print(pep_) #sembra ragionevole

###########################
# DIHEDRALS
###########################

atoms = native_structure.get_atoms()

atom_coord = []
for atoms in native_structure.get_atoms():
    v = atoms.get_vector()
    atom_coord.append(v)

rad_dihedrals = []
for index, row in native_dihedrals.iterrows():

    rad = calc_dihedral(atom_coord[row['ai'] - 1], atom_coord[row['aj'] - 1], atom_coord[row['ak'] - 1], atom_coord[row['al'] - 1])
    
    #rad = calc_dihedral(atom_coord[row['ai'] - 1], atom_coord[row['aj'] - 1], atom_coord[row['ak'] - 1], atom_coord[row['al'] - 1])
    rad_dihedrals.append(rad)


#### GRO DIHEDRALS FIBRILLA NON SERVE, BASTA RICALCOLARSI LE COSE DAL PDB USANDO LO STESSO ELENCO


# AGGIUNGI LE FUNZIONI PER BONDED, ANGLES, IMPROPERS per le exclusion lists, bonds tale e quale, angoli primo e terzo, 
# sui diedri e sugli improri tutte le combinazioni e togli i duplicati keep unique
# nel calcolo delle distanze vengono escluse queste coppie delle exclusion (&) fa parte della stessa catena


print(rad_dihedrals)
print('Radian values : \n', rad_dihedrals)

phi_dihedrals = np.rad2deg(rad_dihedrals)
print('\n Degree values : \n', phi_dihedrals)



native_dihedrals['func'] = 9
native_dihedrals['phi'] = phi_dihedrals
native_dihedrals['kd'] = ''
native_dihedrals['mult'] = ''




print(native_LJ)
print(native_dihedrals)




# mi calcolo i vettori per ogni atomo e poi faccio un append su qualcosa. Poi semplicemente accedo ai valori che corrispondono a cose

# tipo l'atomo 1 sarà il primo elemento dunque 0. 


