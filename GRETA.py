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

atoms = structure.get_atoms()

vectors = []

for atoms in structure.get_atoms():
    v = atoms.get_vector()
    vectors.append(v)

phi_dihedrals = []

for index, row in gro_dihedrals.iterrows():
    phi = calc_dihedral(vectors[row['ai'] - 1], vectors[row['aj'] - 1], vectors[row['ak'] - 1], vectors[row['al'] - 1])
    phi_dihedrals.append(phi)

gro_dihedrals['func'] = 9
gro_dihedrals['phi'] = phi_dihedrals
gro_dihedrals['kd'] = ''
gro_dihedrals['mult'] = ''



print(gro_dihedrals)




# mi calcolo i vettori per ogni atomo e poi faccio un append su qualcosa. Poi semplicemente accedo ai valori che corrispondono a cose

# tipo l'atomo 1 sarà il primo elemento dunque 0. 


