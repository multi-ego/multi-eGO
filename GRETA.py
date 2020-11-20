from Bio.PDB.PDBParser import PDBParser
import pandas as pd

# Creation of PDBParser object

parser = PDBParser(PERMISSIVE=1)

structure_id = 'pep'
filename = 'pep.pdb'
structure = parser.get_structure(structure_id, filename)

LJ_pep = pd.DataFrame(columns = ['atomtype1', 'atomtype2', 'c6', 'c12'])
columns = list(LJ_pep)

data = []




for atom1 in structure.get_atoms():
    for atom2 in structure.get_atoms():
        dist = atom2-atom1
        if dist < 6:
            #residue = atom2.get_parent()
            #print(residue.id[1])
            #print(f'atom 2 {coso} - atom 1 {atom1} = {dist}')
            values = [atom1.name+'_'+str(atom1.get_parent().id[1]), atom2.name+'_'+str(atom2.get_parent().id[1])]
            print(values)
            # ADD C6 AND C12
            

for i in range(4, 10, 3):

  values = [i, i+1, i+2]

  zipped = zip(columns, values)

  a_dictionary = dict(zipped)

  print(a_dictionary)

  data.append(a_dictionary)