from Bio.PDB.PDBParser import PDBParser
import pandas as pd

# Creation of PDBParser object

parser = PDBParser(PERMISSIVE=1)

structure_id = 'pep'
filename = 'pep.pdb'
structure = parser.get_structure(structure_id, filename)

LJ_pep = pd.DataFrame(columns = ['atomtype1', 'atomtype2', 'sigma', 'epsilon'])
columns = list(LJ_pep)

data = []




for atom1 in structure.get_atoms():
    for atom2 in structure.get_atoms():
        if atom2.get_serial_number() < atom1.get_serial_number(): 
            continue
        dist = atom2-atom1
        if dist < 6 and dist > 0:
            #residue = atom2.get_parent()
            #print(residue.id[1])
            #print(f'atom 2 {coso} - atom 1 {atom1} = {dist}')
            atomtype1 = atom1.name+'_'+str(atom1.get_parent().id[1])
            atomtype2 = atom2.name+'_'+str(atom2.get_parent().id[1])

            sigma = (dist/10) / (2**(1/6))
            epsilon = 1

            values = [atomtype1, atomtype2, sigma, epsilon]
            print(values)
            
            zipped = zip(columns, values)
            
            a_dictionary = dict(zipped)

            #print(a_dictionary)

            LJ_pep.append(a_dictionary, ignore_index= True)  

            #print(LJ_pep)