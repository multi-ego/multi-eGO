from Bio.PDB.PDBParser import PDBParser

# Creation of PDBParser object

parser = PDBParser(PERMISSIVE=1)

structure_id = '1fat'
filename = '1fat.pdb'
structure = parser.get_structure(structure_id, filename)
resolution = structure.header["keywords"]


for atom in structure.get_atoms():
    print(atom)