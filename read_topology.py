import gromologist
from read_input import read_top
from gromologist import Top, Pdb




print(gromologist.__file__)
topology = read_top()

protein = topology.molecules[0]
bonds = protein.list_bonds(by_types=True)
print(bonds)

