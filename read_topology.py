import gromologist
from read_input import read_top


topology = read_top()

protein = topology.molecules[0]
protein.list_bonds(by_resid=True)