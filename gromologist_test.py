from read_input import read_top
from gromologist import Top, Pdb
import sys


#GROMOLOGY

native_top = read_top()
print(native_top)
native_top.check_pdb()
native_top.list_molecules()

chains = native_top.molecules[0]

#chains.list_bonds()

original_stdout = sys.stdout

with open('bonds.txt', 'w') as f:
    sys.stdout = f
    chains.list_bonds()
    sys.stdout = original_stdout