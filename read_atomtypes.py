import itertools
import pandas as pd
from read_input import read_top
import numpy as np

# Import the topology informations
topology = read_top()
protein = topology.molecules[0]
atom_topology_num, atom_topology_type, atom_topology_resid, atom_topology_resname, atom_topology_name, atom_topology_mass  = protein.list_atoms()
topology_atoms = pd.DataFrame(np.column_stack([atom_topology_num, atom_topology_type, atom_topology_resid, atom_topology_resname, atom_topology_name, atom_topology_mass]), columns=['nr', 'type','resnr', 'residue', 'atom', 'mass'])