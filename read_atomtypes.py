import itertools
import pandas as pd
from read_input import read_top
import numpy as np

# Import the topology informations
topology = read_top()
protein = topology.molecules[0]
print(protein.list_atoms)