import pandas as pd

def read_gro_atoms():
    # Reading the atoms section from gromos topology
    gro_atoms = pd.read_csv('input/pep_gro_atoms', sep = "\s+", header = None)
    gro_atoms.columns = ["; nr", "type", "resnr", "residue", "atom", "cgnr", 'charge', 'mass']
    gro_atoms["atom_nmr"] = gro_atoms["atom"].apply(str) + '_' + gro_atoms["resnr"].apply(str)
    gro_atoms['res_atom'] = gro_atoms['residue'] + '_' + gro_atoms['atom']
    return gro_atoms