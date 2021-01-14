import pandas as pd
from read_input import read_pdbs


native_pdb, fibril_pdb = read_pdbs()

# Make a dictionary from the pdb of resnum:resid
amino_dict = dict(zip((list(native_pdb.residues.resids)),(list(native_pdb.residues.resnames))))

#print(amino_dict)

pairs_list = pd.read_csv('GRETA/output_TTR/pairs_list.txt', sep='\\s+')

pairs_list[
    
]

print(pairs_list)#.to_string())