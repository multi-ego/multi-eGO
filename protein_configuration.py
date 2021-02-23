# This little module includes the parameters of the protein to study
# In this case there will be TTR and B2m parameters but it is possible to apply them
# theoretically always

#protein = 'TTR'
protein = 'Amylin'
#protein = 'B2m'

# This first version includes all the #flag comments in the script


# The fibril chain length is the number of the residues in the fibril.
# This is used in read_input.py when all the residues numbers are renumbered correctly
# (at the end of one chain, it starts another one instead keeping the same chain)

fibril_chain_length = 11 # 11 TTR
#fibril_chain_length = 64 # 64 B2m



# The fibril offset is used in the case of B2m, where there is not a N or C terminus.
# This way the first residue of the fibril has an offset of N residues compared to the native structure
# In this case is 22: the first residue of the fibril is the 23th of the native

fibril_residue_offset = 0 # TTR peptide has the same length in both fibril and native
#fibril_residue_offset = 22 # 22 B2m residues to add to the fibril



# Only one chain is used for dihedrals selection and the script will map the atom numbers based on the native structure
# Therefore those two values are needed to renumber the atomid of the fibril dihedrals:
# fibril_atom_number is the atom number of the first atom of the chain used for dihedral selection
# fibril_atom_offset is the first atom of the fibril in the native structure (just like the chain_residue offset)

fibril_atom_number = 3061 # TTR is the first atom in dihedrals fibril selection which will be mapped as 
fibril_atom_offset = 1 # the first residue in the case of the native chain
#fibril_atom_number = 540 # B2m 540 is the first residue of the fibril chain used to map the dihedrals
#fibril_atom_offset = 179 # B2m 179 is the same residue in the native structure



# Due to the uneven temperature balance between SMOG and GROMOS FF it is necessary to apply a temperature
# correction based on a ratio. Usually 70 K is the temperature of SMOG and the other one is obtained empirically
# And by empirical i mean to make different simulations and then check the results

temperatura = 300 # TTR
t_ratio = temperatura / 70 # TTR
#temperatura = 296 # B2m
#t_ratio = temperatura / 70 # B2m


########## GRETA
distance_cutoff = 5.5
distance_residue = 3
epsilon_input = 0.285