# This little module includes the parameters of the protein to study
# In this case there will be TTR and B2m parameters but it is possible to apply them
# theoretically always
#protein = 'TTR'
#protein = 'ALA_DP'
#protein = 'PRO_DP'
#protein = 'GLY_DP'
#protein = 'Amylin'
#protein = 'B2m'
#protein = 'harp'
protein = 'ABeta'

distance_cutoff = 5.5
distance_residue = 2
epsilon_input = 0.200

# This option requires a long simulation using explicit solvent and a first run of monomer_pairs.py
# It will make a list of pairs to reweight their epsilon and add to the full force field
# Add monomer_pairs.py in GRETA so that it can do everything automatically
idp = False
#idp = False
ratio_treshold = 0.008

# Does the N_terminal have a protonation?
# Taking off N_1 N_1 and raising c12 of N_1
#N_terminal = True
N_terminal = False

# Settings for LJ 1-4. We introduce some LJ interactions otherwise lost with the removal of explicit H
# The c12 of a LJ 1-4 is too big, therefore we reduce by a factor
lj_reduction = 0.15

# For left alpha we might want to increase the c6 values
multiply_c6 = 1.5

# If you want some tests by using only native or only fibril
greta_to_keep = 'native'
#greta_to_keep = 'fibril'
#greta_to_keep = 'all'

# Acid FFnonbondend it only works on the native pairs
acid_ff = False