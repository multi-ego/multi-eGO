# This little module includes the parameters of the protein to study
# In this case there will be TTR and B2m parameters but it is possible to apply them
# theoretically always

protein = 'TTR'
#protein = 'Amylin'
#protein = 'B2m'
#protein = 'harp'

distance_cutoff = 5.5
distance_residue = 3
epsilon_input = 0.310

# This option requires a long simulation using explicit solvent and a first run of monomer_pairs.py
# It will make a list of pairs to reweight their epsilon and add to the full force field
# Add monomer_pairs.py in GRETA so that it can do everything automatically
idp = True
ratio_treshold = 0.09

# Does the N_terminal have a protonation?
# Taking off N_1 N_1 and raising c12 of N_1
N_terminal = True
 
# NMR or minimum
sigma_method = 'minimum'
#sigma_method = 'NMR'