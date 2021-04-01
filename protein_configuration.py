# This little module includes the parameters of the protein to study
# In this case there will be TTR and B2m parameters but it is possible to apply them
# theoretically always

protein = 'TTR'
#protein = 'Amylin'
#protein = 'B2m'

distance_cutoff = 5.5
distance_residue = 3
epsilon_input = 0.295

# This option requires a long simulation using explicit solvent and a first run of monomer_pairs.py
# It will make a list of pairs to reweight their epsilon and add to the full force field
idp = True
ratio_treshold = 0.09