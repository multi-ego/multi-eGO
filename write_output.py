import os
import datetime

head = "; Made using MAGROS.FF script by Emanuele Scalone at Camilloni Lab"
now = datetime.datetime.now()

def write_greta_atomtypes_atp(atomtypes_atp, protein):
    # This function is used to create the atomtypes.atp.
    directory = 'outputs/output_%s/atomtypes.atp' %(protein)
    file = open(directory, "w")
    file.write(str(head))
    file.write("\n")
    file.write('; ' + f'GRETA {protein}')
    file.write("\n")
    file.write('; ' + str(now))
    file.write("\n")
    file.write("[ atomtypes ]")
    file.write("\n")
    file.write(str(atomtypes_atp.to_string(index = False, header = False)))
    file.close()

def write_greta_topology_atoms(topology_atoms, protein, distance_cutoff, distance_residue):
    directory = 'outputs/output_%s/topology_atoms' %(protein)
    file = open(directory, "w")
    file.write(str(head))
    file.write("\n")
    file.write('; Distance cutoff: ' + str(distance_cutoff))
    file.write("\n")
    file.write('; Residue cutoff: ' + str(distance_residue))
    file.write("\n")
    file.write('; ' + str(now))
    file.write("\n")
    file.write("[ atoms ]")
    file.write("\n")
    file.write(str(topology_atoms.to_string(index = False)))
    file.close()

def write_greta_topology_pairs(pairs_topology, exclusion_topology, protein, distance_cutoff, distance_residue, lj_reduction, greta_to_keep):
    directory = 'outputs/output_%s/topology_pairs' %(protein)
    file = open(directory, "w")
    file.write(str(head))
    file.write("\n")
    file.write('; Distance cutoff: ' + str(distance_cutoff))
    file.write("\n")
    file.write('; Residue cutoff: ' + str(distance_residue))
    file.write("\n")
    file.write('; LJ_reduction: ' + str(lj_reduction))
    file.write("\n")
    file.write('; Protein: ' + str(protein) + '-' + str(greta_to_keep))
    file.write("\n")
    file.write('; ' + str(now))
    file.write("\n\n")
    file.write("[ pairs ]")
    file.write("\n")
    file.write(str(pairs_topology.to_string(index = False)))
    file.write("\n\n")
    file.write("[ exclusions ]")
    file.write("\n")
    file.write(str(exclusion_topology.to_string(index = False)))
    file.close()
    print('- Pairs and Exclusions written')

def write_greta_LJ(atomtypes, greta_LJ, acid_atp, protein, distance_cutoff, distance_residue, greta_to_keep, epsilon_md, epsilon_structure, ratio_treshold, acid_ff):
    if acid_ff == True and acid_atp !=0:
        directory = "outputs/output_%s/acid_ffnonbonded.itp" %(protein)
        file = open(directory, "w")
        file.write(str(head))
        file.write("\n")
        file.write('; Distance cutoff: ' + str(distance_cutoff) + ' ACID')
    else:
        directory = "outputs/output_%s/ffnonbonded.itp" %(protein)
        file = open(directory, "w")
        file.write(str(head))
        file.write("\n")
        file.write(f'; Distance cutoff: {distance_cutoff} - aminoacid exclusion: {distance_residue}')

    file.write("\n")
    file.write(f'; FF parameters {protein}-{greta_to_keep} - epsilon: {epsilon_md, epsilon_structure} - idp contacts treshold: {ratio_treshold}')
    file.write("\n")
    file.write('; ' + str(now))
    file.write("\n")
    file.write("[ atomtypes ]\n")
    file.write(str(atomtypes.to_string(index = False)))
    file.write("\n")
    file.write("\n")
    file.write("[ nonbond_params ]\n")
    if greta_LJ.empty:
        file.write(str('; ai, aj, type, c6, c12, sigma, epsilon'))
    else:
        file.write(str(greta_LJ.to_string(index = False)))
    file.close()

    print('- FF Written.')
    print('\tChange the masses and copy ffnonbonded.itp and atomtypes.atp into the ff folder.')

def write_pairs_list(pairs, ff_name, protein, distance_cutoff, distance_residue):
    pairs.rename(columns = {'; ai':'ai'}, inplace = True)
    pairs_analysis = pairs[['ai', 'aj', 'sigma', 'epsilon']].copy()
    directory = "FF_greta_analysis_%s/%s_%s_pairs_list_c%s_e%s.txt" %(protein, protein, ff_name, distance_cutoff, distance_residue)
    file = open(directory, "w")
    file.write(str(pairs_analysis.to_string(index = False)))
    file.close()
