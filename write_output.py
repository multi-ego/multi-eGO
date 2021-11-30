import os
import datetime
from protein_configuration import protein, distance_residue, distance_cutoff, lj_reduction, acid_ff, greta_to_keep
from topology_definitions import acid_atp


    # Create the folders which will be used by the script
output_folder = 'outputs/output_%s' % (protein)
try:
    os.mkdir(output_folder)
except OSError:
    print (" Creation of the directory %s failed" % output_folder)
else:
    print ("Successfully created the directory %s" % output_folder)


output_folder = 'FF_greta_analysis_%s' % (protein)
try:
    os.mkdir(output_folder)
except OSError:
    print (" Creation of the directory %s failed" % output_folder)
else:
    print ("Successfully created the directory %s" % output_folder)

head = "; Made using MAGROS.FF script by Emanuele Scalone at Camilloni Lab"
now = datetime.datetime.now()

def write_greta_atomtypes_atp(atomtypes_atp):
    # This function is used to create the atomtypes.atp.
    #file = open("../../magros_test/commons/output/atomtypes.atp", "w")
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

def write_greta_topology_atoms(topology_atoms):
    # This function is used to create the atomtypes.atp.
    #file = open("../../magros_test/commons/output/atomtypes.atp", "w")
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

def write_greta_topology_pairs(pairs_topology, exclusion_topology):
    # This function is used to create the atomtypes.atp.
    #file = open("../../magros_test/commons/output/atomtypes.atp", "w")
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
    file.write('; Protein: ' + str(protein) + str(greta_to_keep))
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


def write_greta_LJ(atomtypes, greta_LJ):
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
        file.write('; Distance cutoff: ' + str(distance_cutoff))

    file.write("\n")
    file.write('; Residue cutoff: ' + str(distance_residue) + str(greta_to_keep) + str(protein))
    file.write("\n")
    file.write('; ' + str(now))
    file.write("\n")
    file.write("[ atomtypes ]\n")
    file.write(str(atomtypes.to_string(index = False)))
    file.write("\n")
    file.write("\n")
    file.write("[ nonbond_params ]\n")
    file.write(str(greta_LJ.to_string(index = False)))
    file.close()

def write_pairs_list(pairs, ff_name):
    pairs.rename(columns = {'; ai':'ai'}, inplace = True)
    pairs_analysis = pairs[['ai', 'aj', 'sigma', 'epsilon']].copy()
    directory = "FF_greta_analysis_%s/%s_%s_pairs_list_c%s_e%s.txt" %(protein, protein, ff_name, distance_cutoff, distance_residue)
    file = open(directory, "w")
    file.write(str(pairs_analysis.to_string(index = False)))
    file.close()