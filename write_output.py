import os
import datetime
from protein_configuration import temperatura, protein, distance_residue, distance_cutoff


    # Create the folders which will be used by the script
output_folder = 'output_%s' % (protein)
try:
    os.mkdir(output_folder)
except OSError:
    print (" Creation of the directory %s failed" % output_folder)
else:
    print ("Successfully created the directory %s" % output_folder)



##### GRETA
    # Create the folders which will be used by the script
output_folder = 'GRETA/output_%s' % (protein)
try:
    os.mkdir(output_folder)
except OSError:
    print (" Creation of the directory %s failed" % output_folder)
else:
    print ("Successfully created the directory %s" % output_folder)


head = "; Made using MAGROS.FF script by Emanuele Scalone at Camilloni Lab"
now = datetime.datetime.now()


def write_atomtypes_atp(pep_atomtypes):
    # This function is used to create the atomtypes.atp.
    #file = open("../../magros_test/commons/output/atomtypes.atp", "w")
    directory = 'output_%s/atomtypes.atp' % (protein)
    file = open(directory, "w")
    file.write(str(head))
    file.write("\n")
    file.write('; ' + str(protein))
    file.write("\n")
    file.write('; ' + str(now))
    file.write("\n")
    file.write("[ atomtypes ]")
    file.write("\n")
    file.write(str(pep_atomtypes.to_string(index = False, header = False)))
    file.close()

def write_gromos_topology(gromos_topology):
    # This function creates the atomtypes file to paste into the custom topology
    #file = open("../../magros_test/commons/output/topology_gromos", "w")
    directory = "output_%s/topology_gromos" % (protein)
    file = open(directory, "w")
    file.write(str(head))
    file.write("\n")
    file.write('; ' + str(protein))
    file.write("\n")
    file.write('; ' + str(now))
    file.write("\n")
    file.write("[ atoms ]")
    file.write("\n")
    file.write(str(gromos_topology.to_string(index = False)))
    file.close()


def write_smog_to_gromos_dihedrals(propers_to_gro):
    # This function creates the angles file to paste into the topology
    #file = open("../../magros_test/commons/output/smog_to_gromos_dihedrals", "w")
    directory = "output_%s/smog_to_gromos_dihedrals" % (protein)
    file = open(directory, "w")
    file.write(str(head))
    file.write("\n")
    file.write('; ' + str(protein))
    file.write("\n")
    file.write('; ' + str(now))
    file.write("\n")
    file.write(f'; Temperature Ratio: {temperatura} / 70')
    file.write("\n")
    file.write("[ dihedrals ]")
    file.write("\n")
    file.write(str(propers_to_gro.to_string(index = False)))
    file.close()

def write_merge_ffnonbonded(atomtypes, merge_pairs):
    #file = open("../../magros_test/commons/output/ffnonbonded.itp", "w")
    directory = "output_%s/ffnonbonded.itp" % (protein)
    file = open(directory, "w")
    file.write(str(head))
    file.write("\n")
    file.write('; ' + str(protein))
    file.write("\n")
    file.write('; ' + str(now))
    file.write("\n")
    file.write(f'; Temperature Ratio: {temperatura} / 70')
    file.write("\n")
    file.write("[ atomtypes ]\n")
    file.write(str(atomtypes.to_string(index = False)))
    file.write("\n")
    file.write("\n")
    file.write("[ nonbond_params ]\n")
    file.write(str(merge_pairs.to_string(index = False)))
    file.close()

def write_acid_ffnonbonded(atomtypes, acid_pairs):
    #file = open("../../magros_test/commons/output/ffnonbonded.itp", "w")
    directory = "output_%s/acid_ffnonbonded.itp" % (protein)
    file = open(directory, "w")
    file.write(str(head))
    file.write("\n")
    file.write('; ' + str(protein))
    file.write("\n")
    file.write('; ' + str(now))
    file.write("\n")
    file.write(f'; Temperature Ratio: {temperatura} / 70')
    file.write("\n")
    file.write("; ACID PAIRS")
    file.write("\n")
    file.write("[ atomtypes ]\n")
    file.write(str(atomtypes.to_string(index = False)))
    file.write("\n")
    file.write("\n")
    file.write("[ nonbond_params ]\n")
    file.write(str(acid_pairs.to_string(index = False)))
    file.close()


##################
#### GRETA


def write_greta_atomtypes_atp(atomtypes_atp):
    # This function is used to create the atomtypes.atp.
    #file = open("../../magros_test/commons/output/atomtypes.atp", "w")
    directory = 'GRETA/output_%s/atomtypes.atp' %(protein)
    file = open(directory, "w")
    file.write(str(head))
    file.write("\n")
    file.write('; ' + 'GRETA TTR')
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
    directory = 'GRETA/output_%s/topology_atoms' %(protein)
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

def write_greta_LJ(atomtypes, greta_LJ):
    #file = open("../../magros_test/commons/output/ffnonbonded.itp", "w")
    directory = "GRETA/output_%s/ffnonbonded.itp" %(protein)
    file = open(directory, "w")
    file.write(str(head))
    file.write("\n")
    file.write('; Distance cutoff: ' + str(distance_cutoff))
    file.write("\n")
    file.write('; Residue cutoff: ' + str(distance_residue))
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

def write_acid_greta_LJ(atomtypes, greta_LJ):
    #file = open("../../magros_test/commons/output/ffnonbonded.itp", "w")
    directory = "GRETA/output_%s/acid_ffnonbonded.itp" %(protein)
    file = open(directory, "w")
    file.write(str(head))
    file.write("\n")
    file.write('; Distance cutoff: ' + str(distance_cutoff), 'ACID')
    file.write("\n")
    file.write('; Residue cutoff: ' + str(distance_residue))
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

def write_pairs_list(greta_LJ):
    directory = "GRETA/output_%s/pairs_list.txt" %(protein)
    file = open(directory, "w")
    file.write(str(greta_LJ.to_string(index = False)))
    file.close()