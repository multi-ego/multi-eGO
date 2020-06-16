def write_atomtypes_atp(pep_atomtypes):
    # This function is used to create the atomtypes.atp.
    file = open("../../magros_test/commons/output/atomtypes.atp", "w")
    #file = open("output/atomtypes.atp", "w")
    file.write("[ atomtypes ]")
    file.write("\n")
    file.write(str(pep_atomtypes.to_string(index = False, header = False)))
    file.close()


def write_gromos_topology(gromos_topology):
    # This function creates the atomtypes file to paste into the custom topology
    file = open("../../magros_test/commons/output/topology_gromos", "w")
    #file = open("output/topology_gromos", "w")
    file.write("[ atoms ]")
    file.write("\n")
    file.write(str(gromos_topology.to_string(index = False)))
    file.close()


def write_smog_to_gromos_dihedrals(propers_to_gro):
    # This function creates the angles file to paste into the topology
    file = open("../../magros_test/commons/output/smog_to_gromos_dihedrals", "w")
    #file = open("output/smog_to_gromos_dihedrals", "w")
    file.write("[ dihedrals ]")
    file.write("\n")
    file.write(str(propers_to_gro.to_string(index = False)))
    file.close()

def write_merge_ffnonbonded(atomtypes, merge_pairs):
    file = open("../../magros_test/commons/output/ffnonbonded.itp", "w")
    #file = open("output/ffnonbonded.itp", "w")
    file.write("[ atomtypes ]\n")
    file.write(str(atomtypes.to_string(index = False)))
    file.write("\n")
    file.write("\n")
    file.write("[ nonbond_params ]\n")
    file.write(str(merge_pairs.to_string(index = False)))
    file.close()