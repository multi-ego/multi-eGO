def write_atomtypes_atp(pep_atomtypes):
    # This function is used to create the atomtypes.atp.
    file = open("output/atomtypes.atp", "w")
    file.write("[ atomtypes ]")
    file.write("\n")
    file.write(str(pep_atomtypes.to_string(index = False, header = False)))
    file.close()


def write_gromos_topology(gromos_topology):
    # This function creates the atomtypes file to paste into the custom topology
    file = open("output/topology_gromos", "w")
    file.write("[ atoms ]")
    file.write("\n")
    file.write(str(gromos_topology.to_string(index = False)))
    file.close()


def write_smog_to_gromos_dihedrals(propers_to_gro):
    # This function creates the angles file to paste into the topology
    file = open("output/smog_to_gromos_dihedrals", "w")
    file.write("[ dihedrals ]")
    file.write("\n")
    file.write(str(propers_to_gro.to_string(index = False)))
    file.close()