from time import localtime, strftime
now = strftime("%d-%m-%Y %H:%M", localtime())

def header(parameters):
    header = f'; Multi-eGO force field provided by Emanuele Scalone and Carlo Camilloni at Camilloni Lab \n'
    header = header + f'; Created on the {now}\n'
    header = header + f"; Protein name: {parameters['protein']} \n"
    header = header + f"; The force field type is: {parameters['egos']} \n"
    if parameters['egos'] != 'rc':
        header = header + f"; LJ epsilon: {parameters['epsilon_input']} \n"
    if parameters['ensemble'] == True:
        header = header + f"; LJ potential from a MD/random_coil ratio and threshold: {parameters['ratio_threshold']} \n"
    header = header + f"; Atoms cutoff distance: {parameters['distance_cutoff']} A \n"
    header = header + f"; Skipping contacts within {parameters['distance_residue']} residues \n"
    header = header + f"; Reducing the C12 N-X 1-3 C12 by: {parameters['lj_reduction']} \n"
    header = header + f"; Enhancing C6 for left alpha dihedral by: {parameters['multiply_c6']} \n"
    header = header + "\n"

    return header

def write_atomtypes_atp(atomtypes_atp, parameters):
    # This function is used to create the atomtypes.atp.
    #directory = f"outputs/output_{parameters['protein']}"
    file = open(f'{parameters["output_folder"]}/atomtypes.atp', "w")
    file.write(header(parameters))
    file.write("\n")
    file.write(str(atomtypes_atp.to_string(index = False, header = False)))
    file.close()

def write_topology_atoms(topology_atoms, parameters):
    #directory = f'outputs/output_{parameters["protein"]}'
    file = open(f'{parameters["output_folder"]}/topology_atoms', "w")
    file.write(header(parameters))
    file.write("[ atoms ]")
    file.write("\n")
    file.write(str(topology_atoms.to_string(index = False)))
    file.close()

def write_pairs_exclusion(pairs_topology, exclusion_topology, parameters):
    #directory = f'outputs/output_{parameters["protein"]}'
    file = open(f'{parameters["output_folder"]}/topology_pairs', "w")
    file.write(header(parameters))
    file.write("[ pairs ]")
    file.write("\n")
    file.write(str(pairs_topology.to_string(index = False)))
    file.write("\n\n")
    file.write("[ exclusions ]")
    file.write("\n")
    file.write(str(exclusion_topology.to_string(index = False)))
    file.close()

def write_LJ(atomtypes, greta_LJ, parameters):
    #directory = f"outputs/output_{parameters['protein']}/ffnonbonded.itp"
    file = open(f'{parameters["output_folder"]}/ffnonbonded.itp', "w")
    file.write(header(parameters))

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
