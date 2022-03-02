from time import localtime, strftime
from numpy import column_stack
from topology_parser import topology_moleculetype, topology_atoms, topology_bonds, topology_angles, topology_dihedrals, topology_impropers, topology_system, topology_molecules
import pandas as pd

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


def write_LJ(atomtypes, greta_LJ, parameters):
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

def write_topology(parameters, topology_sections_dict, pairs_topology='', exclusion_topology=''):
    pd.set_option('display.colheader_justify', 'left')
    file = open(f'{parameters["output_folder"]}/topol_GRETA.top', "w")
    file.write(header(parameters))  

    file.write('; Include forcefield parameters\n')  
    file.write('#include "paste/the/ff/folder/here"\n\n')  
    
    file.write('[ moleculetype ]\n')
    df = topology_moleculetype(topology_sections_dict).topology_moleculetype
    file.write(str(df.to_string(index=False)))
    file.write('\n\n')

    top_atoms = topology_atoms(topology_sections_dict).df_topology_atoms
    top_atoms.rename(columns = {'atom_number':'; nr', 'atom_type':'type', 'residue_number':'resnr'}, inplace=True)
    top_atoms['type'] = top_atoms['sb_type']
    top_atoms.drop(columns=['sb_type', 'mass'], inplace=True)
    file.write("[ atoms ]\n")
    file.write(str(top_atoms.to_string(index = False)))
    file.write('\n\n')

    file.write('[ bonds ]\n')
    df = topology_bonds(topology_sections_dict).df_topology_bonds
    df.rename(columns={'ai':'; ai'}, inplace=True)
    file.write(str(df.to_string(index=False)))
    file.write('\n\n')

    file.write('[ angles ]\n')
    df = topology_angles(topology_sections_dict).topology_angles
    df.rename(columns={'ai':'; ai'}, inplace=True)
    file.write(str(df.to_string(index=False)))
    file.write('\n\n')

    file.write('[ dihedrals ]\n')
    df = topology_dihedrals(topology_sections_dict).topology_dihedrals
    df.rename(columns={'ai':'; ai'}, inplace=True)
    file.write(str(df.to_string(index=False)))
    file.write('\n\n')

    file.write('[ dihedrals ]\n')
    df = topology_impropers(topology_sections_dict).topology_impropers
    df.rename(columns={'ai':'; ai'}, inplace=True)
    file.write(str(df.to_string(index=False)))
    file.write('\n\n')
    
    file.write("[ pairs ]")
    file.write("\n")
    file.write(str(pairs_topology.to_string(index = False)))
    file.write("\n\n")
    file.write("[ exclusions ]")
    file.write("\n")
    file.write(str(exclusion_topology.to_string(index = False)))
    file.write('\n\n')

    file.write('; Include Position restraint file\n')
    file.write('#ifdef POSRES\n')
    file.write('#include "posre.itp"\n')
    file.write('#endif\n\n')

    file.write('[ system ]\n')
    df = topology_system(topology_sections_dict).topology_system
    file.write(str(df.to_string(index=False)))
    file.write('\n\n')

    file.write('[ molecules ]\n')
    df = topology_molecules(topology_sections_dict).topology_molecules
    file.write(str(df.to_string(index=False)))
    file.write('\n\n')

    df = ''

    file.close()

