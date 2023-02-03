from time import localtime, strftime
from numpy import column_stack
from topology_parser import topology_parser, topology_ligand_bonds, topology_ligand_angles, topology_ligand_dihedrals
import pandas as pd

now = strftime("%d-%m-%Y %H:%M", localtime())

def header(parameters):
    header = f'; Multi-eGO force field provided by Emanuele Scalone and Carlo Camilloni at Camilloni Lab \n'
    header += f'; Created on the {now}\n'
    header += f"; Protein name: {parameters['protein']} \n"
    header += f"; The force field type is: {parameters['egos']} \n"
    if parameters['egos'] != 'rc':
        header += f"; LJ epsilon: {parameters['epsilon_md']} \n"
        header += f"; LJ amyloid epsilon: {parameters['epsilon_amyl']} \n"
        header += f"; LJ potential from a MD/random_coil ratio and thresholds: {parameters['md_threshold']} {parameters['rc_threshold']}\n"
    header += f"; Atoms cutoff distance: {parameters['distance_cutoff']} A \n"
    header += "\n"

    return header

def write_atomtypes_atp(multiego_ensemble):
    # This function is used to create the atomtypes.atp.
    #directory = f"outputs/output_{parameters['protein']}"
    file = open(f'{multiego_ensemble.parameters["output_folder"]}/atomtypes.atp', "w")
    file.write(header(multiego_ensemble.parameters))
    file.write("\n")
    file.write(str(multiego_ensemble.atomtypes_atp_toWrite))
    file.close()


def write_LJ(ensemble):
    pd.set_option('display.colheader_justify', 'right')

    file = open(f'{ensemble.parameters["output_folder"]}/ffnonbonded.itp', "w")
    file.write(header(ensemble.parameters))
    file.write("[ atomtypes ]\n")
    file.write(str(ensemble.ffnonbonded_atp_toWrite))
    file.write("\n\n")
    file.write("[ nonbond_params ]\n")
    #if ensemble.parameters['egos'] == 'rc':
    #    file.write(str('; ai, aj, type, c6, c12, sigma, epsilon'))
    #else:
    if hasattr(ensemble, 'greta_ffnb_toWrite'):
        file.write(str(ensemble.greta_ffnb_toWrite))

    file.close()

def write_topology(ensemble):
    
    #parameters, topology_sections_dict, ego_topology, **ego_ligand):
    pd.set_option('display.colheader_justify', 'left')
    #top_to_write = topology_parser(topology_sections_dict)

    file = open(f'{ensemble.parameters["output_folder"]}/topol_GRETA.top', "w")
    file.write(header(ensemble.parameters))  

    file.write('; Include forcefield parameters\n')  
    file.write('#include "multi-ego-basic.ff/forcefield.itp"\n\n')  
    
    file.write('[ moleculetype ]\n')
    file.write(str(ensemble.moleculetype_toWrite))
    file.write('\n\n')

    file.write("[ atoms ]\n")
    file.write(str(ensemble.atomtypes_top_toWrite))
    file.write('\n\n')

    file.write('[ bonds ]\n')
    file.write(str(ensemble.bonds_toWrite))
    file.write('\n\n')

    file.write('[ angles ]\n')
    file.write(str(ensemble.angles_toWrite))
    file.write('\n\n')

    file.write('[ dihedrals ]\n')
    file.write(str(ensemble.dihedrals_toWrite))
    file.write('\n\n')

    file.write('[ dihedrals ]\n')
    file.write(str(ensemble.impropers_toWrite))
    file.write('\n\n')

    file.write("[ pairs ]")
    file.write("\n")
    file.write(str(ensemble.pairs_toWrite))
    file.write("\n\n")

    file.write("[ exclusions ]")
    file.write("\n")
    file.write(str(ensemble.exclusions_toWrite))
    file.write('\n\n')

    file.write('; Include Position restraint file\n')
    file.write('#ifdef POSRES\n')
    file.write('#include "posre.itp"\n')
    file.write('#endif\n\n')
    
    if ensemble.parameters['ligand'] == True:
        file.write('; Include ligand topology\n')
        file.write('#include "topol_ligand_GRETA.itp"')
        file.write('\n\n')

    file.write('[ system ]\n')
    file.write(str(ensemble.system_toWrite))
    file.write('\n\n')

    file.write('[ molecules ]\n')
    file.write(str(ensemble.molecules_toWrite))
    if ensemble.parameters['ligand'] == True:
        file.write('\nligand')
    file.write('\n\n')

    file.close()

def write_ligand_topology(ensemble):
    pd.set_option('display.colheader_justify', 'left')
    file = open(f'{ensemble.parameters["output_folder"]}/topol_ligand_GRETA.itp', "w")
    file.write(header(ensemble.parameters))  
    
    file.write('[ moleculetype ]\n')
    file.write(str(ensemble.ligand_moleculetype_toWrite))
    file.write('\n\n')

    file.write("[ atoms ]\n")
    file.write(str(ensemble.ligand_atomtypes_top_toWrite))
    file.write('\n\n')

    file.write('[ bonds ]\n')
    file.write(str(ensemble.ligand_bonds_toWrite))
    file.write('\n\n')

    file.write('[ angles ]\n')
    file.write(str(ensemble.ligand_angles_toWrite))
    file.write('\n\n')

    file.write('[ dihedrals ]\n')
    file.write(str(ensemble.ligand_dihedrals_toWrite))
    file.write('\n\n')

    file.write('[ pairs ]\n')
    file.write(str(ensemble.ligand_pairs_toWrite))
    file.write('\n\n')

    file.write('[ exclusions ]\n')
    file.write(str(ensemble.ligand_exclusions_toWrite))
    file.write('\n\n')

    file.close()
