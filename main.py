#from os import getresgid
import os
from time import process_time_ns
import pandas as pd
import sys, getopt

from pytest import param
from topology_definitions import raw_top 
from read_input import read_pdbs, plainMD_mdmat, random_coil_mdmat
from write_output import write_greta_LJ, write_greta_atomtypes_atp, write_greta_topology_atoms, write_greta_topology_pairs, header
from greta import make_pairs_exclusion_topology, make_pairs, merge_GRETA, make_pdb_atomtypes, make_idp_epsilon
#from mdmat import mdmat 


#TODO unire RC a greta to keep
#TODO nomi delle cartelle in base a greta to keep
#TODO cartelle di input in base a plain, rc


def main(argv):


    parameters = {
        'distance_cutoff':5.5,
        'distance_residue':2,
        'ratio_treshold':0.001,
        #'u_treshold':1.-Self['ratio_treshold'],
        'u_treshold':1.-0.001,
        'N_terminal':False,
        # Settings for LJ 1-4. We introduce some LJ interactions otherwise lost with the removal of explicit H
        # The c12 of a LJ 1-4 is too big, therefore we reduce by a factor
        'lj_reduction':0.15,
        # For left alpha we might want to increase the c6 values
        'multiply_c6':1.5,
        # Acid FFnonbondend it only works on the native pairs
        'acid_ff':False
    }

    print('\n\nMulti-eGO (codename: GRETA)\n')

    make_random_coil = False
    idp = False

    try:
        opts, args = getopt.getopt(argv,"hrip:b:e:",["protein=","build-from=","random-coil","epsilon=","ensemble"])
    except getopt.GetoptError:
        print('main.py -p <protein> -b <native|fibril|all>' )
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('main.py -p <protein> -b <native|fibril|all>')
            sys.exit()
        elif opt in ("-p", "--protein"):
            parameters['protein'] = arg
        elif opt in ("-b", "--build-from"):
            parameters['greta_to_keep'] = arg
        elif opt in ("r", "--random-coil"):
            make_random_coil = True
        elif opt in ("e", "--epsilon"):
            global epsilon_input
            parameters['epsilon_input'] = float(arg)
            parameters['epsilon_structure'] = float(arg)
            parameters['epsilon_md'] = float(arg)
        elif opt in ("i", "--ensemble"):
            idp = True 

    print('Creating a force-field using the following parameters:')
    
    for k,v in parameters.items():
        print(f' -{k,v}')

    # Create the folders which will be used by the script
    output_directory = f"outputs/output_{parameters['protein']}_e{parameters['epsilon_input']}"
    try:
        os.mkdir(output_directory)
    except OSError as error:
        pass

    print('- reading TOPOLOGY')
    raw_topology_atoms, first_resid, acid_atp, topology_bonds, atom_topology_num = raw_top(parameters)
    print('- reading PDB')
    # TODO togli la lettura della fibrilla nel caso fibrilla -> *ARG
    native_pdb, fibril_pdb = read_pdbs(parameters)

    print('- Generating Atomtypes')
    native_atomtypes, fibril_atomtypes, ffnonbonded_atp, atomtypes_atp, topology_atoms, type_c12_dict, proline_n = make_pdb_atomtypes(native_pdb, fibril_pdb, raw_topology_atoms)

    # TODO correggi la massa dell'azoto in prolina
    write_greta_atomtypes_atp(atomtypes_atp, parameters, output_directory)
    write_greta_topology_atoms(topology_atoms, parameters, output_directory)

    print('- Generating LJ Interactions')

    if make_random_coil == True:
        greta_ffnb = pd.DataFrame(columns=['; ai', 'aj', 'type', 'c6', 'c12', '', 'sigma', 'epsilon'])
        write_greta_LJ(ffnonbonded_atp, greta_ffnb, acid_atp, parameters, output_directory)
        topology_pairs, topology_exclusion = make_pairs_exclusion_topology(type_c12_dict, proline_n, raw_topology_atoms, topology_bonds, atom_topology_num, parameters)
        write_greta_topology_pairs(topology_pairs, topology_exclusion, parameters, output_directory)

    else:
        if parameters['greta_to_keep'] == 'native':
            if idp == True:
                atomic_mat_plainMD = plainMD_mdmat(parameters)
                atomic_mat_random_coil = random_coil_mdmat(parameters)
                greta_LJ = make_idp_epsilon(atomic_mat_plainMD, atomic_mat_random_coil, parameters)
                check = set(greta_LJ['ai'].to_list() + greta_LJ['aj'].to_list())
                check = (sorted(check))
            else:
                atomic_mat_random_coil = random_coil_mdmat(parameters)
                greta_LJ = make_pairs(native_pdb, atomic_mat_random_coil, native_atomtypes, parameters)
                if parameters['acid_ff'] == True and acid_atp !=0:
                        greta_LJ = greta_LJ[~greta_LJ.ai.isin(acid_atp)]
                        greta_LJ = greta_LJ[~greta_LJ.aj.isin(acid_atp)]

        elif parameters['greta_to_keep'] == 'fibril':
            atomic_mat_random_coil = random_coil_mdmat(parameters)
            greta_LJ = make_pairs(fibril_pdb, atomic_mat_random_coil, fibril_atomtypes, parameters)

        elif parameters['greta_to_keep'] == 'all':
            if idp == True:
                atomic_mat_plainMD = plainMD_mdmat(parameters)
                atomic_mat_random_coil = random_coil_mdmat(parameters)
                # Contacts are from a plain MD, so at this step we just import the fibril contacts.
                greta_LJ = make_idp_epsilon(atomic_mat_plainMD, atomic_mat_random_coil, parameters).append(make_pairs(structure_pdb=fibril_pdb, atomic_mat_random_coil=atomic_mat_random_coil, atomtypes=fibril_atomtypes, parameters=parameters), sort = False, ignore_index = True)
            else:
                atomic_mat_random_coil = random_coil_mdmat(parameters)
                # Merging native and fibril LJ pairs.
                greta_LJ = make_pairs(native_pdb, atomic_mat_random_coil, native_atomtypes, parameters)
                if parameters['acid_ff'] == True and acid_atp !=0:
                        greta_LJ = greta_LJ[~greta_LJ.ai.isin(acid_atp)]
                        greta_LJ = greta_LJ[~greta_LJ.aj.isin(acid_atp)]
                greta_LJ.append(make_pairs(fibril_pdb, atomic_mat_random_coil, fibril_atomtypes, parameters), sort = False, ignore_index = True)

        else:
            print('ERRORONE')
            exit()

        print('- Merging LJ interactions')
        greta_ffnb = merge_GRETA(greta_LJ, parameters)
        if parameters['N_terminal'] == True:
            greta_ffnb.loc[(greta_ffnb['; ai'] == first_resid) & (greta_ffnb['aj'] == first_resid), '; ai'] = ';'+greta_ffnb['; ai'] 
        write_greta_LJ(ffnonbonded_atp, greta_ffnb, acid_atp, parameters, output_directory)

        print('- Generating Pairs and Exclusions')
        topology_pairs, topology_exclusion = make_pairs_exclusion_topology(type_c12_dict, proline_n, raw_topology_atoms, topology_bonds, atom_topology_num, parameters, greta_ffnb)
        write_greta_topology_pairs(topology_pairs, topology_exclusion, parameters, output_directory)

    print('\nGRETA complete! Carlo is happy\n')


if __name__ == "__main__":
   main(sys.argv[1:])