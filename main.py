#from os import getresgid
import os
import pandas as pd
import sys, getopt
from topology_definitions import raw_top 
from read_input import read_pdbs, plainMD_mdmat, random_coil_mdmat
from write_output import write_greta_LJ, write_greta_atomtypes_atp, write_greta_topology_atoms, write_greta_topology_pairs
from greta import make_pairs_exclusion_topology, make_pairs, merge_GRETA, make_pdb_atomtypes, make_more_atomtypes, make_idp_epsilon


#TODO unire RC a greta to keep


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

    idp = False

    try:
        opts, args = getopt.getopt(argv,"hip:b:e:",["protein=","build-from=","epsilon=","ensemble"])
    except getopt.GetoptError:
        print('main.py -p <protein> -b <single|merge|rc>' )
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('main.py -p <protein> -b <single|merge|rc>')
            sys.exit()
        elif opt in ("-p", "--protein"):
            parameters['protein'] = arg
        elif opt in ("-b", "--build-from"):
            parameters['greta_to_keep'] = arg
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
    output_directory = f"outputs/output_{parameters['protein']}_{parameters['greta_to_keep']}_e{parameters['epsilon_input']}"
    try:
        os.mkdir(output_directory)
    except OSError as error:
        pass

    print('- reading TOPOLOGY')
    raw_topology_atoms, first_resid, acid_atp, topology_bonds, atom_topology_num = raw_top(parameters)
    print('- reading PDB')
    native_pdb = read_pdbs(parameters, False)
    if parameters['greta_to_keep'] == 'merge':
        fibril_pdb = read_pdbs(parameters, True)

    print('- Generating Atomtypes')
    native_atomtypes, ffnonbonded_atp, atomtypes_atp, topology_atoms, type_c12_dict, proline_n = make_pdb_atomtypes(native_pdb, raw_topology_atoms)
    if parameters['greta_to_keep'] == 'merge':
        fibril_atomtypes = make_more_atomtypes(fibril_pdb)

    # TODO correggi la massa dell'azoto in prolina
    write_greta_atomtypes_atp(atomtypes_atp, parameters, output_directory)
    write_greta_topology_atoms(topology_atoms, parameters, output_directory)

    print('- Generating LJ Interactions')

    if parameters['greta_to_keep'] == 'rc':
        greta_ffnb = pd.DataFrame(columns=['; ai', 'aj', 'type', 'c6', 'c12', '', 'sigma', 'epsilon'])
        write_greta_LJ(ffnonbonded_atp, greta_ffnb, acid_atp, parameters, output_directory)
        print('- Generating Pairs and Exclusions')
        topology_pairs, topology_exclusion = make_pairs_exclusion_topology(type_c12_dict, proline_n, raw_topology_atoms, topology_bonds, atom_topology_num, parameters)
        write_greta_topology_pairs(topology_pairs, topology_exclusion, parameters, output_directory)

    elif parameters['greta_to_keep'] == 'single':
        if idp == True:
            atomic_mat_plainMD = plainMD_mdmat(parameters)
            atomic_mat_random_coil = random_coil_mdmat(parameters)
            greta_LJ = make_idp_epsilon(atomic_mat_plainMD, atomic_mat_random_coil, parameters)
        else:
            atomic_mat_random_coil = random_coil_mdmat(parameters)
            greta_LJ = make_pairs(native_pdb, atomic_mat_random_coil, native_atomtypes, parameters)
            if parameters['acid_ff'] == True and acid_atp !=0:
                    greta_LJ = greta_LJ[~greta_LJ.ai.isin(acid_atp)]
                    greta_LJ = greta_LJ[~greta_LJ.aj.isin(acid_atp)]

    elif parameters['greta_to_keep'] == 'merge':
        if idp == True:
            atomic_mat_plainMD = plainMD_mdmat(parameters)
            atomic_mat_random_coil = random_coil_mdmat(parameters)
            greta_LJ = make_idp_epsilon(atomic_mat_plainMD, atomic_mat_random_coil, parameters)
            greta_LJ = greta_LJ.append(make_pairs(fibril_pdb, atomic_mat_random_coil, fibril_atomtypes, parameters), sort = False, ignore_index = True)
        else:
            atomic_mat_random_coil = random_coil_mdmat(parameters)
            greta_LJ = make_pairs(native_pdb, atomic_mat_random_coil, native_atomtypes, parameters)
            if parameters['acid_ff'] == True and acid_atp !=0:
                    greta_LJ = greta_LJ[~greta_LJ.ai.isin(acid_atp)]
                    greta_LJ = greta_LJ[~greta_LJ.aj.isin(acid_atp)]
            greta_LJ = greta_LJ.append(make_pairs(fibril_pdb, atomic_mat_random_coil, fibril_atomtypes, parameters), sort = False, ignore_index = True)

    else:
        print("I dont' understand --build-from=",parameters['greta_to_keep'])
        exit()

    if parameters['greta_to_keep'] != 'rc':
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
