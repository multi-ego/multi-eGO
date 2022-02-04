import os
import pandas as pd
import sys, getopt
from read_input import read_pdbs, plainMD_mdmat, random_coil_mdmat, read_topology
from write_output import write_greta_LJ, write_greta_atomtypes_atp, write_greta_topology_atoms, write_greta_topology_pairs
from greta import make_pairs_exclusion_topology, make_pairs, merge_GRETA, make_pdb_atomtypes, make_more_atomtypes, make_idp_epsilon

def main(argv):

    parameters = {
        'distance_cutoff':5.5,
        'distance_residue':2,
        'ratio_treshold':0.001,
        'N_terminal':False,
        # Settings for LJ 1-4. We introduce some LJ interactions otherwise lost with the removal of explicit H
        # The c12 of a LJ 1-4 is too big, therefore we reduce by a factor
        'lj_reduction':0.15,
        # For left alpha we might want to increase the c6 values
        'multiply_c6':1.5,
        # Acid FFnonbondend it only works on the native pairs
        'acid_ff':False, 
        'idp':True
    }

    print('\n\nMulti-eGO (codename: GRETA)\n')

    readall=0

    try:
        opts, args = getopt.getopt(argv,"",["protein=","egos=","epsilon=","noensemble"])
    except getopt.GetoptError:
        print('multiego.py --protein <protein> --egos <single|merge|rc> --noensemble (optional)' )
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('main.py -p <protein> -b <single|merge|rc>')
            sys.exit()
        elif opt in ("--protein"):
            if not arg:
                print('Provide a protein name')
                sys.exit()
            else:
                parameters['protein'] = arg
                readall +=1
        elif opt in ("--egos"):
            if arg in ('single', 'merge', 'rc'):
                parameters['egos'] = arg
                readall +=1
            else:
                print('--egos accepts <single|merge|rc> options')
                sys.exit()
        elif opt in ("--epsilon"):
            arg = float(arg)
            if arg > 1 or arg < 0:
                print('Epsilon values must be chosen between 0 and 1')
                sys.exit()
            else:
                parameters['epsilon_input'] = float(arg)
                parameters['epsilon_structure'] = float(arg)
                parameters['epsilon_md'] = float(arg)
                readall +=1
        elif opt in ("--noensemble"):
            parameters['idp'] = False 
  
    if readall != 3:
        print('ERROR: missing input argument')
        print('multiego.py --protein <protein> --egos <single|merge|rc> --noensemble (optional)' )
        exit()
 
    print('- Creating a multi-eGO force-field using the following parameters:')
    
    for k,v in parameters.items():
        print(f'\t{k}: {v}')

    # Create the folders which will be used by the script
    output_directory = f"outputs/output_{parameters['protein']}_{parameters['egos']}_e{parameters['epsilon_input']}"
    try:
        os.mkdir(output_directory)
    except OSError as error:
        pass

    print('- reading TOPOLOGY')
    first_resid = read_topology(parameters)[0].first_resid
    acid_atp = read_topology(parameters)[0].acid_atp

    print('- reading PDB')
    native_pdb = read_pdbs(parameters, False)
    if parameters['egos'] == 'merge':
        fibril_pdb = read_pdbs(parameters, True)

    print('- Generating Atomtypes')
    native_atomtypes, ffnonbonded_atp, atomtypes_atp, topology_atoms, type_c12_dict, proline_n = make_pdb_atomtypes(native_pdb, parameters)
    if parameters['egos'] == 'merge':
        fibril_atomtypes = make_more_atomtypes(fibril_pdb)

    # TODO correggi la massa dell'azoto in prolina
    write_greta_atomtypes_atp(atomtypes_atp, parameters, output_directory)
    write_greta_topology_atoms(topology_atoms, parameters, output_directory)

    print('- Generating LJ Interactions')

    if parameters['egos'] == 'rc':
        greta_ffnb = pd.DataFrame(columns=['; ai', 'aj', 'type', 'c6', 'c12', '', 'sigma', 'epsilon'])
        write_greta_LJ(ffnonbonded_atp, greta_ffnb, acid_atp, parameters, output_directory)
        print('- Generating Pairs and Exclusions')
        topology_pairs, topology_exclusion = make_pairs_exclusion_topology(type_c12_dict, proline_n, parameters)
        write_greta_topology_pairs(topology_pairs, topology_exclusion, parameters, output_directory)

    elif parameters['egos'] == 'single':
        if parameters['idp'] == True:
            atomic_mat_plainMD = plainMD_mdmat(parameters)
            atomic_mat_random_coil = random_coil_mdmat(parameters)
            greta_LJ = make_idp_epsilon(atomic_mat_plainMD, atomic_mat_random_coil, parameters)
        else:
            atomic_mat_random_coil = random_coil_mdmat(parameters)
            greta_LJ = make_pairs(native_pdb, atomic_mat_random_coil, native_atomtypes, parameters)
            if parameters['acid_ff'] == True and acid_atp !=0:
                    greta_LJ = greta_LJ[~greta_LJ.ai.isin(acid_atp)]
                    greta_LJ = greta_LJ[~greta_LJ.aj.isin(acid_atp)]

    elif parameters['egos'] == 'merge':
        if parameters['idp'] == True:
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

    else: # Questo serve ancora?
        print("I dont' understand --build-from=",parameters['egos'])
        exit()

    if parameters['egos'] != 'rc':
        print('- Finalising LJ interactions')
        greta_ffnb = merge_GRETA(greta_LJ, parameters)
        if parameters['N_terminal'] == True:
            greta_ffnb.loc[(greta_ffnb['; ai'] == first_resid) & (greta_ffnb['aj'] == first_resid), '; ai'] = ';'+greta_ffnb['; ai'] 
        write_greta_LJ(ffnonbonded_atp, greta_ffnb, acid_atp, parameters, output_directory)

        print('- Generating Pairs and Exclusions')
        topology_pairs, topology_exclusion = make_pairs_exclusion_topology(type_c12_dict, proline_n, parameters, greta_ffnb)
        write_greta_topology_pairs(topology_pairs, topology_exclusion, parameters, output_directory)

    print('\nGRETA completed! Carlo is happy\n')


if __name__ == "__main__":
   main(sys.argv[1:])
