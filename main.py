#from os import getresgid
import os
import pandas as pd
import sys, getopt
from topology_definitions import raw_top 
from read_input import read_pdbs
from write_output import write_greta_LJ, write_greta_atomtypes_atp, write_greta_topology_atoms, write_greta_topology_pairs
from greta import make_pairs_exclusion_topology, make_pairs, merge_GRETA, make_pdb_atomtypes, make_idp_epsilon
from mdmat import mdmat 

# global variables
# modified from input arguments
make_random_coil = False
idp = False
# hardcoded 
distance_cutoff = 5.5
distance_residue = 2
ratio_treshold = 0.001
u_threshold = 1.-ratio_treshold
N_terminal = False
# Settings for LJ 1-4. We introduce some LJ interactions otherwise lost with the removal of explicit H
# The c12 of a LJ 1-4 is too big, therefore we reduce by a factor
lj_reduction = 0.15
# For left alpha we might want to increase the c6 values
multiply_c6 = 1.5
# Acid FFnonbondend it only works on the native pairs
acid_ff = False

def main(argv):
    print('\n\nGRETA:\n')

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
            protein = arg
        elif opt in ("-b", "--build-from"):
            greta_to_keep = arg
        elif opt in ("r", "--random-coil"):
            global make_random_coil 
            make_random_coil = True
        elif opt in ("e", "--epsilon"):
            global epsilon_input
            epsilon_input = arg
        elif opt in ("i", "--ensemble"):
            global idp 
            idp = True 

    epsilon_structure = float(epsilon_input)
    epsilon_md = float(epsilon_input)

    # Create the folders which will be used by the script
    output_folder = 'outputs/output_%s' % (protein)
    try:
        os.mkdir(output_folder)
    except OSError as error:
        pass
        #print(error)

    output_folder = 'FF_greta_analysis_%s' % (protein)
    try:
        os.mkdir(output_folder)
    except OSError as error:
        pass
        #print(error)

    print('- reading TOPOLOGY')
    raw_topology_atoms, first_resid, acid_atp, topology_bonds, atom_topology_num = raw_top(protein)
    print('- reading PDB')
    # TODO togli la lettura della fibrilla nel caso fibrilla -> *ARG
    native_pdb, fibril_pdb = read_pdbs(protein, greta_to_keep)

    print('- Generating Atomtypes')
    native_atomtypes, fibril_atomtypes, ffnonbonded_atp, atomtypes_atp, topology_atoms, type_c12_dict, proline_n = make_pdb_atomtypes(native_pdb, fibril_pdb, raw_topology_atoms)

    # TODO correggi la massa dell'azoto in prolina
    write_greta_atomtypes_atp(atomtypes_atp, protein)
    write_greta_topology_atoms(topology_atoms, protein, distance_cutoff, distance_residue)

    print('- Generating LJ Interactions')

    if make_random_coil == True:
        greta_ffnb = pd.DataFrame(columns=['; ai', 'aj', 'type', 'c6', 'c12', '', 'sigma', 'epsilon'])
        write_greta_LJ(ffnonbonded_atp, greta_ffnb, acid_atp, protein, distance_cutoff, distance_residue, greta_to_keep, epsilon_md, epsilon_structure, ratio_treshold, acid_ff)
        topology_pairs, topology_exclusion = make_pairs_exclusion_topology(type_c12_dict, proline_n, raw_topology_atoms, topology_bonds, atom_topology_num)
        write_greta_topology_pairs(topology_pairs, topology_exclusion, protein, distance_cutoff, distance_residue, lj_reduction, greta_to_keep)

    else:
        if greta_to_keep == 'native':
            if idp == True:
                atomic_mat_plainMD, atomic_mat_random_coil = mdmat(protein, distance_residue, ratio_treshold)
                greta_LJ = make_idp_epsilon(atomic_mat_plainMD, atomic_mat_random_coil, epsilon_md, ratio_treshold)
                check = set(greta_LJ['ai'].to_list() + greta_LJ['aj'].to_list())
                check = (sorted(check))
            else:
                greta_LJ = make_pairs(native_pdb, atomic_mat_random_coil, native_atomtypes, epsilon_structure, ratio_treshold, distance_cutoff, distance_residue, u_threshold)
                if acid_ff == True and acid_atp !=0:
                        greta_LJ = greta_LJ[~greta_LJ.ai.isin(acid_atp)]
                        greta_LJ = greta_LJ[~greta_LJ.aj.isin(acid_atp)]

        elif greta_to_keep == 'fibril':
            greta_LJ = make_pairs(fibril_pdb, atomic_mat_random_coil, fibril_atomtypes, epsilon_structure, ratio_treshold, distance_cutoff, distance_residue, u_threshold)

        elif greta_to_keep == 'all':
            if idp == True:
                atomic_mat_plainMD, atomic_mat_random_coil = mdmat(protein, distance_residue, ratio_treshold)
                # Contacts are from a plain MD, so at this step we just import the fibril contacts.
                greta_LJ = make_idp_epsilon(atomic_mat_plainMD, atomic_mat_random_coil, epsilon_md, ratio_treshold).append(make_pairs(structure_pdb=fibril_pdb, atomic_mat_random_coil=atomic_mat_random_coil, atomtypes=fibril_atomtypes, epsilon_structure=epsilon_structure, ratio_treshold=ratio_treshold, distance_cutoff=distance_cutoff, distance_residue=distance_residue, u_threshold=u_threshold), sort = False, ignore_index = True)
            else:
                # Merging native and fibril LJ pairs.
                greta_LJ = make_pairs(native_pdb, atomic_mat_random_coil, native_atomtypes, epsilon_structure, ratio_treshold, distance_cutoff, distance_residue, u_threshold)
                if acid_ff == True and acid_atp !=0:
                        greta_LJ = greta_LJ[~greta_LJ.ai.isin(acid_atp)]
                        greta_LJ = greta_LJ[~greta_LJ.aj.isin(acid_atp)]
                greta_LJ.append(make_pairs(fibril_pdb, atomic_mat_random_coil, fibril_atomtypes, epsilon_structure, ratio_treshold, distance_cutoff, distance_residue, u_threshold), sort = False, ignore_index = True)

        else:
            print('ERRORONE')
            exit()

        print('- Merging LJ interactions')
        greta_ffnb = merge_GRETA(greta_LJ, epsilon_structure)
        if N_terminal == True:
            greta_ffnb.loc[(greta_ffnb['; ai'] == first_resid) & (greta_ffnb['aj'] == first_resid), '; ai'] = ';'+greta_ffnb['; ai'] 
        write_greta_LJ(ffnonbonded_atp, greta_ffnb, acid_atp, protein, distance_cutoff, distance_residue, greta_to_keep, epsilon_md, epsilon_structure, ratio_treshold, acid_ff)

        print('- Generating Pairs and Exclusions')
        topology_pairs, topology_exclusion = make_pairs_exclusion_topology(type_c12_dict, proline_n, raw_topology_atoms, topology_bonds, atom_topology_num, greta_ffnb)
        write_greta_topology_pairs(topology_pairs, topology_exclusion, protein, distance_cutoff, distance_residue, lj_reduction, greta_to_keep)

    print('\nGRETA complete! Carlo is happy\n')


if __name__ == "__main__":
   main(sys.argv[1:])