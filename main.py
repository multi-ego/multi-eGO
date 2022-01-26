#from os import getresgid
from read_input import read_pdbs
from write_output import write_greta_LJ, write_greta_atomtypes_atp, write_greta_topology_atoms, write_greta_topology_pairs
from greta import make_pairs_exclusion_topology, make_pairs, merge_GRETA, make_pdb_atomtypes, make_idp_epsilon
from topology_definitions import acid_atp, first_resid
from protein_configuration import idp, N_terminal, greta_to_keep, acid_ff, make_random_coil
from mdmat import atomic_mat_plainMD, atomic_mat_random_coil
import pandas as pd

print('\n\n\n\n GRETA\n\n\n\n')
print('GRETA - PDB reading')


# TODO togli la lettura della fibrilla nel caso fibrilla -> *ARG
native_pdb, fibril_pdb = read_pdbs()

print('GRETA - Making Atomtypes')
native_atomtypes, fibril_atomtypes, ffnonbonded_atp, atomtypes_atp, topology_atoms, type_c12_dict, proline_n = make_pdb_atomtypes(native_pdb, fibril_pdb)
print('\n GRETA - Writing outputs')
write_greta_atomtypes_atp(atomtypes_atp)
write_greta_topology_atoms(topology_atoms)

print('\n GRETA - Making native and fibril pairs')

if make_random_coil == True:
    greta_ffnb = pd.DataFrame(columns=['; ai', 'aj', 'type', 'c6', 'c12', '', 'sigma', 'epsilon'])
    topology_pairs, topology_exclusion = make_pairs_exclusion_topology(type_c12_dict, proline_n)


else:
    if greta_to_keep == 'native':
        if idp == True:
            greta_LJ = make_idp_epsilon(atomic_mat_plainMD, atomic_mat_random_coil)
            check = set(greta_LJ['ai'].to_list() + greta_LJ['aj'].to_list())
            check = (sorted(check))
        else:
            greta_LJ = make_pairs(native_pdb, atomic_mat_random_coil, native_atomtypes)
            if acid_ff == True and acid_atp !=0:
                    greta_LJ = greta_LJ[~greta_LJ.ai.isin(acid_atp)]
                    greta_LJ = greta_LJ[~greta_LJ.aj.isin(acid_atp)]

    elif greta_to_keep == 'fibril':
        greta_LJ = make_pairs(fibril_pdb, atomic_mat_random_coil, fibril_atomtypes)

    elif greta_to_keep == 'all':
        if idp == True:
            # Contacts are from a plain MD, so at this step we just import the fibril contacts.
            greta_LJ = make_idp_epsilon(atomic_mat_plainMD, atomic_mat_random_coil).append(make_pairs(structure_pdb=fibril_pdb, atomic_mat_random_coil=atomic_mat_random_coil, atomtypes=fibril_atomtypes), sort = False, ignore_index = True)
        else:
            # Merging native and fibril LJ pairs.
            greta_LJ = make_pairs(native_pdb, atomic_mat_random_coil, native_atomtypes)
            if acid_ff == True and acid_atp !=0:
                    greta_LJ = greta_LJ[~greta_LJ.ai.isin(acid_atp)]
                    greta_LJ = greta_LJ[~greta_LJ.aj.isin(acid_atp)]
            greta_LJ.append(make_pairs(fibril_pdb, atomic_mat_random_coil, fibril_atomtypes), sort = False, ignore_index = True)

    else:
        print('ERRORONE')
        exit()

    print('\n GRETA - native and fibril pairs creation completed')
    print('\n GRETA - Merging native and fibril pairs')

    # TODO correggi la massa dell'azoto in prolina
    greta_ffnb = merge_GRETA(greta_LJ)

    print('\n GRETA - Pairs and Exclusion section preparation')
    topology_pairs, topology_exclusion = make_pairs_exclusion_topology(type_c12_dict, proline_n, greta_ffnb)
    write_greta_topology_pairs(topology_pairs, topology_exclusion)
    print('\n GRETA - Pairs and Exclusion section written')

if N_terminal == True:
    print('Removing N_1 N_1 pair')
    greta_ffnb.loc[(greta_ffnb['; ai'] == first_resid) & (greta_ffnb['aj'] == first_resid), '; ai'] = ';'+greta_ffnb['; ai'] 

write_greta_LJ(ffnonbonded_atp, greta_ffnb)

print('\n GRETA - FF Written. Change the masses and copy ffnonbonded.itp and atomtypes.atp into the ff folder.')

print('GRETA complete! Carlo is happy')
