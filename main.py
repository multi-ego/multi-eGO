from os import getresgid
from read_input import read_pdbs
from write_output import write_greta_LJ, write_greta_atomtypes_atp, write_greta_topology_atoms, write_greta_topology_pairs
from greta import make_pairs_exclusion_topology, make_pairs, merge_GRETA, make_pdb_atomtypes, make_idp_epsilon
from topology_definitions import acid_atp, first_resid
from protein_configuration import idp, N_terminal, greta_to_keep, acid_ff
from mdmat import mdmat_plainMD, mdmat_random_coil



print('\n\n\n\n GRETA\n\n\n\n')
print('GRETA - PDB reading')


# TODO togli la lettura della fibrilla nel caso fibrilla -> *ARG
native_pdb, fibril_pdb = read_pdbs()

print('GRETA - Making Atomtypes')
native_atomtypes, fibril_atomtypes, ffnonbonded_atp, atomtypes_atp, topology_atoms, type_c12_dict, proline_n = make_pdb_atomtypes(native_pdb, fibril_pdb)
native_atomtypes = (sorted(native_atomtypes))

print('\n GRETA - Making native and fibril pairs')

if greta_to_keep == 'native':
    if idp == True:
        #greta_LJ = make_idp_epsilon()
        greta_LJ = make_idp_epsilon(mdmat_plainMD, mdmat_random_coil)
        check = set(greta_LJ['ai'].to_list() + greta_LJ['aj'].to_list())
        check = (sorted(check))
    else:
        greta_LJ = make_pairs(native_pdb, mdmat_random_coil, native_atomtypes)
        if acid_ff == True and acid_atp !=0:
                greta_LJ = greta_LJ[~greta_LJ.ai.isin(acid_atp)]
                greta_LJ = greta_LJ[~greta_LJ.aj.isin(acid_atp)]


elif greta_to_keep == 'fibril':
    greta_LJ = make_pairs(fibril_pdb, fibril_atomtypes)

elif greta_to_keep == 'all':
    if idp == True:
        # Contacts are from a plain MD, so at this step we just import the fibril contacts.
        greta_LJ = make_idp_epsilon().append(make_pairs(fibril_pdb, fibril_atomtypes), sort = False, ignore_index = True)
    else:
        # Merging native and fibril LJ pairs.
        greta_LJ = make_pairs(native_pdb, native_atomtypes)
        if acid_ff == True and acid_atp !=0:
                greta_LJ = greta_LJ[~greta_LJ.ai.isin(acid_atp)]
                greta_LJ = greta_LJ[~greta_LJ.aj.isin(acid_atp)]
        greta_LJ.append(make_pairs(fibril_pdb, fibril_atomtypes), sort = False, ignore_index = True)

else:
    print('ERRORONE')
    exit()


#if native_atomtypes == check:
#    print('same')
#else:
#    print('Check names')
#    #exit()

    

print('\n GRETA - native and fibril pairs creation completed')
print('\n GRETA - Merging native and fibril pairs')

# TODO correggi la massa dell'azoto in prolina
greta_ffnb = merge_GRETA(greta_LJ)
print('\n GRETA - Writing outputs')
write_greta_atomtypes_atp(atomtypes_atp)
write_greta_topology_atoms(topology_atoms)

print('\n GRETA - Pairs and Exclusion section preparation')
topology_pairs, topology_exclusion = make_pairs_exclusion_topology(greta_ffnb, type_c12_dict, proline_n)
write_greta_topology_pairs(topology_pairs, topology_exclusion)
print('\n GRETA - Pairs and Exclusion section written')


if N_terminal == True:
    print('Removing N_1 N_1 pair')
    greta_ffnb.loc[(greta_ffnb['; ai'] == first_resid) & (greta_ffnb['aj'] == first_resid), '; ai'] = ';'+greta_ffnb['; ai'] 

write_greta_LJ(ffnonbonded_atp, greta_ffnb)


print('\n GRETA - FF Written. Change the masses and copy ffnonbonded.itp and atomtypes.atp into the ff folder.')


## Analysis outputs
#print('Writing Force Field Analysis pairs')
#ff_name = 'native'
#if idp == False:
#    write_pairs_list(native_pdb_pairs, ff_name)
#else:
#    write_pairs_list(native_new, ff_name)
#
#ff_name = 'fibril'
#write_pairs_list(fibril_pdb_pairs, ff_name)
#ff_name = 'merge'
#write_pairs_list(greta_ffnb, ff_name)




print('GRETA complete! Carlo is happy')