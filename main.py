from read_input import read_pdbs
from write_output import write_greta_LJ, write_greta_atomtypes_atp, write_greta_topology_atoms, write_greta_topology_pairs, write_pairs_list, write_acid_greta_LJ
from greta import make_pairs_exclusion_topology, make_pairs, merge_GRETA, make_pdb_atomtypes
from topology_definitions import acid_atp, first_resid
from protein_configuration import idp, N_terminal


print('\n\n\n\n GRETA\n\n\n\n')
print('GRETA - PDB reading')
native_pdb, fibril_pdb = read_pdbs()

print('GRETA - Making Atomtypes')
native_atomtypes, fibril_atomtypes, ffnonbonded_atp, atomtypes_atp, topology_atoms, type_c12_dict = make_pdb_atomtypes(native_pdb, fibril_pdb)

print('\n GRETA - Making native and fibril pairs')

native_pdb_pairs = make_pairs(native_pdb, native_atomtypes)
fibril_pdb_pairs = make_pairs(fibril_pdb, fibril_atomtypes)

print('\n GRETA - native and fibril pairs creation completed')
print('\n GRETA - Merging native and fibril pairs')

# TODO correggi la massa dell'azoto in prolina
#print(native_pdb_pairs)


if idp == False:
    greta_ffnb = merge_GRETA(native_pdb_pairs, fibril_pdb_pairs)
else:
    greta_ffnb, native_new = merge_GRETA(native_pdb_pairs, fibril_pdb_pairs)

print('\n GRETA - Writing outputs')
write_greta_atomtypes_atp(atomtypes_atp)
write_greta_topology_atoms(topology_atoms)

print('\n GRETA - Pairs and Exclusion section preparation')
topology_pairs, topology_exclusion = make_pairs_exclusion_topology(greta_ffnb, type_c12_dict)
write_greta_topology_pairs(topology_pairs, topology_exclusion)
print('\n GRETA - Pairs and Exclusion section written')

if len(acid_atp) != 0:
    print('\n GRETA - Acid pH')
    acid_pdb_pairs = native_pdb_pairs.copy()
    acid_pdb_pairs = acid_pdb_pairs[~acid_pdb_pairs.ai.isin(acid_atp)]
    acid_pdb_pairs = acid_pdb_pairs[~acid_pdb_pairs.aj.isin(acid_atp)]
    greta_acid_ffnb = merge_GRETA(acid_pdb_pairs, fibril_pdb_pairs)
    write_acid_greta_LJ(ffnonbonded_atp, greta_acid_ffnb)


if N_terminal == True:
    print('Removing N_1 N_1 pair')
    greta_ffnb.loc[(greta_ffnb['; ai'] == first_resid) & (greta_ffnb['aj'] == first_resid), '; ai'] = ';'+greta_ffnb['; ai'] 
write_greta_LJ(ffnonbonded_atp, greta_ffnb)


print('\n GRETA - FF Written. Change the masses and copy ffnonbonded.itp and atomtypes.atp into the ff folder.')


# Analysis outputs
print('Writing Force Field Analysis pairs')
ff_name = 'native'
if idp == False:
    write_pairs_list(native_pdb_pairs, ff_name)
else:
    write_pairs_list(native_new, ff_name)

ff_name = 'fibril'
write_pairs_list(fibril_pdb_pairs, ff_name)
ff_name = 'merge'
write_pairs_list(greta_ffnb, ff_name)
print('GRETA complete! Carlo is happy')