from read_input import read_pdbs
from write_output import write_greta_LJ, write_greta_atomtypes_atp, write_greta_topology_atoms, write_pairs_list, write_acid_greta_LJ
from greta import make_pairs, merge_GRETA, make_pdb_atomtypes
from topology_definitions import acid_atp


print('\n\n\n\n GRETA\n\n\n\n')
print('GRETA - PDB reading')
native_pdb, fibril_pdb = read_pdbs()

print('GRETA - Making Atomtypes')
native_atomtypes, fibril_atomtypes, ffnonbonded_atp, atomtypes_atp, topology_atoms = make_pdb_atomtypes(native_pdb, fibril_pdb)

print('\n GRETA - Making native and fibril pairs')
# gromologist version, it works same ffnonbonded.itp

# Add idp condition to avoid the useless reading of the native
native_pdb_pairs = make_pairs(native_pdb, native_atomtypes)
fibril_pdb_pairs = make_pairs(fibril_pdb, fibril_atomtypes)

#print(len(native_pdb_pairs))
#print(len(fibril_pdb_pairs))

print('\n GRETA - native and fibril pairs creation completed')
print('\n GRETA - Merging native and fibril pairs')
greta_ffnb = merge_GRETA(native_pdb_pairs, fibril_pdb_pairs)

write_greta_atomtypes_atp(atomtypes_atp)
write_greta_topology_atoms(topology_atoms)
write_greta_LJ(ffnonbonded_atp, greta_ffnb)

if len(acid_atp) != 0:
    print('\n GRETA - Acid pH')
    acid_pdb_pairs = native_pdb_pairs.copy()
    acid_pdb_pairs = acid_pdb_pairs[~acid_pdb_pairs.ai.isin(acid_atp)]
    acid_pdb_pairs = acid_pdb_pairs[~acid_pdb_pairs.aj.isin(acid_atp)]
    greta_acid_ffnb = merge_GRETA(acid_pdb_pairs, fibril_pdb_pairs)
    write_acid_greta_LJ(ffnonbonded_atp, greta_acid_ffnb)

#print(ffnonbonded_atp)
print('\n GRETA - FF Written. Change the masses and copy ffnonbonded.itp and atomtypes.atp into the ff folder.')


# Analysis outputs
print('Writing Force Field Analysis pairs')
ff_name = 'native'
write_pairs_list(native_pdb_pairs, ff_name)
ff_name = 'fibril'
write_pairs_list(fibril_pdb_pairs, ff_name)
ff_name = 'merge'
write_pairs_list(greta_ffnb, ff_name)
print('GRETA complete! Carlo is happy')