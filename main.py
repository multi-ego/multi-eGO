from protein_configuration import protein
from read_input import read_pep_atoms, read_fib_atoms, read_gro_atoms, read_pep_dihedrals, read_fib_dihedrals, read_pep_pairs, read_fib_pairs, read_pdbs, read_gro_bonds, read_gro_angles, read_gro_dihedrals, read_gro_impropers
from functions import make_atomtypes_and_dict, smog_to_gromos_dihedrals, ffnonbonded_merge_pairs, gromos_topology
from write_output import write_atomtypes_atp, write_gromos_topology, write_smog_to_gromos_dihedrals, write_merge_ffnonbonded, write_acid_ffnonbonded, write_greta_LJ, write_greta_atomtypes_atp
from GRETA2 import make_pairs, make_exclusion_list, merge_GRETA, make_pdb_atomtypes


    # Making a dictionary out of it to change the atomnumber to the atomtype
print('')
print('')
print('Creation of peptide and fibril dictionaries')
print('')
    
    # This one is based on the fibril atoms. atp and atomtypes are computed twice but they're the same for peptide and fibril
atp, atomtypes, dict_fib_atomtypes, dict_fib_aminores, fib_smog_to_gro_dict = make_atomtypes_and_dict(read_fib_atoms())

        # dict_fib_atomtypes is used for pairs preparation in FFnonbonded.itp
        # dict_fib_aminores is used for dihedrals
        # pep_smog_to_gro_dict is used for

# IT IS BETTER TO START WITH THE FIBRIL SINCE SOME VARIABLES ARE REPLACED WITH THE MORE COMPLETE INFORMATION FROM THE NATIVE

# The dictionary is based on the peptide atoms it is necessary to create a dictionary used for the peptide part of FFnonbonded.itp
atp, atomtypes, dict_pep_atomtypes, dict_pep_aminores, pep_smog_to_gro_dict = make_atomtypes_and_dict(read_pep_atoms())
    
        # atp is used for the atomtypes.atp
        # atomtypes is used for FFnonbonded.itp
        # dict_pep_atomtypes is used for pairs preparation in FFnonbonded.itp
        # dict_pep_aminores is used for dihedrals
        # pep_smog_to_gro_dict is used for

print('Dictionaries created')
print('')

    # Atomtype.atp

print('Writing atomtype.atp')
print('')

write_atomtypes_atp(atp)

    # Topology section

print('Writing the topology atomtypes section to paste into GROMOS.top')
print('')

gromos_top = gromos_topology(read_gro_atoms())
write_gromos_topology(gromos_top)

print('Writing the proper dihedrals from SMOG to GROMOS')

propers_to_gro = smog_to_gromos_dihedrals(read_pep_dihedrals(), read_fib_dihedrals(), pep_smog_to_gro_dict)
write_smog_to_gromos_dihedrals(propers_to_gro)

print('SMOG to GROMOS topology files ready!')
print('')

    # FFnonbonded section

print('Merge ffnonbonded.itp preparation')

merge_pairs, acid_pairs = ffnonbonded_merge_pairs(read_pep_pairs(), read_fib_pairs(), dict_pep_atomtypes, dict_fib_atomtypes)
write_merge_ffnonbonded(atomtypes, merge_pairs)
write_acid_ffnonbonded(atomtypes, acid_pairs)

print('Merge ffnonbonded.itp created for both neutral and acidic pH')

print(' REMEMBER TO CHANGE THE MASSES IN THE ATOMTYPES.ATP AND FFNONBONDED.ITP, THE H ARE EXPLICIT')





print('\n\n\n\n GRETA TEST \n\n\n\n')
print('GRETA - PDB reading')
native_pdb, fibril_pdb = read_pdbs()

native_atomtypes, fibril_atomtypes, ffnonbonded_atp, atomtypes_atp = make_pdb_atomtypes(native_pdb, fibril_pdb, read_gro_atoms())
print('\n GRETA - Making the exclusion list from bonded')
exclusion_list = make_exclusion_list(native_pdb, read_gro_bonds(), read_gro_angles(), read_gro_dihedrals(), read_gro_impropers())

print('\n GRETA - Making native and fibril pairs')
native_pdb_pairs = make_pairs(native_pdb, exclusion_list, native_atomtypes)
fibril_pdb_pairs = make_pairs(fibril_pdb, exclusion_list, fibril_atomtypes)

print('\n GRETA - native and fibril pairs creation completed')
print('\n GRETA - Merging native and fibril pairs')
greta_ffnb = merge_GRETA(native_pdb_pairs, fibril_pdb_pairs)

#print(ffnonbonded_atp)

write_greta_atomtypes_atp(atomtypes_atp)
write_greta_LJ(ffnonbonded_atp, greta_ffnb)

print('\n GRETA - FF Written. Change the masses and copy ffnonbonded.itp and atomtypes.atp into the ff folder.')

#print(non_bonded.to_string())


# CHECK GRETA WITH SMOG
#print(len(check_SMOG))
#print(len(check_GRETA))

#smog_reverse = check_SMOG[['aj', 'ai']].copy()
#smog_reverse.columns = ['ai', 'aj']
#check_SMOG = check_SMOG.append(smog_reverse, sort = False, ignore_index = True)
#set_check_SMOG = set(check_SMOG['ai'] + '_' + check_SMOG['aj'])

#greta_reverse = check_GRETA[['aj', 'ai']].copy()
#greta_reverse.columns = ['ai', 'aj']
#check_GRETA = check_GRETA.append(greta_reverse, sort = False, ignore_index = True)
#set_check_GRETA = set(check_GRETA['ai'] + '_' + check_GRETA['aj'])

#print(list(set_check_GRETA))
#print(list(set_check_SMOG))

#for g in set_check_GRETA:
#    print('GRETA\t', g)

#for s in set_check_SMOG:
#    print('SMOG\t', s)

#additional_PDB_pairs = set_check_GRETA - set_check_SMOG
#set_SMOG = set(check_SMOG)
#set_GRETA = set(check_GRETA)
#additional_PDB_pairs = [x for x in check_SMOG if x not in set_GRETA]
#print(list(additional_PDB_pairs))
#for a in additional_PDB_pairs:
#    print('diff\t', a)
    
#print(len(additional_PDB_pairs))


# Exclusion list solo se gli atomi sono nella stessa catena
# aggiungi i TER al pdb della fibrilla
#print(native_pdb_pairs)

#native_ex = make_exclusion_list(native_pdb, read_gro_bonds(), read_gro_angles(), read_gro_dihedrals())
#fibril_ex = make_exclusion_list(fibril_pdb, read_gro_bonds(), read_gro_angles(), read_gro_dihedrals())

#print(len(native_ex))
#print(len(fibril_ex))

#if native_ex == fibril_ex:
#    print('same') # SAAAAAMEEEE!!!!
#else:
#    print('diverse')