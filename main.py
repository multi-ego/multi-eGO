from read_input import read_pep_atoms, read_fib_atoms, read_gro_atoms, read_pep_dihedrals, read_fib_dihedrals, read_pep_pairs, read_fib_pairs
from functions import make_atomtypes_and_dict, smog_to_gromos_dihedrals, ffnonbonded_merge_pairs, gromos_topology
from write_output import *



    # Making a dictionary out of it to change the atomnumber to the atomtype
print('')
print('')
print('Creation of peptide and fibril dictionaries')
print('')
    # The dictionary is based on the peptide atoms it is necessary to create a dictionary used for the peptide part of FFnonbonded.itp
atp, atomtypes, dict_pep_atomtypes, dict_pep_aminores, pep_smog_to_gro_dict = make_atomtypes_and_dict(read_pep_atoms())
    
        # atp is used for the atomtypes.atp
        # atomtypes is used for FFnonbonded.itp
        # dict_pep_atomtypes is used for pairs preparation in FFnonbonded.itp
        # dict_pep_aminores is used for dihedrals
        # pep_smog_to_gro_dict is used for

    # This one is based on the fibril atoms. atp and atomtypes are computed twice but they're the same for peptide and fibril
atp, atomtypes, dict_fib_atomtypes, dict_fib_aminores, fib_smog_to_gro_dict = make_atomtypes_and_dict(read_fib_atoms())

        # dict_fib_atomtypes is used for pairs preparation in FFnonbonded.itp
        # dict_fib_aminores is used for dihedrals
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

propers_to_gro = smog_to_gromos_dihedrals(read_pep_dihedrals(), read_fib_dihedrals(), fib_smog_to_gro_dict)
write_smog_to_gromos_dihedrals(propers_to_gro)

print('SMOG to GROMOS topology files ready!')
print('')

    # FFnonbonded section

print('Merge ffnonbonded.itp preparation')

merge_pairs = ffnonbonded_merge_pairs(read_pep_pairs(), read_fib_pairs(), dict_pep_atomtypes, dict_fib_atomtypes)
write_merge_ffnonbonded(atomtypes, merge_pairs)

print('Merge ffnonbonded.itp created')