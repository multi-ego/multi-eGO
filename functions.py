import numpy as np
import pandas as pd
from atomtypes_definitions import gromos_res_atom_dict, gromos_atp, gromos_mass_dict, gromos_resatom_nmr_dict


    # This script includes all the functions used to create a FF

    # Topology section

def gromos_topology(gro_atoms):
    # This function prepares the atomtypes section of gromos topology
    # It will be pasted into the topology along with the dihedrals
    gro_atoms['type'] = gro_atoms['atom_nmr']
    gro_atoms = gro_atoms.drop(['atom_nmr', 'res_atom'], axis=1)
    gro_atoms['mass'] = ''
    gro_atoms['charge'] = ''
    gro_atoms['typeB'] = ''
    gro_atoms['chargeB'] = ''
    gro_atoms['massB'] = ''
    return gro_atoms


def make_atomtypes_and_dict(atomtypes):  # qui si mette l'output di read_*_atoms
    # This function prepare the file for ffnonbonded of the peptide
    # Creation of the atomtypes dictionary    
    # The atomtypes dictionary is used for the FFnonbonded.itp
    dict_atomtypes = atomtypes.set_index("; nr")["type"].to_dict()
    # Handling the information from the topology atomtypes
    # As well atomtypes here is necessary for the creation of FFnonbonded.itp    
    atomtypes['at.group'] = atomtypes['residue'] + '_' + atomtypes['atom']
    
    #print(atomtypes.to_string())
    
    atomtypes['smog_to_gro'] = atomtypes['at.group'] + '_' + atomtypes['resnr'].astype(str)
    smog_to_gro_dict = atomtypes.set_index('; nr')['smog_to_gro'].to_dict()
    # Creation of a dictionary which associates the atom number to the aminoacid and the atom type
    dict_aminores = atomtypes.set_index('; nr')['at.group'].to_dict()
    # Addition of the information from gromos FF (gromos_atp from atomtypes_aa_definitions.py)
    atomtypes['at.group'].replace(gromos_res_atom_dict, inplace = True)
    atomtypes.insert(3, 'at.num', 4)
    atomtypes['at.num'] = atomtypes['at.group'].map(gromos_atp['at.num']) # QUI AD ESEMPIO SI POTREBBE UNIRE CON GROMOS_MASS
    atomtypes.insert(4, 'mass', 5)
    atomtypes['mass'] = atomtypes['at.group'].map(gromos_mass_dict)
    atomtypes["charge"] = '0.000000'
    atomtypes.insert(9, 'ptype', 10)
    atomtypes["ptype"] = 'A'
    atomtypes['c6'] = '0.00000e+00'
    
    print(atomtypes)
    #print(atomtypes.to_string())    
    
    atomtypes['c12'] = atomtypes['at.group'].map(gromos_atp['c12'])
    # Handling the scientific notation
    c12_notation = atomtypes["c12"].map(lambda x:'{:.6e}'.format(x))
    # Removing all the unnecessary columns or duplicated columns in order to prepare the atomtypes for ffnonbonded.itp
    atomtypes = atomtypes.assign(c12 = c12_notation)
    atomtypes.drop(columns = ['; nr', 'resnr', 'residue', 'atom', 'cgnr', 'at.group', 'smog_to_gro'], inplace = True)
    atomtypes.rename(columns = {'type':'; type'}, inplace = True)
    # Since this function is made also for fibrils, a drop duplicate is required, but does not affect the peptide FF
    atomtypes = atomtypes.drop_duplicates(subset = '; type', keep = 'first')
    
    print(atomtypes)
    #print(atomtypes.to_string())
    
    # This last function creates the atomtype for atomtypes.atp
    atp = pd.DataFrame(atomtypes, columns = ['; type', 'mass'])
    return atp, atomtypes, dict_atomtypes, dict_aminores, smog_to_gro_dict


def smog_to_gromos_dihedrals(pep_dihedrals, fib_dihedrals, smog_to_gro_dict): # similar from ffbonded_merge_dihedrals
    # Selection of proper dihedrals from peptide and fibril to create a merged dihedrals,
    # the one which SMOG creates with specific values and to be pasted into Gromacs
    pep_dihedrals = pep_dihedrals.loc[pep_dihedrals['func'] == 1]
    fib_dihedrals = fib_dihedrals.loc[fib_dihedrals['func'] == 1]
    proper_dihedrals = pep_dihedrals.append(fib_dihedrals, sort = False, ignore_index = True)
    proper_dihedrals.loc[:, 'Kd'] = proper_dihedrals.loc[:, 'Kd'].divide(2)
    
    # proper_dihedrals['Kd'] = proper_dihedrals['Kd'] * (300 / 70)
    
    # Actually the thing is on merged dihedrals
    proper_dihedrals[";ai"].replace(smog_to_gro_dict, inplace = True)
    proper_dihedrals["aj"].replace(smog_to_gro_dict, inplace = True)
    proper_dihedrals["ak"].replace(smog_to_gro_dict, inplace = True)
    proper_dihedrals["al"].replace(smog_to_gro_dict, inplace = True)
    # This double dictionary was necessary to map properly the atoms
    proper_dihedrals[";ai"].replace(gromos_resatom_nmr_dict, inplace=True)
    proper_dihedrals["aj"].replace(gromos_resatom_nmr_dict, inplace=True)
    proper_dihedrals["ak"].replace(gromos_resatom_nmr_dict, inplace=True)
    proper_dihedrals["al"].replace(gromos_resatom_nmr_dict, inplace=True)
    proper_dihedrals.to_string(index = False)
    phi0_notation = proper_dihedrals["phi0"].map(lambda x:'{:.9e}'.format(x))
    kd_notation = proper_dihedrals["Kd"].map(lambda x:'{:.9e}'.format(x))
    proper_dihedrals = proper_dihedrals.assign(phi0 = phi0_notation)
    proper_dihedrals = proper_dihedrals.assign(Kd = kd_notation)
    proper_dihedrals["func"] = proper_dihedrals["func"].replace(1, 9)
    proper_dihedrals.columns = ["; ai", "aj", "ak", "al", "func", "phi", "kd", "mult"]
    return proper_dihedrals

    
    # FFnonbonded section


def ffnonbonded_merge_pairs(pep_pairs, fib_pairs, dict_pep_atomtypes, dict_fib_atomtypes):
    # pep_pairs = inp_pep_pairs.copy()
    # fib_pairs = inp_fib_pairs.copy()
    # This script allow to merge the pairs of peptide and fibril.
    # The main difference between the other two pairs function is that peptide C6 and C12 are reweighted.
    # This is because SMOG normalize the LJ potential based on the total number of contacts.
    # Since the peptide has less contacts, the LJ potential is stronger than the fibril
    # Peptide input handling
    pep_pairs[";ai"].replace(dict_pep_atomtypes, inplace = True)
    pep_pairs["aj"].replace(dict_pep_atomtypes, inplace = True)
    pep_pairs.to_string(index = False)
    pep_pairs.columns = ["ai", "aj", "type", "A", "B"]

    # pep_pairs['A'] = pep_pairs['A'] * (300 / 70)
    # pep_pairs['B'] = pep_pairs['B'] * (300 / 70)

    # Fibril input handling
    fib_pairs[';ai'].replace(dict_fib_atomtypes, inplace = True)
    fib_pairs["aj"].replace(dict_fib_atomtypes, inplace = True)
    fib_pairs.to_string(index = False)
    fib_pairs.columns = ["ai", "aj", "type", "A", "B"]

    # fib_pairs['A'] = fib_pairs['A'] * (300 / 70)
    # fib_pairs['B'] = fib_pairs['B'] * (300 / 70)

    # Calcolo di epsilon per peptide e fibrilla
    pep_epsilon = (pep_pairs['A'] ** 2) / (4 * (pep_pairs['B']))
    fib_epsilon = (fib_pairs['A'] ** 2) / (4 * (fib_pairs['B']))
    ratio = pep_epsilon[0] / fib_epsilon[0]
    # QUESTO PRINT CI PIACE MOLTO MA E' SOLO PER PYTHON 3
    print(f'\n'
          f'\tPeptide epsilon: {pep_epsilon[0]}\n'
          f'\tFibril epsilon: {fib_epsilon[0]}\n'
          f'\tRatio: {ratio}'
          f'\n')

    # Reweight peptide LJ
    pep_pairs['A'] = pep_pairs['A'] / ratio
    pep_pairs['B'] = pep_pairs['B'] / ratio

    # From now the function behaves like the others
    A_notation = pep_pairs["A"].map(lambda x:'{:.9e}'.format(x))
    B_notation = pep_pairs["B"].map(lambda x:'{:.9e}'.format(x))
    pep_pairs = pep_pairs.assign(A = A_notation)
    pep_pairs = pep_pairs.assign(B = B_notation)
    A_notation = fib_pairs["A"].map(lambda x:'{:.9e}'.format(x))
    B_notation = fib_pairs["B"].map(lambda x:'{:.9e}'.format(x))
    fib_pairs = fib_pairs.assign(A = A_notation)
    fib_pairs = fib_pairs.assign(B = B_notation)

    # One last step about merging the pairs
    pairs = pep_pairs.append(fib_pairs, sort = False, ignore_index = True)

    # Cleaning the duplicates (the logic has already been explained above
    inv_pairs = pairs[['aj', 'ai', 'type', 'A', 'B']].copy()
    inv_pairs.columns = ['ai', 'aj', 'type', 'A', 'B']
    pairs_full = pairs.append(inv_pairs, sort = False, ignore_index = True)
    n_ai = pairs_full.ai.str.extract('(\d+)')
    n_aj = pairs_full.aj.str.extract('(\d+)')
    pairs_full['n_ai'] = n_ai
    pairs_full['n_aj'] = n_aj
    pairs_full['cond'] = np.where((pairs_full['n_ai'] >= pairs_full['n_aj']), pairs_full['ai'], np.nan)
    pairs_full = pairs_full.dropna()
    pairs_full = pairs_full.drop(['cond', 'n_ai', 'n_aj'], axis = 1)
    # Sorting the pairs
    pairs_full.sort_values(by = ['ai', 'aj', 'A'], inplace = True)
    # Duplicates removal
    # Merging columns in order to drop the duplicates
    # pairs_full['ai'] = pairs_full['ai'].apply(str) + ':' + pairs_full['aj'].apply(str) # questo come funzionava prima
    # Cleaning the remaining duplicates
    pairs_full = pairs_full.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    # Column separation
    ai_aj = pairs_full['ai'].str.split(":", n = 1, expand = True)
    # QUESTI DUE COMANDI SONO DOPO IL SET COPY WARNING
    pairs_full.loc[:, 'ai'] = ai_aj.loc[:, 0]

    # QUESTI DUE COMANDI SONO QUELLI PRIMA DEL COPY WARNING
    # clean_pairs['ai'] = ai_aj[0]
    # clean_pairs['aj'] = ai_aj[1]
    pairs_full.columns = [';ai', 'aj', 'type', 'A', 'B']
    return pairs_full