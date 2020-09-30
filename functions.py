import numpy as np
import pandas as pd
from atomtypes_definitions import gromos_res_atom_dict, gromos_atp, gromos_mass_dict, gromos_resatom_nmr_dict, acid_atp, resnr_pairs
from protein_configuration import t_ratio, temperatura


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
    atomtypes['smog_to_gro'] = atomtypes['at.group'] + '_' + atomtypes['resnr'].astype(str)
    smog_to_gro_dict = atomtypes.set_index('; nr')['smog_to_gro'].to_dict()

    # Creation of a dictionary which associates the atom number to the aminoacid and the atom type
    dict_aminores = atomtypes.set_index('; nr')['at.group'].to_dict()
    # Addition of the information from gromos FF (gromos_atp from atomtypes_aa_definitions.py)
    atomtypes['at.group'].replace(gromos_res_atom_dict, inplace = True)
    atomtypes.insert(3, 'at.num', 4)
    atomtypes['at.num'] = atomtypes['at.group'].map(gromos_atp['at.num']) # QUI AD ESEMPIO SI POTREBBE UNIRE CON GROMOS_MASS
    atomtypes.insert(4, 'mass', 5)
    atomtypes['mass'] = atomtypes['at.group'].map(gromos_atp['mass'])
    atomtypes["charge"] = '0.000000'
    atomtypes.insert(9, 'ptype', 10)
    atomtypes["ptype"] = 'A'
    atomtypes['c6'] = '0.00000e+00'
    atomtypes['c12'] = atomtypes['at.group'].map(gromos_atp['c12'])
    # Handling the scientific notation
    c12_notation = atomtypes["c12"].map(lambda x:'{:.6e}'.format(x))
    # Removing all the unnecessary columns or duplicated columns in order to prepare the atomtypes for ffnonbonded.itp
    atomtypes = atomtypes.assign(c12 = c12_notation)
    atomtypes.drop(columns = ['; nr', 'resnr', 'residue', 'atom', 'cgnr', 'at.group', 'smog_to_gro'], inplace = True)
    atomtypes.rename(columns = {'type':'; type'}, inplace = True)
    # Since this function is made also for fibrils, a drop duplicate is required, but does not affect the peptide FF
    atomtypes = atomtypes.drop_duplicates(subset = '; type', keep = 'first')
    # Change from float to integer the at.num otherwise gromacs does not understand
    atomtypes['at.num'] = atomtypes['at.num'].fillna(0.0).astype(int)
    # This last function creates the atomtype for atomtypes.atp
    atp = pd.DataFrame(atomtypes, columns = ['; type', 'mass'])
    return atp, atomtypes, dict_atomtypes, dict_aminores, smog_to_gro_dict


def smog_to_gromos_dihedrals(pep_dihedrals, fib_dihedrals, smog_to_gro_dict): # similar from ffbonded_merge_dihedrals
    # Selection of proper dihedrals from peptide and fibril to create a merged dihedrals,
    # the one which SMOG creates with specific values and to be pasted into Gromacs
    pep_dihedrals = pep_dihedrals.loc[pep_dihedrals['func'] == 1]
    fib_dihedrals = fib_dihedrals.loc[fib_dihedrals['func'] == 1]

    # The Kds are different between the fibril and the native, therefore here is how to rebalance like the pairs
    # However in the case of dihedrals it is necessary to multiply instead of divide because the way of SMOG
    # computes dihedrals
    
    Kd_pep = pep_dihedrals['Kd']
    Kd_fib = fib_dihedrals['Kd']
    ratio = Kd_pep[0] / Kd_fib[0]
    print(f'\n'
        f'\tDihedral Ratio: {ratio}'
        f'\n')
    pep_dihedrals['Kd'] = pep_dihedrals['Kd'] * ratio
    #pep_dihedrals.loc[:, 'Kd'] = pep_dihedrals.loc[:, 'Kd'].divide(ratio)
    proper_dihedrals = pep_dihedrals.append(fib_dihedrals, sort = False, ignore_index = True)

    # Here all the dihedrals are present (i can read 836 and 715)
    #print(proper_dihedrals.to_string())

    # TUTTI I DIEDRI VENGONO DIVISI PER DUE, CHE SIANO DOPPI (NATIVA E FIBRILLA) O SINGOLI (SOLO NELLA NATIVA)
    proper_dihedrals.loc[:, 'Kd'] = proper_dihedrals.loc[:, 'Kd'].divide(2)
    
    proper_dihedrals['Kd'] = proper_dihedrals['Kd'] * t_ratio

    print(f'\n'
    f'\tTemperature Ratio: {temperatura} / 70'
    f'\n')

    # Actually the thing is on merged dihedrals
    # In this function is necessary to use the native smog_to_gro_dictionary since is the full dictionary
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
    # This script allow to merge the pairs of peptide and fibril.
    # The main difference between the other two pairs function is that peptide C6 and C12 are reweighted.
    # This is because SMOG normalize the LJ potential based on the total number of contacts.
    # Since the peptide has less contacts, the LJ potential is stronger than the fibril
    # Peptide input handling
    pep_pairs[";ai"].replace(dict_pep_atomtypes, inplace = True)
    pep_pairs["aj"].replace(dict_pep_atomtypes, inplace = True)
    pep_pairs.to_string(index = False)
    pep_pairs.columns = ["ai", "aj", "type", "A", "B"]

    pep_pairs['A'] = pep_pairs['A'] * t_ratio
    pep_pairs['B'] = pep_pairs['B'] * t_ratio

    # Fibril input handling
    fib_pairs[';ai'].replace(dict_fib_atomtypes, inplace = True)
    fib_pairs["aj"].replace(dict_fib_atomtypes, inplace = True)
    fib_pairs.to_string(index = False)
    fib_pairs.columns = ["ai", "aj", "type", "A", "B"]

    fib_pairs['A'] = fib_pairs['A'] * t_ratio
    fib_pairs['B'] = fib_pairs['B'] * t_ratio

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

        # If acidic the following pep_pairs will be removed
        # Remove the lines by searching in the two colums 
        # Filter the informations from gromos atomtype and make a list of the atomtypes to remove
        # e.g. OD1_39 -> remove line
        # This step should be done before append the two pairs otherwise some fibril contribution will be lost
        # The list of atoms will be made in atomtypes_definitions.py


    for_acid_pairs = pep_pairs.copy()
    acid_pep_pairs = for_acid_pairs[~for_acid_pairs['ai'].isin(acid_atp)]
    acid_pep_pairs = acid_pep_pairs[~acid_pep_pairs['aj'].isin(acid_atp)]

    # One last step about merging the pairs for both neutral and acid pH
    pairs = pep_pairs.append(fib_pairs, sort = False, ignore_index = True)
    acid_pairs = acid_pep_pairs.append(fib_pairs, sort = False, ignore_index = True)

    # Cleaning the duplicates (the logic has already been explained above)
    inv_pairs = pairs[['aj', 'ai', 'type', 'A', 'B']].copy()
    inv_pairs.columns = ['ai', 'aj', 'type', 'A', 'B']

    inv_acid = acid_pairs[['aj', 'ai', 'type', 'A', 'B']].copy()
    inv_acid.columns = ['ai', 'aj', 'type', 'A', 'B']

    pairs_full = pairs.append(inv_pairs, sort = False, ignore_index = True)
    acid_full = acid_pairs.append(inv_acid, sort = False, ignore_index = True)

    # Sorting the pairs
    pairs_full.sort_values(by = ['ai', 'aj', 'A'], inplace = True)
    # Cleaning the duplicates
    pairs_full = pairs_full.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')

    # Removing the reverse duplicates
    cols = ['ai', 'aj']
    pairs_full[cols] = np.sort(pairs_full[cols].values, axis=1)
    pairs_full = pairs_full.drop_duplicates()

    ### ACID FF
    # Sorting the pairs
    acid_full.sort_values(by = ['ai', 'aj', 'A'], inplace = True)
    # Cleaning the duplicates
    acid_full = acid_full.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')

    # Removing the reverse duplicates
    acid_full[cols] = np.sort(acid_full[cols].values, axis=1)
    acid_full = acid_full.drop_duplicates()


    # Here i should add the duplicates for the N terminal. 

    # Crea lista atomtypes da dizionario


    pairs_full.columns = [';ai', 'aj', 'type', 'A', 'B']
    acid_full.columns = [';ai', 'aj', 'type', 'A', 'B']

    atp_values=list(dict_pep_atomtypes.values())
  
    pairs_full['double'] = ''
    
    for i in atp_values:

        # Questo funziona e riesco a fargli dire quello che voglio.
        # Cioe' flaggare solo i valori che hanno un loro corrispettivo: N_1 N_1, CA_1 CA_1 ...
        pairs_full.loc[(pairs_full[';ai'] == i) & (pairs_full['aj'] == i), 'double'] = 'True'

    doubles = pairs_full.loc[(pairs_full['double'] == 'True')]
    print(doubles)

    #print(pairs_full[pairs_full.values == 'True'].to_string())
    print(len(pairs_full[pairs_full.values == 'True']))
    #print(pairs_full[pairs_full.values == 'False'].to_string())
    #print(len(pairs_full[pairs_full.values == 'False']))    
    print(len(pairs_full))

    # Adesso sia in ai che in aj ci sono gli stessi valori, dunque si puo' vedere la differenza con il dataframe completo
    # e sottrarre questi gia' doppi a quella completa
    # restano quelli da raddoppiare

    
    ##### DA IMPLEMENTARE #####
    
    # https://kanoki.org/2019/07/04/pandas-difference-between-two-dataframes/ 
    ########

    pairs_full = pairs_full.drop(['double'], axis=1)
    return pairs_full, acid_full