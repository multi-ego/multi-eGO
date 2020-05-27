from atomtypes_definitions import *


# This script includes all the functions used to create a FF


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
    atomtypes['at.num'] = atomtypes['at.group'].map(gromos_atp['at.num'])
    atomtypes.insert(4, 'mass', 5)
    atomtypes['mass'] = atomtypes['at.group'].map(gromos_mass_dict)
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
    proper_dihedrals['Kd'] = proper_dihedrals['Kd'] * (300 / 70)
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
    