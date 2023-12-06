import pandas as pd
import numpy as np


def get_bonds(topology):
    """
    Generate bond information DataFrame from the provided topology.

    Args:
    topology: List of bonds in the molecular topology.

    Returns:
    bonds_dataframe: DataFrame containing bond-related information such as atom indices, bond function,
                     equilibrium bond length (req), and force constant (k).
    """
    bonds_dataframe = pd.DataFrame({
        'ai': [bonds.atom1.idx + 1 for bonds in topology],
        'aj': [bonds.atom2.idx + 1 for bonds in topology],
        'funct': [bonds.funct for bonds in topology],
        'req': [bonds.type.req for bonds in topology],
        'k': [bonds.type.k for bonds in topology]
    })
    # Conversion from KCal/mol/A^2 to KJ/mol/nm^2 and from Amber to Gromos
    bonds_dataframe['req'] = bonds_dataframe['req']/10.
    bonds_dataframe['k'] = bonds_dataframe['k']*4.184*100*2
    bonds_dataframe['k'] = bonds_dataframe['k'].map(lambda x: '{:.6e}'.format(x))
    return bonds_dataframe


def get_bond_pairs(topology):
    """
    Generate bond pairs as a list of tuples from the provided topology.

    Args:
    topology: List of bonds in the molecular topology.

    Returns:
    bond_tuple: List of tuples containing pairs of atom indices representing the bonds.
    """
    ai, aj = [], []
    for bonds in topology:
        ai.append(bonds.atom1.idx + 1)
        aj.append(bonds.atom2.idx + 1)
    bond_tuple = list([(str(ai), str(aj)) for ai, aj in zip(ai, aj)])
    return bond_tuple


def get_angles(topology):
    angles_dataframe = pd.DataFrame({
        'ai': [angle.atom1.idx + 1 for angle in topology],
        'aj': [angle.atom2.idx + 1 for angle in topology],
        'ak': [angle.atom3.idx + 1 for angle in topology],
        'funct': [angle.funct for angle in topology],
        'theteq': [angle.type.theteq for angle in topology],
        'k': [angle.type.k for angle in topology]
    })
    angles_dataframe['k'] = angles_dataframe['k']*4.184*2
    angles_dataframe['k'] = angles_dataframe['k'].map(lambda x: '{:.6e}'.format(x))
    return angles_dataframe


def get_dihedrals(topology):
    """
    Extracts dihedral angles information from a molecular topology.

    Args:
    - topology (list): List of dihedral atoms information.

    Returns:
    - dihedrals_dataframe (pandas.DataFrame): DataFrame containing dihedral angles data, including atom indices,
      function type, phase, phi_k, and periodicity.
    """
    dihedrals_dataframe = pd.DataFrame({
        'ai': [dihedral.atom1.idx + 1 for dihedral in topology],
        'aj': [dihedral.atom2.idx + 1 for dihedral in topology],
        'ak': [dihedral.atom3.idx + 1 for dihedral in topology],
        'al': [dihedral.atom4.idx + 1 for dihedral in topology],
        'funct': [dihedral.funct for dihedral in topology],
        'phase': [dihedral.type.phase for dihedral in topology],
        'phi_k': [dihedral.type.phi_k for dihedral in topology],
        'per': [dihedral.type.per for dihedral in topology]
    })
    dihedrals_dataframe['phi_k'] = dihedrals_dataframe['phi_k']*4.184
    return dihedrals_dataframe


def get_impropers(topology):
    """
    Extracts improper torsions information from a molecular topology.

    Args:
    - topology (list): List of improper torsion atoms information.

    Returns:
    - impropers_dataframe (pandas.DataFrame): DataFrame containing improper torsion data, including atom indices,
      function type, psi_eq, and psi_k.
    """
    impropers_dataframe = pd.DataFrame({
        'ai': [improper.atom1.idx + 1 for improper in topology],
        'aj': [improper.atom2.idx + 1 for improper in topology],
        'ak': [improper.atom3.idx + 1 for improper in topology],
        'al': [improper.atom4.idx + 1 for improper in topology],
        'funct': [improper.funct for improper in topology],
        'psi_eq': [improper.type.psi_eq for improper in topology],
        'psi_k': [improper.type.psi_k for improper in topology]
    })
    impropers_dataframe['psi_k'] = impropers_dataframe['psi_k']*4.184*2
    return impropers_dataframe


def get_pairs(topology):
    """
    Extracts pair information from a molecular topology.

    Args:
    - topology (list): List of pair atoms information.

    Returns:
    - pairs_dataframe (pandas.DataFrame): DataFrame containing pair data, including atom indices, function type, and pair type.
    """
    pairs_dataframe = pd.DataFrame({
        'ai': [pair.atom1.idx + 1 for pair in topology],
        'aj': [pair.atom2.idx + 1 for pair in topology],
        'funct': [pair.funct for pair in topology],
        'type': [pair.type for pair in topology],
    })
    return pairs_dataframe


def get_14_interaction_list(reduced_topology, bond_pair):
    """
    Creates lists containing a atoms involved in 1-4 interactions.

    Parameters
    ----------
    reduced_topology: pandas.DataFrame
        function
    bond_bair: i don't know man
        function

    Returns
    -------
    exclusion_bonds: list
        Contains interaction with distance up to 3 bonds
    p14: list
        Contains interactions with distance exactly 3 bonds
    """
    # Building the exclusion bonded list
    # exclusion_bonds are all the interactions within 3 bonds
    # p14 are specifically the interactions at exactly 3 bonds
    ex, ex14, p14, exclusion_bonds = [], [], [], []
    for atom in reduced_topology['number'].to_list():
        for t in bond_pair:
            if t[0] == atom:
                first = t[1]
                ex.append(t[1])
            elif t[1] == atom:
                first = t[0]
                ex.append(t[0])
            else:
                continue
            for tt in bond_pair:
                if (tt[0] == first) & (tt[1] != atom):
                    second = tt[1]
                    ex.append(tt[1])
                elif (tt[1] == first) & (tt[0] != atom):
                    second = tt[0]
                    ex.append(tt[0])
                else:
                    continue
                for ttt in bond_pair:
                    if (ttt[0] == second) & (ttt[1] != first):
                        ex.append(ttt[1])
                        ex14.append(ttt[1])
                    elif (ttt[1] == second) & (ttt[0] != first):
                        ex.append(ttt[0])
                        ex14.append(ttt[0])

        for e in ex:
            exclusion_bonds.append((str(str(atom) + '_' + str(e))))
            exclusion_bonds.append((str(str(e) + '_' + str(atom))))

        ex = []
        for e in ex14:
            p14.append((str(str(atom) + '_' + str(e))))
            p14.append((str(str(e) + '_' + str(atom))))
        ex14 = []

    return exclusion_bonds, p14


def create_pairs_14_dataframe(atomtype1, atomtype2, c6=0.0, shift=0, prefactor=None, constant=None):
    '''
    Used to create additional or modified, multi-eGO-specific 1-4 (like) interactions. Two sets of atomtypes with
    specific shifts in the residue index can be fed to the function to obtain a new set of 1-4 interaction pairs.

    Parameters
    ----------
    atomtype1: list or list-like
        Contains the first set of atomtypes
    atomtype2: list or list-like
        Contains the second set of atomtypes
    c6: float
        Sets a fixed c6 LJ parameters for the specificied type of interaction (default = 0.0)
    shift: int
        Defines the shift in residue index in which to apply the shift. Positive shifts apply the function
        to the atom of the next residue. Negative shifts apply the function to the atom of the previous residue
    prefactor: float
        Factor which to multiply the c12 with after using the combination rule for LJ parameters
    constant: float
        A constant c12 value to use for LJ c12 parameters

    Returns
    -------
    pairs_14: pd.DataFrame
        A DataFrame containing output containing the additional atom indices and LJ parameters
    '''
    if prefactor is not None and constant is not None:
        raise ValueError("Either prefactor or constant has to be set.")
    if prefactor is None and constant is None:
        raise ValueError("Neither prefactor nor constant has been set.")
    pairs_14_ai, pairs_14_aj, pairs_14_c6, pairs_14_c12 = [], [], [], []

    for index, line_atomtype1 in atomtype1.iterrows():
        line_atomtype2 = atomtype2.loc[(atomtype2['resnum'] == line_atomtype1['resnum']+shift)].squeeze(axis=None)
        if not line_atomtype2.empty:
            pairs_14_ai.append(line_atomtype1['number'])
            pairs_14_aj.append(line_atomtype2['number'])
            pairs_14_c6.append(c6)
            if constant is not None:
                pairs_14_c12.append(constant)
            if prefactor is not None:
                pairs_14_c12.append(prefactor*np.sqrt(line_atomtype1['c12']*line_atomtype2['c12']))

    pairs_14 = pd.DataFrame(columns=['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source'])
    pairs_14['ai'] = pairs_14_ai
    pairs_14['aj'] = pairs_14_aj
    pairs_14['func'] = 1
    pairs_14['c6'] = pairs_14_c6
    pairs_14['c12'] = pairs_14_c12
    pairs_14['source'] = '1-4'
    pairs_14['probability'] = 1.0
    pairs_14['rc_probability'] = 1.0

    return pairs_14


def protein_LJ14(reduced_topology):
    """
    Generates Lennard-Jones 14 (LJ14) pairs specific to protein structure.

    Args:
    - reduced_topology (pd.DataFrame): DataFrame containing reduced topology information.

    Returns:
    - pairs (pd.DataFrame): DataFrame with LJ14 pairs for protein interactions.
    """
    # Here we make a dictionary of the atoms used for local geometry
    first_backbone_nitrogen = reduced_topology.loc[(reduced_topology['name'] == 'N')&(reduced_topology['type'] == 'NL')]
    backbone_nitrogen = reduced_topology.loc[(reduced_topology['name'] == 'N')&(reduced_topology['type'] != 'NL')]
    backbone_carbonyl = reduced_topology.loc[reduced_topology['name'] == 'C']
    backbone_oxygen = reduced_topology.loc[reduced_topology['name']=='O']
    ct_oxygen = reduced_topology.loc[(reduced_topology['name']=='O1')|(reduced_topology['name']=='O2')]
    sidechain_cb = reduced_topology.loc[reduced_topology['name'] == 'CB']
    sidechain_cgs = reduced_topology.loc[(reduced_topology['name'] == 'CG')|(reduced_topology['name'] == 'CG1')|(reduced_topology['name'] == 'CG2')|(reduced_topology['name'] == 'SG')|(reduced_topology['name'] == 'OG')|(reduced_topology['name'] == 'OG1')&(reduced_topology['resname'] != 'PRO')]
    pairs = pd.DataFrame()

    # For backbone carbonyl take the CB of the next residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=backbone_carbonyl, atomtype2=sidechain_cb, prefactor=0.275, shift=+1)], axis=0, sort=False, ignore_index=True)
    # For backbone oxygen take the CB of the same residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=backbone_oxygen, atomtype2=sidechain_cb, prefactor=0.1)], axis=0, sort=False, ignore_index=True)
    # now we add the pair between the last CB and the two OCT ones
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=ct_oxygen, atomtype2=sidechain_cb, prefactor=0.1)], axis=0, sort=False, ignore_index=True)
    # For each backbone nitrogen take the CB of the previuos residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=backbone_nitrogen, atomtype2=sidechain_cb, prefactor=0.65, shift=-1)], axis=0, sort=False, ignore_index=True)
    # For the first backbone nitrogen take the N of the next residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=first_backbone_nitrogen, atomtype2=backbone_nitrogen, constant=4.e-06, shift=+1)], axis=0, sort=False, ignore_index=True)
    # For each backbone nitrogen take the N of the next residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=backbone_nitrogen, atomtype2=backbone_nitrogen, prefactor=0.343, shift=+1)], axis=0, sort=False, ignore_index=True)
    # For each backbone carbonyl take the carbonyl of the next residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=backbone_carbonyl, atomtype2=backbone_carbonyl, prefactor=0.5, shift=-1)], axis=0, sort=False, ignore_index=True)
    # For each backbone carbonyl take the CGs of the same residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=sidechain_cgs, atomtype2=backbone_carbonyl, prefactor=0.078)], axis=0, sort=False, ignore_index=True)
    # For each backbone nitrogen take the CGs of the same residue and save in a pairs tuple
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=sidechain_cgs, atomtype2=backbone_nitrogen, prefactor=0.087)], axis=0, sort=False, ignore_index=True)
    pairs = pd.concat([pairs, create_pairs_14_dataframe(atomtype1=sidechain_cgs, atomtype2=first_backbone_nitrogen, prefactor=0.087)], axis=0, sort=False, ignore_index=True)

    # make it symmetric
    inv_LJ = pairs[['aj', 'ai', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source']].copy()
    inv_LJ.columns = ['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source']

    pairs = pd.concat([pairs, inv_LJ], axis=0, sort=False, ignore_index=True)
    pairs['ai'] = pairs['ai'].astype(str)
    pairs['aj'] = pairs['aj'].astype(str)

    return pairs
