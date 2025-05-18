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
    bonds_data = []

    for bond in topology:
        ai = bond.atom1.idx + 1
        aj = bond.atom2.idx + 1
        funct = bond.funct
        if bond.type is not None:
            req = bond.type.req
            k = bond.type.k
        else:
            if bond.atom1.name == "SG" and bond.atom2.name == "SG":
                req = 0.204 * 10.0
                k = 2.5e05 / (4.184 * 100 * 2)
            else:
                req = None
                k = None
                print("WARNING: there is an unparametrized bond in your reference topology: ", bond)

        bonds_data.append({"ai": ai, "aj": aj, "funct": funct, "req": req, "k": k})

    if bonds_data:
        bonds_dataframe = pd.DataFrame(bonds_data)
        # Conversion from KCal/mol/A^2 to KJ/mol/nm^2 and from Amber to Gromos
        bonds_dataframe["req"] = bonds_dataframe["req"] / 10.0
        bonds_dataframe["k"] = bonds_dataframe["k"] * 4.184 * 100 * 2
        bonds_dataframe["k"] = bonds_dataframe["k"].map(lambda x: "{:.6e}".format(x))
    else:
        bonds_dataframe = pd.DataFrame(columns=["ai", "aj", "funct", "req", "k"])

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
    """
    Extracts angle data from a topology object and constructs a pandas DataFrame.

    Parameters:
    - topology (parmed.Structure): The ParmEd topology object containing angle information.

    Returns:
    - pandas.DataFrame: A DataFrame containing angle data, including atom indices (ai, aj, ak),
                        function type (funct), equilibrium angles (theteq), and force constants (k).
    """
    # List to store angle data dictionaries
    angles_data = []

    # Iterate over each angle in the topology
    for angle in topology:
        # Extract indices of atoms involved in the angle
        ai = angle.atom1.idx + 1
        aj = angle.atom2.idx + 1
        ak = angle.atom3.idx + 1

        # Extract angle function type
        funct = angle.funct

        # Check if the angle type is not None
        if angle.type is not None:
            # If angle type is not None, extract parameters theteq and k
            theteq = angle.type.theteq
            k = angle.type.k
        else:
            # If angle type is None, handle the unparametrized angle
            if angle.atom2.name == "SG" and (angle.atom1.name == "SG" or angle.atom3.name == "SG"):
                # Handle special case for cysteine residues
                theteq = 104.0
                k = 4.613345e02 / (4.184 * 2)  # Converting kcal/(mol*rad^2) to kJ/(mol*rad^2)
            else:
                # Set the parameters to None and print a warning
                theteq = None
                k = None
                print("WARNING: there is an unparametrized angle in your reference topology:", angle)

        # Append the angle data to the angles_data list as a dictionary
        angles_data.append({"ai": ai, "aj": aj, "ak": ak, "funct": funct, "theteq": theteq, "k": k})

    if angles_data:
        # Create a pandas DataFrame from the angles_data list
        angles_dataframe = pd.DataFrame(angles_data)

        # Convert k values from kcal/(mol*rad^2) to kJ/(mol*rad^2)
        angles_dataframe["k"] = angles_dataframe["k"] * 4.184 * 2

        # Format k values in scientific notation with 6 decimal places
        angles_dataframe["k"] = angles_dataframe["k"].map(lambda x: "{:.6e}".format(x))
    else:
        angles_dataframe = pd.DataFrame(columns=["ai", "aj", "ak", "funct", "theteq", "k"])
    # Return the DataFrame containing angle data
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
    dihedrals_data = []

    # Iterate over each dihedral in the topology
    for dihedral in topology:
        # Extract indices of atoms involved in the dihedral
        ai = dihedral.atom1.idx + 1
        aj = dihedral.atom2.idx + 1
        ak = dihedral.atom3.idx + 1
        al = dihedral.atom4.idx + 1

        # Extract dihedral function type
        funct = dihedral.funct

        # Check if the dihedral type is not None
        if dihedral.type is not None:
            # If dihedral type is not None, extract parameters theteq and k
            phase = dihedral.type.phase
            phi_k = dihedral.type.phi_k
            per = dihedral.type.per
        else:
            # If dihedral type is None, handle the unparametrized dihedral
            if dihedral.atom2.name == "SG" and dihedral.atom3.name == "SG":
                # Handle special case for cysteine residues
                phase = 0.0
                phi_k = 16.7 / (4.184)  # Converting kcal/(mol*rad^2) to kJ/(mol*rad^2)
                per = 2
            elif dihedral.atom3.name == "SG" and dihedral.atom4.name == "SG":
                # Handle special case for cysteine residues
                phase = 0.0
                phi_k = 2.93 / (4.184)  # Converting kcal/(mol*rad^2) to kJ/(mol*rad^2)
                per = 3
            else:
                # Set the parameters to None and print a warning
                phase = None
                phi_k = None
                per = None
                print("WARNING: there is an unparametrized dihedral in your reference topology:", dihedral)

        # Append the dihedral data to the dihedrals_data list as a dictionary
        dihedrals_data.append(
            {"ai": ai, "aj": aj, "ak": ak, "al": al, "funct": funct, "phase": phase, "phi_k": phi_k, "per": per}
        )
    if dihedrals_data:
        # Create a pandas DataFrame from the dihedrals_data list
        dihedrals_dataframe = pd.DataFrame(dihedrals_data)

        # Convert k values from kcal/(mol*rad^2) to kJ/(mol*rad^2)
        dihedrals_dataframe["phi_k"] = dihedrals_dataframe["phi_k"] * 4.184
    else:
        dihedrals_dataframe = pd.DataFrame(columns=["ai", "aj", "ak", "al", "funct", "phase", "phi_k", "per"])

    # Return the DataFrame containing dihedral data
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
    impropers_data = []

    # Iterate over each improper in the topology
    for improper in topology:
        # Extract indices of atoms involved in the improper
        ai = improper.atom1.idx + 1
        aj = improper.atom2.idx + 1
        ak = improper.atom3.idx + 1
        al = improper.atom4.idx + 1

        # Extract improper function type
        funct = improper.funct

        # Check if the improper type is not None
        if improper.type is not None:
            # If improper type is not None, extract parameters theteq and k
            psi_eq = improper.type.psi_eq
            psi_k = improper.type.psi_k
        else:
            # If improper type is None, handle the unparametrized improper
            # Set the parameters to None and print a warning
            psi_eq = None
            psi_k = None
            print("WARNING: there is an unparametrized improper in your reference topology:", improper)

        # Append the improper data to the impropers_data list as a dictionary
        impropers_data.append({"ai": ai, "aj": aj, "ak": ak, "al": al, "funct": funct, "psi_eq": psi_eq, "psi_k": psi_k})

    if impropers_data:
        # Create a pandas DataFrame from the impropers_data list
        impropers_dataframe = pd.DataFrame(impropers_data)

        # Convert k values from kcal/(mol*rad^2) to kJ/(mol*rad^2)
        impropers_dataframe["psi_k"] = impropers_dataframe["psi_k"] * 4.184 * 2
    else:
        impropers_dataframe = pd.DataFrame(columns=["ai", "aj", "ak", "al", "funct", "psi_eq", "psi_k"])

    # Return the DataFrame containing improper data
    return impropers_dataframe


def get_pairs(topology):
    """
    Extracts pair information from a molecular topology.

    Args:
    - topology (list): List of pair atoms information.

    Returns:
    - pairs_dataframe (pandas.DataFrame): DataFrame containing pair data, including atom indices, function type, and pair type.
    """
    pairs_dataframe = pd.DataFrame(
        {
            "ai": [pair.atom1.idx + 1 for pair in topology],
            "aj": [pair.atom2.idx + 1 for pair in topology],
            "funct": [pair.funct for pair in topology],
            "type": [pair.type for pair in topology],
        }
    )
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
    for atom in reduced_topology["number"].to_list():
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
            exclusion_bonds.append((str(str(atom) + "_" + str(e))))
            exclusion_bonds.append((str(str(e) + "_" + str(atom))))

        ex = []
        for e in ex14:
            p14.append((str(str(atom) + "_" + str(e))))
            p14.append((str(str(e) + "_" + str(atom))))
        ex14 = []

    return exclusion_bonds, p14


def get_lj_params(topology):
    lj_params = pd.DataFrame(columns=["ai", "c6", "c12"], index=np.arange(len(topology.atoms)))
    for i, atom in enumerate(topology.atoms):
        c6, c12 = atom.sigma * 0.1, atom.epsilon * 4.184
        lj_params.loc[i] = [atom.atom_type, c6, c12]

    return lj_params


def get_lj_pairs(topology):
    """
    Extracts Lennard-Jones pair information from a molecular topology.

    Parameters
    ----------
    topology: parmed.topology object
        Contains the molecular topology information

    Returns
    -------
    pairs_dataframe: pd.DataFrame
        DataFrame containing Lennard-Jones pair information
    """
    lj_pairs = pd.DataFrame(columns=["ai", "aj", "epsilon", "sigma"], index=np.arange(len(topology.parameterset.nbfix_types)))
    for i, (sbtype_i, sbtype_j) in enumerate(topology.parameterset.nbfix_types):
        key = (sbtype_i, sbtype_j)
        # This is read as rmin not as sigma --> must be scaled by 1/2**(1/6)
        # Any contact present more then once is overwritten by the last one in the nonbond_params
        c12, c6 = topology.parameterset.nbfix_types[key][0] * 4.184, topology.parameterset.nbfix_types[key][1] * 0.1 / (
            2 ** (1 / 6)
        )
        epsilon = c6**2 / (4 * c12) if c6 > 0 else -c12
        sigma = (c12 / c6) ** (1 / 6) if c6 > 0 else c12 ** (1 / 12) / (2.0 ** (1.0 / 6.0))
        lj_pairs.loc[i] = [sbtype_i, sbtype_j, epsilon, sigma]

    return lj_pairs


def get_lj14_pairs(topology):
    """
    Extracts Lennard-Jones pair information from a molecular topology.

    Parameters
    ----------
    topology: parmed.topology object
        Contains the molecular topology information

    Returns
    -------
    pairs_dataframe: pd.DataFrame
        DataFrame containing Lennard-Jones pair information
    """
    lj14_pairs = pd.DataFrame()
    for mol, top in topology.molecules.items():
        pair14 = [
            {
                "ai": pair.atom1.type,
                "aj": pair.atom2.type,
                "c6": pair.type.sigma * 0.1,
                "c12": pair.type.epsilon * 4.184,
            }
            for pair in top[0].adjusts
        ]
        df = pd.DataFrame(pair14)
        lj14_pairs = pd.concat([lj14_pairs, df])

    lj14_pairs = lj14_pairs.reset_index()

    # Calculate "epsilon" using a vectorized conditional expression
    lj14_pairs["epsilon"] = np.where(lj14_pairs["c6"] > 0, lj14_pairs["c6"] ** 2 / (4 * lj14_pairs["c12"]), -lj14_pairs["c12"])

    # Calculate "sigma" using a vectorized conditional expression
    lj14_pairs["sigma"] = np.where(
        lj14_pairs["c6"] > 0,
        (lj14_pairs["c12"] / lj14_pairs["c6"]) ** (1 / 6),
        lj14_pairs["c12"] ** (1 / 12) / (2.0 ** (1.0 / 6.0)),
    )

    lj14_pairs.drop(columns=["c6", "c12"], inplace=True)

    return lj14_pairs


def create_pairs_14_dataframe(atomtype1, atomtype2, c6=0.0, shift=0, prefactor=None, constant=None):
    """
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
    """
    # if prefactor is not None and constant is not None:
    #    raise ValueError("Either prefactor or constant has to be set.")
    if prefactor is None and constant is None:
        raise ValueError("Neither prefactor nor constant has been set.")
    pairs_14_ai, pairs_14_aj, pairs_14_c6, pairs_14_c12 = [], [], [], []

    for index, line_atomtype1 in atomtype1.iterrows():
        line_atomtype2 = atomtype2.loc[(atomtype2["resnum"] == line_atomtype1["resnum"] + shift)].squeeze(axis=None)
        if not line_atomtype2.empty:
            pairs_14_ai.append(line_atomtype1["number"])
            pairs_14_aj.append(line_atomtype2["number"])
            pairs_14_c6.append(c6)
            c12 = 1
            if constant is not None:
                c12 = constant
            if prefactor is not None:
                mixed_c12 = prefactor * np.sqrt(line_atomtype1["c12"] * line_atomtype2["c12"])
                c12 = min(c12, mixed_c12)

            pairs_14_c12.append(c12)

    pairs_14 = pd.DataFrame(
        columns=[
            "ai",
            "aj",
            "func",
            "c6",
            "c12",
            "probability",
            "rc_probability",
            "source",
        ]
    )
    pairs_14["ai"] = pairs_14_ai
    pairs_14["aj"] = pairs_14_aj
    pairs_14["func"] = 1
    pairs_14["c6"] = pairs_14_c6
    pairs_14["c12"] = pairs_14_c12
    pairs_14["source"] = "1-4"
    pairs_14["probability"] = 1.0
    pairs_14["rc_probability"] = 1.0

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
    first_backbone_nitrogen = reduced_topology.loc[(reduced_topology["name"] == "N") & (reduced_topology["type"] == "NL")]
    backbone_nitrogen = reduced_topology.loc[(reduced_topology["name"] == "N") & (reduced_topology["type"] != "NL")]
    backbone_carbonyl = reduced_topology.loc[reduced_topology["name"] == "C"]
    backbone_calpha = reduced_topology.loc[reduced_topology["name"] == "CA"]
    backbone_oxygen = reduced_topology.loc[reduced_topology["name"] == "O"]
    ct_oxygen = reduced_topology.loc[(reduced_topology["name"] == "O1") | (reduced_topology["name"] == "O2")]
    sidechain_cb = reduced_topology.loc[reduced_topology["name"] == "CB"]
    sidechain_cgs = reduced_topology.loc[
        (reduced_topology["name"] == "CG")
        | (reduced_topology["name"] == "CG1")
        | (reduced_topology["name"] == "CG2")
        | (reduced_topology["name"] == "SG")
        | (reduced_topology["name"] == "OG")
        | (reduced_topology["name"] == "OG1") & (reduced_topology["resname"] != "PRO")
    ]
    sidechain_cds = reduced_topology.loc[
        (reduced_topology["name"] == "CD")
        | (reduced_topology["name"] == "CD1")
        | (reduced_topology["name"] == "CD2")
        | (reduced_topology["name"] == "SD")
        | (reduced_topology["name"] == "OD")
        | (reduced_topology["name"] == "OD1")
        | (reduced_topology["name"] == "OD2")
        | (reduced_topology["name"] == "ND1")
        | (reduced_topology["name"] == "ND2") & (reduced_topology["resname"] != "PRO")
    ]
    pairs = pd.DataFrame()

    # For backbone oxygen take the CB of the same residue and save in a pairs tuple
    pairs = pd.concat(
        [
            pairs,
            create_pairs_14_dataframe(atomtype1=backbone_oxygen, atomtype2=sidechain_cb, constant=1.5e-6, prefactor=1),
            # create_pairs_14_dataframe(atomtype1=backbone_oxygen, atomtype2=sidechain_cb, prefactor=1),
        ],
        axis=0,
        sort=False,
        ignore_index=True,
    )
    # now we add the pair between the last CB and the two OCT ones
    pairs = pd.concat(
        [
            pairs,
            create_pairs_14_dataframe(atomtype1=ct_oxygen, atomtype2=sidechain_cb, constant=1.5e-6, prefactor=1),
            # create_pairs_14_dataframe(atomtype1=ct_oxygen, atomtype2=sidechain_cb, prefactor=1),
        ],
        axis=0,
        sort=False,
        ignore_index=True,
    )
    # For each backbone nitrogen take the CB of the previuos residue and save in a pairs tuple
    pairs = pd.concat(
        [
            pairs,
            create_pairs_14_dataframe(
                atomtype1=backbone_nitrogen,
                atomtype2=sidechain_cb,
                constant=2.7e-6,
                prefactor=1.00,
                shift=-1,
            ),
        ],
        axis=0,
        sort=False,
        ignore_index=True,
    )
    # For the first backbone nitrogen take the N of the next residue and save in a pairs tuple
    pairs = pd.concat(
        [
            pairs,
            create_pairs_14_dataframe(
                atomtype1=first_backbone_nitrogen,
                atomtype2=backbone_nitrogen,
                constant=4.0e-06,
                shift=+1,
            ),
        ],
        axis=0,
        sort=False,
        ignore_index=True,
    )
    # For each backbone carbonyl take the CGs of the same residue and save in a pairs tuple
    pairs = pd.concat(
        [
            pairs,
            create_pairs_14_dataframe(atomtype1=sidechain_cgs, atomtype2=backbone_carbonyl, constant=1.2e-6, prefactor=0.25),
        ],
        axis=0,
        sort=False,
        ignore_index=True,
    )
    # For each backbone nitrogen take the CGs of the same residue and save in a pairs tuple
    pairs = pd.concat(
        [
            pairs,
            create_pairs_14_dataframe(atomtype1=sidechain_cgs, atomtype2=backbone_nitrogen, constant=5.5e-7, prefactor=0.2),
        ],
        axis=0,
        sort=False,
        ignore_index=True,
    )
    pairs = pd.concat(
        [
            pairs,
            create_pairs_14_dataframe(
                atomtype1=sidechain_cgs,
                atomtype2=first_backbone_nitrogen,
                constant=5.5e-7,
                prefactor=0.2,
            ),
        ],
        axis=0,
        sort=False,
        ignore_index=True,
    )
    pairs = pd.concat(
        [
            pairs,
            create_pairs_14_dataframe(atomtype1=sidechain_cds, atomtype2=backbone_calpha, constant=5e-7, prefactor=0.100),
        ],
        axis=0,
        sort=False,
        ignore_index=True,
    )

    # make it symmetric
    inv_LJ = pairs[["aj", "ai", "func", "c6", "c12", "probability", "rc_probability", "source"]].copy()
    inv_LJ.columns = [
        "ai",
        "aj",
        "func",
        "c6",
        "c12",
        "probability",
        "rc_probability",
        "source",
    ]

    pairs = pd.concat([pairs, inv_LJ], axis=0, sort=False, ignore_index=True)
    pairs["ai"] = pairs["ai"].astype(str)
    pairs["aj"] = pairs["aj"].astype(str)

    return pairs
