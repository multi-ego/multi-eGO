import pandas as pd


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
        impropers_data.append(
            {"ai": ai, "aj": aj, "ak": ak, "al": al, "funct": funct, "psi_eq": psi_eq, "psi_k": psi_k}
        )

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
