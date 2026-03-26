from . import type_definitions
from . import io

import numpy as np
import pandas as pd


def assign_molecule_type(molecule_type_dict, molecule_name, molecule_topology):
    """
    Decides if the molecule type of the system is 'protein', 'nucleic_acid' or 'other'
    and writes said information into molecule_type_dict before returning it

    Parameters
    ----------
    molecule_type_dict : dict
        Contains the molecule type information per system
    molecule_name : str
        The name of system
    molecule_topology : parmed.Topology
        The topology of the molecule, which will be used to figure out the molecule_type

    Returns
    -------
    molecule_type_dict : dict
        Updated molecule_type_dict with the added new system name
    """
    first_aminoacid = molecule_topology.residues[0].name

    if first_aminoacid in type_definitions.aminoacids_list:
        molecule_type_dict[molecule_name] = "protein"
    elif first_aminoacid in type_definitions.nucleic_acid_list:
        molecule_type_dict[molecule_name] = "nucleic_acid"
    else:
        molecule_type_dict[molecule_name] = "other"

    return molecule_type_dict


def initialize_topology(topology, custom_dict, args):
    """
    Initializes a topology DataFrame using provided molecule information.
    """

    (
        ensemble_topology_dataframe,
        new_number,
        col_molecule,
        new_resnum,
        ensemble_molecules_idx_sbtype_dictionary,
    ) = (pd.DataFrame(), [], [], [], {})

    molecule_type_dict = {}

    for molecule_number, (molecule_name, molecule_topology) in enumerate(topology.molecules.items(), 1):
        molecule_type_dict = assign_molecule_type(molecule_type_dict, molecule_name, molecule_topology[0])
        ensemble_molecules_idx_sbtype_dictionary[f"{str(molecule_number)}_{molecule_name}"] = {}
        ensemble_topology_dataframe = pd.concat([ensemble_topology_dataframe, molecule_topology[0].to_dataframe()], axis=0)
        for atom in molecule_topology[0].atoms:
            new_number.append(str(atom.idx + 1))
            col_molecule.append(f"{molecule_number}_{molecule_name}")
            new_resnum.append(str(atom.residue.number))

    ensemble_topology_dataframe["number"] = new_number
    ensemble_topology_dataframe["molecule"] = col_molecule
    ensemble_topology_dataframe["molecule_number"] = col_molecule
    ensemble_topology_dataframe[["molecule_number", "molecule_name"]] = ensemble_topology_dataframe.molecule.str.split(
        "_", expand=True, n=1
    )
    ensemble_topology_dataframe["resnum"] = new_resnum
    ensemble_topology_dataframe["cgnr"] = ensemble_topology_dataframe["resnum"]
    ensemble_topology_dataframe["ptype"] = "A"

    # Extending the from_ff_to_multiego dictionary to include the custom dictionary for special molecules
    # (if none is present the extended dictionary is equivalent to the standard one)
    from_ff_to_multiego_extended = type_definitions.from_ff_to_multiego.copy()
    from_ff_to_multiego_extended.update(custom_dict)

    ensemble_topology_dataframe = ensemble_topology_dataframe.replace({"name": from_ff_to_multiego_extended})
    ensemble_topology_dataframe["sb_type"] = (
        ensemble_topology_dataframe["name"]
        + "_"
        + ensemble_topology_dataframe["molecule_name"]
        + "_"
        + ensemble_topology_dataframe["resnum"].astype(str)
    )

    atp_c12_map = {k: v for k, v in zip(type_definitions.gromos_atp["name"], type_definitions.gromos_atp["rc_c12"])}
    atp_mg_c6_map = {k: v for k, v in zip(type_definitions.gromos_atp["name"], type_definitions.gromos_atp["mg_c6"])}
    atp_mg_c12_map = {k: v for k, v in zip(type_definitions.gromos_atp["name"], type_definitions.gromos_atp["mg_c12"])}

    # TODO: we will need to extend this to mg c6/c12 stuff?
    if args.custom_c12 is not None:
        custom_c12_dict = io.read_custom_c12_parameters(args.custom_c12)
        name_to_c12_appo = {key: val for key, val in zip(custom_c12_dict.name, custom_c12_dict.c12)}
        atp_c12_map.update(name_to_c12_appo)

    ensemble_topology_dataframe["charge"] = 0.0
    ensemble_topology_dataframe["rc_c6"] = 0.0
    ensemble_topology_dataframe["rc_c12"] = ensemble_topology_dataframe["type"].map(atp_c12_map)
    ensemble_topology_dataframe["mg_c6"] = ensemble_topology_dataframe["type"].map(atp_mg_c6_map)
    ensemble_topology_dataframe["mg_c12"] = ensemble_topology_dataframe["type"].map(atp_mg_c12_map)
    ensemble_topology_dataframe["molecule_type"] = ensemble_topology_dataframe["molecule_name"].map(molecule_type_dict)

    for molecule in ensemble_molecules_idx_sbtype_dictionary.keys():
        temp_topology_dataframe = ensemble_topology_dataframe.loc[ensemble_topology_dataframe["molecule"] == molecule]
        number_sbtype_dict = temp_topology_dataframe[["number", "sb_type"]].set_index("number")["sb_type"].to_dict()
        ensemble_molecules_idx_sbtype_dictionary[molecule] = number_sbtype_dict

    sbtype_c12_dict = ensemble_topology_dataframe[["sb_type", "rc_c12"]].set_index("sb_type")["rc_c12"].to_dict()
    sbtype_mg_c12_dict = ensemble_topology_dataframe[["sb_type", "mg_c12"]].set_index("sb_type")["mg_c12"].to_dict()
    sbtype_mg_c6_dict = ensemble_topology_dataframe[["sb_type", "mg_c6"]].set_index("sb_type")["mg_c6"].to_dict()
    sbtype_name_dict = ensemble_topology_dataframe[["sb_type", "name"]].set_index("sb_type")["name"].to_dict()
    sbtype_moltype_dict = (
        ensemble_topology_dataframe[["sb_type", "molecule_type"]].set_index("sb_type")["molecule_type"].to_dict()
    )

    return (
        ensemble_topology_dataframe,
        ensemble_molecules_idx_sbtype_dictionary,
        sbtype_c12_dict,
        sbtype_mg_c12_dict,
        sbtype_mg_c6_dict,
        sbtype_name_dict,
        sbtype_moltype_dict,
        molecule_type_dict,
    )


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
