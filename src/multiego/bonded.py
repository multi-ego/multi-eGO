from .model_config import config
from . import type_definitions

import pandas as pd
import networkx as nx
import numpy as np
from itertools import combinations
import sys


def compute_bond_distances(reduced_topology, bond_pair, max_distance=config.max_bond_separation):
    """
    Computes the shortest bond-path distance between every pair of atoms in a
    molecule, up to a configurable maximum.

    A graph is built from the provided bond pairs and shortest paths are
    computed using NetworkX. Pairs separated by more than ``max_distance``
    bonds — or not connected at all — are assigned a distance of
    ``max_distance + 1``. Self-distances are always 0. The result is symmetric
    (both (ai, aj) and (aj, ai) are present).

    Parameters
    ----------
    reduced_topology : pd.DataFrame
        Per-atom topology slice containing at least the columns ``sb_type``
        and ``number`` (atom index as string).
    bond_pair : list of tuple
        List of (atom_number_i, atom_number_j) pairs defining covalent bonds,
        as returned by ``topology.get_bond_pairs``.
    max_distance : int, optional
        Maximum bond separation to resolve; pairs beyond this distance are
        capped at ``max_distance + 1``. Defaults to
        ``config.max_bond_separation``.

    Returns
    -------
    pd.DataFrame
        Symmetric DataFrame with columns ``ai`` (sb_type), ``aj`` (sb_type),
        and ``bond_distance`` (int). Contains one row per ordered pair
        including self-pairs.
    """
    # Build atom number ↔ sb_type mapping
    sbtype_to_atnum = reduced_topology.set_index("sb_type")["number"].to_dict()

    sbtypes = list(sbtype_to_atnum.keys())

    # Build graph from bond pairs
    G = nx.Graph()
    G.add_edges_from(bond_pair)

    # Compute shortest path lengths (up to max_distance)
    all_lengths = dict(nx.all_pairs_shortest_path_length(G, cutoff=max_distance))

    # Collect all pairs with distances
    data = []

    for ai, aj in combinations(sbtypes, 2):
        ai_num = sbtype_to_atnum[ai]
        aj_num = sbtype_to_atnum[aj]

        dist = all_lengths.get(ai_num, {}).get(aj_num, config.max_bond_separation + 1)
        if dist > max_distance:
            dist = config.max_bond_separation + 1

        data.append((ai, aj, dist))
        data.append((aj, ai, dist))  # symmetric

    # Also include self-distances (0)
    for ai in sbtypes:
        data.append((ai, ai, 0))

    df = pd.DataFrame(data, columns=["ai", "aj", "bond_distance"])
    return df


def create_pairs_14_dataframe(atomtype1, atomtype2, c6=0.0, shift=0, prefactor=None, constant=None):
    """
    Builds a DataFrame of multi-eGO-specific 1-4 (or 1-4-like) LJ interaction
    pairs between two sets of atom types, optionally offset by a residue index
    shift.

    For each atom in ``atomtype1``, the matching atom in ``atomtype2`` is found
    at residue index ``resnum + shift``. When a match exists, the c12 parameter
    is determined by ``constant`` and/or ``prefactor``:

    - If only ``constant`` is given, c12 is set to that value.
    - If only ``prefactor`` is given, c12 is the combination-rule geometric mean
      scaled by ``prefactor``: ``prefactor * sqrt(c12_i * c12_j)``.
    - If both are given, c12 is the minimum of ``constant`` and the scaled
      geometric mean, allowing ``constant`` to act as an upper bound.

    At least one of ``prefactor`` or ``constant`` must be provided.

    Parameters
    ----------
    atomtype1 : pd.DataFrame
        Topology slice for the first atom group. Must contain columns
        ``number``, ``resnum``, and ``c12``.
    atomtype2 : pd.DataFrame
        Topology slice for the second atom group. Must contain columns
        ``number``, ``resnum``, and ``c12``.
    c6 : float, optional
        Fixed LJ c6 value applied to all pairs (default 0.0, i.e. purely
        repulsive).
    shift : int, optional
        Residue index offset applied when matching atoms in ``atomtype2``
        to those in ``atomtype1``. A positive shift matches the atom in the
        *next* residue; a negative shift matches the *previous* residue
        (default 0, same residue).
    prefactor : float or None, optional
        Multiplicative factor applied to the geometric-mean c12. If ``None``,
        the geometric mean is not computed.
    constant : float or None, optional
        Fixed c12 value. Acts as a ceiling when ``prefactor`` is also given.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``ai``, ``aj``, ``func``, ``c6``, ``c12``,
        ``probability``, ``rc_probability``, and ``source``. Contains one row
        per matched pair; empty if no residue-shifted matches are found.

    Raises
    ------
    ValueError
        If neither ``prefactor`` nor ``constant`` is provided.
    """
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
    Generates multi-eGO-specific 1-4 LJ pairs for a protein molecule.

    Protein backbone and sidechain atom groups are extracted from the
    topology and paired according to the rules defined in
    ``type_definitions.atom_type_combinations``. Each rule specifies two
    atom-group names, an optional prefactor applied to the combination-rule
    c12, a constant c12 ceiling, and a residue-index shift. The resulting
    pairs encode local geometric constraints specific to the multi-eGO
    representation (e.g. N–CB repulsion across the peptide bond).

    The output is symmetrised so that both (ai, aj) and (aj, ai) directions
    are present.

    Parameters
    ----------
    reduced_topology : pd.DataFrame
        Per-atom topology slice for a single protein molecule. Must contain
        columns ``number``, ``sb_type``, ``resnum``, ``name``, ``type``,
        ``resname``, ``molecule_type``, and ``c12``.

    Returns
    -------
    pd.DataFrame
        Symmetric DataFrame with columns ``ai``, ``aj``, ``func``, ``c6``,
        ``c12``, ``probability``, ``rc_probability``, and ``source``,
        where ``source`` is always ``"1-4"``. Atom indices (``ai``, ``aj``)
        are returned as strings matching the ``number`` column.
    """
    # Build named selections for the atom groups referenced in atom_type_combinations
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

    loc_geom_atoms = {
        "first_backbone_nitrogen": first_backbone_nitrogen,
        "backbone_nitrogen": backbone_nitrogen,
        "backbone_carbonyl": backbone_carbonyl,
        "backbone_calpha": backbone_calpha,
        "backbone_oxygen": backbone_oxygen,
        "ct_oxygen": ct_oxygen,
        "sidechain_cb": sidechain_cb,
        "sidechain_cgs": sidechain_cgs,
        "sidechain_cds": sidechain_cds,
    }

    pairs = pd.DataFrame()
    # iterate over all atom type combinations defined in the type_definitions
    # and create pairs for each combination
    for element in type_definitions.atom_type_combinations:
        pairs = pd.concat(
            [
                pairs,
                create_pairs_14_dataframe(
                    atomtype1=loc_geom_atoms[element[0]],
                    atomtype2=loc_geom_atoms[element[1]],
                    constant=element[3],
                    prefactor=element[2],
                    shift=element[4],
                ),
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


def generate_14_data(meGO_ensemble):
    """
    Generates 1-4 LJ interaction pairs for all molecules in the ensemble.

    For protein molecules, pairs are generated by ``protein_LJ14`` using the
    multi-eGO-specific local geometry rules. For non-protein molecules, pairs
    are read directly from the ``[pairs]`` section of the GROMACS topology
    (``meGO_ensemble.user_pairs``); those pairs must carry explicit C6/C12
    values or the function exits with an error.

    Parameters
    ----------
    meGO_ensemble : MeGOEnsemble
        Fully initialised ensemble. The following attributes are accessed:
        ``bond_pairs``, ``topology_dataframe``, ``sbtype_c12_dict``,
        ``molecule_type_dict``, and ``user_pairs``.

    Returns
    -------
    pairs14 : pd.DataFrame
        Combined 1-4 pairs for all molecules with columns ``ai``, ``aj``
        (sb_type as category), ``func``, ``c6``, ``c12``, ``probability``,
        ``rc_probability``, ``source``, ``rep``, ``same_chain``,
        ``molecule_name_ai``, and ``molecule_name_aj``.
    """
    pairs14 = pd.DataFrame()
    for idx, (molecule, bond_pair) in enumerate(meGO_ensemble.bond_pairs.items(), start=1):
        if not bond_pair:
            continue
        reduced_topology = meGO_ensemble.topology_dataframe.loc[meGO_ensemble.topology_dataframe["molecule_name"] == molecule][
            [
                "number",
                "sb_type",
                "resnum",
                "name",
                "type",
                "resname",
                "molecule_type",
            ]
        ].copy()

        reduced_topology["number"] = reduced_topology["number"].astype(str)
        reduced_topology["resnum"] = reduced_topology["resnum"].astype(int)
        type_atnum_dict = reduced_topology.set_index("number")["sb_type"].to_dict()
        reduced_topology["c12"] = reduced_topology["sb_type"].map(meGO_ensemble.sbtype_c12_dict)

        pairs = pd.DataFrame()
        if meGO_ensemble.molecule_type_dict[molecule] == "protein":
            pairs = protein_LJ14(reduced_topology)
            pairs["ai"] = pairs["ai"].map(type_atnum_dict).astype("category")
            pairs["aj"] = pairs["aj"].map(type_atnum_dict).astype("category")
            pairs["rep"] = pairs["c12"]
            pairs["source"] = pairs["source"].astype("category")
            pairs["same_chain"] = True
        else:
            pairs["ai"] = meGO_ensemble.user_pairs[molecule].ai.astype(str)
            pairs["aj"] = meGO_ensemble.user_pairs[molecule].aj.astype(str)
            pairs["ai"] = pairs["ai"].map(type_atnum_dict).astype("category")
            pairs["aj"] = pairs["aj"].map(type_atnum_dict).astype("category")

            nonprotein_c12 = []
            for test in meGO_ensemble.user_pairs[molecule].type:
                if test is None:
                    print(
                        "\nERROR: you have 1-4 pairs defined in your reference topology without the associated C6/C12 values"
                    )
                    print("       user provided 1-4 pairs need to define also the C6/C12\n")
                    sys.exit()
                nonprotein_c12.append(float(test.epsilon) * 4.184)

            pairs["func"] = 1
            pairs["c6"] = 0.0
            pairs["c12"] = nonprotein_c12
            pairs["probability"] = 1.0
            pairs["rc_probability"] = 1.0
            pairs["source"] = pd.Series(["1-4"] * len(pairs), dtype="category")
            pairs["rep"] = pairs["c12"]
            pairs["same_chain"] = True
            tmp = pairs.copy()
            tmp["ai"], tmp["aj"] = tmp["aj"], tmp["ai"]
            pairs = pd.concat([pairs, tmp], axis=0, sort=False, ignore_index=True)

        mol_ai = f"{idx}_{molecule}"
        pairs["molecule_name_ai"] = mol_ai
        pairs["molecule_name_aj"] = mol_ai

        if not pairs.empty:
            # this collect only the interactions that are exactly 1-4
            pairs14 = pd.concat([pairs14, pairs], axis=0, sort=False, ignore_index=True)

    return pairs14


def generate_bond_distance_data(meGO_ensemble):
    """
    Computes pairwise bond-path distances for all molecules in the ensemble.

    Iterates over every molecule in ``meGO_ensemble.bond_pairs`` and calls
    ``compute_bond_distances`` to build a symmetric distance table up to
    ``config.max_bond_separation``. Results from all molecules are
    concatenated into a single DataFrame.

    Parameters
    ----------
    meGO_ensemble : MeGOEnsemble
        Fully initialised ensemble. The following attributes are accessed:
        ``bond_pairs`` and ``topology_dataframe``.

    Returns
    -------
    pd.DataFrame
        Concatenated bond-distance table for all molecules with columns
        ``ai`` (sb_type), ``aj`` (sb_type), and ``bond_distance`` (int).
        Pairs beyond ``config.max_bond_separation`` bonds are assigned
        ``max_bond_separation + 1``.
    """
    all_bd = pd.DataFrame()
    for idx, (molecule, bond_pair) in enumerate(meGO_ensemble.bond_pairs.items(), start=1):
        if not bond_pair:
            continue
        reduced_topology = meGO_ensemble.topology_dataframe.loc[meGO_ensemble.topology_dataframe["molecule_name"] == molecule][
            [
                "number",
                "sb_type",
            ]
        ].copy()

        reduced_topology["number"] = reduced_topology["number"].astype(str)

        # this assigns a bond distance to any pairs of atom in a molecule
        mol_bd = compute_bond_distances(reduced_topology, bond_pair)
        all_bd = pd.concat([all_bd, mol_bd], axis=0, sort=False, ignore_index=True)

    return all_bd
