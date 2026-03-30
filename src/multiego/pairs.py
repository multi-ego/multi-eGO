from . import bonded, type_definitions
from .model_config import config

from collections import defaultdict, deque
import numpy as np
import pandas as pd


def generate_bond_exclusions(reduced_topology, bond_pair):
    """
    Enumerate nth-bond pairs for a molecule.

    Performs a BFS from every atom up to ``config.max_bond_separation`` bonds.
    Atoms between ``config.bond14_separation + 1`` and
    ``config.max_bond_separation`` bonds form the nth-bond list.

    Each entry in the returned list is a ``"<atom_i>_<atom_j>"`` string
    identifier, and every pair is recorded in both directions.

    Parameters
    ----------
    reduced_topology : pd.DataFrame
        Per-atom topology slice for a single molecule.  Must contain the
        column ``number`` (atom index as string).
    bond_pair : list of tuple
        List of ``(atom_number_i, atom_number_j)`` pairs defining covalent
        bonds, as integers.

    Returns
    -------
    nth_bonds : list of str
        Bidirectional ``"i_j"`` identifiers for pairs between
        ``config.bond14_separation + 1`` and ``config.max_bond_separation``
        bonds apart.
    """
    # Build the connectivity graph
    graph = defaultdict(set)
    for a, b in bond_pair:
        graph[a].add(b)
        graph[b].add(a)

    nth_bonds = []

    for atom in reduced_topology["number"].to_list():
        visited = set([atom])
        nth_tmp = set()
        queue = deque([(atom, 0)])

        while queue:
            current_atom, depth = queue.popleft()

            if depth <= config.bond14_separation:
                pass  # skip the origin atom
            elif config.bond14_separation < depth <= config.max_bond_separation:
                nth_tmp.add(current_atom)

            # Stop traversal after depth 5
            if depth < config.max_bond_separation:
                for neighbor in graph[current_atom]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append((neighbor, depth + 1))

        # Add bidirectional string identifiers
        for e in nth_tmp:
            nth_bonds.append(f"{atom}_{e}")
            nth_bonds.append(f"{e}_{atom}")

    return nth_bonds


def make_pairs_exclusion_topology(meGO_ensemble, args, meGO_LJ_14=None):
    """
    Prepares the [ pairs ] and [ exclusions ] sections for output to topol_mego.top.

    The 1-4 pairs are always generated internally via
    ``bonded.generate_14_data``.  For ``egos == "production"`` the caller may
    additionally supply *meGO_LJ_14* (the trained LJ interactions), which are
    merged into the pair list.

    Parameters
    ----------
    meGO_ensemble : dict
        The meGO_ensemble object contains all the relevant system information.
    args : argparse.Namespace
        Parsed command-line arguments.
    meGO_LJ_14 : pd.DataFrame, optional
        Trained LJ interactions for production runs.  Not used for
        ``egos == "mg"``.

    Returns
    -------
    pairs_molecule_dict : dict
        Write-ready pairs/exclusions interactions keyed by molecule name.
    """
    # Generate immutable 1-4 pairs from the bonded topology
    pairs14_all = bonded.generate_14_data(meGO_ensemble)

    pairs_molecule_dict = {}
    for idx, (molecule, bond_pair) in enumerate(meGO_ensemble.bond_pairs.items(), start=1):
        reduced_topology = meGO_ensemble.topology_dataframe.loc[
            meGO_ensemble.topology_dataframe["molecule_name"] == molecule
        ][
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
        mol_ai = f"{idx}_{molecule}"

        type_atnum_dict = reduced_topology.astype({"number": int}).set_index("number")["sb_type"].to_dict()
        atnum_type_dict = reduced_topology.set_index("sb_type")["number"].to_dict()

        nthbond = generate_bond_exclusions(reduced_topology, bond_pair)

        pairs = pd.DataFrame()

        # in this case, on top of the already prepared 1-4 interactions, we need to set as repulsive
        # all the interactions up to nth-bonds. Some are already repulsive, but other not.
        # So here we generate them iterating over the list of relevant atom pairs. We need to be carefull
        # to assign the correct value for the repulsive interaction. So every pair that has some special repulsive
        # value in type_definitions need to be included here.
        # Build nth-bond pair list, filtering out H-X pairs where X is
        # not an allowed partner (e.g. H-CH2 is skipped).
        filtered_combinations = []
        for pair in nthbond:
            a_str, b_str = pair.split("_")
            a, b = int(a_str), int(b_str)

            if a not in type_atnum_dict or b not in type_atnum_dict:
                continue

            ai = type_atnum_dict[a]
            aj = type_atnum_dict[b]
            ti = meGO_ensemble.sbtype_type_dict[ai]
            tj = meGO_ensemble.sbtype_type_dict[aj]

            if (ti == "H" and tj not in type_definitions.H_ALLOWED_PARTNERS) or (
                tj == "H" and ti not in type_definitions.H_ALLOWED_PARTNERS
            ):
                continue

            filtered_combinations.append((ai, aj))

        pairs = pd.DataFrame(filtered_combinations, columns=["ai", "aj"])
        pairs["c6"] = 0.0
        pairs["c12"] = np.sqrt(
            pairs["ai"].map(meGO_ensemble.sbtype_c12_dict) * pairs["aj"].map(meGO_ensemble.sbtype_c12_dict)
        )

        # Apply special c12 overrides from the table in type_definitions.
        # Map atom types once, then loop over rules.
        type_ai = pairs["ai"].map(meGO_ensemble.sbtype_type_dict)
        type_aj = pairs["aj"].map(meGO_ensemble.sbtype_type_dict)

        for types_i, types_j, c12_val in type_definitions.NTHBOND_C12_OVERRIDES:
            mask = type_ai.isin(types_i) & type_aj.isin(types_j)
            if types_i != types_j:
                mask |= type_ai.isin(types_j) & type_aj.isin(types_i)
            pairs.loc[mask, "c12"] = c12_val

        pairs["same_chain"] = True
        pairs["probability"] = 1.0
        pairs["rc_probability"] = 1.0
        pairs["source"] = "mg"
        pairs["rep"] = pairs["c12"]
        pairs["func"] = 1

        # then we add intramolecular intereactions paired to an intermolecular one
        # and we regenerate prior intramolecular interactions
        if args.egos == "production" and meGO_LJ_14 is not None and not meGO_LJ_14.empty:
            meGO_pairs = meGO_LJ_14[meGO_LJ_14["molecule_name_ai"] == mol_ai][
                [
                    "ai",
                    "aj",
                    "c6",
                    "c12",
                    "probability",
                    "rc_probability",
                    "source",
                ]
            ].copy()

            meGO_pairs["func"] = 1
            meGO_pairs = meGO_pairs[meGO_pairs["c12"] > 0.0]
            pairs = pd.concat([pairs, meGO_pairs], ignore_index=True)

        # Append the immutable 1-4 pairs for this molecule
        if not pairs14_all.empty:
            pairs14 = pairs14_all[pairs14_all["molecule_name_ai"] == mol_ai][
                [
                    "ai",
                    "aj",
                    "func",
                    "c6",
                    "c12",
                    "probability",
                    "rc_probability",
                    "source",
                ]
            ].copy()
            pairs = pd.concat([pairs, pairs14], ignore_index=True)

        if not pairs.empty:
            pairs["ai"] = pairs["ai"].map(atnum_type_dict)
            pairs["aj"] = pairs["aj"].map(atnum_type_dict)
            pairs = pairs[pairs["c12"] > 0.0]
            pairs = pairs[
                [
                    "ai",
                    "aj",
                    "func",
                    "c6",
                    "c12",
                    "probability",
                    "rc_probability",
                    "source",
                ]
            ]
            pairs.dropna(inplace=True)
            pairs["ai"] = pairs["ai"].astype(int)
            pairs["aj"] = pairs["aj"].astype(int)

            inv_pairs = pairs[["aj", "ai", "func", "c6", "c12", "probability", "rc_probability", "source"]].copy()
            inv_pairs.columns = ["ai", "aj", "func", "c6", "c12", "probability", "rc_probability", "source"]
            pairs = pd.concat([pairs, inv_pairs], axis=0, sort=False, ignore_index=True)
            pairs = pairs[pairs["ai"] < pairs["aj"]]

            # Deduplicate (ai, aj) pairs, keeping the highest-priority source.
            # Priority: "1-4" (highest) > any other source > "mg" (lowest).
            pairs["_priority"] = pairs["source"].map(lambda s: 0 if s == "1-4" else 2 if s == "mg" else 1)
            pairs.sort_values(by=["ai", "aj", "_priority"], inplace=True)
            pairs.drop_duplicates(subset=["ai", "aj"], keep="first", inplace=True, ignore_index=True)
            pairs.drop(columns="_priority", inplace=True)
            pairs.sort_values(by=["ai", "aj"], inplace=True, ignore_index=True)

        pairs_molecule_dict[molecule] = pairs

    return pairs_molecule_dict
