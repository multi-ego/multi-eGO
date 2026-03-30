from . import type_definitions
from .model_config import config

from collections import defaultdict, deque
import numpy as np
import pandas as pd


def generate_bond_exclusions(reduced_topology, bond_pair):
    """
    Enumerate bond exclusions, 1-4 pairs, and nth-bond pairs for a molecule.

    Performs a BFS from every atom up to ``config.max_bond_separation`` bonds.
    Atoms within ``config.bond14_separation`` bonds are collected as exclusions;
    atoms at exactly ``config.bond14_separation`` bonds form the 1-4 list; atoms
    between ``config.bond14_separation + 1`` and ``config.max_bond_separation``
    bonds form the nth-bond list.

    Each entry in the returned lists is a ``"<atom_i>_<atom_j>"`` string
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
    exclusion_bonds : list of str
        Bidirectional ``"i_j"`` identifiers for all pairs within
        ``config.bond14_separation`` bonds (1-2, 1-3, 1-4).
    p14 : list of str
        Bidirectional ``"i_j"`` identifiers for pairs at exactly
        ``config.bond14_separation`` bonds (1-4 only).
    nth_bonds : list of str
        Bidirectional ``"i_j"`` identifiers for all pairs within
        ``config.max_bond_separation`` bonds (superset of ``exclusion_bonds``).
    """
    # Build the connectivity graph
    graph = defaultdict(set)
    for a, b in bond_pair:
        graph[a].add(b)
        graph[b].add(a)

    exclusion_bonds, p14, nth_bonds = [], [], []

    for atom in reduced_topology["number"].to_list():
        visited = set([atom])
        ex_tmp = set()
        nth_tmp = set()
        p14_tmp = set()
        queue = deque([(atom, 0)])

        while queue:
            current_atom, depth = queue.popleft()

            if depth == 0:
                pass  # skip the origin atom
            elif 0 < depth <= config.bond14_separation:
                ex_tmp.add(current_atom)
                if depth == config.bond14_separation:
                    p14_tmp.add(current_atom)
            elif config.bond14_separation < depth <= config.max_bond_separation:
                nth_tmp.add(current_atom)

            # Stop traversal after depth 5
            if depth < config.max_bond_separation:
                for neighbor in graph[current_atom]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append((neighbor, depth + 1))

        # Add bidirectional string identifiers
        for e in ex_tmp:
            exclusion_bonds.append(f"{atom}_{e}")
            exclusion_bonds.append(f"{e}_{atom}")
        for e in nth_tmp:
            nth_bonds.append(f"{atom}_{e}")
            nth_bonds.append(f"{e}_{atom}")
        for e in p14_tmp:
            p14.append(f"{atom}_{e}")
            p14.append(f"{e}_{atom}")

    return exclusion_bonds, p14, nth_bonds


def make_pairs_exclusion_topology(meGO_ensemble, meGO_LJ_14, args):
    """
    Prepares the [ pairs ] and [ exclusions ] sections for output to topol_mego.top.

    Parameters
    ----------
    meGO_ensemble : dict
        The meGO_ensemble object contains all the relevant system information.
    meGO_LJ_14 : pd.DataFrame
        Contains the contact information for the 1-4 interactions.
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    pairs_molecule_dict : dict
        Write-ready pairs/exclusions interactions keyed by molecule name.
    """
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

        exclusion_bonds, p14, nthbond = generate_bond_exclusions(reduced_topology, bond_pair)

        pairs = pd.DataFrame()

        # in this case, on top of the already prepared 1-4 interactions, we need to set as repulsive
        # all the interactions up to nth-bonds. Some are already repulsive, but other not.
        # So here we generate them iterating over the list of relevant atom pairs. We need to be carefull
        # to assign the correct value for the repulsive interaction. So every pair that has some special repulsive
        # value in type_definitions need to be included here.
        if args.egos == "mg":
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

            df = pd.DataFrame(filtered_combinations, columns=["ai", "aj"])
            df["c6"] = 0.0
            df["c12"] = np.sqrt(
                df["ai"].map(meGO_ensemble.sbtype_c12_dict) * df["aj"].map(meGO_ensemble.sbtype_c12_dict)
            )

            # Apply special c12 overrides from the table in type_definitions.
            # Map atom types once, then loop over rules.
            type_ai = df["ai"].map(meGO_ensemble.sbtype_type_dict)
            type_aj = df["aj"].map(meGO_ensemble.sbtype_type_dict)

            for types_i, types_j, c12_val in type_definitions.NTHBOND_C12_OVERRIDES:
                mask = type_ai.isin(types_i) & type_aj.isin(types_j)
                if types_i != types_j:
                    mask |= type_ai.isin(types_j) & type_aj.isin(types_i)
                df.loc[mask, "c12"] = c12_val

            df["same_chain"] = True
            df["probability"] = 1.0
            df["rc_probability"] = 1.0
            df["source"] = "mg"
            df["rep"] = df["c12"]
            df["func"] = 1
            df["molecule_name_ai"] = mol_ai
            df["molecule_name_aj"] = mol_ai
            pairs = pd.concat([meGO_LJ_14, df], axis=0, sort=False, ignore_index=True)

        elif args.egos == "production" and not meGO_LJ_14.empty:
            pairs = meGO_LJ_14[meGO_LJ_14["molecule_name_ai"] == mol_ai][
                [
                    "ai",
                    "aj",
                    "c6",
                    "c12",
                    "same_chain",
                    "probability",
                    "rc_probability",
                    "mg_sigma",
                    "mg_epsilon",
                    "source",
                    "rep",
                    "bond_distance",
                ]
            ].copy()

            if not pairs.empty:
                # Within 5 bonds: use default repulsion
                pairs.loc[(~pairs["same_chain"]) & (pairs["bond_distance"] <= config.max_bond_separation), "c6"] = 0.0
                pairs.loc[(~pairs["same_chain"]) & (pairs["bond_distance"] <= config.max_bond_separation), "c12"] = (
                    pairs["rep"]
                )
                # Beyond 5 bonds: use MG prior
                pairs.loc[
                    (~pairs["same_chain"])
                    & (pairs["bond_distance"] > config.max_bond_separation)
                    & (pairs["mg_epsilon"] < 0.0),
                    "c6",
                ] = 0.0
                pairs.loc[
                    (~pairs["same_chain"])
                    & (pairs["bond_distance"] > config.max_bond_separation)
                    & (pairs["mg_epsilon"] < 0.0),
                    "c12",
                ] = -pairs["mg_epsilon"]
                pairs.loc[
                    (~pairs["same_chain"])
                    & (pairs["bond_distance"] > config.max_bond_separation)
                    & (pairs["mg_epsilon"] > 0.0),
                    "c6",
                ] = (
                    4 * pairs["mg_epsilon"] * (pairs["mg_sigma"] ** 6)
                )
                pairs.loc[
                    (~pairs["same_chain"])
                    & (pairs["bond_distance"] > config.max_bond_separation)
                    & (pairs["mg_epsilon"] > 0.0),
                    "c12",
                ] = (
                    4 * pairs["mg_epsilon"] * (pairs["mg_sigma"] ** 12)
                )

        if not pairs.empty:
            pairs["ai"] = pairs["ai"].map(atnum_type_dict)
            pairs["aj"] = pairs["aj"].map(atnum_type_dict)
            pairs["check"] = pairs["ai"].astype(str) + "_" + pairs["aj"].astype(str)
            pairs["remove"] = ""
            pairs.loc[(pairs["check"].isin(exclusion_bonds)), "remove"] = "Yes"
            pairs.loc[(pairs["check"].isin(p14) & (pairs["same_chain"])), "remove"] = "No"
            mask = pairs.remove == "Yes"
            pairs = pairs[~mask]
            pairs["func"] = 1
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
            pairs.drop_duplicates(inplace=True, ignore_index=True)
            pairs.sort_values(by=["ai", "aj"], inplace=True)

        pairs_molecule_dict[molecule] = pairs

    return pairs_molecule_dict
