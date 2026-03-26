from . import topology

import pandas as pd
import sys


def generate_bonded_interactions(meGO_ensemble):
    """
    Generates the bonded interactions and stores them in meGO_ensemble

    Parameters
    ----------
    meGO_ensemble : dict
        The meGO_ensemble object containing all the relevant system information

    Returns
    -------
    meGO_ensemble : dict
        The updated meGO_ensemble object with updated/added bonded parameters
    """
    meGO_ensemble.meGO_bonded_interactions = {}
    meGO_ensemble.bond_pairs = {}
    meGO_ensemble.user_pairs = {}

    for molecule, topol in meGO_ensemble.topology.molecules.items():
        meGO_ensemble.meGO_bonded_interactions[molecule] = {
            "bonds": topology.get_bonds(topol[0].bonds),
            "angles": topology.get_angles(topol[0].angles),
            "dihedrals": topology.get_dihedrals(topol[0].dihedrals),
            "impropers": topology.get_impropers(topol[0].impropers),
            "pairs": topology.get_pairs(topol[0].adjusts),
        }
        # The following bonds are used in the parametrization of LJ 1-4
        meGO_ensemble.bond_pairs[molecule] = topology.get_bond_pairs(topol[0].bonds)
        meGO_ensemble.user_pairs[molecule] = topology.get_pairs(topol[0].adjusts)

    return meGO_ensemble


def generate_14_data(meGO_ensemble):
    """
    Generates data for 1-4 interactions within a molecular ensemble.

    Parameters
    ----------
    meGO_ensemble : dict
        A dictionary containing information about the molecular ensemble.

    Returns
    -------
    pairs14 : pd.DataFrame
        DataFrame containing information about 1-4 interactions.
    exclusion_bonds14 : pd.DataFrame
        DataFrame containing exclusion bonded interactions.
    """
    pairs14 = pd.DataFrame()
    exclusion_bonds14 = pd.DataFrame()
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

        bd = topology.compute_bond_distances(reduced_topology, bond_pair)
        exclusion_bonds14 = pd.concat([exclusion_bonds14, bd], axis=0, sort=False, ignore_index=True)

        reduced_topology["c12"] = reduced_topology["sb_type"].map(meGO_ensemble.sbtype_c12_dict)

        pairs = pd.DataFrame()
        if meGO_ensemble.molecule_type_dict[molecule] == "protein":
            pairs = topology.protein_LJ14(reduced_topology)
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
            pairs14 = pd.concat([pairs14, pairs], axis=0, sort=False, ignore_index=True)

    return pairs14, exclusion_bonds14
