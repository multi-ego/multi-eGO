from . import type_definitions
from .model_config import config

import numpy as np
import pandas as pd
import itertools


def _make_combinations(sbtype1, sbtype2):
    """Return the full ordered list of (ai, aj) pairs for the given sb_type sets."""
    if sbtype1 == sbtype2:
        return list(itertools.product(sbtype1, repeat=2))
    return list(set(itertools.product(sbtype1, sbtype2)) | set(itertools.product(sbtype2, sbtype1)))


def generate_MG_LJ_pairs_rep(sbtype1, sbtype2, sbtype_c12_dict, c12_rep=None, factor=1.0):
    """
    Generates repulsive LJ pairs for the MG prior force field.

    Parameters
    ----------
    sbtype1, sbtype2 : list of str
        Lists of sb_type identifiers for the two interacting groups.
    sbtype_c12_dict : dict
        Mapping from sb_type to rc_c12 value, used when c12_rep is not specified.
    c12_rep : float or None
        Fixed repulsive c12 value. If None, computed as geometric mean from dictionary.
    factor : float
        Divisor applied to c12_rep (default 1.0).

    Returns
    -------
    pd.DataFrame
        Pairs DataFrame with columns ai, aj, c6, c12, epsilon, sigma, mg_sigma, mg_epsilon.
    """
    combinations = _make_combinations(sbtype1, sbtype2)

    if c12_rep is None:
        c12_rep = np.array([np.sqrt(sbtype_c12_dict[ai] * sbtype_c12_dict[aj]) for ai, aj in combinations])

    pairs_LJ = pd.DataFrame(combinations, columns=["ai", "aj"])
    pairs_LJ["c12"] = c12_rep / factor
    pairs_LJ["c6"] = 0.0
    pairs_LJ["epsilon"] = -c12_rep / factor
    pairs_LJ["sigma"] = (c12_rep / factor) ** (1.0 / 12.0) / 2.0 ** (1.0 / 6.0)
    pairs_LJ["mg_sigma"] = pairs_LJ["sigma"]
    pairs_LJ["mg_epsilon"] = -c12_rep / factor

    return pairs_LJ.reset_index(drop=True)


def generate_MG_LJ_pairs_attr(sbtype1, sbtype2, sbtype_mg_c12_dict, sbtype_mg_c6_dict, epsilon=None, sigma=None):
    """
    Generates attractive LJ pairs for the MG prior force field.

    Parameters
    ----------
    sbtype1, sbtype2 : list of str
        Lists of sb_type identifiers for the two interacting groups.
    sbtype_mg_c12_dict : dict
        Mapping from sb_type to mg_c12, used to derive sigma when not specified.
    sbtype_mg_c6_dict : dict
        Mapping from sb_type to mg_c6, used to derive sigma when not specified.
    epsilon : float or None
        Interaction well depth.
    sigma : float or None
        Interaction distance parameter. Cannot be provided without epsilon.

    Returns
    -------
    pd.DataFrame
        Pairs DataFrame with columns ai, aj, c6, c12, epsilon, sigma, mg_sigma, mg_epsilon.
    """
    combinations = _make_combinations(sbtype1, sbtype2)

    if epsilon is not None and sigma is None:
        c6 = np.array([np.sqrt(sbtype_mg_c6_dict[ai] * sbtype_mg_c6_dict[aj]) for ai, aj in combinations])
        c12_rep = np.array([np.sqrt(sbtype_mg_c12_dict[ai] * sbtype_mg_c12_dict[aj]) for ai, aj in combinations])
        sigma = (c12_rep / c6) ** (1.0 / 6.0)
    elif epsilon is None and sigma is not None:
        raise ValueError("You can provide a custom value for epsilon but not for sigma alone")

    pairs_LJ = pd.DataFrame(combinations, columns=["ai", "aj"])
    pairs_LJ["c12"] = 4.0 * epsilon * sigma**12.0
    pairs_LJ["c6"] = 4.0 * epsilon * sigma**6.0
    pairs_LJ["epsilon"] = epsilon
    pairs_LJ["sigma"] = sigma
    pairs_LJ["mg_sigma"] = sigma
    pairs_LJ["mg_epsilon"] = epsilon

    return pairs_LJ.reset_index(drop=True)


def generate_MG_LJ(meGO_ensemble):
    """
    Generates the molten-globule prior non-bonded interactions.

    The MG force field includes special repulsive and attractive interaction pairs
    (e.g. O-O, H-H, O-H) defined in type_definitions.special_non_local.

    Parameters
    ----------
    meGO_ensemble : dict
        The initialized meGO ensemble.

    Returns
    -------
    rc_LJ : pd.DataFrame
        DataFrame of MG prior LJ interactions.
    """
    chunks = []
    for special in type_definitions.special_non_local:
        sbtype_a = [
            sbtype for sbtype, atomtype in meGO_ensemble.sbtype_type_dict.items() if atomtype in special["atomtypes"][0]
        ]
        sbtype_b = [
            sbtype for sbtype, atomtype in meGO_ensemble.sbtype_type_dict.items() if atomtype in special["atomtypes"][1]
        ]
        if special["interaction"] == "rep":
            temp_LJ = generate_MG_LJ_pairs_rep(
                sbtype_a,
                sbtype_b,
                meGO_ensemble.sbtype_c12_dict,
                c12_rep=special["epsilon"],
                factor=special.get("factor", 1.0),
            )
        elif special["interaction"] == "att":
            temp_LJ = generate_MG_LJ_pairs_attr(
                sbtype_a,
                sbtype_b,
                meGO_ensemble.sbtype_mg_c12_dict,
                meGO_ensemble.sbtype_mg_c6_dict,
                epsilon=special["epsilon"],
                sigma=special["sigma"],
            )
        else:
            raise ValueError(f"Unknown interaction type {special['interaction']} in special_non_local definition")

        chunks.append(temp_LJ)

    rc_LJ = pd.concat(chunks, axis=0, ignore_index=True)

    rc_LJ["type"] = 1
    rc_LJ["same_chain"] = False
    rc_LJ["source"] = "mg"
    rc_LJ["reference"] = "mg"
    rc_LJ["rep"] = rc_LJ["c12"]
    rc_LJ["probability"] = 1.0
    rc_LJ["rc_probability"] = 1.0
    rc_LJ["rc_threshold"] = 1.0
    rc_LJ["md_threshold"] = 1.0
    rc_LJ["learned"] = 0
    rc_LJ["bond_distance"] = config.max_bond_separation + 1
    molecule_names_dictionary = {name.split("_", 1)[1]: name for name in meGO_ensemble.molecules_idx_sbtype_dictionary}
    rc_LJ["molecule_name_ai"] = rc_LJ["ai"].apply(lambda x: "_".join(x.split("_")[1:-1])).map(molecule_names_dictionary)
    rc_LJ["molecule_name_aj"] = rc_LJ["aj"].apply(lambda x: "_".join(x.split("_")[1:-1])).map(molecule_names_dictionary)
    rc_LJ["ai"] = rc_LJ["ai"].astype("category")
    rc_LJ["aj"] = rc_LJ["aj"].astype("category")
    rc_LJ["molecule_name_ai"] = rc_LJ["molecule_name_ai"].astype("category")
    rc_LJ["molecule_name_aj"] = rc_LJ["molecule_name_aj"].astype("category")

    return rc_LJ
