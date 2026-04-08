from . import type_definitions
from . import mg
from . import bonded
from . import _term
from .model_config import config
from . import io

import numpy as np
import pandas as pd
import itertools
import sys
import time


def create_linearized_mask(
    set1,
    set2,
    types,
    symmetrize=False,
    inner_op=lambda x, y: x == y,
    outer_op=lambda x, y: x * y,
):
    """
    Creates a linearized boolean mask based on comparison operations between two sets.

    For each ``(type1, type2)`` pair in ``types``, ``inner_op`` is applied
    element-wise to each set and the results are combined with ``outer_op``.
    The final mask is the logical OR across all type pairs.  When
    ``symmetrize=True`` each ``(type1, type2)`` entry is automatically
    complemented by ``(type2, type1)``.

    Parameters
    ----------
    set1 : numpy.ndarray
        First 1-D array of values (e.g. atom-type labels).
    set2 : numpy.ndarray
        Second 1-D array of values, same length as ``set1``.
    types : list of tuple
        Pairs ``(type1, type2)`` used for element-wise comparisons.
    symmetrize : bool, optional
        If ``True``, automatically add the reversed pair ``(type2, type1)``
        for every entry in ``types`` (default ``False``).
    inner_op : callable, optional
        Element-wise comparison applied to each array against a type value
        (default ``lambda x, y: x == y``).
    outer_op : callable, optional
        Combines the two per-set results into a single boolean array
        (default ``lambda x, y: x * y``, i.e. logical AND).

    Returns
    -------
    numpy.ndarray
        1-D boolean mask of length ``set1.shape[0]``.
    """
    if symmetrize:
        types = list(set(types + [(t[1], t[0]) for t in types]))

    mask = np.full(set1.shape[0], False)
    for type1, type2 in types:
        mask |= outer_op(inner_op(set1, type1), inner_op(set2, type2))

    return mask


def set_sig_epsilon(meGO_LJ, parameters):
    """
    Set sigma and epsilon for all LJ interactions based on probability and distance.

    Attractive interactions are assigned when the training probability exceeds the
    adaptive threshold relative to the RC distribution. Repulsive interactions are
    assigned when training probability is above MD threshold but below the attractive
    threshold. 1-4 interactions are always kept repulsive.

    Parameters
    ----------
    meGO_LJ : pd.DataFrame
        DataFrame containing LJ parameters including probability thresholds.
    parameters : argparse.Namespace
        Parsed arguments (``learn_tolerance`` is now read from ``config``).

    Returns
    -------
    pd.DataFrame
        Updated DataFrame with epsilon and sigma columns set.
    """
    # first: all contacts are set as for the prior model
    # these contacts are not considered as learned so can be overriden
    meGO_LJ["learned"] = 0
    meGO_LJ["epsilon"] = meGO_LJ["epsilon_prior"]
    meGO_LJ["sigma"] = meGO_LJ["sigma_prior"]

    # Attractive interactions
    # These are defined only if the training probability is greater than MD_threshold and
    # by comparing them with RC_probabilities so that the resulting epsilon is between eps_min and eps_0
    condition = (
        meGO_LJ["probability"]
        > meGO_LJ["limit_rc_att"] * np.maximum(meGO_LJ["rc_probability"], meGO_LJ["rc_threshold"])
    ) & (meGO_LJ["probability"] > meGO_LJ["md_threshold"])
    with np.errstate(divide="ignore"):
        meGO_LJ.loc[condition, "epsilon"] = np.maximum(0.0, meGO_LJ["epsilon_prior"]) - (
            (meGO_LJ["epsilon_0"] - np.maximum(0.0, meGO_LJ["epsilon_prior"])) / np.log(meGO_LJ["rc_threshold"])
        ) * (np.log(meGO_LJ["probability"] / (np.maximum(meGO_LJ["rc_probability"], meGO_LJ["rc_threshold"]))))
    meGO_LJ.loc[condition, "learned"] = 1
    meGO_LJ.loc[condition, "sigma"] = meGO_LJ["distance"] / 2.0 ** (1.0 / 6.0)

    # Repulsive interactions (probability above MD threshold but below attractive threshold)
    # this is used only when MD_th < MD_p < limit_rc_att*RC_p
    condition = (
        (
            meGO_LJ["probability"]
            <= meGO_LJ["limit_rc_att"] * np.maximum(meGO_LJ["rc_probability"], meGO_LJ["rc_threshold"])
        )
        & (meGO_LJ["probability"] > meGO_LJ["md_threshold"])
        & (meGO_LJ["rc_probability"] > meGO_LJ["md_threshold"])
    )
    meGO_LJ.loc[condition, "epsilon"] = (-meGO_LJ["rep"] * (meGO_LJ["distance"] / meGO_LJ["rc_distance"]) ** 12).clip(
        lower=-20 * meGO_LJ["rep"], upper=-0.05 * meGO_LJ["rep"]
    )[condition]
    meGO_LJ.loc[condition, "learned"] = 1

    # for repulsive interactions reset sigma to its effective value for consistent merging
    meGO_LJ.loc[(meGO_LJ["epsilon"] < 0.0), "sigma"] = (-meGO_LJ["epsilon"]) ** (1.0 / 12.0) / (2.0 ** (1.0 / 6.0))

    meGO_LJ.dropna(subset=["epsilon"], inplace=True)
    meGO_LJ = meGO_LJ[meGO_LJ.epsilon != 0]

    return meGO_LJ


def apply_symmetries(meGO_ensemble, meGO_input, symmetry):
    """
    Apply symmetries to the molecular ensemble.

    Parameters
    ----------
    meGO_ensemble : dict
        A dictionary containing relevant meGO data.
    meGO_input : pd.DataFrame
        Input DataFrame containing molecular ensemble data.
    symmetry : list
        List of symmetry definitions parsed from arguments.

    Returns
    -------
    pd.DataFrame
        DataFrame with applied symmetries.
    """
    if meGO_input.empty:
        return pd.DataFrame(columns=meGO_input.columns)

    if not symmetry:
        return pd.DataFrame(columns=meGO_input.columns)

    topology_df = meGO_ensemble.topology_dataframe
    dict_sbtype_to_resname = topology_df.set_index("sb_type")["resname"].to_dict()

    # Pre-compute sb_types that belong to C-terminal (CTER) and N-terminal (NTER)
    # residues — the last and first residue of each molecule respectively.
    cter_sbtypes: set = set()
    nter_sbtypes: set = set()
    topology_df["_resnum_int"] = topology_df["resnum"].astype(int)
    for _, mol_group in topology_df.groupby("molecule_name"):
        resnums = mol_group["_resnum_int"]
        cter_sbtypes |= set(mol_group.loc[resnums == resnums.max(), "sb_type"])
        nter_sbtypes |= set(mol_group.loc[resnums == resnums.min(), "sb_type"])
    topology_df.drop(columns=["_resnum_int"], inplace=True)

    def _filter_and_masks(df, ai_col, aj_col, sym0):
        """Return (filtered_df, ai_mask, aj_mask) for a given residue specifier."""
        if sym0 == "CTER":
            terminal = cter_sbtypes
            ai_mask = df[ai_col].isin(terminal)
            aj_mask = df[aj_col].isin(terminal)
        elif sym0 == "NTER":
            terminal = nter_sbtypes
            ai_mask = df[ai_col].isin(terminal)
            aj_mask = df[aj_col].isin(terminal)
        else:
            resn_ai = df[ai_col].map(dict_sbtype_to_resname)
            resn_aj = df[aj_col].map(dict_sbtype_to_resname)
            ai_mask = resn_ai == sym0
            aj_mask = resn_aj == sym0
        return df[ai_mask | aj_mask], ai_mask, aj_mask

    df_list = []

    for sym in symmetry:
        if not sym:
            continue
        meGO_filtered, mgf_mask_ai, mgf_mask_aj = _filter_and_masks(meGO_input, "ai", "aj", sym[0])
        # Reindex masks to match the filtered frame
        mgf_mask_ai = mgf_mask_ai.reindex(meGO_filtered.index)
        mgf_mask_aj = mgf_mask_aj.reindex(meGO_filtered.index)
        for atypes in itertools.permutations(sym[1:]):
            t_df_ai = meGO_filtered[meGO_filtered["ai"].str.startswith(f"{atypes[0]}_") & mgf_mask_ai]
            t_df_aj = meGO_filtered[meGO_filtered["aj"].str.startswith(f"{atypes[0]}_") & mgf_mask_aj]
            t_df_ai.loc[:, "ai"] = t_df_ai["ai"].str.replace(r"^(.*?)_", atypes[1] + "_", regex=True)
            t_df_aj.loc[:, "aj"] = t_df_aj["aj"].str.replace(r"^(.*?)_", atypes[1] + "_", regex=True)
            df_list.extend([t_df_ai, t_df_aj])

    df_tmp = pd.concat(df_list, ignore_index=True)
    df_list = []
    df_tmp.drop_duplicates(inplace=True)

    for sym in symmetry:
        if not sym:
            continue
        df_tmp_filt, df_mask_ai, df_mask_aj = _filter_and_masks(df_tmp, "ai", "aj", sym[0])
        df_mask_ai = df_mask_ai.reindex(df_tmp_filt.index)
        df_mask_aj = df_mask_aj.reindex(df_tmp_filt.index)
        for atypes in itertools.permutations(sym[1:]):
            t_df_ai = df_tmp_filt[df_tmp_filt["ai"].str.startswith(f"{atypes[0]}_") & df_mask_ai]
            t_df_aj = df_tmp_filt[df_tmp_filt["aj"].str.startswith(f"{atypes[0]}_") & df_mask_aj]
            t_df_ai.loc[:, "ai"] = t_df_ai["ai"].str.replace(r"^(.*?)_", atypes[1] + "_", regex=True)
            t_df_aj.loc[:, "aj"] = t_df_aj["aj"].str.replace(r"^(.*?)_", atypes[1] + "_", regex=True)
            df_list.extend([t_df_ai, t_df_aj])

    tmp_df = pd.concat(df_list + [df_tmp], ignore_index=True)
    tmp_df.drop_duplicates(inplace=True)

    return tmp_df


def init_LJ_datasets(meGO_ensemble, matrices, args):
    """
    Assembles the full training dataset by merging train/reference contact matrices
    with 1-4 pair data and computing default repulsive and MG sigma/epsilon values.

    Parameters
    ----------
    meGO_ensemble : dict
        The initialized meGO ensemble.
    matrices : dict
        Contains 'reference_matrices' and 'train_matrices'.
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    train_dataset : pd.DataFrame
        Fully assembled and annotated training dataset ready for generate_LJ.
    """
    td_fields = [
        "molecule_name_ai",
        "ai",
        "molecule_name_aj",
        "aj",
        "distance",
        "probability",
        "cutoff",
        "same_chain",
        "source",
        "reference",
        "epsilon_0",
        "epsilon_prior",
        "sigma_prior",
        "md_threshold",
        "rc_threshold",
        "limit_rc_att",
        "rc_distance",
        "rc_probability",
        "learned",
    ]

    chunks = []
    for name, ref_name in meGO_ensemble.train_matrix_tuples:
        if ref_name not in matrices["reference_matrices"].keys():
            raise RuntimeError(
                f'Encountered error while trying to find {ref_name} in reference matrices {matrices["reference_matrices"].keys()}'
            )

        temp_merged = pd.merge(
            matrices["train_matrices"][name],
            matrices["reference_matrices"][ref_name],
            left_index=True,
            right_index=True,
            how="outer",
        )

        if (np.abs(temp_merged["rc_cutoff"] - temp_merged["cutoff"])).max() > 0:
            print(
                temp_merged[["ai", "aj", "rc_ai", "rc_aj", "source", "rc_source", "cutoff", "rc_cutoff", "same_chain"]]
                .loc[(np.abs(temp_merged["rc_cutoff"] - temp_merged["cutoff"]) > 0)]
                .to_string()
            )
            sys.exit("ERROR: Inconsistent cutoff values between the TRAINING and corresponding REFERENCE input data")

        # Restrict the check to contacts that appear in both matrices: the outer
        # merge leaves rc_same_chain = NaN for training contacts that have no
        # counterpart in the reference, and those rows must not trigger the error.
        # Cast to bool before comparing: after an outer join, boolean columns that
        # receive any NaN are upcast to object dtype, and Series.equals() is
        # dtype-sensitive, so object([True]) != bool([True]) even though the
        # values are identical. A plain != on bool-cast values avoids this.
        both_present = temp_merged["same_chain"].notna() & temp_merged["rc_same_chain"].notna()
        mismatched = both_present & (
            temp_merged["same_chain"].astype(bool) != temp_merged["rc_same_chain"].astype(bool)
        )
        if mismatched.any():
            diff_indices = temp_merged.index[mismatched].tolist()
            print(f"Difference found at indices: {diff_indices}")
            sys.exit("ERROR: You are pairing intra and inter molecular training and reference data")

        chunks.append(temp_merged[td_fields])

    train_dataset = (
        pd.concat(chunks, axis=0, sort=False, ignore_index=True) if chunks else pd.DataFrame(columns=td_fields)
    )

    train_dataset["molecule_name_ai"] = train_dataset["molecule_name_ai"].astype("category")
    train_dataset["molecule_name_aj"] = train_dataset["molecule_name_aj"].astype("category")
    train_dataset["source"] = train_dataset["source"].astype("category")

    all_bd = bonded.generate_bond_distance_data(meGO_ensemble)

    train_dataset = pd.merge(
        train_dataset,
        all_bd[["ai", "aj", "bond_distance"]],
        how="left",
        on=["ai", "aj"],
    )
    train_dataset["bond_distance"] = train_dataset["bond_distance"].fillna(config.max_bond_separation + 1).astype(int)

    train_dataset["ai"] = train_dataset["ai"].astype("category")
    train_dataset["aj"] = train_dataset["aj"].astype("category")

    # Remove 0_1_2_3 intramolecular interactions
    train_dataset = train_dataset[
        ~((train_dataset["same_chain"]) & (train_dataset["bond_distance"] <= config.bond14_separation))
    ]
    train_dataset.reset_index(inplace=True)

    # default repulsive C12
    train_dataset["rep"] = np.sqrt(
        train_dataset["ai"].map(meGO_ensemble.sbtype_c12_dict) * train_dataset["aj"].map(meGO_ensemble.sbtype_c12_dict)
    )

    # default (mg) sigma
    pairwise_mg_sigma = (
        train_dataset["ai"].map(meGO_ensemble.sbtype_mg_c12_dict)
        * train_dataset["aj"].map(meGO_ensemble.sbtype_mg_c12_dict)
        / (
            train_dataset["ai"].map(meGO_ensemble.sbtype_mg_c6_dict)
            * train_dataset["aj"].map(meGO_ensemble.sbtype_mg_c6_dict)
        )
    ) ** (1 / 12)
    train_dataset["mg_sigma"] = pd.Series(pairwise_mg_sigma)

    # default (mg) epsilon
    pairwise_mg_epsilon = (
        train_dataset["ai"].map(meGO_ensemble.sbtype_mg_c6_dict)
        * train_dataset["aj"].map(meGO_ensemble.sbtype_mg_c6_dict)
    ) / (
        4
        * np.sqrt(
            train_dataset["ai"].map(meGO_ensemble.sbtype_mg_c12_dict)
            * train_dataset["aj"].map(meGO_ensemble.sbtype_mg_c12_dict)
        )
    )
    train_dataset["mg_epsilon"] = pd.Series(pairwise_mg_epsilon)

    # Apply special non-local interaction rules from type_definitions.
    #
    # The naive approach iterates over every rule, computes a boolean mask by
    # scanning type_ai/type_aj for each (type_i, type_j) pair (186 scans total),
    # then writes results with .loc — slow for large datasets.
    #
    # Instead we build two tiny 2D lookup tables indexed by integer atom-type
    # codes, do a single vectorised fancy-index to get the per-row rule outcome,
    # then apply all updates with np.where on bare numpy arrays.

    # --- Step 1: assign a small integer code to every atom type in the rules ---
    # Types absent from the rules get code 0 ("no rule applies").
    all_rule_types = sorted(
        {atomtype for rule in type_definitions.special_non_local for side in rule["atomtypes"] for atomtype in side}
    )
    type_to_code = {atomtype: code for code, atomtype in enumerate(all_rule_types, start=1)}
    n_codes = len(all_rule_types) + 1  # +1 for the 0 sentinel

    # Rule-outcome codes (int8 to keep the tables tiny):
    #   1 = repulsive with a fixed epsilon value
    #   2 = repulsive, derive mg parameters from the existing rep (epsilon=None)
    #   3 = attractive
    REP_FIXED, REP_NONE, ATT = np.int8(1), np.int8(2), np.int8(3)

    # --- Step 2: fill the lookup tables (one entry per ordered type pair) ---
    lut_outcome = np.zeros((n_codes, n_codes), dtype=np.int8)
    lut_epsilon = np.full((n_codes, n_codes), np.nan)
    lut_sigma = np.full((n_codes, n_codes), np.nan)

    for rule in type_definitions.special_non_local:
        types_i, types_j = rule["atomtypes"]
        epsilon = rule["epsilon"]
        sigma = rule.get("sigma")
        if rule["interaction"] == "rep":
            outcome = REP_FIXED if epsilon is not None else REP_NONE
        else:
            outcome = ATT
        for ti in types_i:
            for tj in types_j:
                for ci, cj in ((type_to_code[ti], type_to_code[tj]), (type_to_code[tj], type_to_code[ti])):
                    lut_outcome[ci, cj] = outcome
                    lut_epsilon[ci, cj] = epsilon if epsilon is not None else np.nan
                    lut_sigma[ci, cj] = sigma if sigma is not None else np.nan

    # --- Step 3: encode every row's ai/aj type string as an integer code ---
    # Two chained pandas .map() calls; no Python loop over rows.
    ai_codes = (
        train_dataset["ai"].map(meGO_ensemble.sbtype_type_dict).map(type_to_code).fillna(0).to_numpy(dtype=np.intp)
    )
    aj_codes = (
        train_dataset["aj"].map(meGO_ensemble.sbtype_type_dict).map(type_to_code).fillna(0).to_numpy(dtype=np.intp)
    )

    # --- Step 4: single vectorised lookup — one numpy fancy-index per table ---
    row_outcome = lut_outcome[ai_codes, aj_codes]
    row_epsilon = lut_epsilon[ai_codes, aj_codes]
    row_sigma = lut_sigma[ai_codes, aj_codes]

    # --- Step 5: apply all rule outcomes in one pass with np.where on numpy arrays ---
    # Working on raw arrays avoids the per-call pandas alignment overhead of .loc.
    rep = train_dataset["rep"].to_numpy(dtype=float)
    mg_epsilon = train_dataset["mg_epsilon"].to_numpy(dtype=float)
    mg_sigma = train_dataset["mg_sigma"].to_numpy(dtype=float)

    is_rep_fixed = row_outcome == REP_FIXED
    is_rep_none = row_outcome == REP_NONE
    is_att = row_outcome == ATT

    # repulsive, fixed epsilon: override rep, mg_epsilon, and mg_sigma
    rep = np.where(is_rep_fixed, row_epsilon, rep)
    mg_epsilon = np.where(is_rep_fixed, -row_epsilon, mg_epsilon)
    mg_sigma = np.where(is_rep_fixed, row_epsilon ** (1 / 12) / 2 ** (1 / 6), mg_sigma)

    # repulsive, epsilon=None: derive mg_epsilon and mg_sigma from the existing rep
    mg_epsilon = np.where(is_rep_none, -rep, mg_epsilon)
    mg_sigma = np.where(is_rep_none, rep ** (1 / 12) / 2 ** (1 / 6), mg_sigma)

    # attractive: override mg_epsilon and/or mg_sigma where the rule specifies them
    mg_epsilon = np.where(is_att & ~np.isnan(row_epsilon), row_epsilon, mg_epsilon)
    mg_sigma = np.where(is_att & ~np.isnan(row_sigma), row_sigma, mg_sigma)

    train_dataset["rep"] = rep
    train_dataset["mg_epsilon"] = mg_epsilon
    train_dataset["mg_sigma"] = mg_sigma

    train_dataset.dropna(subset=["mg_sigma"], inplace=True)
    train_dataset = train_dataset.loc[train_dataset["rep"] > 0.0]

    if (np.abs(train_dataset["cutoff"] - 1.45 * train_dataset["rep"] ** (1 / 12))).max() > 10e-6:
        print(
            train_dataset[["ai", "aj", "source", "same_chain", "cutoff", "rep"]]
            .loc[(np.abs(train_dataset["cutoff"] - 1.45 * train_dataset["rep"] ** (1 / 12)) > 10e-6)]
            .to_string()
        )
        sys.exit("ERROR: Inconsistent cutoff and C12 repulsive values")

    return train_dataset


def generate_LJ(meGO_ensemble, train_dataset, parameters):
    """
    Generates production LJ interactions from the training dataset.

    Parameters
    ----------
    meGO_ensemble : dict
        The initialized meGO ensemble.
    train_dataset : pd.DataFrame
        Output of init_LJ_datasets.
    parameters : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    meGO_LJ : pd.DataFrame
        Non-bonded LJ interactions for ffnonbonded.itp.
    meGO_LJ_14 : pd.DataFrame
        1-4 LJ interactions for the topology pairs section.
    stat_str : str
        Statistics string for the header.
    """
    st = time.time()
    _term.sub("Set sigma and epsilon")
    meGO_LJ = train_dataset[train_dataset["learned"] == 1].copy()

    needed_fields = [
        "molecule_name_ai",
        "ai",
        "molecule_name_aj",
        "aj",
        "probability",
        "same_chain",
        "source",
        "reference",
        "rc_probability",
        "sigma",
        "epsilon",
        "bond_distance",
        "rep",
        "sigma_prior",
        "epsilon_prior",
        "mg_sigma",
        "mg_epsilon",
        "md_threshold",
        "rc_threshold",
        "learned",
    ]

    meGO_LJ = set_sig_epsilon(meGO_LJ, parameters)[needed_fields]

    et = time.time()
    _term.timing(et - st)
    st = et

    if parameters.symmetry:
        _term.sub("Apply the defined atomic symmetries")
        meGO_LJ_sym = apply_symmetries(meGO_ensemble, meGO_LJ, parameters.symmetry)
        meGO_LJ = pd.concat([meGO_LJ, meGO_LJ_sym])
        meGO_LJ.reset_index(inplace=True)
        et = time.time()
        _term.timing(et - st)
        st = et

    _term.sub("Merging multiple states (training, symmetries, inter/intra)")

    # Merging priority: learned > not learned, attractive > repulsive, shorter > longer, stronger attractive, weaker repulsive
    meGO_LJ["type"] = np.sign(meGO_LJ["epsilon"])
    meGO_LJ.sort_values(
        by=["ai", "aj", "same_chain", "learned", "type", "sigma", "epsilon"],
        ascending=[True, True, True, False, False, True, False],
        inplace=True,
    )
    meGO_LJ = meGO_LJ.drop_duplicates(subset=["ai", "aj", "same_chain"], keep="first")

    # Remove contacts that are effectively identical to the default MG prior
    meGO_LJ = meGO_LJ.loc[
        ~(
            (meGO_LJ["epsilon"] > 0)
            & (meGO_LJ["mg_epsilon"] > 0)
            & ((abs(meGO_LJ["epsilon"] - meGO_LJ["mg_epsilon"]) / meGO_LJ["mg_epsilon"]) < config.learn_tolerance)
            & ((abs(meGO_LJ["sigma"] - meGO_LJ["mg_sigma"]) / meGO_LJ["mg_sigma"]) < config.learn_tolerance)
            & ((meGO_LJ["bond_distance"] > config.bond14_separation) | (~meGO_LJ["same_chain"]))
        )
    ]
    meGO_LJ = meGO_LJ.loc[
        ~(
            (meGO_LJ["epsilon"] < 0)
            & (meGO_LJ["mg_epsilon"] < 0)
            & ((abs(meGO_LJ["epsilon"] - meGO_LJ["mg_epsilon"]) / abs(meGO_LJ["mg_epsilon"])) < config.learn_tolerance)
            & ((meGO_LJ["bond_distance"] > config.bond14_separation) | (~meGO_LJ["same_chain"]))
            & ~((meGO_LJ["bond_distance"] <= config.max_bond_separation) & (meGO_LJ["same_chain"]))
        )
    ]

    meGO_LJ = meGO_LJ[needed_fields]

    stat_str = io.print_stats(meGO_LJ)

    # --- identify duplicates ON ORIGINAL DATA ---
    dup_mask = meGO_LJ.duplicated(subset=["ai", "aj"], keep=False)
    # --- pairs section (ONLY intramolecular interactions from duplicated pairs) ---
    dup_df = meGO_LJ.loc[dup_mask].copy()
    meGO_LJ_14 = (
        dup_df.assign(_priority=~dup_df["same_chain"])
        .sort_values(["ai", "aj", "_priority"])
        .drop_duplicates(subset=["ai", "aj"], keep="first")
        .drop(columns="_priority")
    )
    # ffnonbonded: prioritise intermolecular interactions on duplicate ai/aj
    meGO_LJ.sort_values(by=["ai", "aj", "same_chain"], ascending=[True, True, True], inplace=True)
    meGO_LJ = meGO_LJ.drop_duplicates(subset=["ai", "aj"], keep="first")

    # no cross interactions
    common = meGO_LJ_14["molecule_name_ai"].cat.categories.union(meGO_LJ_14["molecule_name_aj"].cat.categories)
    meGO_LJ_14["molecule_name_ai"] = meGO_LJ_14["molecule_name_ai"].cat.set_categories(common)
    meGO_LJ_14["molecule_name_aj"] = meGO_LJ_14["molecule_name_aj"].cat.set_categories(common)
    meGO_LJ_14 = meGO_LJ_14[meGO_LJ_14["molecule_name_ai"] == meGO_LJ_14["molecule_name_aj"]]
    # intramolecular interactions within few bonds should be move in pairs
    copy_intra = meGO_LJ.loc[(meGO_LJ["same_chain"]) & (meGO_LJ["bond_distance"] <= config.max_bond_separation)]
    meGO_LJ_14 = pd.concat([meGO_LJ_14, copy_intra], axis=0, sort=False, ignore_index=True)
    meGO_LJ = meGO_LJ.loc[~((meGO_LJ["same_chain"]) & (meGO_LJ["bond_distance"] <= config.max_bond_separation))]

    if parameters.force_split:
        split_ii = meGO_LJ.loc[(meGO_LJ["same_chain"])]
        meGO_LJ_14 = pd.concat([meGO_LJ_14, split_ii], axis=0, sort=False, ignore_index=True)
        meGO_LJ = meGO_LJ.loc[(~meGO_LJ["same_chain"])]

    # Add MG prior interactions for pairs not learned in any other way
    needed_fields_mg = [
        "molecule_name_ai",
        "ai",
        "molecule_name_aj",
        "aj",
        "probability",
        "same_chain",
        "source",
        "reference",
        "rc_probability",
        "sigma",
        "epsilon",
        "bond_distance",
        "rep",
        "mg_sigma",
        "mg_epsilon",
        "md_threshold",
        "rc_threshold",
        "learned",
    ]
    basic_LJ = mg.generate_MG_LJ(meGO_ensemble)[needed_fields_mg]
    meGO_LJ = pd.concat([meGO_LJ, basic_LJ])

    # Symmetrize
    inverse_meGO_LJ = meGO_LJ.rename(
        columns={"ai": "aj", "aj": "ai", "molecule_name_ai": "molecule_name_aj", "molecule_name_aj": "molecule_name_ai"}
    ).copy()
    meGO_LJ = pd.concat([meGO_LJ, inverse_meGO_LJ], axis=0, sort=False, ignore_index=True)

    meGO_LJ["ai"] = meGO_LJ["ai"].astype("category")
    meGO_LJ["aj"] = meGO_LJ["aj"].astype("category")
    meGO_LJ["molecule_name_ai"] = meGO_LJ["molecule_name_ai"].astype("category")
    meGO_LJ["molecule_name_aj"] = meGO_LJ["molecule_name_aj"].astype("category")

    meGO_LJ.sort_values(by=["ai", "aj", "learned", "sigma"], ascending=[True, True, False, True], inplace=True)
    meGO_LJ = meGO_LJ.drop_duplicates(subset=["ai", "aj"], keep="first")

    meGO_LJ["c6"] = np.where(meGO_LJ["epsilon"] < 0.0, 0.0, 4 * meGO_LJ["epsilon"] * (meGO_LJ["sigma"] ** 6))
    meGO_LJ["c12"] = np.where(
        meGO_LJ["epsilon"] < 0.0, -meGO_LJ["epsilon"], 4 * meGO_LJ["epsilon"] * (meGO_LJ["sigma"] ** 12)
    )

    meGO_LJ_14["c6"] = np.where(
        meGO_LJ_14["epsilon"] < 0.0, 0.0, 4 * meGO_LJ_14["epsilon"] * (meGO_LJ_14["sigma"] ** 6)
    )
    meGO_LJ_14["c12"] = np.where(
        meGO_LJ_14["epsilon"] < 0.0, -meGO_LJ_14["epsilon"], 4 * meGO_LJ_14["epsilon"] * (meGO_LJ_14["sigma"] ** 12)
    )

    et = time.time()
    _term.timing(et - st)

    return meGO_LJ, meGO_LJ_14, stat_str
