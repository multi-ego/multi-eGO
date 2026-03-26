from . import type_definitions
from .mg import generate_MG_LJ
from .bonded import generate_bond_distance_data
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

    Args:
    - set1 (numpy.ndarray): First set of values.
    - set2 (numpy.ndarray): Second set of values.
    - types (list): List of tuples representing types for comparison.
    - symmetrize (bool): Flag to determine whether to symmetrize type selection.
    - inner_op (function): Operation for element-wise comparison within each set.
    - outer_op (function): Operation for the combination of comparisons between sets.

    Returns:
    - numpy.ndarray: A linearized boolean mask reflecting the comparison operations.

    Note:
    - This function does not provide a flattened mask array but operates on a 1D array.
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
        Parsed arguments, used for epsilon_min.

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
        (meGO_LJ["probability"] > meGO_LJ["limit_rc_att"] * np.maximum(meGO_LJ["rc_probability"], meGO_LJ["rc_threshold"]))
        & (meGO_LJ["probability"] > meGO_LJ["md_threshold"])
        # & (meGO_LJ["mg_epsilon"] > 0.0)
    )
    meGO_LJ.loc[condition, "epsilon"] = np.maximum(0.0, meGO_LJ["epsilon_prior"]) - (
        (meGO_LJ["epsilon_0"] - np.maximum(0.0, meGO_LJ["epsilon_prior"])) / np.log(meGO_LJ["rc_threshold"])
    ) * (np.log(meGO_LJ["probability"] / (np.maximum(meGO_LJ["rc_probability"], meGO_LJ["rc_threshold"]))))
    meGO_LJ.loc[condition, "learned"] = 1
    meGO_LJ.loc[condition, "sigma"] = meGO_LJ["distance"] / 2.0 ** (1.0 / 6.0)

    # this is used only when MD_th < MD_p < limit_rc_att*RC_p
    # condition = (
    #    (meGO_LJ["probability"] <= meGO_LJ["limit_rc_att"] * np.maximum(meGO_LJ["rc_probability"], meGO_LJ["rc_threshold"]))
    #    & (meGO_LJ["probability"] > meGO_LJ["md_threshold"])
    #    & (meGO_LJ["rc_probability"] > meGO_LJ["md_threshold"])
    #    & (meGO_LJ["mg_epsilon"] > 0.0)
    # )
    # meGO_LJ.loc[condition, "epsilon"] = parameters.epsilon_min
    # meGO_LJ.loc[condition, "learned"] = 1
    # meGO_LJ.loc[condition, "sigma"] = meGO_LJ["distance"] / 2.0 ** (1.0 / 6.0)

    # Repulsive interactions
    # negative epsilon are used to identify non-attractive interactions
    # These are defined only if the training probability is greater than MD_threshold and
    # by comparing them with RC_probabilities so that the resulting epsilon is between eps_min and eps_0
    # condition = (
    #    (meGO_LJ["probability"] > meGO_LJ["limit_rc_att"] * np.maximum(meGO_LJ["rc_probability"], meGO_LJ["rc_threshold"]))
    #    & (meGO_LJ["probability"] > meGO_LJ["md_threshold"])
    #    & (meGO_LJ["mg_epsilon"] < 0.0)
    #    & (meGO_LJ["rc_probability"] > meGO_LJ["md_threshold"])
    # )
    # meGO_LJ.loc[condition, "epsilon"] = (-meGO_LJ["rep"] * (meGO_LJ["distance"] / meGO_LJ["rc_distance"]) ** 12).clip(
    #    lower=-20 * meGO_LJ["rep"], upper=-0.05 * meGO_LJ["rep"]
    # )[condition]
    # meGO_LJ.loc[condition, "learned"] = 1

    # Repulsive interactions (probability above MD threshold but below attractive threshold)
    # this is used only when MD_th < MD_p < limit_rc_att*RC_p
    condition = (
        (meGO_LJ["probability"] <= meGO_LJ["limit_rc_att"] * np.maximum(meGO_LJ["rc_probability"], meGO_LJ["rc_threshold"]))
        & (meGO_LJ["probability"] > meGO_LJ["md_threshold"])
        & (meGO_LJ["rc_probability"] > meGO_LJ["md_threshold"])
        # & (meGO_LJ["mg_epsilon"] < 0.0)
    )
    meGO_LJ.loc[condition, "epsilon"] = (-meGO_LJ["rep"] * (meGO_LJ["distance"] / meGO_LJ["rc_distance"]) ** 12).clip(
        lower=-20 * meGO_LJ["rep"], upper=-0.05 * meGO_LJ["rep"]
    )[condition]
    meGO_LJ.loc[condition, "learned"] = 1

    # 1-4 interactions are part of bonded interactions and cannot become attractive
    condition = (meGO_LJ["bond_distance"] <= config.bond14_separation) & (meGO_LJ["same_chain"])
    meGO_LJ.loc[condition, "epsilon"] = -meGO_LJ["rep"]
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
    dict_sbtype_to_resname = meGO_ensemble.topology_dataframe.set_index("sb_type")["resname"].to_dict()
    mglj_resn_ai = meGO_input["ai"].map(dict_sbtype_to_resname)
    mglj_resn_aj = meGO_input["aj"].map(dict_sbtype_to_resname)
    df_list = []

    for sym in symmetry:
        if not sym:
            continue
        meGO_filtered = meGO_input[(mglj_resn_ai == sym[0]) | (mglj_resn_aj == sym[0])]
        mgf_resn_ai = meGO_filtered["ai"].map(dict_sbtype_to_resname)
        mgf_resn_aj = meGO_filtered["aj"].map(dict_sbtype_to_resname)
        for atypes in itertools.permutations(sym[1:]):
            t_df_ai = meGO_filtered[meGO_filtered["ai"].str.startswith(f"{atypes[0]}_") & (mgf_resn_ai == sym[0])]
            t_df_aj = meGO_filtered[meGO_filtered["aj"].str.startswith(f"{atypes[0]}_") & (mgf_resn_aj == sym[0])]
            t_df_ai.loc[:, "ai"] = t_df_ai["ai"].str.replace(r"^(.*?)_", atypes[1] + "_", regex=True)
            t_df_aj.loc[:, "aj"] = t_df_aj["aj"].str.replace(r"^(.*?)_", atypes[1] + "_", regex=True)
            df_list.extend([t_df_ai, t_df_aj])

    df_tmp = pd.concat(df_list, ignore_index=True)
    df_list = []
    df_tmp.drop_duplicates(inplace=True)

    df_resn_ai = df_tmp["ai"].map(dict_sbtype_to_resname)
    df_resn_aj = df_tmp["aj"].map(dict_sbtype_to_resname)

    for sym in symmetry:
        if not sym:
            continue
        df_tmp_filt = df_tmp[(df_resn_ai == sym[0]) | (df_resn_aj == sym[0])]
        df_resn_ai_f = df_tmp_filt["ai"].map(dict_sbtype_to_resname)
        df_resn_aj_f = df_tmp_filt["aj"].map(dict_sbtype_to_resname)
        for atypes in itertools.permutations(sym[1:]):
            t_df_ai = df_tmp_filt[df_tmp_filt["ai"].str.startswith(f"{atypes[0]}_") & (df_resn_ai_f == sym[0])]
            t_df_aj = df_tmp_filt[df_tmp_filt["aj"].str.startswith(f"{atypes[0]}_") & (df_resn_aj_f == sym[0])]
            t_df_ai.loc[:, "ai"] = t_df_ai["ai"].str.replace(r"^(.*?)_", atypes[1] + "_", regex=True)
            t_df_aj.loc[:, "aj"] = t_df_aj["aj"].str.replace(r"^(.*?)_", atypes[1] + "_", regex=True)
            df_list.extend([t_df_ai, t_df_aj])

    tmp_df = pd.concat(df_list + [df_tmp], ignore_index=True)
    tmp_df.drop_duplicates(inplace=True)

    return tmp_df


def init_LJ_datasets(meGO_ensemble, matrices, pairs14, args):
    """
    Assembles the full training dataset by merging train/reference contact matrices
    with 1-4 pair data and computing default repulsive and MG sigma/epsilon values.

    Parameters
    ----------
    meGO_ensemble : dict
        The initialized meGO ensemble.
    matrices : dict
        Contains 'reference_matrices' and 'train_matrices'.
    pairs14 : pd.DataFrame
        1-4 pair interactions from generate_14_data.
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    train_dataset : pd.DataFrame
        Fully assembled and annotated training dataset ready for generate_LJ.
    """
    train_dataset = pd.DataFrame()

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
            sys.exit(
                "HERE SOMETHING BAD HAPPEND: There are inconsistent cutoff values between the TRAINING and corresponding REFERENCE input data"
            )

        if not temp_merged["rc_same_chain"].equals(temp_merged["rc_same_chain"]):
            diff_indices = temp_merged.index[temp_merged["same_chain"] != temp_merged["rc_same_chain"]].tolist()
            print(f"Difference found at indices: {diff_indices}")
            sys.exit("HERE SOMETHING BAD HAPPEND: You are pairing intra and inter molecular training and reference data")

        temp_merged = temp_merged[td_fields]
        train_dataset = pd.concat([train_dataset, temp_merged], axis=0, sort=False, ignore_index=True)

    train_dataset["molecule_name_ai"] = train_dataset["molecule_name_ai"].astype("category")
    train_dataset["molecule_name_aj"] = train_dataset["molecule_name_aj"].astype("category")
    train_dataset["source"] = train_dataset["source"].astype("category")

    all_bd = generate_bond_distance_data(meGO_ensemble)

    train_dataset = pd.merge(
        pd.merge(
            train_dataset,
            pairs14[["ai", "aj", "same_chain", "rep"]],
            how="left",
            on=["ai", "aj", "same_chain"],
        ),
        all_bd[["ai", "aj", "bond_distance"]],
        how="left",
        on=["ai", "aj"],
    )
    train_dataset["bond_distance"] = train_dataset["bond_distance"].fillna(config.max_bond_separation + 1).astype(int)

    train_dataset["ai"] = train_dataset["ai"].astype("category")
    train_dataset["aj"] = train_dataset["aj"].astype("category")

    # Remove 0_1_2_3 intramolecular interactions
    train_dataset = train_dataset[
        ~((train_dataset["same_chain"]) & (train_dataset["bond_distance"] < config.bond14_separation))
    ]
    train_dataset.reset_index(inplace=True)

    train_dataset.loc[
        (train_dataset["bond_distance"] == config.bond14_separation)
        & (train_dataset["same_chain"])
        & (train_dataset["rep"].isnull()),
        "rep",
    ] = 0.0

    type_to_c12 = {key: val for key, val in zip(type_definitions.gromos_atp.name, type_definitions.gromos_atp.rc_c12)}

    if args.custom_c12 is not None:
        from . import io as _io

        custom_c12_dict = _io.read_custom_c12_parameters(args.custom_c12)
        type_to_c12_appo = {key: val for key, val in zip(custom_c12_dict.name, custom_c12_dict.c12)}
        type_to_c12.update(type_to_c12_appo)

    # default repulsive C12
    pairwise_c12 = np.sqrt(
        train_dataset["ai"].map(meGO_ensemble.sbtype_c12_dict) * train_dataset["aj"].map(meGO_ensemble.sbtype_c12_dict)
    )
    train_dataset["rep"] = train_dataset["rep"].fillna(pd.Series(pairwise_c12))

    # default (mg) sigma
    pairwise_mg_sigma = (
        train_dataset["ai"].map(meGO_ensemble.sbtype_mg_c12_dict)
        * train_dataset["aj"].map(meGO_ensemble.sbtype_mg_c12_dict)
        / (train_dataset["ai"].map(meGO_ensemble.sbtype_mg_c6_dict) * train_dataset["aj"].map(meGO_ensemble.sbtype_mg_c6_dict))
    ) ** (1 / 12)
    train_dataset["mg_sigma"] = pd.Series(pairwise_mg_sigma)

    # default (mg) epsilon
    pairwise_mg_epsilon = (
        train_dataset["ai"].map(meGO_ensemble.sbtype_mg_c6_dict) * train_dataset["aj"].map(meGO_ensemble.sbtype_mg_c6_dict)
    ) / (
        4
        * np.sqrt(
            train_dataset["ai"].map(meGO_ensemble.sbtype_mg_c12_dict)
            * train_dataset["aj"].map(meGO_ensemble.sbtype_mg_c12_dict)
        )
    )
    train_dataset["mg_epsilon"] = pd.Series(pairwise_mg_epsilon)

    # special cases from type_definitions
    type_ai = train_dataset["ai"].map(meGO_ensemble.sbtype_type_dict).to_numpy()
    type_aj = train_dataset["aj"].map(meGO_ensemble.sbtype_type_dict).to_numpy()

    for rule in type_definitions.special_non_local:
        types_i, types_j = rule["atomtypes"]
        pair_list = [(i, j) for i in types_i for j in types_j]
        mask = create_linearized_mask(type_ai, type_aj, pair_list, symmetrize=True)

        if rule["interaction"] == "rep":
            if rule["epsilon"] is not None:
                train_dataset.loc[mask, "rep"] = rule["epsilon"]
                train_dataset.loc[mask, "mg_epsilon"] = -rule["epsilon"]
                train_dataset.loc[mask, "mg_sigma"] = rule["epsilon"] ** (1 / 12) / 2 ** (1 / 6)
            else:
                rep_vals = train_dataset.loc[mask, "rep"]
                train_dataset.loc[mask, "mg_epsilon"] = -rep_vals
                train_dataset.loc[mask, "mg_sigma"] = rep_vals ** (1 / 12) / 2 ** (1 / 6)
        elif rule["interaction"] == "att":
            if rule["epsilon"] is not None:
                train_dataset.loc[mask, "mg_epsilon"] = rule["epsilon"]
            if rule["sigma"] is not None:
                train_dataset.loc[mask, "mg_sigma"] = rule["sigma"]

    train_dataset.dropna(subset=["mg_sigma"], inplace=True)
    train_dataset = train_dataset.loc[train_dataset["rep"] > 0.0]

    if (np.abs(train_dataset["cutoff"] - 1.45 * train_dataset["rep"] ** (1 / 12))).max() > 10e-6:
        print(
            train_dataset[["ai", "aj", "source", "same_chain", "cutoff", "rep"]]
            .loc[(np.abs(train_dataset["cutoff"] - 1.45 * train_dataset["rep"] ** (1 / 12)) > 10e-6)]
            .to_string()
        )
        sys.exit("HERE SOMETHING BAD HAPPEND: There are inconsistent cutoff and C12 repulsive values")

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
    print("\t- Set sigma and epsilon")
    meGO_LJ = train_dataset[train_dataset["learned"]].copy()

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
    elapsed_time = et - st
    st = et
    print("\t- Done in:", elapsed_time, "seconds")

    if parameters.symmetry:
        print("\t- Apply the defined atomic symmetries")
        meGO_LJ_sym = apply_symmetries(meGO_ensemble, meGO_LJ, parameters.symmetry)
        meGO_LJ = pd.concat([meGO_LJ, meGO_LJ_sym])
        meGO_LJ.reset_index(inplace=True)
        et = time.time()
        elapsed_time = et - st
        st = et
        print("\t- Done in:", elapsed_time, "seconds")

    print("\t- Merging multiple states (training, symmetries, inter/intra)")

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
            & ((abs(meGO_LJ["epsilon"] - meGO_LJ["mg_epsilon"]) / meGO_LJ["mg_epsilon"]) < parameters.relative_c12d)
            & ((abs(meGO_LJ["sigma"] - meGO_LJ["mg_sigma"]) / meGO_LJ["mg_sigma"]) < parameters.relative_c12d)
            & ((meGO_LJ["bond_distance"] > config.bond14_separation) | (~meGO_LJ["same_chain"]))
        )
    ]
    meGO_LJ = meGO_LJ.loc[
        ~(
            (meGO_LJ["epsilon"] < 0)
            & (meGO_LJ["mg_epsilon"] < 0)
            & ((abs(meGO_LJ["epsilon"] - meGO_LJ["mg_epsilon"]) / abs(meGO_LJ["mg_epsilon"])) < parameters.relative_c12d)
            & ((meGO_LJ["bond_distance"] > config.bond14_separation) | (~meGO_LJ["same_chain"]))
            & ~((meGO_LJ["bond_distance"] <= config.max_bond_separation) & (meGO_LJ["same_chain"]))
        )
    ]

    meGO_LJ = meGO_LJ[needed_fields]

    stat_str = io.print_stats(meGO_LJ)

    meGO_LJ_14 = meGO_LJ.copy()

    # ffnonbonded: prioritise intermolecular interactions on duplicate ai/aj
    meGO_LJ.sort_values(by=["ai", "aj", "same_chain"], ascending=[True, True, True], inplace=True)
    meGO_LJ = meGO_LJ.drop_duplicates(subset=["ai", "aj"], keep="first")

    # pairs section: prioritise intramolecular interactions on duplicate ai/aj
    meGO_LJ_14.sort_values(by=["ai", "aj", "same_chain"], ascending=[True, True, False], inplace=True)
    meGO_LJ_14 = meGO_LJ_14.drop_duplicates(subset=["ai", "aj"], keep="first")

    test = pd.merge(meGO_LJ_14, meGO_LJ, how="right", on=["ai", "aj"])
    meGO_LJ_14 = test.loc[(~test["same_chain_x"]) | ((test["same_chain_x"]) & (~test["same_chain_y"]))]
    meGO_LJ_14 = meGO_LJ_14.loc[:, ~meGO_LJ_14.columns.str.endswith("_y")]
    meGO_LJ_14.columns = meGO_LJ_14.columns.str.rstrip("_x")

    meGO_LJ_14 = meGO_LJ_14[meGO_LJ_14["molecule_name_ai"] == meGO_LJ_14["molecule_name_aj"]]

    copy14 = meGO_LJ.loc[(meGO_LJ["bond_distance"] == config.bond14_separation) & (meGO_LJ["same_chain"])]
    meGO_LJ_14 = pd.concat([meGO_LJ_14, copy14], axis=0, sort=False, ignore_index=True)
    meGO_LJ = meGO_LJ.loc[~((meGO_LJ["bond_distance"] == config.bond14_separation) & (meGO_LJ["same_chain"]))]

    if not parameters.single_molecule:
        copy_intra = meGO_LJ.loc[(meGO_LJ["same_chain"]) & (meGO_LJ["bond_distance"] <= config.max_bond_separation)]
        meGO_LJ_14 = pd.concat([meGO_LJ_14, copy_intra], axis=0, sort=False, ignore_index=True)
        meGO_LJ = meGO_LJ.loc[~((meGO_LJ["same_chain"]) & (meGO_LJ["bond_distance"] <= config.max_bond_separation))]

    if not parameters.force_split:
        meGO_LJ_14 = meGO_LJ_14.loc[
            ~(
                (~meGO_LJ_14["same_chain"])
                & (meGO_LJ_14["molecule_name_ai"] == meGO_LJ_14["molecule_name_aj"])
                & (meGO_LJ_14["epsilon"] > 0.0)
                & (meGO_LJ_14["bond_distance"] > config.max_bond_separation)
            )
        ]
    else:
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
    basic_LJ = generate_MG_LJ(meGO_ensemble)[needed_fields_mg]
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
    meGO_LJ["c12"] = np.where(meGO_LJ["epsilon"] < 0.0, -meGO_LJ["epsilon"], 4 * meGO_LJ["epsilon"] * (meGO_LJ["sigma"] ** 12))

    meGO_LJ_14["c6"] = np.where(meGO_LJ_14["epsilon"] < 0.0, 0.0, 4 * meGO_LJ_14["epsilon"] * (meGO_LJ_14["sigma"] ** 6))
    meGO_LJ_14["c12"] = np.where(
        meGO_LJ_14["epsilon"] < 0.0, -meGO_LJ_14["epsilon"], 4 * meGO_LJ_14["epsilon"] * (meGO_LJ_14["sigma"] ** 12)
    )

    et = time.time()
    elapsed_time = et - st
    print("\t- Done in:", elapsed_time, "seconds")

    return meGO_LJ, meGO_LJ_14, stat_str
