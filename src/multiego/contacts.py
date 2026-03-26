from . import io
from . import contacts_init

import numpy as np
import pandas as pd
import parmed
import warnings
import os
import sys
import time


def check_intra_domain_complementarity(matrices):
    """
    Checks that each non-single reference matrix (double intramat_1_1 from two different
    references) have non overlapping learning flags.
    """
    mats_names = []
    for key, _ in matrices.items():
        mats_names.append("_".join(key.split("_")[-3:]))
    to_check_names = list(set([a for a in mats_names if mats_names.count(a) > 1]))
    to_check = [[k for k in matrices.keys() if "_".join(k.split("_")[-3:]) in check] for check in to_check_names]
    for check in to_check:
        intra_flags = []
        for key in check:
            intra_flags.append(matrices[key]["rc_learned"].to_numpy())
        if np.any(np.sum(intra_flags, axis=0) > 1):
            raise ValueError(f"Learning flag complementarity not satisfied for {check} (e.g. intra-inter domain splitting)")


def initialize_molecular_contacts(contact_matrix, prior_matrix, args, reference):
    """
    Initializes a contact matrix for a given simulation.
    """
    # remove un-learned contacts (intra-inter domain)
    # contact_matrix["learned"] = prior_matrix["rc_learned"].to_numpy()
    contact_matrix["reference"] = reference["reference"]
    # calculate adaptive rc/md threshold
    p_sort = np.sort(contact_matrix["probability"].loc[(contact_matrix["learned"])].to_numpy())[::-1]
    norm = np.sum(p_sort)
    if norm == 0:
        p_sort_normalized = 0
        md_threshold = 1
    else:
        p_sort_normalized = np.cumsum(p_sort) / norm
        md_threshold = p_sort[np.min(np.where(p_sort_normalized > args.p_to_learn)[0])]

    contact_matrix["epsilon_0"] = reference["epsilon"]
    contact_matrix["md_threshold"] = md_threshold
    contact_matrix["rc_threshold"] = contact_matrix["md_threshold"] ** (
        (contact_matrix["epsilon_0"] - np.maximum(0, prior_matrix["epsilon_prior"]))
        / (contact_matrix["epsilon_0"] - args.epsilon_min)
    )
    contact_matrix["limit_rc_att"] = contact_matrix["rc_threshold"] ** (
        (np.maximum(0, prior_matrix["epsilon_prior"]) - args.epsilon_min)
        / (contact_matrix["epsilon_0"] - np.maximum(0, prior_matrix["epsilon_prior"]))
    )
    # this is for 0 : + eps
    # contact_matrix["limit_rc_rep"] = contact_matrix["rc_threshold"] ** (
    #    (np.maximum(0, prior_matrix["epsilon_prior"]))
    #    / (contact_matrix["epsilon_0"] - np.maximum(0, prior_matrix["epsilon_prior"]))
    # )
    # this is for -eps : + eps
    # contact_matrix["limit_rc_rep"] = contact_matrix["rc_threshold"] ** (
    #     (np.maximum(0, prior_matrix["epsilon_prior"]) + args.epsilon_min)
    #     / (contact_matrix["epsilon_0"] - np.maximum(0, prior_matrix["epsilon_prior"]))
    # )

    # modify limit_rc_att in the cases where epsilon_prior is negative and limit_rc_att is below 1
    contact_matrix.loc[(contact_matrix["limit_rc_att"] < 1) & (prior_matrix["epsilon_prior"] < 0), "limit_rc_att"] = 1

    return contact_matrix


def _path_to_matrix_name(path, root_dir):
    """Converts a full matrix file path to a flat unique name used as dict key."""
    name = path.replace(f"{root_dir}/inputs/", "")
    name = name.replace("/", "_")
    name = name.replace(".ndx", "")
    name = name.replace(".gz", "")
    name = name.replace(".h5", "")
    return name


def init_meGO_matrices(ensemble, args, custom_dict):
    """
    Initializes meGO contact matrices.

    Reads reference and training contact matrices for all input references,
    computes prior sigma/epsilon values, and pairs each training matrix with
    its corresponding reference matrix.

    Parameters
    ----------
    ensemble : dict
        The initialized meGO ensemble (from topology_init.init_meGO_ensemble).
    args : argparse.Namespace
        Parsed command-line / config-file arguments.
    custom_dict : dict
        Custom atom-name mapping dictionary.

    Returns
    -------
    ensemble : dict
        Updated ensemble with train_matrix_tuples populated.
    matrices : dict
        Contains 'reference_matrices' and 'train_matrices' keyed by unique names.
    """
    st = time.time()
    reference_contact_matrices = {}
    matrices = {}

    for reference in args.input_refs:
        print("\t-", f"Initializing {reference['reference']} ensemble data")
        reference_path = f"{args.root_dir}/inputs/{args.system}/{reference['reference']}"
        topol_files = [f for f in os.listdir(reference_path) if ".top" in f]
        if len(topol_files) > 1:
            raise RuntimeError(f"More than 1 topology file found in {reference_path}. Only one should be used")

        topology_path = f"{reference_path}/{topol_files[0]}"
        if not os.path.isfile(topology_path):
            raise FileNotFoundError(f"{topology_path} not found.")

        print("\t\t-", f"Reading {topology_path}")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            topol = parmed.load_file(topology_path)

        lj_data = contacts_init.get_lj_params(topol)
        lj_pairs = contacts_init.get_lj_pairs(topol)
        reversed_lj_pairs = lj_pairs.rename(columns={"ai": "aj", "aj": "ai"})
        symmetric_lj_pairs = pd.concat([lj_pairs, reversed_lj_pairs])
        symmetric_lj_pairs = symmetric_lj_pairs.drop_duplicates(subset=["ai", "aj"]).reset_index(drop=True)

        lj14_pairs = contacts_init.get_lj14_pairs(topol)
        reversed_lj14_pairs = lj14_pairs.rename(columns={"ai": "aj", "aj": "ai"})
        symmetric_lj14_pairs = pd.concat([lj14_pairs, reversed_lj14_pairs])
        symmetric_lj14_pairs = symmetric_lj14_pairs.drop_duplicates(subset=["ai", "aj"]).reset_index(drop=True)

        lj_data_dict = {str(key): val for key, val in zip(lj_data["ai"], lj_data[["c6", "c12"]].values)}

        ensemble["topology_dataframe"]["c6"] = lj_data["c6"].to_numpy()
        ensemble["topology_dataframe"]["c12"] = lj_data["c12"].to_numpy()

        matrix_paths = [f"{reference_path}/{a}" for a in os.listdir(reference_path) if reference["matrix"] in a]

        if len(matrix_paths) > 1:
            raise ValueError(f"More than 1 matrix found in {reference_path}: {matrix_paths}")
        if matrix_paths == []:
            raise FileNotFoundError(
                f"Contact matrix file(s) must be named as intramat_X_X.ndx(.gz/.h5) or intermat_X_Y.ndx(.gz/.h5). Found instead: {reference['matrix']}"
            )

        path = matrix_paths[0]
        name = _path_to_matrix_name(path, args.root_dir)
        reference_contact_matrices[name] = io.read_molecular_contacts(
            path, ensemble["molecules_idx_sbtype_dictionary"], reference["reference"], path.endswith(".h5")
        )

        reference_contact_matrices[name] = reference_contact_matrices[name].add_prefix("rc_")
        reference_contact_matrices[name]["c6_i"] = [lj_data_dict[x][0] for x in reference_contact_matrices[name]["rc_ai"]]
        reference_contact_matrices[name]["c6_j"] = [lj_data_dict[x][0] for x in reference_contact_matrices[name]["rc_aj"]]
        reference_contact_matrices[name]["c6"] = np.sqrt(
            reference_contact_matrices[name]["c6_i"] * reference_contact_matrices[name]["c6_j"]
        )
        reference_contact_matrices[name]["c12_i"] = [lj_data_dict[x][1] for x in reference_contact_matrices[name]["rc_ai"]]
        reference_contact_matrices[name]["c12_j"] = [lj_data_dict[x][1] for x in reference_contact_matrices[name]["rc_aj"]]
        reference_contact_matrices[name]["c12"] = np.sqrt(
            reference_contact_matrices[name]["c12_i"] * reference_contact_matrices[name]["c12_j"]
        )
        reference_contact_matrices[name]["sigma_prior"] = np.where(
            reference_contact_matrices[name]["c6"] > 0,
            (reference_contact_matrices[name]["c12"] / reference_contact_matrices[name]["c6"]) ** (1 / 6),
            reference_contact_matrices[name]["c12"] ** (1 / 12) / (2.0 ** (1.0 / 6.0)),
        )
        reference_contact_matrices[name]["epsilon_prior"] = np.where(
            reference_contact_matrices[name]["c6"] > 0,
            reference_contact_matrices[name]["c6"] ** 2 / (4 * reference_contact_matrices[name]["c12"]),
            -reference_contact_matrices[name]["c12"],
        )

        lj_sigma_map = symmetric_lj_pairs.set_index(["ai", "aj"])["sigma"]
        lj_epsilon_map = symmetric_lj_pairs.set_index(["ai", "aj"])["epsilon"]
        lj14_sigma_map = symmetric_lj14_pairs.set_index(["ai", "aj"])["sigma"]
        lj14_epsilon_map = symmetric_lj14_pairs.set_index(["ai", "aj"])["epsilon"]

        common_indices = lj_sigma_map.index.intersection(reference_contact_matrices[name].set_index(["rc_ai", "rc_aj"]).index)
        common_indices_14 = lj14_sigma_map.index.intersection(
            reference_contact_matrices[name][reference_contact_matrices[name]["rc_same_chain"]]
            .set_index(["rc_ai", "rc_aj"])
            .index
        )

        reference_contact_matrices[name].loc[common_indices, "sigma_prior"] = lj_sigma_map.astype("float64")
        if not common_indices_14.empty:
            reference_contact_matrices[name].loc[common_indices_14, "sigma_prior"] = lj14_sigma_map.astype("float64")
        reference_contact_matrices[name].loc[common_indices, "epsilon_prior"] = lj_epsilon_map.astype("float64")
        if not common_indices_14.empty:
            reference_contact_matrices[name].loc[common_indices_14, "epsilon_prior"] = lj14_epsilon_map.astype("float64")

        reference_contact_matrices[name].drop(columns=["c6_i", "c6_j", "c12_i", "c12_j", "c6", "c12"], inplace=True)

        et = time.time()
        elapsed_time = et - st
        st = et
        print("\t- Done in:", elapsed_time, "seconds")

    matrices["reference_matrices"] = reference_contact_matrices
    reference_set = set(ensemble["topology_dataframe"]["name"].to_list())

    # TODO check intra domain complementarity
    # check_intra_domain_complementarity(matrices["reference_matrices"])

    computed_contact_matrices = []
    train_contact_matrices_general = {}
    train_contact_matrices = {}
    train_topology_dataframe = pd.DataFrame()

    for reference in args.input_refs:
        trainings = reference["train"]
        for simulation in trainings:
            print("\t-", f"Initializing {simulation} ensemble data")
            simulation_path = f"{args.root_dir}/inputs/{args.system}/{simulation}"
            topology_path = f"{simulation_path}/topol.top"
            if not os.path.isfile(topology_path):
                raise FileNotFoundError(f"{topology_path} not found.")

            print("\t\t-", f"Reading {topology_path}")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                topol = parmed.load_file(topology_path)

            (
                temp_topology_dataframe,
                ensemble["molecules_idx_sbtype_dictionary"],
                _,
                _,
                _,
                _,
                _,
                _,
            ) = contacts_init.initialize_topology(topol, custom_dict, args)

            train_topology_dataframe = pd.concat(
                [train_topology_dataframe, temp_topology_dataframe],
                axis=0,
                ignore_index=True,
            )

            matrix_paths = [f"{simulation_path}/{a}" for a in os.listdir(simulation_path) if reference["matrix"] in a]
            if len(matrix_paths) > 1:
                raise ValueError(f"More than 1 matrix found in {simulation_path}: {matrix_paths}")
            if matrix_paths == []:
                raise FileNotFoundError(
                    f"Contact matrix file(s) must be named as intramat_X_X.ndx(.gz/.h5) or intermat_X_Y.ndx(.gz/.h5). Found instead {reference['matrix']}"
                )

            path = matrix_paths[0]
            train_name = _path_to_matrix_name(path, args.root_dir)
            name = f"{args.system}/{reference['reference']}/{simulation}/{reference['matrix']}".replace("/", "_")

            if train_name not in computed_contact_matrices:
                train_contact_matrices_general[train_name] = io.read_molecular_contacts(
                    path, ensemble["molecules_idx_sbtype_dictionary"], simulation, path.endswith(".h5")
                )
                computed_contact_matrices.append(train_name)
                train_contact_matrices[name] = train_contact_matrices_general[train_name]
            else:
                train_contact_matrices[name] = train_contact_matrices_general[train_name].copy()

            ref_name = f"{args.system}/{reference['reference']}/{reference['matrix']}".replace("/", "_")

            if ref_name == []:
                raise FileNotFoundError(f"No corresponding reference matrix found for {path}")
            ensemble["train_matrix_tuples"].append((name, ref_name))
            train_contact_matrices[name] = initialize_molecular_contacts(
                train_contact_matrices[name],
                reference_contact_matrices[ref_name],
                args,
                reference,
            )

            et = time.time()
            elapsed_time = et - st
            st = et
            print("\t- Done in:", elapsed_time, "seconds")

    del train_contact_matrices_general
    matrices["train_matrices"] = train_contact_matrices

    comparison_set = set()
    for number, molecule in enumerate(ensemble["topology"].molecules, 1):
        comparison_dataframe = train_topology_dataframe.loc[train_topology_dataframe["molecule"] == f"{number}_{molecule}"]
        if not comparison_dataframe.empty:
            comparison_set = set(
                comparison_dataframe[
                    # TODO use a nicer way to do this (Use a list of possible "H" or dictionary of names external to multiego.py and import it)
                    (~comparison_dataframe["name"].str.startswith("H"))
                    | (comparison_dataframe["name"].str == "H")
                ]["name"].to_list()
            )
        else:
            raise RuntimeError("the molecule names in the training topologies do not match those in the reference")

    difference_set = comparison_set.difference(reference_set)
    if difference_set:
        print(
            f'The following atomtypes are not converted:\n{difference_set} \nYou MUST add them in "from_ff_to_multiego" dictionary to properly merge all the contacts.'
        )
        sys.exit()

    return ensemble, matrices


def init_meGO_ensemble(args, custom_dict):
    print("\t-", "Initializing system topology")
    base_topology_path = f"{args.root_dir}/inputs/{args.system}/topol.top"

    if not os.path.isfile(base_topology_path):
        raise FileNotFoundError(f"{base_topology_path} not found.")

    print("\t\t-", f"Reading {base_topology_path}")
    # ignore the dihedral type overriding in parmed
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        defines = {"DISULFIDE": 1}
        base_reference_topology = parmed.load_file(base_topology_path, defines)
    (
        topology_dataframe,
        molecules_idx_sbtype_dictionary,
        sbtype_c12_dict,
        sbtype_mg_c12_dict,
        sbtype_mg_c6_dict,
        sbtype_name_dict,
        sbtype_moltype_dict,
        molecule_type_dict,
    ) = contacts_init.initialize_topology(base_reference_topology, custom_dict, args)

    ensemble = {}
    ensemble["topology"] = base_reference_topology
    ensemble["topology_dataframe"] = topology_dataframe
    ensemble["molecules_idx_sbtype_dictionary"] = (
        molecules_idx_sbtype_dictionary  # molecule, {index, mego_type} -> 1: N_mol_resnum
    )
    ensemble["sbtype_c12_dict"] = sbtype_c12_dict  # {mego_type: c12} N_mol_resnum: c12
    ensemble["sbtype_mg_c12_dict"] = sbtype_mg_c12_dict  # {mego_type: c12} N_mol_resnum: c12
    ensemble["sbtype_mg_c6_dict"] = sbtype_mg_c6_dict  # {mego_type: c12} N_mol_resnum: c12
    ensemble["sbtype_name_dict"] = sbtype_name_dict  # {mego_type: atom_name} N_mol_resnum: N
    ensemble["sbtype_moltype_dict"] = sbtype_moltype_dict  # {mego_type: moltype} N_mol_resnum: protein
    ensemble["sbtype_number_dict"] = (
        ensemble["topology_dataframe"][["sb_type", "number"]].set_index("sb_type")["number"].to_dict()
    )  # {mego_type: atomnumer} N_mol_resnum: 1
    # {mego_type: atomtype} N_mol_resnum: NL
    ensemble["sbtype_type_dict"] = {key: name for key, name in ensemble["topology_dataframe"][["sb_type", "type"]].values}
    # {molecule: moltype} mol: protein
    ensemble["molecule_type_dict"] = molecule_type_dict
    ensemble["train_matrix_tuples"] = []

    return ensemble
