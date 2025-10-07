from .resources import type_definitions
from . import io
from . import topology
from .util import masking

# import glob
import numpy as np
import pandas as pd
import parmed
import os
import warnings
import itertools
import time
import networkx as nx



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
    from_ff_to_multiego_extended = type_definitions.from_ff_to_multiego
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


def initialize_molecular_contacts(contact_matrix, prior_matrix, args, reference):
    """
    This function initializes a contact matrix for a given simulation.
    """

    # remove un-learned contacts (intra-inter domain)
    contact_matrix["learned"] = prior_matrix["rc_learned"].to_numpy()
    contact_matrix["reference"] = reference["reference"]
    # calculate adaptive rc/md threshold
    # sort probabilities, and calculate the normalized cumulative distribution
    p_sort = np.sort(contact_matrix["probability"].loc[(contact_matrix["learned"])].to_numpy())[::-1]
    norm = np.sum(p_sort)
    if norm == 0:
        p_sort_normalized = 0
        md_threshold = 1
    else:
        # find md threshold
        p_sort_normalized = np.cumsum(p_sort) / norm
        md_threshold = p_sort[np.min(np.where(p_sort_normalized > args.p_to_learn)[0])]

    contact_matrix["epsilon_0"] = reference["epsilon"]
    # add the columns for rc, md threshold
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

    # modify limit_rc_att in the cases where epsilon_prior is negative and limit_rc_att is below 1 == epsilon_0 < epsilon_min)
    contact_matrix.loc[(contact_matrix["limit_rc_att"] < 1) & (prior_matrix["epsilon_prior"] < 0), "limit_rc_att"] = 1

    return contact_matrix


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
    ) = initialize_topology(base_reference_topology, custom_dict, args)

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


def check_intra_domain_complementarity(matrices):
    """
    This function checks that each non-single reference matrix (double intramat_1_1 from two different references) have non overlapping learning flag
    """

    mats_names = []
    for key, _ in matrices.items():
        # get matrices names e.g. intramat_1_1, intermat_1_2
        mats_names.append("_".join(key.split("_")[-3:]))
    # check all matrices which are present in more than 1 reference
    to_check_names = list(set([a for a in mats_names if mats_names.count(a) > 1]))
    to_check = [[k for k in matrices.keys() if "_".join(k.split("_")[-3:]) in check] for check in to_check_names]
    for check in to_check:
        intra_flags = []
        for key in check:
            intra_flags.append(matrices[key]["rc_learned"].to_numpy())
        if np.any(np.sum(intra_flags, axis=0) > 1):
            raise ValueError(f"Learning flag complementarity not satisfied for {check} (e.g. intra-inter domain splitting)")


# TODO this hole function should iterate over references and than internally over the trainings keeping stored the already processed training by path name
# Even though in this way the check consinstency between reference matrices is faster
def init_meGO_matrices(ensemble, args, custom_dict):
    """
    Initializes meGO.

    Args:
    - args (object): Object containing arguments for initializing the ensemble.

    Returns:
    - ensemble (dict): A dictionary containing the initialized ensemble with various molecular attributes and contact matrices.

    This function sets up meGO by initializing the reference topology and processing train and check contact matrices based
    on the provided arguments. It reads topology files, loads molecular information, and sets up dictionaries and data frames
    to organize molecular data and contact matrices.

    The function initializes the reference topology and extracts essential molecular details such as topological data frames,
    subtype dictionaries, c12 values, names, molecule types, and contact matrices for the reference ensemble. It then processes
    train and check contact matrices, aligning them with the reference ensemble to detect any differences in atom types.

    If atom type differences are found between ensembles, the function prints a warning message and exits, indicating
    the need to add missing atom types to the conversion dictionary for proper contact merging.

    Note:
    - This function assumes the availability of various directories, files, and modules (e.g., 'parmed', 'io').
    - The 'args' object should contain necessary arguments for setting up the ensemble.
    - The returned 'ensemble' dictionary encapsulates crucial details of the initialized ensemble for further analysis or processing.
    """

    st = time.time()
    reference_contact_matrices = {}
    matrices = {}

    # if there are more than 1 reference associated to the same
    # check if reference are associated to the same molecule pair
    # if intramat> check for intra domain complementarity
    for reference in args.input_refs:  # reference_paths:
        print("\t-", f"Initializing {reference['reference']} ensemble data")
        reference_path = f"{args.root_dir}/inputs/{args.system}/{reference['reference']}"
        topol_files = [f for f in os.listdir(reference_path) if ".top" in f]
        if len(topol_files) > 1:
            raise RuntimeError(f"More than 1 topology file found in {reference_path}. Only one should be used")

        topology_path = f"{reference_path}/{topol_files[0]}"
        if not os.path.isfile(topology_path):
            raise FileNotFoundError(f"{topology_path} not found.")

        print("\t\t-", f"Reading {topology_path}")
        # ignore the dihedral type overriding in parmed
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            topol = parmed.load_file(topology_path)

        # these are the atom type c6_i,c12_j
        lj_data = topology.get_lj_params(topol)
        # these are the combined cases (c6_ij, c12_ij)
        lj_pairs = topology.get_lj_pairs(topol)
        # Create reversed pairs
        reversed_lj_pairs = lj_pairs.rename(columns={"ai": "aj", "aj": "ai"})
        # Combine original and reversed
        symmetric_lj_pairs = pd.concat([lj_pairs, reversed_lj_pairs])
        # Remove duplicates to avoid duplication of symmetric pairs
        # This step ensures that if a pair already exists in both directions, it's not duplicated
        symmetric_lj_pairs = symmetric_lj_pairs.drop_duplicates(subset=["ai", "aj"]).reset_index(drop=True)

        # these are the combined cases in the [pairs] section (c6_ij, c12_ij)
        lj14_pairs = topology.get_lj14_pairs(topol)
        # Create reversed pairs
        reversed_lj14_pairs = lj14_pairs.rename(columns={"ai": "aj", "aj": "ai"})
        # Combine original and reversed
        symmetric_lj14_pairs = pd.concat([lj14_pairs, reversed_lj14_pairs])
        # Remove duplicates to avoid duplication of symmetric pairs
        # This step ensures that if a pair already exists in both directions, it's not duplicated
        symmetric_lj14_pairs = symmetric_lj14_pairs.drop_duplicates(subset=["ai", "aj"]).reset_index(drop=True)

        lj_data_dict = {str(key): val for key, val in zip(lj_data["ai"], lj_data[["c6", "c12"]].values)}

        ensemble["topology_dataframe"]["c6"] = lj_data["c6"].to_numpy()
        ensemble["topology_dataframe"]["c12"] = lj_data["c12"].to_numpy()

        matrix_paths = [f"{reference_path}/{a}" for a in os.listdir(reference_path) if reference["matrix"] in a]

        # if matrix path is more than 1 raise error
        if len(matrix_paths) > 1:
            raise ValueError(f"More than 1 matrix found in {reference_path}: {matrix_paths}")

        if matrix_paths == []:
            raise FileNotFoundError(
                f"Contact matrix file(s) must be named as intramat_X_X.ndx(.gz/.h5) or intermat_X_Y.ndx(.gz/.h5). Found instead: {reference['matrix']}"
            )

        path = matrix_paths[0]
        name = path.replace(f"{args.root_dir}/inputs/", "")
        name = name.replace("/", "_")
        name = name.replace(".ndx", "")
        name = name.replace(".gz", "")
        name = name.replace(".h5", "")
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

        # Create a mapping from lj_pairs for sigma and epsilon
        lj_sigma_map = symmetric_lj_pairs.set_index(["ai", "aj"])["sigma"]
        lj_epsilon_map = symmetric_lj_pairs.set_index(["ai", "aj"])["epsilon"]
        lj14_sigma_map = symmetric_lj14_pairs.set_index(["ai", "aj"])["sigma"]
        lj14_epsilon_map = symmetric_lj14_pairs.set_index(["ai", "aj"])["epsilon"]

        # Filter lj_sigma_map to include only indices that exist in reference_contact_matrices[name]
        common_indices = lj_sigma_map.index.intersection(reference_contact_matrices[name].set_index(["rc_ai", "rc_aj"]).index)
        # in this case we want to apply it only for intramolecular contacts
        common_indices_14 = lj14_sigma_map.index.intersection(
            reference_contact_matrices[name][reference_contact_matrices[name]["rc_same_chain"]]
            .set_index(["rc_ai", "rc_aj"])
            .index
        )

        # Update sigma values where they exist in lj_pairs
        reference_contact_matrices[name].loc[common_indices, "sigma_prior"] = lj_sigma_map.astype("float64")
        reference_contact_matrices[name].loc[common_indices_14, "sigma_prior"] = lj14_sigma_map.astype("float64")

        # Update epsilon values where they exist in lj_pairs
        reference_contact_matrices[name].loc[common_indices, "epsilon_prior"] = lj_epsilon_map.astype("float64")
        reference_contact_matrices[name].loc[common_indices_14, "epsilon_prior"] = lj14_epsilon_map.astype("float64")

        reference_contact_matrices[name].drop(columns=["c6_i", "c6_j", "c12_i", "c12_j", "c6", "c12"], inplace=True)

        et = time.time()
        elapsed_time = et - st
        st = et
        print("\t- Done in:", elapsed_time, "seconds")

    matrices["reference_matrices"] = reference_contact_matrices
    reference_set = set(ensemble["topology_dataframe"]["name"].to_list())

    # check intra domain complementarity
    # check_intra_domain_complementarity(matrices["reference_matrices"])

    # now we process the train contact matrices
    # keep track of the training matrices already processed
    computed_contact_matrices = []
    # Store un-processed train matrices (is multiple reference input use the same train matrix use this instead of re-reading)
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
            # ignore the dihedral type overriding in parmed
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
            ) = initialize_topology(topol, custom_dict, args)

            train_topology_dataframe = pd.concat(
                [train_topology_dataframe, temp_topology_dataframe],
                axis=0,
                ignore_index=True,
            )
            matrix_paths = [f"{simulation_path}/{a}" for a in os.listdir(simulation_path) if reference["matrix"] in a]
            # if matrix path is more than 1 raise error
            if len(matrix_paths) > 1:
                raise ValueError(f"More than 1 matrix found in {reference_path}: {matrix_paths}")

            if matrix_paths == []:
                raise FileNotFoundError(
                    f"Contact matrix file(s) must be named as intramat_X_X.ndx(.gz/.h5) or intermat_X_Y.ndx(.gz/.h5). Found instead {reference['matrix']}"
                )

            path = matrix_paths[0]
            # needed to check if training wa already read to avoid reading it multiple times
            train_name = path.replace(f"{args.root_dir}/inputs/", "")
            train_name = train_name.replace("/", "_")
            train_name = train_name.replace(".ndx", "")
            train_name = train_name.replace(".gz", "")
            train_name = train_name.replace(".h5", "")
            # Use the name containing both the reference and the training in order to have a unique training name for each training and reference
            name = f"{args.system}/{reference['reference']}/{simulation}/{reference['matrix']}"
            name = name.replace("/", "_")
            # if the training was already read just copy it instead of re-reading it
            if train_name not in computed_contact_matrices:
                train_contact_matrices_general[train_name] = io.read_molecular_contacts(
                    path, ensemble["molecules_idx_sbtype_dictionary"], simulation, path.endswith(".h5")
                )
                computed_contact_matrices.append(train_name)
                train_contact_matrices[name] = train_contact_matrices_general[train_name]
            else:
                train_contact_matrices[name] = train_contact_matrices_general[train_name].copy()

            # reference name is already uniquely associated to the training
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

    # force memory cleaning to decrease footprint in case of large dataset
    del train_contact_matrices_general
    matrices["train_matrices"] = train_contact_matrices

    comparison_set = set()
    for number, molecule in enumerate(ensemble["topology"].molecules, 1):
        comparison_dataframe = train_topology_dataframe.loc[train_topology_dataframe["molecule"] == f"{number}_{molecule}"]
        if not comparison_dataframe.empty:
            comparison_set = set(
                comparison_dataframe[
                    # TODO use a nicer way to to this (Use a list of possible "H" or dictionary of names external to multiego.py and import it)
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
        exit()

    return ensemble, matrices


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
    if "meGO_bonded_interactions" not in meGO_ensemble.keys():
        meGO_ensemble["meGO_bonded_interactions"] = {}
    if "bond_pairs" not in meGO_ensemble.keys():
        meGO_ensemble["bond_pairs"] = {}
    if "user_pairs" not in meGO_ensemble.keys():
        meGO_ensemble["user_pairs"] = {}

    for molecule, topol in meGO_ensemble["topology"].molecules.items():
        meGO_ensemble["meGO_bonded_interactions"][molecule] = {
            "bonds": topology.get_bonds(topol[0].bonds),
            "angles": topology.get_angles(topol[0].angles),
            "dihedrals": topology.get_dihedrals(topol[0].dihedrals),
            "impropers": topology.get_impropers(topol[0].impropers),
            "pairs": topology.get_pairs(topol[0].adjusts),
        }
        # The following bonds are used in the parametrization of LJ 1-4
        meGO_ensemble["bond_pairs"][molecule] = topology.get_bond_pairs(topol[0].bonds)
        meGO_ensemble["user_pairs"][molecule] = topology.get_pairs(topol[0].adjusts)

    return meGO_ensemble


def generate_14_data(meGO_ensemble):
    """
    Generates data for 1-4 interactions within a molecular ensemble.

    Args:
    - meGO_ensemble (dict): A dictionary containing information about the molecular ensemble.

    Returns:
    - pairs14 (DataFrame): DataFrame containing information about 1-4 interactions.
    - exclusion_bonds14 (DataFrame): DataFrame containing exclusion bonded interactions.

    This function generates data for 1-4 interactions within a molecular ensemble.
    It iterates through each molecule in the ensemble, processes the topology, and computes
    exclusion bonded interactions and specific 1-4 interactions.

    The function creates DataFrames 'pairs14' and 'exclusion_bonds14' containing information
    about 1-4 interactions and exclusion bonded interactions, respectively.
    It extracts details such as atom numbers, subtypes, residue numbers, names, types, residue names,
    molecule types, and interaction characteristics.

    Note:
    - The 'meGO_ensemble' dictionary is expected to contain necessary details regarding the molecular ensemble.
    - The returned DataFrames provide comprehensive information about 1-4 interactions and exclusion bonded
      interactions within the ensemble for further analysis or processing.
    """
    # First of all we generate the random-coil 1-4 interactions:
    pairs14 = pd.DataFrame()
    exclusion_bonds14 = pd.DataFrame()
    for idx, (molecule, bond_pair) in enumerate(meGO_ensemble["bond_pairs"].items(), start=1):
        if not bond_pair:
            continue
        reduced_topology = (
            meGO_ensemble["topology_dataframe"]
            .loc[meGO_ensemble["topology_dataframe"]["molecule_name"] == molecule][
                [
                    "number",
                    "sb_type",
                    "resnum",
                    "name",
                    "type",
                    "resname",
                    "molecule_type",
                ]
            ]
            .copy()
        )

        reduced_topology["number"] = reduced_topology["number"].astype(str)
        reduced_topology["resnum"] = reduced_topology["resnum"].astype(int)
        # Dictionaries definitions to map values
        type_atnum_dict = reduced_topology.set_index("number")["sb_type"].to_dict()

        # Building the exclusion bonded list
        # exclusion_bonds are all the interactions within 3 bonds
        # p14 are specifically the interactions at exactly 3 bonds
        # b6 are the interactions up to 6 bonds
        # exclusion_bonds, tmp_p14, b6 = topology.generate_bond_exclusions(reduced_topology, bond_pair)
        bd = topology.compute_bond_distances(reduced_topology, bond_pair)

        # split->convert->remerge:
        # tmp_ex = pd.DataFrame(columns=["ai", "aj", "exclusion_bonds"])
        # tmp_ex["exclusion_bonds"] = pd.Series(exclusion_bonds).astype("category")
        # tmp_ex[["ai", "aj"]] = tmp_ex["exclusion_bonds"].str.split("_", expand=True)
        # tmp_ex["ai"] = tmp_ex["ai"].map(type_atnum_dict).astype("category")
        # tmp_ex["aj"] = tmp_ex["aj"].map(type_atnum_dict).astype("category")
        # tmp_ex["same_chain"] = True
        # tmp_ex["1-4"] = "1_2_3"
        # tmp_ex.loc[(tmp_ex["exclusion_bonds"].isin(tmp_p14)), "1-4"] = "1_4"
        # tmp_ex["1-4"] = tmp_ex["1-4"].astype("category")
        # exclusion_bonds14 = pd.concat([exclusion_bonds14, tmp_ex], axis=0, sort=False, ignore_index=True)
        exclusion_bonds14 = pd.concat([exclusion_bonds14, bd], axis=0, sort=False, ignore_index=True)

        # Adding the c12 for 1-4 interactions
        reduced_topology["c12"] = reduced_topology["sb_type"].map(meGO_ensemble["sbtype_c12_dict"])

        pairs = pd.DataFrame()
        if meGO_ensemble["molecule_type_dict"][molecule] == "protein":
            pairs = topology.protein_LJ14(reduced_topology)
            pairs["ai"] = pairs["ai"].map(type_atnum_dict).astype("category")
            pairs["aj"] = pairs["aj"].map(type_atnum_dict).astype("category")
            pairs["rep"] = pairs["c12"]
            pairs["source"] = pairs["source"].astype("category")
            pairs["same_chain"] = True
            # pairs["1-4"] = "1_4"
        else:
            pairs["ai"] = meGO_ensemble["user_pairs"][molecule].ai.astype(str)
            pairs["aj"] = meGO_ensemble["user_pairs"][molecule].aj.astype(str)
            pairs["ai"] = pairs["ai"].map(type_atnum_dict).astype("category")
            pairs["aj"] = pairs["aj"].map(type_atnum_dict).astype("category")

            nonprotein_c12 = []
            for test in meGO_ensemble["user_pairs"][molecule].type:
                if test is None:
                    print(
                        "\nERROR: you have 1-4 pairs defined in your reference topology without the associated C6/C12 values"
                    )
                    print("       user provided 1-4 pairs need to define also the C6/C12\n")
                    exit()
                nonprotein_c12.append(float(test.epsilon) * 4.184)

            pairs["func"] = 1
            pairs["c6"] = 0.0
            pairs["c12"] = nonprotein_c12
            pairs["probability"] = 1.0
            pairs["rc_probability"] = 1.0
            pairs["source"] = pd.Series(["1-4"] * len(pairs), dtype="category")
            pairs["rep"] = pairs["c12"]
            pairs["same_chain"] = True
            # copy and symmetrize
            tmp = pairs.copy()
            tmp["ai"], tmp["aj"] = tmp["aj"], tmp["ai"]
            pairs = pd.concat([pairs, tmp], axis=0, sort=False, ignore_index=True)

        mol_ai = f"{idx}_{molecule}"
        pairs["molecule_name_ai"] = mol_ai
        pairs["molecule_name_aj"] = mol_ai

        if not pairs.empty:
            pairs14 = pd.concat([pairs14, pairs], axis=0, sort=False, ignore_index=True)

    return pairs14, exclusion_bonds14


def annotate_bond_distances(reduced_topology, bond_pair, pair_df, max_distance=6):
    """
    Given a topology and a DataFrame with sb_type pairs (ai, aj),
    return the same DataFrame with an extra column `bond_distance` (0–6, or -1).
    """

    # Map sb_type → atom number
    sbtype_to_atnum = reduced_topology.set_index("sb_type")["number"].to_dict()

    # Build bond graph (undirected)
    G = nx.Graph()
    G.add_edges_from(bond_pair)

    # Precompute all shortest paths up to max_distance
    all_lengths = dict(nx.all_pairs_shortest_path_length(G, cutoff=max_distance))

    # For each row in the pair_df, determine bond distance
    bond_distances = []
    for _, row in pair_df.iterrows():
        try:
            ai = sbtype_to_atnum[row["ai"]]
            aj = sbtype_to_atnum[row["aj"]]
        except KeyError:
            bond_distances.append(-1)
            continue

        dist = all_lengths.get(ai, {}).get(aj, -1)
        if dist > max_distance:
            dist = -1

        bond_distances.append(dist)

    # Add column to dataframe
    pair_df = pair_df.copy()
    pair_df["bond_distance"] = bond_distances
    return pair_df


def init_LJ_datasets(meGO_ensemble, matrices, pairs14, exclusion_bonds14, args):
    # we cycle over train matrices to pair them with reference matrices and
    # then we add 1-4 assignments and defaults c12s and concatenate everything
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

    for name, ref_name in meGO_ensemble["train_matrix_tuples"]:
        # sysname_train_intramat_1_1 <-> sysname_reference_intramat_1_1
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

        # This is a debug check to avoid data inconsistencies
        if (np.abs(temp_merged["rc_cutoff"] - temp_merged["cutoff"])).max() > 0:
            print(
                temp_merged[["ai", "aj", "rc_ai", "rc_aj", "source", "rc_source", "cutoff", "rc_cutoff"]]
                .loc[(np.abs(temp_merged["rc_cutoff"] - temp_merged["cutoff"]) > 0)]
                .to_string()
            )
            exit(
                "HERE SOMETHING BAD HAPPEND: There are inconsistent cutoff values between the MD and corresponding RC input data"
            )

        # This is a debug check to avoid data inconsistencies
        if not temp_merged["rc_same_chain"].equals(temp_merged["rc_same_chain"]):
            diff_indices = temp_merged.index[temp_merged["same_chain"] != temp_merged["rc_same_chain"]].tolist()
            print(f"Difference found at indices: {diff_indices}")
            exit("HERE SOMETHING BAD HAPPEND: You are pairing intra and inter molecular training and reference data")

        temp_merged = temp_merged[td_fields]
        train_dataset = pd.concat([train_dataset, temp_merged], axis=0, sort=False, ignore_index=True)

    train_dataset["molecule_name_ai"] = train_dataset["molecule_name_ai"].astype("category")
    train_dataset["molecule_name_aj"] = train_dataset["molecule_name_aj"].astype("category")
    train_dataset["source"] = train_dataset["source"].astype("category")

    # now we need to set the bond distance field and the default repulsive C12 "rep"

    train_dataset = pd.merge(
        pd.merge(
            train_dataset,
            pairs14[["ai", "aj", "same_chain", "rep"]],
            how="left",
            on=["ai", "aj", "same_chain"],
        ),
        # exclusion_bonds14[["ai", "aj", "same_chain", "1-4"]],
        exclusion_bonds14[["ai", "aj", "bond_distance"]],
        how="left",
        on=["ai", "aj"],
    )
    train_dataset["bond_distance"] = train_dataset["bond_distance"].fillna(7).astype(int)

    train_dataset["ai"] = train_dataset["ai"].astype("category")
    train_dataset["aj"] = train_dataset["aj"].astype("category")

    # We remove from train the 0_1_2_3 intramolecolar interactions
    train_dataset = train_dataset[
        # ~(((train_dataset["ai"] == train_dataset["aj"]) & train_dataset["same_chain"]) | (train_dataset["1-4"] == "1_2_3"))
        ~((train_dataset["same_chain"]) & (train_dataset["bond_distance"] < 3))
    ]
    train_dataset.reset_index(inplace=True)

    # train_dataset["1-4"] = train_dataset["1-4"].cat.add_categories(["1>4"])
    # train_dataset["1-4"] = train_dataset["1-4"].fillna("1>4").astype("category")
    # train_dataset.loc[(train_dataset["1-4"] == "1_4") & (train_dataset["rep"].isnull()), "rep"] = 0.0
    train_dataset.loc[
        (train_dataset["bond_distance"] == 3) & (train_dataset["same_chain"]) & (train_dataset["rep"].isnull()), "rep"
    ] = 0.0

    type_to_c12 = {key: val for key, val in zip(type_definitions.gromos_atp.name, type_definitions.gromos_atp.rc_c12)}

    if args.custom_c12 is not None:
        custom_c12_dict = io.read_custom_c12_parameters(args.custom_c12)
        type_to_c12_appo = {key: val for key, val in zip(custom_c12_dict.name, custom_c12_dict.c12)}
        type_to_c12.update(type_to_c12_appo)

    type_ai_mapped = train_dataset["ai"].map(meGO_ensemble["sbtype_type_dict"])
    type_aj_mapped = train_dataset["aj"].map(meGO_ensemble["sbtype_type_dict"])

    # set of interactions with parameters not resulting from the combination rule
    # oxygen-oxygen repulsion
    OO_mask = masking.create_linearized_mask(
        type_ai_mapped.to_numpy(),
        type_aj_mapped.to_numpy(),
        [("O", "O"), ("OM", "OM"), ("O", "OM")],
        symmetrize=True,
    )

    # hydrongen-oxygen attraction
    HO_mask = masking.create_linearized_mask(
        type_ai_mapped.to_numpy(),
        type_aj_mapped.to_numpy(),
        [("H", "O"), ("H", "OM"), ("H", "OA")],
        symmetrize=True,
    )

    # oxygen-nitrogen repulsion (when not attractive)
    ON_mask = masking.create_linearized_mask(
        type_ai_mapped.to_numpy(),
        type_aj_mapped.to_numpy(),
        [
            ("O", "N"),
            ("OM", "N"),
        ],
        symmetrize=True,
    )

    # hydrogen-hydrogen repulsion
    # Define condition where only ai or aj (but not both) starts with "H"
    H_mask = train_dataset["ai"].str.startswith("H") ^ train_dataset["aj"].str.startswith("H")
    HH_mask = train_dataset["ai"].str.startswith("H") & train_dataset["aj"].str.startswith("H")

    # TODO
    # here we should iterate over the pairs in type definition amd generate the repulsions
    # NL-NZ repulsion
    NN_mask = masking.create_linearized_mask(
        type_ai_mapped.to_numpy(),
        type_aj_mapped.to_numpy(),
        [("NL", "NL"), ("NZ", "NZ"), ("NL", "NZ"), ("N", "N")],
        symmetrize=True,
    )

    CC_mask = masking.create_linearized_mask(
        type_ai_mapped.to_numpy(),
        type_aj_mapped.to_numpy(),
        [("CH3", "CH3"), ("C", "C"), ("C", "CH3"), ("C", "CAH")],
        symmetrize=True,
    )

    # default repulsive C12 (rep)
    pairwise_c12 = np.sqrt(
        train_dataset["ai"].map(meGO_ensemble["sbtype_c12_dict"]) * train_dataset["aj"].map(meGO_ensemble["sbtype_c12_dict"])
    )
    train_dataset["rep"] = train_dataset["rep"].fillna(pd.Series(pairwise_c12))
    # special REP cases:
    # train_dataset.loc[OO_mask & (train_dataset["1-4"] != "1_4"), "rep"] = type_definitions.mg_OO_c12_rep
    # train_dataset.loc[ON_mask & (train_dataset["1-4"] != "1_4"), "rep"] = type_definitions.mg_ON_c12_rep
    # train_dataset.loc[HH_mask & (train_dataset["1-4"] != "1_4"), "rep"] = type_definitions.mg_HH_c12_rep
    # train_dataset.loc[NN_mask & (train_dataset["1-4"] != "1_4"), "rep"] = type_definitions.mg_NN_c12_rep
    train_dataset.loc[OO_mask & ((train_dataset["bond_distance"] != 3) | (~train_dataset["same_chain"])), "rep"] = (
        type_definitions.mg_OO_c12_rep
    )
    train_dataset.loc[ON_mask & ((train_dataset["bond_distance"] != 3) | (~train_dataset["same_chain"])), "rep"] = (
        type_definitions.mg_ON_c12_rep
    )
    train_dataset.loc[HH_mask & ((train_dataset["bond_distance"] != 3) | (~train_dataset["same_chain"])), "rep"] = (
        type_definitions.mg_HH_c12_rep
    )
    train_dataset.loc[NN_mask & ((train_dataset["bond_distance"] != 3) | (~train_dataset["same_chain"])), "rep"] = (
        type_definitions.mg_NN_c12_rep
    )
    # train_dataset.loc[CC_mask & ((train_dataset["bond_distance"] != 3) | (~train_dataset["same_chain"])), "rep"] = type_definitions.mg_CC_c12_rep

    # default (mg) sigma
    pairwise_mg_sigma = (
        train_dataset["ai"].map(meGO_ensemble["sbtype_mg_c12_dict"])
        * train_dataset["aj"].map(meGO_ensemble["sbtype_mg_c12_dict"])
        / (
            train_dataset["ai"].map(meGO_ensemble["sbtype_mg_c6_dict"])
            * train_dataset["aj"].map(meGO_ensemble["sbtype_mg_c6_dict"])
        )
    ) ** (1 / 12)
    train_dataset["mg_sigma"] = pd.Series(pairwise_mg_sigma)
    # special mg sigma cases:
    train_dataset.loc[OO_mask, "mg_sigma"] = (type_definitions.mg_OO_c12_rep) ** (1 / 12) / 2 ** (1 / 6)
    train_dataset.loc[HH_mask, "mg_sigma"] = type_definitions.mg_HH_c12_rep ** (1 / 12) / 2 ** (1 / 6)
    train_dataset.loc[NN_mask, "mg_sigma"] = (type_definitions.mg_NN_c12_rep) ** (1 / 12) / 2 ** (1 / 6)
    train_dataset.loc[HO_mask, "mg_sigma"] = type_definitions.mg_HO_sigma
    # train_dataset.loc[CC_mask, "mg_sigma"] = type_definitions.mg_CC_c12_rep ** (1 / 12) / 2 ** (1 / 6)

    # default (mg) epsilon
    pairwise_mg_epsilon = (
        train_dataset["ai"].map(meGO_ensemble["sbtype_mg_c6_dict"])
        * train_dataset["aj"].map(meGO_ensemble["sbtype_mg_c6_dict"])
    ) / (
        4
        * np.sqrt(
            train_dataset["ai"].map(meGO_ensemble["sbtype_mg_c12_dict"])
            * train_dataset["aj"].map(meGO_ensemble["sbtype_mg_c12_dict"])
        )
    )
    train_dataset["mg_epsilon"] = pd.Series(pairwise_mg_epsilon)
    # special cases:
    train_dataset.loc[OO_mask, "mg_epsilon"] = -type_definitions.mg_OO_c12_rep
    train_dataset.loc[H_mask, "mg_epsilon"] = 0.0
    train_dataset.loc[HH_mask, "mg_epsilon"] = -type_definitions.mg_HH_c12_rep
    train_dataset.loc[NN_mask, "mg_epsilon"] = -type_definitions.mg_NN_c12_rep
    train_dataset.loc[HO_mask, "mg_epsilon"] = type_definitions.mg_eps_HO
    # train_dataset.loc[CC_mask, "mg_epsilon"] = -type_definitions.mg_CC_c12_rep

    # final cleaning
    train_dataset.dropna(subset=["mg_sigma"], inplace=True)
    train_dataset = train_dataset.loc[train_dataset["rep"] > 0.0]

    # This is a debug check to avoid data inconsistencies
    if (np.abs(train_dataset["cutoff"] - 1.45 * train_dataset["rep"] ** (1 / 12))).max() > 10e-6:
        print(
            train_dataset[["ai", "aj", "source", "same_chain", "cutoff", "rep"]]
            .loc[(np.abs(train_dataset["cutoff"] - 1.45 * train_dataset["rep"] ** (1 / 12)) > 10e-6)]
            .to_string()
        )
        # exit(
        #    "HERE SOMETHING BAD HAPPEND: There are inconsistent cutoff and C12 repulsive values"
        # )

    return train_dataset


def generate_MG_LJ_pairs_rep(sbtype1, sbtype2, dictionary_name_rc_c12, c12_rep=None, factor=1.0):

    # if sbtype1==sbtype2: use repeat
    if sbtype1 == sbtype2:
        combinations = list(itertools.product(sbtype1, repeat=2))
    else:
        combinations = list(itertools.product(sbtype1, sbtype2)) + list(itertools.product(sbtype2, sbtype1))
    # print("Combinations:", combinations)
    if c12_rep is None:
        c12_rep = np.array([np.sqrt(dictionary_name_rc_c12[ai] * dictionary_name_rc_c12[aj]) for ai, aj in combinations])

    pairs_LJ = pd.DataFrame(combinations, columns=["ai", "aj"])
    pairs_LJ["c12"] = c12_rep / factor
    pairs_LJ["c6"] = 0.0
    pairs_LJ["epsilon"] = -c12_rep / factor
    pairs_LJ["sigma"] = c12_rep / factor ** (1.0 / 12.0) / 2.0 ** (1.0 / 6.0)
    pairs_LJ["mg_sigma"] = pairs_LJ["sigma"]
    pairs_LJ["mg_epsilon"] = -c12_rep / factor

    return pairs_LJ.reset_index(drop=True)


def generate_MG_LJ_pairs_attr(sbtype1, sbtype2, dictionary_name_mg_c12, dictionary_name_mg_c6, epsilon=None, sigma=None):

    if sbtype1 == sbtype2:
        combinations = list(itertools.product(sbtype1, repeat=2))
    else:
        combinations = list(itertools.product(sbtype1, sbtype2)) + list(itertools.product(sbtype2, sbtype1))

    # TODO remove this case: is usless (already handled in gromacs combination rule)
    # if epsilon is None and sigma is None:
    #     c6 = np.array([np.sqrt(dictionary_name_mg_c6[ai]*dictionary_name_mg_c6[aj]) for ai, aj in combinations])
    #     c12_rep = np.array([np.sqrt(dictionary_name_mg_c12[ai]*dictionary_name_mg_c12[aj]) for ai, aj in combinations])
    #     epsilon = c6 ** 2.0 / (4.0 * c12_rep)
    #     sigma = (c12_rep/c6) ** (1.0 / 6.0)
    if epsilon is not None and sigma is None:
        # define sigma as combination rule of all combinations pairs
        c6 = np.array([np.sqrt(dictionary_name_mg_c6[ai] * dictionary_name_mg_c6[aj]) for ai, aj in combinations])
        c12_rep = np.array([np.sqrt(dictionary_name_mg_c12[ai] * dictionary_name_mg_c12[aj]) for ai, aj in combinations])
        sigma = (c12_rep / c6) ** (1.0 / 6.0)

    elif epsilon is None and sigma is not None:
        # raise error
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
    The multi-eGO molten-globule force-field includes special repulsive and attractive interaction pairs like O-O, H-H, and O-H.
    TODO: define them by means of an external dictionary instead of hardcoding them. This dictionary should be used also from make_mat
    these are generate in the following
    """
    # reconstruct mapping dictionaries for c12 and c6 which are already mapped in topology_dataframe
    dictionary_name_rc_c12 = {
        name: rc12
        for name, rc12 in zip(meGO_ensemble["topology_dataframe"]["sb_type"], meGO_ensemble["topology_dataframe"]["rc_c12"])
    }
    dictionary_name_mg_c12 = {
        name: mg12
        for name, mg12 in zip(meGO_ensemble["topology_dataframe"]["sb_type"], meGO_ensemble["topology_dataframe"]["mg_c12"])
    }
    dictionary_name_mg_c6 = {
        name: mg6
        for name, mg6 in zip(meGO_ensemble["topology_dataframe"]["sb_type"], meGO_ensemble["topology_dataframe"]["mg_c6"])
    }

    # OO in MG are repulsive (Ramachandran and negatively charged sidechains)
    O_OM_sbtype = [sbtype for sbtype, atomtype in meGO_ensemble["sbtype_type_dict"].items() if atomtype in ["O", "OM"]]
    OO_LJ = generate_MG_LJ_pairs_rep(O_OM_sbtype, O_OM_sbtype, dictionary_name_rc_c12, type_definitions.mg_OO_c12_rep)

    # HH in MG are repulsive (Ramachandran)
    H_H_sbtype = [sbtype for sbtype, atomtype in meGO_ensemble["sbtype_type_dict"].items() if atomtype == "H"]
    HH_LJ = generate_MG_LJ_pairs_rep(H_H_sbtype, H_H_sbtype, dictionary_name_rc_c12, type_definitions.mg_HH_c12_rep)

    # HO in MG are attractive (H-bonds)
    O_OM_OA_sbtype = [
        sbtype for sbtype, atomtype in meGO_ensemble["sbtype_type_dict"].items() if atomtype in ["O", "OM", "OA"]
    ]
    HO_LJ = generate_MG_LJ_pairs_attr(
        O_OM_OA_sbtype,
        H_H_sbtype,
        dictionary_name_mg_c12,
        dictionary_name_mg_c6,
        epsilon=type_definitions.mg_eps_HO,
        sigma=type_definitions.mg_HO_sigma,
    )

    pol_sbtype = [
        sbtype
        for sbtype, atomtype in meGO_ensemble["sbtype_type_dict"].items()
        if atomtype in ["O", "OM", "OA", "N", "NT", "NL", "NR", "NZ", "NE", "CAH", "C", "S", "P", "OE", "CR1"]
    ]
    # CAH in pol --> deve diventare sidechain-backbone weak attr
    hyd_sbtype = [
        sbtype
        for sbtype, atomtype in meGO_ensemble["sbtype_type_dict"].items()
        if atomtype in ["CH", "CH3", "CH3p", "CH2", "CH2r", "CH1"]
    ]
    # pol_hyd_LJ = generate_MG_LJ_pairs_rep(
    #     pol_sbtype, hyd_sbtype, dictionary_name_rc_c12, factor=10)
    pol_hyd_LJ = generate_MG_LJ_pairs_attr(
        pol_sbtype, hyd_sbtype, dictionary_name_mg_c12, dictionary_name_mg_c6, epsilon=type_definitions.mg_eps_pol
    )

    # NL/NZ in MG are repulsive (positevely charged sidechains and N-terminus)
    NL_NZ_sbtype = [
        sbtype
        for sbtype, atomtype in meGO_ensemble["sbtype_type_dict"].items()
        if atomtype in ["NL", "NZ", "N"]  # TODO add also the other
    ]
    NN_LJ = generate_MG_LJ_pairs_rep(NL_NZ_sbtype, NL_NZ_sbtype, dictionary_name_rc_c12, type_definitions.mg_NN_c12_rep)

    # CC_sbtype = [
    #     sbtype for sbtype, atomtype in meGO_ensemble["sbtype_type_dict"].items() if atomtype in ["CH3", "C", "CAH"]
    # ]
    # CC_LJ = generate_MG_LJ_pairs_rep(CC_sbtype, CC_sbtype, dictionary_name_rc_c12, type_definitions.mg_CC_c12_rep)
    # combine them:
    rc_LJ = pd.concat([OO_LJ, HH_LJ, HO_LJ, NN_LJ, pol_hyd_LJ], axis=0)
    # rc_LJ = pd.concat([OO_LJ, HH_LJ, HO_LJ, NN_LJ, CC_LJ], axis=0)
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
    # rc_LJ["1-4"] = "1>4"
    rc_LJ["bond_distance"] = 7
    molecule_names_dictionary = {name.split("_", 1)[1]: name for name in meGO_ensemble["molecules_idx_sbtype_dictionary"]}
    rc_LJ["molecule_name_ai"] = rc_LJ["ai"].apply(lambda x: "_".join(x.split("_")[1:-1])).map(molecule_names_dictionary)
    rc_LJ["molecule_name_aj"] = rc_LJ["aj"].apply(lambda x: "_".join(x.split("_")[1:-1])).map(molecule_names_dictionary)
    rc_LJ["ai"] = rc_LJ["ai"].astype("category")
    rc_LJ["aj"] = rc_LJ["aj"].astype("category")
    rc_LJ["molecule_name_ai"] = rc_LJ["molecule_name_ai"].astype("category")
    rc_LJ["molecule_name_aj"] = rc_LJ["molecule_name_aj"].astype("category")

    return rc_LJ


def set_sig_epsilon(meGO_LJ, parameters):
    """
    Set the epsilon parameter for LJ interactions based on probability and distance.

    This function sets the epsilon parameter for LJ interactions based on the probability and distance of interactions.
    It adjusts epsilon values to represent the strength of LJ interactions, considering both attractive and repulsive forces.

    Parameters
    ----------
    meGO_LJ : pd.DataFrame
        DataFrame containing LJ parameters.

    Returns
    -------
    pd.DataFrame
        DataFrame containing LJ parameters with updated epsilon values.

    Notes
    -----
    This function calculates epsilon values for LJ interactions based on the probability and distance of interactions,
    adjusting them to represent the strength of attractive and repulsive forces. It ensures that LJ parameters are
    consistent with the given probability and distance thresholds, maintaining the accuracy of simulations or calculations.
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
        meGO_LJ["probability"] > meGO_LJ["limit_rc_att"] * np.maximum(meGO_LJ["rc_probability"], meGO_LJ["rc_threshold"])
    ) & (meGO_LJ["probability"] > meGO_LJ["md_threshold"])

    meGO_LJ.loc[condition, "epsilon"] = np.maximum(0.0, meGO_LJ["epsilon_prior"]) - (
        (meGO_LJ["epsilon_0"] - np.maximum(0.0, meGO_LJ["epsilon_prior"])) / np.log(meGO_LJ["rc_threshold"])
    ) * (np.log(meGO_LJ["probability"] / (np.maximum(meGO_LJ["rc_probability"], meGO_LJ["rc_threshold"]))))
    meGO_LJ.loc[condition, "learned"] = 1
    meGO_LJ.loc[condition, "sigma"] = meGO_LJ["distance"] / 2.0 ** (1.0 / 6.0)

    # Not-attractive interactions
    # this is used only when MD_th < MD_p < limit_rc_att*RC_p
    # negative epsilon are used to identify non-attractive interactions
    condition = (
        meGO_LJ["probability"] <= meGO_LJ["limit_rc_att"] * np.maximum(meGO_LJ["rc_probability"], meGO_LJ["rc_threshold"])
    ) & (meGO_LJ["probability"] > meGO_LJ["md_threshold"])
    meGO_LJ.loc[condition, "epsilon"] = -meGO_LJ["rep"] * (
        1.0 + (np.maximum(meGO_LJ["rc_probability"], meGO_LJ["rc_threshold"]) - meGO_LJ["probability"])
    )
    meGO_LJ.loc[condition, "learned"] = 1
    # for repulsive interaction we reset sigma to its effective value
    # this because when merging repulsive contacts from different sources what will matters
    # will be the repulsive strength that in this way is consistent
    meGO_LJ.loc[(meGO_LJ["epsilon"] < 0.0), "sigma"] = (-meGO_LJ["epsilon"]) ** (1.0 / 12.0) / (2.0 ** (1.0 / 6.0))

    # # 1-4 restored to the default values
    # meGO_LJ.loc[(meGO_LJ["bond_distance"] <7) & (meGO_LJ["same_chain"]), "sigma"]   =    meGO_LJ["mg_sigma"]
    # meGO_LJ.loc[(meGO_LJ["bond_distance"] <7) & (meGO_LJ["same_chain"]), "epsilon"] =  - meGO_LJ["rep"]

    # clean NaN and zeros
    meGO_LJ.dropna(subset=["epsilon"], inplace=True)
    meGO_LJ = meGO_LJ[meGO_LJ.epsilon != 0]

    return meGO_LJ


def apply_symmetries(meGO_ensemble, meGO_input, symmetry):
    """
    Apply symmetries to the molecular ensemble.

    This function applies symmetries to the input molecular ensemble based on the provided symmetry parameters.

    Parameters
    ----------
    meGO_ensemble : dict
        A dictionary containing relevant meGO data such as interactions and statistics within the molecular ensemble.
    meGO_input : pd.DataFrame
        Input DataFrame containing molecular ensemble data.
    parameters : dict
        A dictionary containing parameters parsed from the command-line, including symmetry information.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the molecular ensemble data with applied symmetries.
    """

    # Step 1: Initialize variables
    dict_sbtype_to_resname = meGO_ensemble["topology_dataframe"].set_index("sb_type")["resname"].to_dict()
    mglj_resn_ai = meGO_input["ai"].map(dict_sbtype_to_resname)
    mglj_resn_aj = meGO_input["aj"].map(dict_sbtype_to_resname)
    df_list = []

    # Step 2: Loop through symmetries and permutations
    for sym in symmetry:
        if not sym:
            continue
        # Pre-filter the DataFrame to speed up when there are multiple equivalent atoms
        meGO_filtered = meGO_input[(mglj_resn_ai == sym[0]) | (mglj_resn_aj == sym[0])]
        mgf_resn_ai = meGO_filtered["ai"].map(dict_sbtype_to_resname)
        mgf_resn_aj = meGO_filtered["aj"].map(dict_sbtype_to_resname)
        for atypes in itertools.permutations(sym[1:]):
            t_df_ai = meGO_filtered[meGO_filtered["ai"].str.startswith(f"{atypes[0]}_") & (mgf_resn_ai == sym[0])]
            t_df_aj = meGO_filtered[meGO_filtered["aj"].str.startswith(f"{atypes[0]}_") & (mgf_resn_aj == sym[0])]
            t_df_ai.loc[:, "ai"] = t_df_ai["ai"].str.replace(r"^(.*?)_", atypes[1] + "_", regex=True)
            t_df_aj.loc[:, "aj"] = t_df_aj["aj"].str.replace(r"^(.*?)_", atypes[1] + "_", regex=True)
            df_list.extend([t_df_ai, t_df_aj])

    # Step 3: Concatenate DataFrames
    df_tmp = pd.concat(df_list, ignore_index=True)
    df_list = []
    df_tmp.drop_duplicates(inplace=True)

    # Step 4: Filter and concatenate again
    df_resn_ai = df_tmp["ai"].map(dict_sbtype_to_resname)
    df_resn_aj = df_tmp["aj"].map(dict_sbtype_to_resname)

    for sym in symmetry:
        if not sym:
            continue
        # Pre-filter the DataFrame to speed up when there are multiple equivalent atoms
        df_tmp_filt = df_tmp[(df_resn_ai == sym[0]) | (df_resn_aj == sym[0])]
        df_resn_ai_f = df_tmp_filt["ai"].map(dict_sbtype_to_resname)
        df_resn_aj_f = df_tmp_filt["aj"].map(dict_sbtype_to_resname)
        for atypes in itertools.permutations(sym[1:]):
            t_df_ai = df_tmp_filt[df_tmp_filt["ai"].str.startswith(f"{atypes[0]}_") & (df_resn_ai_f == sym[0])]
            t_df_aj = df_tmp_filt[df_tmp_filt["aj"].str.startswith(f"{atypes[0]}_") & (df_resn_aj_f == sym[0])]
            t_df_ai.loc[:, "ai"] = t_df_ai["ai"].str.replace(r"^(.*?)_", atypes[1] + "_", regex=True)
            t_df_aj.loc[:, "aj"] = t_df_aj["aj"].str.replace(r"^(.*?)_", atypes[1] + "_", regex=True)
            df_list.extend([t_df_ai, t_df_aj])

    # Step 5: Concatenate and remove duplicates
    tmp_df = pd.concat(df_list + [df_tmp], ignore_index=True)
    tmp_df.drop_duplicates(inplace=True)

    return tmp_df


def generate_LJ(meGO_ensemble, train_dataset, parameters):
    """
    Generates LJ (Lennard-Jones) interactions and associated atomic contacts within a molecular ensemble.

    Parameters
    ----------
    meGO_ensemble : dict
        Contains relevant meGO data such as interactions and statistics within the molecular ensemble.
    train_dataset : pd.DataFrame
        DataFrame containing training dataset information for LJ interactions.
    parameters : dict
        Contains parameters parsed from the command-line.

    Returns
    -------
    meGO_LJ : pd.DataFrame
        Contains non-bonded atomic contacts associated with LJ parameters and statistics.
    meGO_LJ_14 : pd.DataFrame
        Contains 1-4 atomic contacts associated with LJ parameters and statistics.
    """

    st = time.time()
    print("\t- Set sigma and epsilon")
    # copy only learned contacts
    meGO_LJ = train_dataset[train_dataset["learned"]].copy()
    # meGO needed fields
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
        # "1-4",
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

    # generate attractive and repulsive interactions
    meGO_LJ = set_sig_epsilon(meGO_LJ, parameters)[needed_fields]

    et = time.time()
    elapsed_time = et - st
    st = et
    print("\t- Done in:", elapsed_time, "seconds")

    # apply symmetries for equivalent atoms
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

    # Merging of multiple simulations:
    # 1. learned over not learned
    # 2. attractive over repulsive
    # 3. shorter over longer
    # 4. stronger over weaker attractive
    # 5. wearker over stronger repulsive
    meGO_LJ["type"] = np.sign(meGO_LJ["epsilon"])
    meGO_LJ.sort_values(
        by=["ai", "aj", "same_chain", "learned", "type", "sigma", "epsilon"],
        ascending=[True, True, True, False, False, True, False],
        inplace=True,
    )
    # Cleaning the duplicates
    meGO_LJ = meGO_LJ.drop_duplicates(subset=["ai", "aj", "same_chain"], keep="first")

    # now we can remove contacts with default c6/c12 becasue these
    # are uninformative and predefined. This also allow to replace them with contact learned
    # by either intra/inter training. We cannot remove 1-4 interactions.
    # we should not remove default interactions in the window of 2 neighor AA to
    # avoid replacing them with unwanted interactions

    # this removes attractive/repulsive contacts that are default
    meGO_LJ = meGO_LJ.loc[
        ~(
            (meGO_LJ["epsilon"] > 0)
            & (meGO_LJ["mg_epsilon"] > 0)
            & ((abs(meGO_LJ["epsilon"] - meGO_LJ["mg_epsilon"]) / meGO_LJ["mg_epsilon"]) < parameters.relative_c12d)
            & ((abs(meGO_LJ["sigma"] - meGO_LJ["mg_sigma"]) / meGO_LJ["mg_sigma"]) < parameters.relative_c12d)
            & ((meGO_LJ["bond_distance"] > 3) | (~meGO_LJ["same_chain"]))
            # & (meGO_LJ["1-4"] == "1>4")
        )
    ]
    meGO_LJ = meGO_LJ.loc[
        ~(
            (meGO_LJ["epsilon"] < 0)
            & (meGO_LJ["mg_epsilon"] < 0)
            & ((abs(meGO_LJ["epsilon"] - meGO_LJ["mg_epsilon"]) / abs(meGO_LJ["mg_epsilon"])) < parameters.relative_c12d)
            # & (meGO_LJ["1-4"] == "1>4")
            & ((meGO_LJ["bond_distance"] > 3) | (~meGO_LJ["same_chain"]))
            & ~(
                (meGO_LJ["bond_distance"] < 7)
                # (abs(meGO_LJ["ai"].apply(get_residue_number) - meGO_LJ["aj"].apply(get_residue_number)) < 3)
                & (meGO_LJ["same_chain"])
            )
        )
    ]

    meGO_LJ = meGO_LJ[needed_fields]

    # now is a good time to acquire statistics on the parameters
    # this should be done per interaction pair (cycling over all molecules combinations) and inter/intra/intra_d
    stat_str = io.print_stats(meGO_LJ)

    # Here we create a copy of contacts to be added in pairs-exclusion section in topol.top.
    meGO_LJ_14 = meGO_LJ.copy()

    # Sorting the pairs prioritising intermolecular interactions
    meGO_LJ.sort_values(by=["ai", "aj", "same_chain"], ascending=[True, True, True], inplace=True)
    meGO_LJ = meGO_LJ.drop_duplicates(subset=["ai", "aj"], keep="first")

    # Pairs prioritise intramolecular interactions
    meGO_LJ_14.sort_values(by=["ai", "aj", "same_chain"], ascending=[True, True, False], inplace=True)
    meGO_LJ_14 = meGO_LJ_14.drop_duplicates(subset=["ai", "aj"], keep="first")

    # where meGO_LJ_14 is the same of meGO_LJ and same_chain is yes that the line can be dropped
    # that is I want to keep lines with same_chain no or lines with same chain yes that have same_chain no in meGO_LJ
    test = pd.merge(meGO_LJ_14, meGO_LJ, how="right", on=["ai", "aj"])
    meGO_LJ_14 = test.loc[(~test["same_chain_x"]) | ((test["same_chain_x"]) & (~test["same_chain_y"]))]
    # removes columns ending with _y
    meGO_LJ_14 = meGO_LJ_14.loc[:, ~meGO_LJ_14.columns.str.endswith("_y")]
    # rename the columns _x
    meGO_LJ_14.columns = meGO_LJ_14.columns.str.rstrip("_x")

    # remove intermolecular interactions across molecules from meGO_LJ_14
    meGO_LJ_14 = meGO_LJ_14[meGO_LJ_14["molecule_name_ai"] == meGO_LJ_14["molecule_name_aj"]]

    # copy 1-4 interactions into meGO_LJ_14
    # copy14 = meGO_LJ.loc[(meGO_LJ["1-4"] == "1_4")]
    copy14 = meGO_LJ.loc[(meGO_LJ["bond_distance"] == 3) & (meGO_LJ["same_chain"])]
    meGO_LJ_14 = pd.concat([meGO_LJ_14, copy14], axis=0, sort=False, ignore_index=True)
    # remove them from the default force-field
    # meGO_LJ = meGO_LJ.loc[(meGO_LJ["1-4"] != "1_4")]
    meGO_LJ = meGO_LJ.loc[~((meGO_LJ["bond_distance"] == 3) & (meGO_LJ["same_chain"]))]

    if not parameters.single_molecule:
        # neighbour intramolecular interactions are not used as intermolecular
        copy_intra = meGO_LJ.loc[
            (meGO_LJ["same_chain"])
            & (meGO_LJ["bond_distance"] < 7)
            # & (abs(meGO_LJ["ai"].apply(get_residue_number) - meGO_LJ["aj"].apply(get_residue_number)) < 3)
        ]
        meGO_LJ_14 = pd.concat([meGO_LJ_14, copy_intra], axis=0, sort=False, ignore_index=True)
        # remove them from the default force-field
        meGO_LJ = meGO_LJ.loc[
            ~(
                (meGO_LJ["same_chain"])
                & (meGO_LJ["bond_distance"] < 7)
                # & (abs(meGO_LJ["ai"].apply(get_residue_number) - meGO_LJ["aj"].apply(get_residue_number)) < 3)
            )
        ]

    # now we can decide to keep intermolecular interactions as intramolecular ones
    # to do this is enough to remove it from meGO_LJ_14, in this way the value used for the contact is the one meGO_LJ
    if not parameters.force_split:
        # Filter rows in meGO_LJ_14 that meet the condition
        meGO_LJ_14 = meGO_LJ_14.loc[
            ~(
                (~meGO_LJ_14["same_chain"])
                & (meGO_LJ_14["molecule_name_ai"] == meGO_LJ_14["molecule_name_aj"])
                & (meGO_LJ_14["epsilon"] > 0.0)
                & (meGO_LJ_14["bond_distance"] > 6)
                # & (abs(meGO_LJ_14["ai"].apply(get_residue_number) - meGO_LJ_14["aj"].apply(get_residue_number)) > 2)
            )
        ]
    else:
        split_ii = meGO_LJ.loc[(meGO_LJ["same_chain"])]
        # move the intramolecular interaction in the topology
        meGO_LJ_14 = pd.concat([meGO_LJ_14, split_ii], axis=0, sort=False, ignore_index=True)
        # remove them from the default force-field
        meGO_LJ = meGO_LJ.loc[(~meGO_LJ["same_chain"])]

    # Now is time to add masked default interactions for pairs
    # that have not been learned in any other way
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
        # "1-4",
        "bond_distance",
        "rep",
        "mg_sigma",
        "mg_epsilon",
        "md_threshold",
        "rc_threshold",
        "learned",
    ]
    basic_LJ = generate_MG_LJ(meGO_ensemble)[needed_fields]
    meGO_LJ = pd.concat([meGO_LJ, basic_LJ])

    # make meGO_LJ fully symmetric
    # Create inverse DataFrame
    inverse_meGO_LJ = meGO_LJ.rename(
        columns={"ai": "aj", "aj": "ai", "molecule_name_ai": "molecule_name_aj", "molecule_name_aj": "molecule_name_ai"}
    ).copy()
    # Concatenate original and inverse DataFrames
    # Here we have a fully symmetric matrix for both intra/intersame/intercross
    meGO_LJ = pd.concat([meGO_LJ, inverse_meGO_LJ], axis=0, sort=False, ignore_index=True)

    meGO_LJ["ai"] = meGO_LJ["ai"].astype("category")
    meGO_LJ["aj"] = meGO_LJ["aj"].astype("category")
    meGO_LJ["molecule_name_ai"] = meGO_LJ["molecule_name_ai"].astype("category")
    meGO_LJ["molecule_name_aj"] = meGO_LJ["molecule_name_aj"].astype("category")

    # Sorting the pairs prioritising learned interactions
    meGO_LJ.sort_values(by=["ai", "aj", "learned", "sigma"], ascending=[True, True, False, True], inplace=True)
    # Cleaning the duplicates, that is that we retained a not learned interaction only if it is unique
    # first we remove duplicated masked interactions
    meGO_LJ = meGO_LJ.drop_duplicates(subset=["ai", "aj"], keep="first")
    # meGO_LJ = meGO_LJ.loc[(~(meGO_LJ.duplicated(subset=["ai", "aj"], keep="first")))]

    # we are ready to finalize the setup
    # Calculate c6 and c12 for meGO_LJ
    meGO_LJ["c6"] = np.where(meGO_LJ["epsilon"] < 0.0, 0.0, 4 * meGO_LJ["epsilon"] * (meGO_LJ["sigma"] ** 6))

    meGO_LJ["c12"] = np.where(meGO_LJ["epsilon"] < 0.0, -meGO_LJ["epsilon"], 4 * meGO_LJ["epsilon"] * (meGO_LJ["sigma"] ** 12))

    # Calculate c6 and c12 for meGO_LJ_14
    meGO_LJ_14["c6"] = np.where(meGO_LJ_14["epsilon"] < 0.0, 0.0, 4 * meGO_LJ_14["epsilon"] * (meGO_LJ_14["sigma"] ** 6))

    meGO_LJ_14["c12"] = np.where(
        meGO_LJ_14["epsilon"] < 0.0, -meGO_LJ_14["epsilon"], 4 * meGO_LJ_14["epsilon"] * (meGO_LJ_14["sigma"] ** 12)
    )

    et = time.time()
    elapsed_time = et - st
    print("\t- Done in:", elapsed_time, "seconds")

    return meGO_LJ, meGO_LJ_14, stat_str


def sort_LJ(meGO_ensemble, meGO_LJ):
    # Add or modify columns in the original DataFrame
    meGO_LJ["type"] = 1
    meGO_LJ["number_ai"] = meGO_LJ["ai"].map(meGO_ensemble["sbtype_number_dict"]).astype(int)
    meGO_LJ["number_aj"] = meGO_LJ["aj"].map(meGO_ensemble["sbtype_number_dict"]).astype(int)

    # Filter and explicitly create a copy to avoid the warning
    meGO_LJ = meGO_LJ[(meGO_LJ["ai"].cat.codes <= meGO_LJ["aj"].cat.codes)].copy()

    # across molecules use molecule_ai<=molecule_aj
    (
        meGO_LJ["ai"],
        meGO_LJ["aj"],
        meGO_LJ["molecule_name_ai"],
        meGO_LJ["molecule_name_aj"],
        meGO_LJ["number_ai"],
        meGO_LJ["number_aj"],
    ) = np.where(
        (meGO_LJ["molecule_name_ai"].astype(str) <= meGO_LJ["molecule_name_aj"].astype(str)),
        [
            meGO_LJ["ai"],
            meGO_LJ["aj"],
            meGO_LJ["molecule_name_ai"],
            meGO_LJ["molecule_name_aj"],
            meGO_LJ["number_ai"],
            meGO_LJ["number_aj"],
        ],
        [
            meGO_LJ["aj"],
            meGO_LJ["ai"],
            meGO_LJ["molecule_name_aj"],
            meGO_LJ["molecule_name_ai"],
            meGO_LJ["number_aj"],
            meGO_LJ["number_ai"],
        ],
    )

    # in the same molecule use ai<=aj
    # Apply np.where to swap values only when molecule_name_ai == molecule_name_aj
    (
        meGO_LJ["ai"],
        meGO_LJ["aj"],
        meGO_LJ["molecule_name_ai"],
        meGO_LJ["molecule_name_aj"],
        meGO_LJ["number_ai"],
        meGO_LJ["number_aj"],
    ) = np.where(
        meGO_LJ["molecule_name_ai"] == meGO_LJ["molecule_name_aj"],  # Only apply when names are equal
        np.where(
            meGO_LJ["number_ai"] <= meGO_LJ["number_aj"],  # Condition to check number_ai vs number_aj
            [
                meGO_LJ["ai"],
                meGO_LJ["aj"],
                meGO_LJ["molecule_name_ai"],
                meGO_LJ["molecule_name_aj"],
                meGO_LJ["number_ai"],
                meGO_LJ["number_aj"],
            ],
            [
                meGO_LJ["aj"],
                meGO_LJ["ai"],
                meGO_LJ["molecule_name_aj"],
                meGO_LJ["molecule_name_ai"],
                meGO_LJ["number_aj"],
                meGO_LJ["number_ai"],
            ],
        ),
        # If molecule_name_ai != molecule_name_aj, keep the values unchanged
        [
            meGO_LJ["ai"],
            meGO_LJ["aj"],
            meGO_LJ["molecule_name_ai"],
            meGO_LJ["molecule_name_aj"],
            meGO_LJ["number_ai"],
            meGO_LJ["number_aj"],
        ],
    )

    meGO_LJ.sort_values(by=["molecule_name_ai", "molecule_name_aj", "number_ai", "number_aj"], inplace=True)

    final_fields = [
        "ai",
        "aj",
        "type",
        "c6",
        "c12",
        "sigma",
        "epsilon",
        "probability",
        "rc_probability",
        "md_threshold",
        "rc_threshold",
        "same_chain",
        "source",
        "bond_distance",
        "number_ai",
        "number_aj",
    ]

    meGO_LJ = meGO_LJ[final_fields]

    return meGO_LJ


def get_residue_number(s):
    return int(s.split("_")[-1])


def make_pairs_exclusion_topology(meGO_ensemble, meGO_LJ_14, args):
    """
    This function prepares the [ exclusion ] and [ pairs ] section to output to topology.top

    Parameters
    ----------
    meGO_ensemlbe : dict
        The meGO_ensemble object contains all the relevant system information
    meGO_LJ_14 : pd.DataFrame
        Contains the contact information for the 1-4 interactions

    Returns
    -------
    pairs_molecule_dict : dict
        Contains the "write out"-ready pairs-exclusions interactions for each molecule
    """
    # pairs and exclusions are built per molecule type and saved in a dictionary
    pairs_molecule_dict = {}
    for idx, (molecule, bond_pair) in enumerate(meGO_ensemble["bond_pairs"].items(), start=1):
        reduced_topology = (
            meGO_ensemble["topology_dataframe"]
            .loc[meGO_ensemble["topology_dataframe"]["molecule_name"] == molecule][
                [
                    "number",
                    "sb_type",
                    "resnum",
                    "name",
                    "type",
                    "resname",
                    "molecule_type",
                ]
            ]
            .copy()
        )

        reduced_topology["number"] = reduced_topology["number"].astype(str)
        # reduced_topology["resnum"] = reduced_topology["resnum"].astype(int)

        type_atnum_dict = reduced_topology.astype({"number": int}).set_index("number")["sb_type"].to_dict()
        atnum_type_dict = reduced_topology.set_index("sb_type")["number"].to_dict()
        # resnum_type_dict = reduced_topology.set_index("sb_type")["resnum"].to_dict()

        # Building the exclusion bonded list
        # exclusion_bonds are all the interactions within 3 bonds
        # p14 are specifically the interactions at exactly 3 bonds
        # b6 are all the interactions within 6 bonds
        exclusion_bonds, p14, b6 = topology.generate_bond_exclusions(reduced_topology, bond_pair)

        pairs = pd.DataFrame()
        # in the case of the MG prior we need to remove interactions in a window of 2 residues
        if args.egos == "mg":
            filtered_combinations = []
            for pair in b6:
                a_str, b_str = pair.split("_")
                a, b = int(a_str), int(b_str)

                if a not in type_atnum_dict or b not in type_atnum_dict:
                    continue

                ai = type_atnum_dict[a]
                aj = type_atnum_dict[b]

                # Hydrogen filtering
                if (
                    meGO_ensemble["sbtype_type_dict"][ai] == "H"
                    and meGO_ensemble["sbtype_type_dict"][aj] not in {"H", "O", "OM", "OA"}
                ) or (
                    meGO_ensemble["sbtype_type_dict"][aj] == "H"
                    and meGO_ensemble["sbtype_type_dict"][ai] not in {"H", "O", "OM", "OA"}
                ):
                    continue

                filtered_combinations.append((ai, aj))

            ## Create a list of tuples (sbtype, residue_number)
            # sbtype_with_residue = [(sbtype, resnum_type_dict[sbtype]) for sbtype in reduced_topology["sb_type"]]
            ## Sort the list by residue numbers
            # sbtype_with_residue.sort(key=lambda x: x[1])
            ## Initialize a list to hold the filtered combinations
            # filtered_combinations = []
            ## Use two pointers to find valid pairs
            # n = len(sbtype_with_residue)
            # for i in range(n):
            #    j = i + 1  # Start with the current sbtype
            #    # Find the range of valid sbtypes
            #    while j < n and abs(sbtype_with_residue[j][1] - sbtype_with_residue[i][1]) <= 2:
            #        filtered_combinations.append((sbtype_with_residue[i][0], sbtype_with_residue[j][0]))
            #        j += 1

            ## Filter out invalid combinations
            # valid_combinations = [
            #    (ai, aj)
            #    for ai, aj in filtered_combinations
            #    # this is to remove all interaction of H with the rest exept for O, OM, and OA
            #    if not (
            #        (
            #            meGO_ensemble["sbtype_type_dict"][ai] == "H"
            #            and meGO_ensemble["sbtype_type_dict"][aj] not in {"H", "O", "OM", "OA"}
            #        )
            #        or (
            #            meGO_ensemble["sbtype_type_dict"][aj] == "H"
            #            and meGO_ensemble["sbtype_type_dict"][ai] not in {"H", "O", "OM", "OA"}
            #        )
            #    )
            # ]

            # Create a DataFrame from the filtered combinations
            # df = pd.DataFrame(valid_combinations, columns=["ai", "aj"])
            df = pd.DataFrame(filtered_combinations, columns=["ai", "aj"])
            df["c6"] = 0.0
            df["c12"] = np.sqrt(
                df["ai"].map(meGO_ensemble["sbtype_c12_dict"]) * df["aj"].map(meGO_ensemble["sbtype_c12_dict"])
            )

            df.loc[
                (
                    (df["ai"].map(meGO_ensemble["sbtype_type_dict"]) == "OM")
                    | (df["ai"].map(meGO_ensemble["sbtype_type_dict"]) == "O")
                )
                & (
                    (df["aj"].map(meGO_ensemble["sbtype_type_dict"]) == "OM")
                    | (df["aj"].map(meGO_ensemble["sbtype_type_dict"]) == "O")
                ),
                "c12",
            ] = type_definitions.mg_OO_c12_rep

            df.loc[
                ((df["ai"].map(meGO_ensemble["sbtype_type_dict"]) == "H"))
                & ((df["aj"].map(meGO_ensemble["sbtype_type_dict"]) == "H")),
                "c12",
            ] = type_definitions.mg_HH_c12_rep

            df.loc[
                (
                    (df["ai"].map(meGO_ensemble["sbtype_type_dict"]) == "NL")
                    | (df["ai"].map(meGO_ensemble["sbtype_type_dict"]) == "NZ")
                )
                & (
                    (df["aj"].map(meGO_ensemble["sbtype_type_dict"]) == "NL")
                    | (df["aj"].map(meGO_ensemble["sbtype_type_dict"]) == "NZ")
                ),
                "c12",
            ] = type_definitions.mg_NN_c12_rep

            df.loc[
                (
                    (
                        (df["ai"].map(meGO_ensemble["sbtype_type_dict"]) == "OM")
                        | (df["ai"].map(meGO_ensemble["sbtype_type_dict"]) == "O")
                    )
                    & ((df["aj"].map(meGO_ensemble["sbtype_type_dict"]) == "N"))
                )
                | (
                    ((df["ai"].map(meGO_ensemble["sbtype_type_dict"]) == "N"))
                    & (
                        (df["aj"].map(meGO_ensemble["sbtype_type_dict"]) == "OM")
                        | (df["aj"].map(meGO_ensemble["sbtype_type_dict"]) == "O")
                    )
                ),
                "c12",
            ] = type_definitions.mg_ON_c12_rep

            df["same_chain"] = True
            df["probability"] = 1.0
            df["rc_probability"] = 1.0
            df["source"] = "mg"
            df["rep"] = df["c12"]
            df["1-4"] = "1>4"
            # The exclusion list was made based on the atom number
            df["check"] = df["ai"].map(atnum_type_dict).astype(str) + "_" + df["aj"].map(atnum_type_dict).astype(str)
            # Here the drop the contacts which are already defined by GROMACS, including the eventual 1-4 exclusion defined in the LJ_df
            df["remove"] = ""
            df.loc[(df["check"].isin(exclusion_bonds)), "remove"] = "Yes"
            df.loc[(df["check"].isin(p14) & (df["same_chain"])), "remove"] = "Yes"
            mask = df.remove == "Yes"
            df = df[~mask]
            df.drop(columns=["check", "remove"], inplace=True)
            pairs = pd.concat([meGO_LJ_14, df], axis=0, sort=False, ignore_index=True)

        elif args.egos == "production" and not meGO_LJ_14.empty:
            mol_ai = f"{idx}_{molecule}"
            # pairs do not have duplicates because these have been cleaned before
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
            # Intermolecular interactions are excluded
            # this need to be the default repulsion if within two residue
            if not pairs.empty:
                pairs.loc[
                    (~pairs["same_chain"]) & (pairs["bond_distance"] < 7),
                    # & (abs(pairs["ai"].apply(get_residue_number) - pairs["aj"].apply(get_residue_number)) < 3),
                    "c6",
                ] = 0.0
                pairs.loc[
                    (~pairs["same_chain"]) & (pairs["bond_distance"] < 7),
                    # & (abs(pairs["ai"].apply(get_residue_number) - pairs["aj"].apply(get_residue_number)) < 3),
                    "c12",
                ] = pairs["rep"]
                # else it should be default mg
                pairs.loc[
                    (~pairs["same_chain"]) & (pairs["bond_distance"] > 6)
                    # & (abs(pairs["ai"].apply(get_residue_number) - pairs["aj"].apply(get_residue_number)) > 2)
                    & (pairs["mg_epsilon"] < 0.0),
                    "c6",
                ] = 0.0
                pairs.loc[
                    (~pairs["same_chain"]) & (pairs["bond_distance"] > 6)
                    # & (abs(pairs["ai"].apply(get_residue_number) - pairs["aj"].apply(get_residue_number)) > 2)
                    & (pairs["mg_epsilon"] < 0.0),
                    "c12",
                ] = -pairs["mg_epsilon"]
                pairs.loc[
                    (~pairs["same_chain"]) & (pairs["bond_distance"] > 6)
                    # & (abs(pairs["ai"].apply(get_residue_number) - pairs["aj"].apply(get_residue_number)) > 2)
                    & (pairs["mg_epsilon"] > 0.0),
                    "c6",
                ] = (
                    4 * pairs["mg_epsilon"] * (pairs["mg_sigma"] ** 6)
                )
                pairs.loc[
                    (~pairs["same_chain"]) & (pairs["bond_distance"] > 6)
                    # & (abs(pairs["ai"].apply(get_residue_number) - pairs["aj"].apply(get_residue_number)) > 2)
                    & (pairs["mg_epsilon"] > 0.0),
                    "c12",
                ] = (
                    4 * pairs["mg_epsilon"] * (pairs["mg_sigma"] ** 12)
                )

        # now we are ready to finalize
        if not pairs.empty:
            # The exclusion list was made based on the atom number
            pairs["ai"] = pairs["ai"].map(atnum_type_dict)
            pairs["aj"] = pairs["aj"].map(atnum_type_dict)
            pairs["check"] = pairs["ai"].astype(str) + "_" + pairs["aj"].astype(str)
            # Here the drop the contacts which are already defined by GROMACS, including the eventual 1-4 exclusion defined in the LJ_pairs
            pairs["remove"] = ""
            pairs.loc[(pairs["check"].isin(exclusion_bonds)), "remove"] = "Yes"
            pairs.loc[(pairs["check"].isin(p14) & (pairs["same_chain"])), "remove"] = "No"
            mask = pairs.remove == "Yes"
            pairs = pairs[~mask]
            # finalize
            pairs["func"] = 1
            # this is a safety check
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

            # Here we want to sort so that ai is smaller than aj
            inv_pairs = pairs[
                [
                    "aj",
                    "ai",
                    "func",
                    "c6",
                    "c12",
                    "probability",
                    "rc_probability",
                    "source",
                ]
            ].copy()
            inv_pairs.columns = [
                "ai",
                "aj",
                "func",
                "c6",
                "c12",
                "probability",
                "rc_probability",
                "source",
            ]
            pairs = pd.concat([pairs, inv_pairs], axis=0, sort=False, ignore_index=True)
            pairs = pairs[pairs["ai"] < pairs["aj"]]
            pairs.drop_duplicates(inplace=True, ignore_index=True)
            pairs.sort_values(by=["ai", "aj"], inplace=True)

        pairs_molecule_dict[molecule] = pairs

    return pairs_molecule_dict
