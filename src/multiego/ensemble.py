from .resources import type_definitions
from . import io
from . import topology
from .util import masking

import glob
import numpy as np
import pandas as pd
import parmed
import os
import warnings
import itertools
import time


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

    Args:
    - topology (object): An object containing information about the molecules.

    Returns:
    - ensemble_topology_dataframe (DataFrame): DataFrame containing ensemble topology information.
    - ensemble_molecules_idx_sbtype_dictionary (dict): Dictionary mapping molecule indexes to their respective subtypes.
    - sbtype_c12_dict (dict): Dictionary mapping subtype to their c12 values.
    - sbtype_name_dict (dict): Dictionary mapping subtype to their names.
    - sbtype_moltype_dict (dict): Dictionary mapping subtype to their molecule types.
    - molecule_type_dict (dict): Dictionary mapping molecule names to their types.

    This function initializes a topology DataFrame by extracting information about molecules and their atoms.
    It creates a DataFrame containing details about atoms, molecules, their types, and assigns specific values based on the provided information.
    The function also generates dictionaries mapping different atom subtypes to their respective characteristics and molecule types.

    Note:
    - The 'topology' object is expected to contain molecule information.
    - The returned DataFrame and dictionaries provide comprehensive details about the molecular structure and characteristics.
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
    ensemble_topology_dataframe.rename(columns={"epsilon": "c12"}, inplace=True)

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
    ensemble_topology_dataframe["c6"] = ensemble_topology_dataframe["type"].map(atp_mg_c6_map)
    ensemble_topology_dataframe["c12"] = ensemble_topology_dataframe["type"].map(atp_mg_c12_map)
    ensemble_topology_dataframe["molecule_type"] = ensemble_topology_dataframe["molecule_name"].map(molecule_type_dict)

    for molecule in ensemble_molecules_idx_sbtype_dictionary.keys():
        temp_topology_dataframe = ensemble_topology_dataframe.loc[ensemble_topology_dataframe["molecule"] == molecule]
        number_sbtype_dict = temp_topology_dataframe[["number", "sb_type"]].set_index("number")["sb_type"].to_dict()
        ensemble_molecules_idx_sbtype_dictionary[molecule] = number_sbtype_dict

    sbtype_c12_dict = ensemble_topology_dataframe[["sb_type", "rc_c12"]].set_index("sb_type")["rc_c12"].to_dict()
    sbtype_name_dict = ensemble_topology_dataframe[["sb_type", "name"]].set_index("sb_type")["name"].to_dict()
    sbtype_moltype_dict = (
        ensemble_topology_dataframe[["sb_type", "molecule_type"]].set_index("sb_type")["molecule_type"].to_dict()
    )

    return (
        ensemble_topology_dataframe,
        ensemble_molecules_idx_sbtype_dictionary,
        sbtype_c12_dict,
        sbtype_name_dict,
        sbtype_moltype_dict,
        molecule_type_dict,
    )


def initialize_molecular_contacts(contact_matrix, args):
    """
    This function initializes a contact matrix for a given simulation.

    Parameters
    ----------
    contact_matrix : pd.DataFrame
        Contains contact information read from intra-/intermat
    path : str
        Path to the simulation folder
    ensemble_molecules_idx_sbtype_dictionary : dict
        Associates atom indices to atoms named according to multi-eGO conventions
    simulation : str
        The simulation classified equivalent to the input folder
    args : argparse.Namespace
        Parsed arguments

    Returns
    -------
    contact_matrix : pd.DataFrame
        A contact matrix containing contact data for each of the different simulations
    """
    # calculate adaptive rc/md threshold
    # sort probabilities, and calculate the normalized cumulative distribution
    p_sort = np.sort(contact_matrix["probability"].loc[(contact_matrix["intra_domain"])].to_numpy())[::-1]
    p_sort_id = np.sort(contact_matrix["probability"].loc[~(contact_matrix["intra_domain"])].to_numpy())[::-1]
    norm = np.sum(p_sort)
    norm_id = np.sum(p_sort_id)
    if norm == 0:
        p_sort_normalized = 0
        md_threshold = 1
    else:
        # find md threshold
        p_sort_normalized = np.cumsum(p_sort) / norm
        md_threshold = p_sort[np.min(np.where(p_sort_normalized > args.p_to_learn)[0])]

    if norm_id == 0:
        p_sort_normalized = 0
        md_threshold_id = 1
    else:
        # find md threshold
        p_sort_normalized = np.cumsum(p_sort_id) / norm_id
        md_threshold_id = p_sort_id[np.min(np.where(p_sort_normalized > args.p_to_learn)[0])]

    # set zf this simplify a lot of the following code
    # for intra-domain
    contact_matrix.loc[(contact_matrix["same_chain"]) & (contact_matrix["intra_domain"]), "zf"] = args.f
    # for inter-domain
    contact_matrix.loc[(contact_matrix["same_chain"]) & (~contact_matrix["intra_domain"]), "zf"] = args.inter_domain_f
    # for inter-molecular
    contact_matrix.loc[(~contact_matrix["same_chain"]), "zf"] = args.inter_f

    # Extract the second part of the strings after splitting
    molecule_ai = contact_matrix["molecule_name_ai"].str.split("_", expand=True, n=1)[1]
    molecule_aj = contact_matrix["molecule_name_aj"].str.split("_", expand=True, n=1)[1]

    # Combine both Series into a single Series and find unique values
    molecules = pd.Series(np.concatenate([molecule_ai, molecule_aj])).unique()

    # Now determine mol_1 and mol_2
    mol_1 = molecules[0]
    mol_2 = mol_1 if len(molecules) == 1 else molecules[1]

    contact_matrix.loc[(contact_matrix["same_chain"]) & (contact_matrix["intra_domain"]), "epsilon_0"] = (
        args.multi_epsilon_intra[mol_1]
    )
    contact_matrix.loc[(contact_matrix["same_chain"]) & (~contact_matrix["intra_domain"]), "epsilon_0"] = (
        args.multi_epsilon_inter_domain[mol_1]
    )
    contact_matrix.loc[(~contact_matrix["same_chain"]), "epsilon_0"] = args.multi_epsilon_inter[mol_1][mol_2]

    # add the columns for rc, md threshold
    contact_matrix["md_threshold"] = md_threshold
    contact_matrix.loc[~(contact_matrix["intra_domain"]), "md_threshold"] = md_threshold_id
    contact_matrix["rc_threshold"] = contact_matrix["md_threshold"] ** (
        1.0 / (1.0 - (args.epsilon_min / contact_matrix["epsilon_0"]))
    )
    contact_matrix["limit_rc"] = (
        1.0
        / contact_matrix["rc_threshold"] ** (args.epsilon_min / contact_matrix["epsilon_0"])
        * contact_matrix["zf"] ** (1 - (args.epsilon_min / contact_matrix["epsilon_0"]))
    )

    # TODO think on the limits of f (should be those for which all repulsive/attractive interactions are removed)
    f_min = md_threshold

    if args.f != 1:
        tmp_f_max = contact_matrix["rc_threshold"].loc[(contact_matrix["same_chain"]) & (contact_matrix["intra_domain"])]
        if not tmp_f_max.empty:
            f_max = 1.0 / tmp_f_max.iloc[0]
            if args.f > f_max:
                print(f"f is not in the correct range:\n f_max={f_max} > f={args.f} > f_min={f_min}. Choose a proper value")
                exit()

    if args.inter_f != 1:
        tmp_f_max = contact_matrix["rc_threshold"].loc[(~contact_matrix["same_chain"])]
        if not tmp_f_max.empty:
            f_max = 1.0 / tmp_f_max.iloc[0]
            if args.inter_f > f_max:
                print(
                    f"f is not in the correct range:\n f_max={f_max} > f={args.inter_f} > f_min={f_min}. Choose a proper value"
                )
                exit()

    if args.inter_domain_f != 1:
        tmp_f_max = contact_matrix["rc_threshold"].loc[(contact_matrix["same_chain"]) & (~contact_matrix["intra_domain"])]
        if not tmp_f_max.empty:
            f_max = 1.0 / tmp_f_max.iloc[0]
            if args.inter_domain_f > f_max:
                print(
                    f"f is not in the correct range:\n f_max={f_max} > f={args.inter_domain_f} > f_min={f_min}. Choose a proper value"
                )
                exit()

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
    for reference in args.reference:  # reference_paths:
        print("\t-", f"Initializing {reference} ensemble data")
        reference_path = f"{args.root_dir}/inputs/{args.system}/{reference}"
        # path = f"{args.root_dir}/inputs/{args.system}/{reference_path}"
        topology_path = f"{reference_path}/topol.top"
        matrix_paths = glob.glob(f"{reference_path}/int??mat_?_?.ndx")
        matrix_paths = matrix_paths + glob.glob(f"{reference_path}/int??mat_?_?.ndx.gz")
        matrix_paths = matrix_paths + glob.glob(f"{reference_path}/int??mat_?_?.ndx.h5")
        if matrix_paths == []:
            raise FileNotFoundError(
                "Contact matrix file(s) must be named as intramat_X_X.ndx(.gz/.h5) or intermat_X_Y.ndx(.gz/.h5)"
            )
        for path in matrix_paths:
            name = path.replace(f"{args.root_dir}/inputs/", "")
            name = name.replace("/", "_")
            name = name.replace(".ndx", "")
            name = name.replace(".gz", "")
            name = name.replace(".h5", "")
            reference_contact_matrices[name] = io.read_molecular_contacts(
                path, ensemble["molecules_idx_sbtype_dictionary"], reference, path.endswith(".h5")
            )
            reference_contact_matrices[name] = reference_contact_matrices[name].add_prefix("rc_")

        et = time.time()
        elapsed_time = et - st
        st = et
        print("\t- Done in:", elapsed_time, "seconds")

    matrices["reference_matrices"] = reference_contact_matrices
    reference_set = set(ensemble["topology_dataframe"]["name"].to_list())

    # now we process the train contact matrices
    train_contact_matrices = {}
    train_topology_dataframe = pd.DataFrame()
    for simulation in args.train:
        print("\t-", f"Initializing {simulation} ensemble data")
        simulation_path = f"{args.root_dir}/inputs/{args.system}/{simulation}"
        topology_path = f"{simulation_path}/topol.top"
        if not os.path.isfile(topology_path):
            raise FileNotFoundError(f"{topology_path} not found.")

        print("\t\t-", f"Reading {topology_path}")
        # ignore the dihedral type overriding in parmed
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            topology = parmed.load_file(topology_path)

        (
            temp_topology_dataframe,
            ensemble["molecules_idx_sbtype_dictionary"],
            _,
            _,
            _,
            _,
        ) = initialize_topology(topology, custom_dict, args)

        train_topology_dataframe = pd.concat(
            [train_topology_dataframe, temp_topology_dataframe],
            axis=0,
            ignore_index=True,
        )
        matrix_paths = glob.glob(f"{simulation_path}/int??mat_?_?.ndx")
        matrix_paths = matrix_paths + glob.glob(f"{simulation_path}/int??mat_?_?.ndx.gz")
        matrix_paths = matrix_paths + glob.glob(f"{simulation_path}/int??mat_?_?.ndx.h5")
        if matrix_paths == []:
            raise FileNotFoundError(
                "Contact matrix file(s) must be named as intramat_X_X.ndx(.gz/.h5) or intermat_X_Y.ndx(.gz/.h5)"
            )
        for path in matrix_paths:
            name = path.replace(f"{args.root_dir}/inputs/", "")
            name = name.replace("/", "_")
            name = name.replace(".ndx", "")
            name = name.replace(".gz", "")
            name = name.replace(".h5", "")
            train_contact_matrices[name] = io.read_molecular_contacts(
                path, ensemble["molecules_idx_sbtype_dictionary"], simulation, path.endswith(".h5")
            )
            train_contact_matrices[name] = initialize_molecular_contacts(
                train_contact_matrices[name],
                args,
            )
            # ref_name = reference_path + "_" + path.split("/")[-1]
            # find corresponding reference matrix (given by the number_number at the end of the name)
            # using reference contact matrices
            identifier = (
                f'_{("_").join(path.split("/")[-1].replace(".ndx", "").replace(".gz", "").replace(".h5", "").split("_")[-3:])}'
            )
            ref_name = [key for key in reference_contact_matrices.keys() if key.endswith(identifier)]
            if ref_name == []:
                raise FileNotFoundError(f"No corresponding reference matrix found for {path}")
            ref_name = ref_name[0]
            ensemble["train_matrix_tuples"].append((name, ref_name))

        et = time.time()
        elapsed_time = et - st
        st = et
        print("\t- Done in:", elapsed_time, "seconds")

    matrices["train_matrices"] = train_contact_matrices

    comparison_set = set()
    for number, molecule in enumerate(ensemble["topology"].molecules, 1):
        comparison_dataframe = train_topology_dataframe.loc[train_topology_dataframe["molecule"] == f"{number}_{molecule}"]
        if not comparison_dataframe.empty:
            comparison_set = set(comparison_dataframe[~comparison_dataframe["name"].str.startswith("H")]["name"].to_list())
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
    for molecule, bond_pair in meGO_ensemble["bond_pairs"].items():
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
        exclusion_bonds, tmp_p14 = topology.get_14_interaction_list(reduced_topology, bond_pair)
        # split->convert->remerge:
        tmp_ex = pd.DataFrame(columns=["ai", "aj", "exclusion_bonds"])
        tmp_ex["exclusion_bonds"] = pd.Series(exclusion_bonds).astype("category")
        tmp_ex[["ai", "aj"]] = tmp_ex["exclusion_bonds"].str.split("_", expand=True)
        tmp_ex["ai"] = tmp_ex["ai"].map(type_atnum_dict).astype("category")
        tmp_ex["aj"] = tmp_ex["aj"].map(type_atnum_dict).astype("category")
        tmp_ex["same_chain"] = True
        tmp_ex["1-4"] = "1_2_3"
        tmp_ex.loc[(tmp_ex["exclusion_bonds"].isin(tmp_p14)), "1-4"] = "1_4"
        tmp_ex["1-4"] = tmp_ex["1-4"].astype("category")
        exclusion_bonds14 = pd.concat([exclusion_bonds14, tmp_ex], axis=0, sort=False, ignore_index=True)

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
            pairs["1-4"] = "1_4"
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

        if not pairs.empty:
            pairs14 = pd.concat([pairs14, pairs], axis=0, sort=False, ignore_index=True)

    return pairs14, exclusion_bonds14


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
        "zf",
        "epsilon_0",
        "md_threshold",
        "rc_threshold",
        "limit_rc",
        "rc_distance",
        "rc_probability",
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

    train_dataset = pd.merge(
        pd.merge(
            train_dataset,
            pairs14[["ai", "aj", "same_chain", "rep"]],
            how="left",
            on=["ai", "aj", "same_chain"],
        ),
        exclusion_bonds14[["ai", "aj", "same_chain", "1-4"]],
        how="left",
        on=["ai", "aj", "same_chain"],
    )

    train_dataset["ai"] = train_dataset["ai"].astype("category")
    train_dataset["aj"] = train_dataset["aj"].astype("category")

    # We remove from train the 0_1_2_3 intramolecolar interactions
    train_dataset = train_dataset[
        ~(((train_dataset["ai"] == train_dataset["aj"]) & train_dataset["same_chain"]) | (train_dataset["1-4"] == "1_2_3"))
    ]
    train_dataset.reset_index(inplace=True)

    train_dataset["1-4"] = train_dataset["1-4"].cat.add_categories(["1>4"])
    train_dataset["1-4"] = train_dataset["1-4"].fillna("1>4").astype("category")
    train_dataset.loc[(train_dataset["1-4"] == "1_4") & (train_dataset["rep"].isnull()), "rep"] = 0.0

    type_to_c12 = {key: val for key, val in zip(type_definitions.gromos_atp.name, type_definitions.gromos_atp.rc_c12)}
    if args.custom_c12 is not None:
        custom_c12_dict = io.read_custom_c12_parameters(args.custom_c12)
        type_to_c12_appo = {key: val for key, val in zip(custom_c12_dict.name, custom_c12_dict.c12)}
        type_to_c12.update(type_to_c12_appo)

    type_ai_mapped = train_dataset["ai"].map(meGO_ensemble["sbtype_type_dict"])
    type_aj_mapped = train_dataset["aj"].map(meGO_ensemble["sbtype_type_dict"])

    oxygen_mask = masking.create_linearized_mask(
        type_ai_mapped.to_numpy(),
        type_aj_mapped.to_numpy(),
        [("O", "O"), ("OM", "OM"), ("O", "OM")],
        symmetrize=True,
    )

    pairwise_c12 = np.where(
        oxygen_mask,
        11.4 * np.sqrt(type_ai_mapped.map(type_to_c12) * type_aj_mapped.map(type_to_c12)),
        np.sqrt(
            train_dataset["ai"].map(meGO_ensemble["sbtype_c12_dict"])
            * train_dataset["aj"].map(meGO_ensemble["sbtype_c12_dict"])
        ),
    )
    train_dataset["rep"] = train_dataset["rep"].fillna(pd.Series(pairwise_c12))

    return train_dataset


def get_residue_number(s):
    return int(s.split("_")[-1])


def generate_rc_LJ(meGO_ensemble):
    """
    The multi-eGO random coil force-field includes special repulsive interaction only for oxygen-oxygen pairs
    these are generate in the following
    """
    O_OM_sbtype = [
        sbtype for sbtype, atomtype in meGO_ensemble["sbtype_type_dict"].items() if atomtype == "O" or atomtype == "OM"
    ]

    # Generate all possible combinations
    combinations = list(itertools.combinations_with_replacement(O_OM_sbtype, 2))

    # Create a DataFrame from the combinations
    rc_LJ = pd.DataFrame(combinations, columns=["ai", "aj"])
    rc_LJ["type"] = 1
    rc_LJ["c6"] = 0.0
    rc_LJ["c12"] = 11.4 * np.sqrt(
        rc_LJ["ai"].map(meGO_ensemble["sbtype_c12_dict"]) * rc_LJ["aj"].map(meGO_ensemble["sbtype_c12_dict"])
    )

    return rc_LJ


def generate_mg_LJ(meGO_ensemble):
    """
    The multi-eGO random coil force-field includes special repulsive interaction for oxygen-oxygen pairs and for the
    ch1a with all other atoms, these are generate in the following
    """
    O_OM_sbtype = [
        sbtype for sbtype, atomtype in meGO_ensemble["sbtype_type_dict"].items() if atomtype == "O" or atomtype == "OM"
    ]

    # Generate all possible combinations
    combinations = list(itertools.combinations_with_replacement(O_OM_sbtype, 2))

    # Create a DataFrame from the combinations
    rc_LJ = pd.DataFrame(combinations, columns=["ai", "aj"])
    rc_LJ["type"] = 1
    rc_LJ["c6"] = 0.0
    rc_LJ["c12"] = 11.4 * np.sqrt(
        rc_LJ["ai"].map(meGO_ensemble["sbtype_c12_dict"]) * rc_LJ["aj"].map(meGO_ensemble["sbtype_c12_dict"])
    )

    CH1a_sbtype = [sbtype for sbtype, atomtype in meGO_ensemble["sbtype_type_dict"].items() if atomtype == "CH1a"]

    all_sbtypes = list(meGO_ensemble["sbtype_type_dict"].keys())
    # Step 3: Create a list for 'aj' by excluding 'CH1a_sbtype'
    aj_sbtype = [sbtype for sbtype in all_sbtypes if sbtype not in CH1a_sbtype]
    # Step 4: Build the DataFrame with all combinations of 'ai' and 'aj'
    ai_aj_combinations = [(ai, aj) for ai in CH1a_sbtype for aj in aj_sbtype]
    # Create the DataFrame
    df = pd.DataFrame(ai_aj_combinations, columns=["ai", "aj"])

    df["type"] = 1
    df["c6"] = 0.0
    df["c12"] = np.sqrt(df["ai"].map(meGO_ensemble["sbtype_c12_dict"]) * df["aj"].map(meGO_ensemble["sbtype_c12_dict"]))

    mg_LJ = pd.concat([rc_LJ, df])

    return mg_LJ


def generate_basic_LJ(meGO_ensemble, args, matrices=None):
    """
    Generates basic LJ (Lennard-Jones) interactions DataFrame within a molecular ensemble.

    Args:
    - meGO_ensemble (dict): A dictionary containing information about the molecular ensemble.

    Returns:
    - basic_LJ (DataFrame): DataFrame containing basic LJ interactions.

    This function generates a DataFrame 'basic_LJ' containing basic LJ interactions within a molecular ensemble.
    It calculates LJ interactions based on atom types, molecules, and reference matrices present in the ensemble.

    Note:
    - The 'meGO_ensemble' dictionary is expected to contain necessary details regarding the molecular ensemble.
    - The returned DataFrame 'basic_LJ' includes columns defining LJ interaction properties such as atom indices,
      types, c6, c12, sigma, epsilon, probability, rc_probability, molecule names, source, and thresholds.
    - The generated DataFrame provides basic LJ interactions for further analysis or processing within the ensemble.
    """
    columns = [
        "ai",
        "aj",
        "type",
        "c6",
        "c12",
        "sigma",
        "epsilon",
        "probability",
        "rc_probability",
        "distance",
        "molecule_name_ai",
        "molecule_name_aj",
        "same_chain",
        "source",
        "md_threshold",
        "rc_threshold",
        "number_ai",
        "number_aj",
        "cutoff",
    ]

    basic_LJ = pd.DataFrame()
    topol_df = meGO_ensemble["topology_dataframe"]

    name_to_c12 = {key: val for key, val in zip(type_definitions.gromos_atp.name, type_definitions.gromos_atp.rc_c12)}
    if args.custom_c12 is not None:
        custom_c12_dict = io.read_custom_c12_parameters(args.custom_c12)
        name_to_c12_appo = {key: val for key, val in zip(custom_c12_dict.name, custom_c12_dict.c12)}
        name_to_c12.update(name_to_c12_appo)

    for name in matrices["reference_matrices"].keys():
        temp_basic_LJ = pd.DataFrame(columns=columns)
        mol_num_i = str(name.split("_")[-2])
        mol_num_j = str(name.split("_")[-1])
        ensemble = matrices["reference_matrices"][name]
        temp_basic_LJ["ai"] = ensemble["rc_ai"]
        temp_basic_LJ["aj"] = ensemble["rc_aj"]
        temp_basic_LJ["type"] = 1
        temp_basic_LJ["c6"] = 0.0
        temp_basic_LJ["c12"] = 0.0
        temp_basic_LJ["same_chain"] = ensemble["rc_same_chain"]
        temp_basic_LJ["molecule_name_ai"] = ensemble["rc_molecule_name_ai"]
        temp_basic_LJ["molecule_name_aj"] = ensemble["rc_molecule_name_aj"]
        temp_basic_LJ["source"] = "basic"

        atom_set_i = topol_df[topol_df["molecule_number"] == mol_num_i]["type"]
        atom_set_j = topol_df[topol_df["molecule_number"] == mol_num_j]["type"]
        c12_list_i = atom_set_i.map(name_to_c12).to_numpy(dtype=np.float64)
        c12_list_j = atom_set_j.map(name_to_c12).to_numpy(dtype=np.float64)
        ai_name = atom_set_i.to_numpy(dtype=str)
        aj_name = atom_set_j.to_numpy(dtype=str)
        oxygen_mask = masking.create_array_mask(ai_name, aj_name, [("O", "OM"), ("O", "O"), ("OM", "OM")], symmetrize=True)
        temp_basic_LJ["c12"] = 11.4 * np.sqrt(c12_list_i * c12_list_j[:, np.newaxis]).flatten()
        temp_basic_LJ["rep"] = temp_basic_LJ["c12"]
        temp_basic_LJ = temp_basic_LJ[oxygen_mask]
        temp_basic_LJ = temp_basic_LJ.dropna(axis=1, how="all")
        temp_basic_LJ = temp_basic_LJ.drop_duplicates(subset=["ai", "aj", "same_chain"], keep="first")

        basic_LJ = pd.concat([basic_LJ, temp_basic_LJ])

    basic_LJ["probability"] = 1.0
    basic_LJ["rc_probability"] = 1.0
    basic_LJ["rc_threshold"] = 1.0
    basic_LJ["md_threshold"] = 1.0
    basic_LJ["epsilon"] = -basic_LJ["c12"]
    basic_LJ["cutoff"] = 1.45 * basic_LJ["c12"] ** (1.0 / 12.0)
    basic_LJ["sigma"] = basic_LJ["cutoff"] / (2.0 ** (1.0 / 6.0))
    basic_LJ["distance"] = basic_LJ["cutoff"]
    basic_LJ["learned"] = 0
    basic_LJ["1-4"] = "1>4"
    # Sorting the pairs prioritising intermolecular interactions
    basic_LJ.sort_values(by=["ai", "aj", "same_chain"], ascending=[True, True, True], inplace=True)
    # Cleaning the duplicates
    basic_LJ = basic_LJ.drop_duplicates(subset=["ai", "aj"], keep="first")

    return basic_LJ


def set_sig_epsilon(meGO_LJ, needed_fields):
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
    # when distance estimates are poor we use the cutoff value
    meGO_LJ.loc[(meGO_LJ["probability"] <= meGO_LJ["md_threshold"]), "distance"] = (meGO_LJ["rep"] / meGO_LJ["epsilon_0"]) ** (
        1.0 / 12.0
    )
    meGO_LJ.loc[(meGO_LJ["rc_probability"] <= meGO_LJ["md_threshold"]), "rc_distance"] = (
        meGO_LJ["rep"] / meGO_LJ["epsilon_0"]
    ) ** (1.0 / 12.0)

    # Epsilon is initialised to a rescaled C12
    # This is always correct becasue distance is always well defined by either training data
    # or using default C12 values
    # negative epsilon are used to identify non-attractive interactions
    meGO_LJ["epsilon"] = -meGO_LJ["rep"] * (meGO_LJ["distance"] / meGO_LJ["rc_distance"]) ** 12

    # Attractive interactions
    # These are defined only if the training probability is greater than MD_threshold and
    # by comparing them with RC_probabilities
    meGO_LJ.loc[
        (meGO_LJ["probability"] > meGO_LJ["limit_rc"] * np.maximum(meGO_LJ["rc_probability"], meGO_LJ["rc_threshold"]))
        & (meGO_LJ["probability"] > meGO_LJ["md_threshold"]),
        "epsilon",
    ] = -(meGO_LJ["epsilon_0"] / np.log(meGO_LJ["zf"] * meGO_LJ["rc_threshold"])) * (
        np.log(meGO_LJ["probability"] / (meGO_LJ["zf"] * np.maximum(meGO_LJ["rc_probability"], meGO_LJ["rc_threshold"])))
    )

    # General repulsive term
    # this is used only when MD_p < RC_p eventually corrected by the ZF
    # negative epsilon are used to identify non-attractive interactions
    meGO_LJ.loc[
        (
            np.maximum(meGO_LJ["probability"], meGO_LJ["rc_threshold"])
            < meGO_LJ["zf"] * np.maximum(meGO_LJ["rc_probability"], meGO_LJ["rc_threshold"])
        ),
        "epsilon",
    ] = -(meGO_LJ["epsilon_0"] / (np.log(meGO_LJ["zf"] * meGO_LJ["rc_threshold"]))) * meGO_LJ["distance"] ** 12 * (
        np.log(
            np.maximum(meGO_LJ["probability"], meGO_LJ["rc_threshold"])
            / (meGO_LJ["zf"] * np.maximum(meGO_LJ["rc_probability"], meGO_LJ["rc_threshold"]))
        )
    ) - (
        meGO_LJ["rep"] * (meGO_LJ["distance"] / meGO_LJ["rc_distance"]) ** 12
    )

    # update the c12 1-4 interactions
    # in principle this is not needed but we redefine them to be safe
    meGO_LJ.loc[(meGO_LJ["1-4"] == "1_4"), "epsilon"] = -meGO_LJ["rep"] * (meGO_LJ["distance"] / meGO_LJ["rc_distance"]) ** 12

    # clean NaN and zeros
    meGO_LJ.dropna(subset=["epsilon"], inplace=True)
    meGO_LJ = meGO_LJ[meGO_LJ.epsilon != 0]

    # lower value for repulsion
    meGO_LJ.loc[
        (meGO_LJ["1-4"] != "1_4") & (meGO_LJ["epsilon"] < 0.0) & (-meGO_LJ["epsilon"] < 0.1 * meGO_LJ["rep"]),
        "epsilon",
    ] = (
        -0.1 * meGO_LJ["rep"]
    )
    # higher value for repulsion
    meGO_LJ.loc[
        (meGO_LJ["1-4"] != "1_4") & (meGO_LJ["epsilon"] < 0.0) & (-meGO_LJ["epsilon"] > 20.0 * meGO_LJ["rep"]),
        "epsilon",
    ] = (
        -20.0 * meGO_LJ["rep"]
    )

    # but within a lower
    meGO_LJ.loc[
        (meGO_LJ["1-4"] == "1_4") & (-meGO_LJ["epsilon"] < 0.1 * meGO_LJ["rep"]),
        "epsilon",
    ] = (
        -0.1 * meGO_LJ["rep"]
    )
    # and an upper value
    meGO_LJ.loc[
        (meGO_LJ["1-4"] == "1_4") & (-meGO_LJ["epsilon"] > meGO_LJ["rep"]),
        "epsilon",
    ] = -meGO_LJ["rep"]

    # Sigma is set from the estimated interaction length
    meGO_LJ = meGO_LJ.assign(sigma=meGO_LJ["distance"] / (2.0 ** (1.0 / 6.0)))

    # for repulsive interaction we reset sigma to its effective value
    # this because when merging repulsive contacts from different sources what will matters
    # will be the repulsive strength
    meGO_LJ.loc[(meGO_LJ["epsilon"] < 0.0), "sigma"] = (-meGO_LJ["epsilon"]) ** (1.0 / 12.0)

    # add a flag to identify learned contacts
    meGO_LJ.loc[:, "learned"] = 1

    # keep only needed fields
    meGO_LJ = meGO_LJ[needed_fields]

    return meGO_LJ


def consistency_checks(meGO_LJ):
    """
    Perform consistency checks on LJ parameters.

    This function performs consistency checks on LJ (Lennard-Jones) parameters to avoid data inconsistencies.

    Parameters
    ----------
    meGO_LJ : pd.DataFrame
        DataFrame containing LJ parameters such as repulsive term (rep) and cutoff distance (cutoff).

    Raises
    ------
    RuntimeError
        If inconsistencies are found between the calculated cutoff distance and the provided cutoff values.

    Notes
    -----
    This function is primarily used for debugging purposes to ensure data integrity and consistency.
    """
    # This is a debug check to avoid data inconsistencies
    if (np.abs(1.45 * meGO_LJ["rep"] ** (1 / 12) - meGO_LJ["cutoff"])).max() > 10e-6:
        print(
            meGO_LJ[["ai", "aj", "same_chain", "source", "rep", "cutoff"]]
            .loc[(np.abs(1.45 * meGO_LJ["rep"] ** (1 / 12) - meGO_LJ["cutoff"]) > 10e-6)]
            .to_string()
        )
        exit("HERE SOMETHING BAD HAPPEND: There are inconsistent cutoff/c12 values")


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


def generate_LJ(meGO_ensemble, train_dataset, basic_LJ, parameters):
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
    meGO_LJ = train_dataset.copy()
    # meGO needed fields
    needed_fields = [
        "molecule_name_ai",
        "ai",
        "molecule_name_aj",
        "aj",
        "probability",
        "same_chain",
        "source",
        "rc_probability",
        "sigma",
        "epsilon",
        "1-4",
        "distance",
        "cutoff",
        "rep",
        "md_threshold",
        "rc_threshold",
        "learned",
    ]

    # generate attractive and repulsive interactions
    meGO_LJ = set_sig_epsilon(meGO_LJ, needed_fields)

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
    # meGO consistency checks
    consistency_checks(meGO_LJ)

    # Merging of multiple simulations:
    # Here we sort all the atom pairs based on the distance and the probability.
    # among attractive we keep the shortest the same among repulsive.
    meGO_LJ["type"] = np.sign(meGO_LJ["epsilon"])
    meGO_LJ.sort_values(
        by=["ai", "aj", "same_chain", "type", "sigma", "epsilon"],
        ascending=[True, True, True, False, True, False],
        inplace=True,
    )
    # Cleaning the duplicates
    meGO_LJ = meGO_LJ.drop_duplicates(subset=["ai", "aj", "same_chain"], keep="first")

    # now we can remove all repulsive contacts with default (i.e., rep) c12 becasue these
    # are uninformative and predefined. This also allow to replace them with contact learned
    # by either intra/inter training. We cannot remove 1-4 interactions.
    meGO_LJ = meGO_LJ.loc[
        ~(
            (meGO_LJ["epsilon"] < 0)
            & ((abs(-meGO_LJ["epsilon"] - meGO_LJ["rep"]) / meGO_LJ["rep"]) < parameters.relative_c12d)
            & (meGO_LJ["1-4"] == "1>4")
        )
    ]

    # transfer rule for inter/intra contacts:
    # 1) only attractive contacts can be transferd
    # 2) attractive contacts that can be transferd are those non affected by their random coil (prc <= rc_threshold)
    # 3) an attractive contacts can only take the place of a trivial repulsive contact (i.e. a repulsive contact with prc <= rc_threshold)
    meGO_LJ["trivial"] = False
    meGO_LJ.loc[(meGO_LJ["epsilon"] > 0) & (meGO_LJ["rc_probability"] > meGO_LJ["rc_threshold"]), "trivial"] = True
    meGO_LJ.loc[(meGO_LJ["epsilon"] < 0) & (meGO_LJ["rc_probability"] <= meGO_LJ["rc_threshold"]), "trivial"] = True
    # Identify rows where "trivial repulsive" is True and there exists another duplicate row
    duplicated_at_least_one_trivial = (
        (meGO_LJ.duplicated(subset=["ai", "aj", "1-4"], keep=False)) & (meGO_LJ["trivial"]) & (meGO_LJ["type"] == -1)
    )
    # Identify rows where both are trivial/not trivial
    duplicated_same_trivial = meGO_LJ.duplicated(subset=["ai", "aj", "trivial"], keep=False)
    # Identify rows where both are attractive or repulsive
    duplicated_same_type = meGO_LJ.duplicated(subset=["ai", "aj", "type"], keep=False)
    # Combine the conditions to remove only the rows that are trivial but duplicated with a non-trivial counterpart
    remove_duplicates_mask = duplicated_at_least_one_trivial & ~duplicated_same_trivial & ~duplicated_same_type
    # Remove rows where "trivial" is True and there exists another duplicate row with "trivial" as False with the not Trivial attractive and the Trivial repulsive
    meGO_LJ = meGO_LJ[~remove_duplicates_mask]
    meGO_LJ = meGO_LJ[needed_fields]

    # now is a good time to acquire statistics on the parameters
    # this should be done per interaction pair (cycling over all molecules combinations) and inter/intra/intra_d
    io.print_stats(meGO_LJ)

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

    # copy 1-4 interactions into meGO_LJ_14
    copy14 = meGO_LJ.loc[(meGO_LJ["1-4"] == "1_4")]
    meGO_LJ_14 = pd.concat([meGO_LJ_14, copy14], axis=0, sort=False, ignore_index=True)
    # remove them from the default force-field
    meGO_LJ = meGO_LJ.loc[(meGO_LJ["1-4"] != "1_4")]

    if not parameters.single_molecule:
        # if an intramolecular interactions is associated with a large rc_probability then it is moved to meGO_LJ_14 to
        # avoid its use as intermolecular, this includes all repulsive
        copy_intra = meGO_LJ.loc[(meGO_LJ["same_chain"]) & (meGO_LJ["rc_probability"] > meGO_LJ["rc_threshold"])]
        meGO_LJ_14 = pd.concat([meGO_LJ_14, copy_intra], axis=0, sort=False, ignore_index=True)
        # remove them from the default force-field
        meGO_LJ = meGO_LJ.loc[~((meGO_LJ["same_chain"]) & (meGO_LJ["rc_probability"] > meGO_LJ["rc_threshold"]))]

    # now we can decide to keep intermolecular interactions as intramolecular ones
    # to do this is enough to remove it from meGO_LJ_14, in this way the value used for the contact is the one meGO_LJ
    if not parameters.force_split:
        # Reset index of meGO_LJ_14
        meGO_LJ_14_reset_index = meGO_LJ_14.reset_index()
        # Filter rows in meGO_LJ_14 that meet the condition
        filtered_meGO_LJ_14 = meGO_LJ_14_reset_index.loc[
            (~meGO_LJ_14_reset_index["same_chain"])
            & (meGO_LJ_14_reset_index["molecule_name_ai"] == meGO_LJ_14_reset_index["molecule_name_aj"])
            & (meGO_LJ_14_reset_index["epsilon"] > 0.0)
            & (meGO_LJ_14_reset_index["rc_probability"] <= meGO_LJ_14_reset_index["rc_threshold"])
        ]
        filtered_train_dataset = train_dataset.loc[train_dataset["same_chain"]].copy()
        filtered_train_dataset.sort_values(by=["ai", "aj", "rc_threshold"], ascending=[True, True, True], inplace=True)
        filtered_train_dataset = filtered_train_dataset.drop_duplicates(subset=["ai", "aj"], keep="first")[
            ["ai", "aj", "rc_probability", "rc_threshold"]
        ]
        # Merge filtered_meGO_LJ_14 with train_dataset based on ai and aj columns
        merged = pd.merge(filtered_meGO_LJ_14, filtered_train_dataset, on=["ai", "aj"], how="inner")
        # Filter rows where rc_probability is less than or equal to rc_threshold
        to_remove_indices = merged[merged["rc_probability_y"] <= merged["rc_threshold_y"]]["index"]
        meGO_LJ_14.drop(to_remove_indices, inplace=True)
    else:
        split_ii = meGO_LJ.loc[(meGO_LJ["same_chain"])]
        # move the intramolecular interaction in the topology
        meGO_LJ_14 = pd.concat([meGO_LJ_14, split_ii], axis=0, sort=False, ignore_index=True)
        # remove them from the default force-field
        meGO_LJ = meGO_LJ.loc[(~meGO_LJ["same_chain"])]

    # Now is time to add masked default interactions for pairs
    # that have not been learned in any other way
    basic_LJ = basic_LJ[needed_fields]
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
    meGO_LJ.sort_values(by=["ai", "aj", "same_chain", "learned"], ascending=[True, True, True, False], inplace=True)
    # Cleaning the duplicates, that is that we retained a not learned interaction only if it is unique
    # first we remove duplicated masked interactions
    meGO_LJ = meGO_LJ.drop_duplicates(subset=["ai", "aj", "same_chain", "learned"], keep="first")
    meGO_LJ = meGO_LJ.loc[(~(meGO_LJ.duplicated(subset=["ai", "aj"], keep=False)) | (meGO_LJ["learned"] == 1))]

    # we are ready to finalize the setup
    # Calculate c6 and c12 for meGO_LJ
    meGO_LJ["c6"] = np.where(meGO_LJ["epsilon"] < 0.0, 0.0, 4 * meGO_LJ["epsilon"] * (meGO_LJ["sigma"] ** 6))

    meGO_LJ["c12"] = np.where(
        meGO_LJ["epsilon"] < 0.0, -meGO_LJ["epsilon"], abs(4 * meGO_LJ["epsilon"] * (meGO_LJ["sigma"] ** 12))
    )

    # Calculate c6 and c12 for meGO_LJ_14
    meGO_LJ_14["c6"] = np.where(meGO_LJ_14["epsilon"] < 0.0, 0.0, 4 * meGO_LJ_14["epsilon"] * (meGO_LJ_14["sigma"] ** 6))

    meGO_LJ_14["c12"] = np.where(
        meGO_LJ_14["epsilon"] < 0.0, -meGO_LJ_14["epsilon"], abs(4 * meGO_LJ_14["epsilon"] * (meGO_LJ_14["sigma"] ** 12))
    )

    meGO_LJ["type"] = 1
    meGO_LJ["number_ai"] = meGO_LJ["ai"].map(meGO_ensemble["sbtype_number_dict"]).astype(int)
    meGO_LJ["number_aj"] = meGO_LJ["aj"].map(meGO_ensemble["sbtype_number_dict"]).astype(int)

    # Here we want to sort so that ai is smaller than aj
    meGO_LJ = meGO_LJ[(meGO_LJ["ai"].cat.codes <= meGO_LJ["aj"].cat.codes)]
    (
        meGO_LJ["ai"],
        meGO_LJ["aj"],
        meGO_LJ["molecule_name_ai"],
        meGO_LJ["molecule_name_aj"],
        meGO_LJ["number_ai"],
        meGO_LJ["number_aj"],
    ) = np.where(
        meGO_LJ["number_ai"] <= meGO_LJ["number_aj"],
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

    meGO_LJ.sort_values(by=["number_ai", "number_aj"], inplace=True)

    # meGO consistency checks
    consistency_checks(meGO_LJ)
    consistency_checks(meGO_LJ_14)

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
        "rep",
        "cutoff",
        "same_chain",
        "source",
        "number_ai",
        "number_aj",
    ]

    meGO_LJ = meGO_LJ[final_fields]

    et = time.time()
    elapsed_time = et - st
    print("\t- Done in:", elapsed_time, "seconds")

    return meGO_LJ, meGO_LJ_14


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
    pairs_molecule_dict = {}
    for molecule, bond_pair in meGO_ensemble["bond_pairs"].items():
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

        atnum_type_dict = reduced_topology.set_index("sb_type")["number"].to_dict()
        resnum_type_dict = reduced_topology.set_index("sb_type")["resnum"].to_dict()

        # Building the exclusion bonded list
        # exclusion_bonds are all the interactions within 3 bonds
        # p14 are specifically the interactions at exactly 3 bonds
        exclusion_bonds, p14 = topology.get_14_interaction_list(reduced_topology, bond_pair)

        # in the case of the MG prior we need to remove interactions in a window of 2 residues
        if args.egos == "mg":
            # Create a list of tuples (sbtype, residue_number)
            sbtype_with_residue = [
                (sbtype, resnum_type_dict[sbtype])
                for sbtype in reduced_topology["sb_type"]
                if meGO_ensemble["sbtype_type_dict"][sbtype] != "CH1a"
            ]
            # Sort the list by residue numbers
            sbtype_with_residue.sort(key=lambda x: x[1])
            # Initialize a list to hold the filtered combinations
            filtered_combinations = []
            # Use two pointers to find valid pairs
            n = len(sbtype_with_residue)
            for i in range(n):
                j = i + 1  # Start with the current sbtype
                # Find the range of valid sbtypes
                while j < n and abs(sbtype_with_residue[j][1] - sbtype_with_residue[i][1]) <= 2:
                    filtered_combinations.append((sbtype_with_residue[i][0], sbtype_with_residue[j][0]))
                    j += 1

            # Create a DataFrame from the filtered combinations
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
            ] *= 11.4
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
            meGO_LJ_14 = pd.concat([meGO_LJ_14, df], axis=0, sort=False, ignore_index=True)

        pairs = pd.DataFrame()
        if not meGO_LJ_14.empty:
            # pairs do not have duplicates because these have been cleaned before
            pairs = meGO_LJ_14[
                [
                    "ai",
                    "aj",
                    "c6",
                    "c12",
                    "same_chain",
                    "probability",
                    "rc_probability",
                    "source",
                    "rep",
                ]
            ].copy()

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
            pairs["func"] = 1
            # Intermolecular interactions are excluded
            pairs.loc[(~pairs["same_chain"]), "c6"] = 0.0
            pairs.loc[(~pairs["same_chain"]), "c12"] = pairs["rep"]
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
