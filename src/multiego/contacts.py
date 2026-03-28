"""
contacts.py — contact matrix processing for multi-eGO.

This module covers two related stages of the pipeline:

1. **Training topology indexing** (``_index_training_topology``):
   Lightweight helper that builds only the per-atom DataFrame and the
   atom-index → sb_type mapping needed when loading training contact matrices.
   This is intentionally narrower than the full initialisation in
   ``ensemble_data._initialize_topology`` — LJ parameter lookup and all
   ``sbtype_*`` dict extractions are omitted because they are not needed for
   training topologies.

2. **Contact matrix loading** (``init_meGO_matrices`` and its private helpers):
   Reads reference and training contact matrices, annotates them with prior LJ
   parameters from the force-field topology, and computes the adaptive
   probability thresholds used by ``lj.set_sig_epsilon``.

Private helpers are prefixed with ``_`` and are not part of the public API.
"""

from . import io
from . import type_definitions

import numpy as np
import os
import pandas as pd
import parmed
import sys
import time
import warnings

# ---------------------------------------------------------------------------
# Training topology indexing (private — used only by _load_train_matrix)
# ---------------------------------------------------------------------------


def _index_training_topology(topology, custom_dict):
    """
    Build a per-atom DataFrame and the atom-index → sb_type mapping for a
    training topology.

    This is a lightweight subset of the full topology initialisation performed
    by ``MeGOEnsemble.from_topology``. It only constructs the two data
    structures that ``_load_train_matrix`` actually needs:

    - the topology DataFrame (used later in ``init_meGO_matrices`` to check
      that all non-hydrogen atom names have been mapped to multi-eGO sb_types)
    - the ``molecules_idx_sbtype_dictionary`` (used by
      ``io.read_molecular_contacts`` to remap raw atom indices to sb_type
      labels)

    LJ parameter lookup and all ``sbtype_*`` dict extractions are intentionally
    omitted — those already live on the ``MeGOEnsemble`` built from the base
    topology and are not needed here.

    Parameters
    ----------
    topology : parmed.Structure
        Loaded GROMACS training topology.
    custom_dict : dict
        Additional atom-name remappings for non-standard residues, merged on
        top of ``type_definitions.from_ff_to_multiego``.

    Returns
    -------
    topology_dataframe : pd.DataFrame
        Per-atom DataFrame with columns ``number``, ``molecule``,
        ``molecule_number``, ``molecule_name``, ``resnum``, ``name``,
        and ``sb_type``.
    molecules_idx_sbtype_dictionary : dict
        ``{molecule_key: {atom_number_str: sb_type}}`` mapping atom indices
        to sb_types per molecule, where ``molecule_key`` is
        ``"{number}_{molecule_name}"``.
    """
    frames = []
    new_number, col_molecule, new_resnum = [], [], []
    molecules_idx_sbtype_dictionary = {}

    for molecule_number, (molecule_name, molecule_topology) in enumerate(topology.molecules.items(), 1):
        molecules_idx_sbtype_dictionary[f"{molecule_number}_{molecule_name}"] = {}
        frames.append(molecule_topology[0].to_dataframe())
        for atom in molecule_topology[0].atoms:
            new_number.append(str(atom.idx + 1))
            col_molecule.append(f"{molecule_number}_{molecule_name}")
            new_resnum.append(str(atom.residue.number))

    topology_dataframe = pd.concat(frames, axis=0, ignore_index=True)
    topology_dataframe["number"] = new_number
    topology_dataframe["molecule"] = col_molecule
    topology_dataframe["molecule_number"] = col_molecule
    topology_dataframe[["molecule_number", "molecule_name"]] = topology_dataframe.molecule.str.split(
        "_", expand=True, n=1
    )
    topology_dataframe["resnum"] = new_resnum

    from_ff_to_multiego_extended = type_definitions.from_ff_to_multiego.copy()
    from_ff_to_multiego_extended.update(custom_dict)
    topology_dataframe = topology_dataframe.replace({"name": from_ff_to_multiego_extended})

    topology_dataframe["sb_type"] = (
        topology_dataframe["name"]
        + "_"
        + topology_dataframe["molecule_name"]
        + "_"
        + topology_dataframe["resnum"].astype(str)
    )

    for molecule in molecules_idx_sbtype_dictionary:
        tmp = topology_dataframe.loc[topology_dataframe["molecule"] == molecule]
        molecules_idx_sbtype_dictionary[molecule] = tmp[["number", "sb_type"]].set_index("number")["sb_type"].to_dict()

    return topology_dataframe, molecules_idx_sbtype_dictionary


# ---------------------------------------------------------------------------
# LJ parameter extraction (private — used only by _load_reference_matrix)
# ---------------------------------------------------------------------------


def _get_lj_params(topology):
    """
    Extract per-atom LJ parameters from a parmed topology.

    Parameters are read from each atom's ``sigma`` (Angstrom to nm via x0.1)
    and ``epsilon`` (kcal/mol to kJ/mol via x4.184) attributes and returned as
    geometric-mean-ready c6/c12 values.

    Parameters
    ----------
    topology : parmed.Structure
        Loaded GROMACS topology.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``ai`` (atom type string), ``c6``, and ``c12``,
        indexed 0...N-1 for each atom in the topology.
    """
    lj_params = pd.DataFrame(
        {
            "ai": [atom.atom_type for atom in topology.atoms],
            "c6": [atom.sigma * 0.1 for atom in topology.atoms],
            "c12": [atom.epsilon * 4.184 for atom in topology.atoms],
        }
    )
    return lj_params


def _get_lj_pairs(topology):
    """
    Extract explicit nonbonded pair parameters (``[nonbond_params]``) from a
    parmed topology as LJ sigma/epsilon values.

    Values are stored internally as rmin (not sigma), so the c6 column is
    scaled by ``1 / 2^(1/6)`` on read. If a pair type appears more than once,
    the last entry is used (consistent with GROMACS behaviour).

    Parameters
    ----------
    topology : parmed.Structure
        Loaded GROMACS topology.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``ai``, ``aj``, ``epsilon``, and ``sigma``.
    """
    nbfix = topology.parameterset.nbfix_types
    lj_pairs = pd.DataFrame(columns=["ai", "aj", "epsilon", "sigma"], index=np.arange(len(nbfix)))
    for i, (sbtype_i, sbtype_j) in enumerate(nbfix):
        c12 = nbfix[(sbtype_i, sbtype_j)][0] * 4.184
        c6 = nbfix[(sbtype_i, sbtype_j)][1] * 0.1 / (2 ** (1 / 6))
        epsilon = c6**2 / (4 * c12) if c6 > 0 else -c12
        sigma = (c12 / c6) ** (1 / 6) if c6 > 0 else c12 ** (1 / 12) / (2.0 ** (1.0 / 6.0))
        lj_pairs.loc[i] = [sbtype_i, sbtype_j, epsilon, sigma]
    return lj_pairs


def _get_lj14_pairs(topology):
    """
    Extract 1-4 pair parameters (``[pairs]``) from all molecules in a parmed
    topology as LJ sigma/epsilon values.

    The parmed ``adjusts`` list stores sigma (Angstrom to nm via x0.1) and
    epsilon (kcal/mol to kJ/mol via x4.184). Epsilon and sigma are derived from
    c6/c12 using standard LJ formulae; when c6 = 0 the interaction is purely
    repulsive and epsilon is set to -c12.

    Parameters
    ----------
    topology : parmed.Structure
        Loaded GROMACS topology containing one or more molecules.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``ai``, ``aj``, ``epsilon``, and ``sigma``,
        concatenated across all molecules.
    """
    lj14_pairs = pd.DataFrame()
    for mol, top in topology.molecules.items():
        pair14 = [
            {
                "ai": pair.atom1.type,
                "aj": pair.atom2.type,
                "c6": pair.type.sigma * 0.1,
                "c12": pair.type.epsilon * 4.184,
            }
            for pair in top[0].adjusts
        ]
        lj14_pairs = pd.concat([lj14_pairs, pd.DataFrame(pair14)])

    lj14_pairs = lj14_pairs.reset_index(drop=True)
    lj14_pairs["epsilon"] = np.where(
        lj14_pairs["c6"] > 0,
        lj14_pairs["c6"] ** 2 / (4 * lj14_pairs["c12"]),
        -lj14_pairs["c12"],
    )
    lj14_pairs["sigma"] = np.where(
        lj14_pairs["c6"] > 0,
        (lj14_pairs["c12"] / lj14_pairs["c6"]) ** (1 / 6),
        lj14_pairs["c12"] ** (1 / 12) / (2.0 ** (1.0 / 6.0)),
    )
    lj14_pairs.drop(columns=["c6", "c12"], inplace=True)
    return lj14_pairs


def _make_symmetric(df):
    """
    Return a symmetrised copy of a pairwise DataFrame by concatenating the
    original with its (ai<->aj)-swapped mirror and deduplicating.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with at least columns ``ai`` and ``aj``.

    Returns
    -------
    pd.DataFrame
        Deduplicated symmetric DataFrame with reset index.
    """
    return (
        pd.concat([df, df.rename(columns={"ai": "aj", "aj": "ai"})])
        .drop_duplicates(subset=["ai", "aj"])
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# Contact matrix processing
# ---------------------------------------------------------------------------


def check_intra_domain_complementarity(matrices):
    """
    Verify that reference matrices covering the same molecule pair (e.g. two
    separate intramat_1_1 entries from different references) have
    non-overlapping ``rc_learned`` flags.

    If any atom-pair position is marked as learned in more than one reference,
    the domain splitting is inconsistent and a ``ValueError`` is raised.

    Parameters
    ----------
    matrices : dict
        Mapping of flat matrix name to contact matrix DataFrame. Each
        DataFrame must contain a boolean ``rc_learned`` column.

    Raises
    ------
    ValueError
        If any position has ``rc_learned == True`` in more than one reference
        matrix for the same molecule pair.
    """
    mats_names = ["_".join(key.split("_")[-3:]) for key in matrices]
    to_check_names = [n for n in set(mats_names) if mats_names.count(n) > 1]
    to_check = [[k for k in matrices if "_".join(k.split("_")[-3:]) == name] for name in to_check_names]
    for group in to_check:
        flags = [matrices[k]["rc_learned"].to_numpy() for k in group]
        if np.any(np.sum(flags, axis=0) > 1):
            raise ValueError(
                f"Learning flag complementarity not satisfied for {group} " "(e.g. intra-inter domain splitting)"
            )


def initialize_molecular_contacts(contact_matrix, prior_matrix, args, reference):
    """
    Annotate a training contact matrix with adaptive probability thresholds
    and prior LJ parameters derived from the corresponding reference matrix.

    The adaptive MD threshold ``md_threshold`` is the smallest probability
    value such that the cumulative sum of sorted (descending) learned
    probabilities covers ``args.p_to_learn`` of the total mass. When no
    contacts are learned (norm = 0), the threshold defaults to 1.

    Two per-contact threshold columns are added:

    - ``rc_threshold``: the RC probability equivalent to ``md_threshold``,
      scaled by the ratio of epsilon ranges.
    - ``limit_rc_att``: the RC probability above which the MD contact is
      considered genuinely attractive. Values below 1 arising from negative
      ``epsilon_prior`` are clamped to 1.

    Parameters
    ----------
    contact_matrix : pd.DataFrame
        Training contact matrix with at least columns ``probability`` and
        ``learned`` (boolean flag derived from the reference ``rc_learned``).
    prior_matrix : pd.DataFrame
        Corresponding reference contact matrix with column ``epsilon_prior``.
    args : argparse.Namespace
        Parsed arguments; ``args.p_to_learn`` and ``args.epsilon_min`` are
        used.
    reference : dict
        Single input_refs entry; ``reference["reference"]`` (name string) and
        ``reference["epsilon"]`` (= epsilon_0, the maximum interaction energy)
        are used.

    Returns
    -------
    pd.DataFrame
        The input ``contact_matrix`` with added columns: ``reference``,
        ``epsilon_0``, ``md_threshold``, ``rc_threshold``, and
        ``limit_rc_att``.
    """
    contact_matrix["reference"] = reference["reference"]

    p_sort = np.sort(contact_matrix["probability"].loc[contact_matrix["learned"]].to_numpy())[::-1]
    norm = np.sum(p_sort)
    if norm == 0:
        md_threshold = 1
    else:
        md_threshold = p_sort[np.min(np.where(np.cumsum(p_sort) / norm > args.p_to_learn)[0])]

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

    # When epsilon_prior is negative the formula can produce limit_rc_att < 1,
    # which would incorrectly suppress learning. Clamp to 1 in those cases.
    contact_matrix.loc[
        (contact_matrix["limit_rc_att"] < 1) & (prior_matrix["epsilon_prior"] < 0),
        "limit_rc_att",
    ] = 1

    return contact_matrix


def _path_to_matrix_name(path, root_dir):
    """
    Convert a full matrix file path to a flat unique dict key by stripping
    the root prefix, replacing path separators with underscores, and removing
    all recognised file extensions (``.ndx``, ``.gz``, ``.h5``).

    Parameters
    ----------
    path : str
        Absolute path to the contact matrix file.
    root_dir : str
        Absolute path to the multi-eGO root directory.

    Returns
    -------
    str
        Flat name such as ``"GB1_reference_intramat_1_1"``.
    """
    name = path.replace(f"{root_dir}/inputs/", "").replace("/", "_")
    for ext in (".ndx", ".gz", ".h5"):
        name = name.replace(ext, "")
    return name


def _load_reference_matrix(reference, ensemble, args):
    """
    Load a single reference topology and contact matrix, then annotate the
    matrix with per-contact prior sigma and epsilon values.

    Prior parameters are derived in two steps:

    1. **Combination-rule defaults** — geometric means of per-atom c6/c12
       values read from the reference topology via ``_get_lj_params``.
    2. **Explicit overrides** — where the topology defines explicit nonbonded
       pairs (``[nonbond_params]``) or 1-4 pairs (``[pairs]``), the
       combination-rule values are replaced by the explicit ones.

    As a side effect, ``ensemble.topology_dataframe["c6"]`` and ``["c12"]``
    are updated in-place with the reference topology per-atom values.

    Parameters
    ----------
    reference : dict
        Single ``input_refs`` entry with keys ``'reference'``, ``'matrix'``,
        ``'epsilon'``, and ``'train'``.
    ensemble : MeGOEnsemble
        Initialised ensemble; ``topology_dataframe`` and
        ``molecules_idx_sbtype_dictionary`` are accessed.
    args : argparse.Namespace
        Parsed arguments providing ``root_dir`` and ``system``.

    Returns
    -------
    name : str
        Unique flat key for this reference matrix.
    contact_matrix : pd.DataFrame
        Enriched reference contact matrix with ``rc_`` prefix on all original
        columns plus new columns ``sigma_prior`` and ``epsilon_prior``.

    Raises
    ------
    RuntimeError
        If more than one ``.top`` file is found in the reference directory.
    FileNotFoundError
        If the topology or contact matrix file cannot be located.
    """
    reference_path = f"{args.root_dir}/inputs/{args.system}/{reference['reference']}"
    topol_files = [f for f in os.listdir(reference_path) if ".top" in f]
    if not topol_files:
        raise FileNotFoundError(f"No topology file (.top) found in {reference_path}")
    if len(topol_files) > 1:
        raise RuntimeError(f"More than 1 topology file found in {reference_path}. Only one should be used")

    topology_path = f"{reference_path}/{topol_files[0]}"
    if not os.path.isfile(topology_path):
        raise FileNotFoundError(f"{topology_path} not found.")

    print("\t\t-", f"Reading {topology_path}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        topol = parmed.load_file(topology_path)

    lj_data = _get_lj_params(topol)
    lj_data_dict = {str(key): val for key, val in zip(lj_data["ai"], lj_data[["c6", "c12"]].values)}
    ensemble.topology_dataframe["c6"] = lj_data["c6"].to_numpy()
    ensemble.topology_dataframe["c12"] = lj_data["c12"].to_numpy()

    symmetric_lj_pairs = _make_symmetric(_get_lj_pairs(topol))
    symmetric_lj14_pairs = _make_symmetric(_get_lj14_pairs(topol))

    matrix_paths = [f"{reference_path}/{a}" for a in os.listdir(reference_path) if reference["matrix"] in a]
    if len(matrix_paths) > 1:
        raise ValueError(f"More than 1 matrix found in {reference_path}: {matrix_paths}")
    if not matrix_paths:
        raise FileNotFoundError(
            f"Contact matrix file(s) must be named as intramat_X_X.ndx(.gz/.h5) or "
            f"intermat_X_Y.ndx(.gz/.h5). Found instead: {reference['matrix']}"
        )

    path = matrix_paths[0]
    name = _path_to_matrix_name(path, args.root_dir)
    contact_matrix = io.read_molecular_contacts(
        path, ensemble.molecules_idx_sbtype_dictionary, reference["reference"], path.endswith(".h5")
    )

    contact_matrix = contact_matrix.add_prefix("rc_")
    contact_matrix["c6_i"] = [lj_data_dict[x][0] for x in contact_matrix["rc_ai"]]
    contact_matrix["c6_j"] = [lj_data_dict[x][0] for x in contact_matrix["rc_aj"]]
    contact_matrix["c6"] = np.sqrt(contact_matrix["c6_i"] * contact_matrix["c6_j"])
    contact_matrix["c12_i"] = [lj_data_dict[x][1] for x in contact_matrix["rc_ai"]]
    contact_matrix["c12_j"] = [lj_data_dict[x][1] for x in contact_matrix["rc_aj"]]
    contact_matrix["c12"] = np.sqrt(contact_matrix["c12_i"] * contact_matrix["c12_j"])

    # Combination-rule defaults; overridden below where explicit entries exist.
    contact_matrix["sigma_prior"] = np.where(
        contact_matrix["c6"] > 0,
        (contact_matrix["c12"] / contact_matrix["c6"]) ** (1 / 6),
        contact_matrix["c12"] ** (1 / 12) / (2.0 ** (1.0 / 6.0)),
    )
    contact_matrix["epsilon_prior"] = np.where(
        contact_matrix["c6"] > 0,
        contact_matrix["c6"] ** 2 / (4 * contact_matrix["c12"]),
        -contact_matrix["c12"],
    )

    lj_sigma_map = symmetric_lj_pairs.set_index(["ai", "aj"])["sigma"]
    lj_epsilon_map = symmetric_lj_pairs.set_index(["ai", "aj"])["epsilon"]
    lj14_sigma_map = symmetric_lj14_pairs.set_index(["ai", "aj"])["sigma"]
    lj14_epsilon_map = symmetric_lj14_pairs.set_index(["ai", "aj"])["epsilon"]

    common_indices = lj_sigma_map.index.intersection(contact_matrix.set_index(["rc_ai", "rc_aj"]).index)
    common_indices_14 = lj14_sigma_map.index.intersection(
        contact_matrix[contact_matrix["rc_same_chain"]].set_index(["rc_ai", "rc_aj"]).index
    )

    contact_matrix.loc[common_indices, "sigma_prior"] = lj_sigma_map.astype("float64")
    contact_matrix.loc[common_indices, "epsilon_prior"] = lj_epsilon_map.astype("float64")
    if not common_indices_14.empty:
        contact_matrix.loc[common_indices_14, "sigma_prior"] = lj14_sigma_map.astype("float64")
        contact_matrix.loc[common_indices_14, "epsilon_prior"] = lj14_epsilon_map.astype("float64")

    contact_matrix.drop(columns=["c6_i", "c6_j", "c12_i", "c12_j", "c6", "c12"], inplace=True)

    return name, contact_matrix


def _load_train_matrix(
    simulation,
    reference,
    ensemble,
    args,
    custom_dict,
    reference_contact_matrices,
    computed_contact_matrices,
    train_contact_matrices_general,
    computed_train_topologies,
):
    """
    Load a single training contact matrix, pair it with its reference, and
    initialise the adaptive probability thresholds.

    Training matrices are cached in ``train_contact_matrices_general`` so that
    the same trajectory data is not read from disk more than once when multiple
    references share the same training simulation.

    Parameters
    ----------
    simulation : str
        Name of the training simulation subfolder under ``inputs/SYSTEM/``.
    reference : dict
        Single ``input_refs`` entry owning this simulation.
    ensemble : MeGOEnsemble
        The meGO ensemble; ``molecules_idx_sbtype_dictionary`` is updated
        in-place with the training topology atom-index mapping.
    args : argparse.Namespace
        Parsed arguments providing ``root_dir`` and ``system``.
    custom_dict : dict
        Custom atom-name mapping dictionary forwarded to
        ``initialize_topology``.
    reference_contact_matrices : dict
        Already-populated reference matrices keyed by flat name.
    computed_contact_matrices : set
        Tracks which raw training matrices have already been read; modified
        in-place to avoid re-reading the same file.
    train_contact_matrices_general : dict
        Cache of raw (un-initialised) training matrices; modified in-place.
    computed_train_topologies : dict
        Cache mapping simulation path → ``(topology_dataframe, sbtype_dict)``
        so that the same training topology is not re-parsed when it appears in
        more than one YAML input block; modified in-place.

    Returns
    -------
    name : str
        Unique flat key for this (reference, simulation, matrix) combination.
    ref_name : str
        Flat key of the corresponding reference matrix.
    contact_matrix : pd.DataFrame
        Initialised training contact matrix with threshold columns added by
        ``initialize_molecular_contacts``.
    topology_dataframe : pd.DataFrame
        Topology DataFrame for this training simulation, used for atom-type
        validation in ``init_meGO_matrices``.

    Raises
    ------
    FileNotFoundError
        If the training topology or contact matrix file cannot be located.
    ValueError
        If more than one contact matrix file matches the requested matrix name.
    """
    simulation_path = f"{args.root_dir}/inputs/{args.system}/{simulation}"
    topology_path = f"{simulation_path}/topol.top"
    if not os.path.isfile(topology_path):
        raise FileNotFoundError(f"{topology_path} not found.")

    if simulation_path in computed_train_topologies:
        temp_topology_dataframe, ensemble.molecules_idx_sbtype_dictionary = computed_train_topologies[simulation_path]
    else:
        print("\t\t-", f"Reading {topology_path}")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            topol = parmed.load_file(topology_path)
        temp_topology_dataframe, ensemble.molecules_idx_sbtype_dictionary = _index_training_topology(topol, custom_dict)
        computed_train_topologies[simulation_path] = (temp_topology_dataframe, ensemble.molecules_idx_sbtype_dictionary)

    matrix_paths = [f"{simulation_path}/{a}" for a in os.listdir(simulation_path) if reference["matrix"] in a]
    if len(matrix_paths) > 1:
        raise ValueError(f"More than 1 matrix found in {simulation_path}: {matrix_paths}")
    if not matrix_paths:
        raise FileNotFoundError(
            f"Contact matrix file(s) must be named as intramat_X_X.ndx(.gz/.h5) or "
            f"intermat_X_Y.ndx(.gz/.h5). Found instead {reference['matrix']}"
        )

    path = matrix_paths[0]
    train_name = _path_to_matrix_name(path, args.root_dir)
    name = f"{args.system}/{reference['reference']}/{simulation}/{reference['matrix']}".replace("/", "_")

    if train_name not in computed_contact_matrices:
        train_contact_matrices_general[train_name] = io.read_molecular_contacts(
            path, ensemble.molecules_idx_sbtype_dictionary, simulation, path.endswith(".h5")
        )
        computed_contact_matrices.add(train_name)

    contact_matrix = train_contact_matrices_general[train_name].copy()

    ref_name = f"{args.system}/{reference['reference']}/{reference['matrix']}".replace("/", "_")

    contact_matrix = initialize_molecular_contacts(
        contact_matrix,
        reference_contact_matrices[ref_name],
        args,
        reference,
    )

    return name, ref_name, contact_matrix, temp_topology_dataframe


def init_meGO_matrices(ensemble, args, custom_dict):
    """
    Load and initialise all reference and training contact matrices for the
    system.

    For each entry in ``args.input_refs`` the corresponding reference topology
    and contact matrix are loaded and enriched with prior LJ parameters
    (``_load_reference_matrix``). Training matrices are then loaded, paired
    with their reference, and annotated with adaptive probability thresholds
    (``_load_train_matrix``).

    After loading, atom names in the training topologies are compared against
    those in the base topology. Any names not mapped to multi-eGO sb_types are
    reported and the program exits, requiring the user to add the missing
    entries to ``from_ff_to_multiego`` or ``custom_dict``.

    Parameters
    ----------
    ensemble : MeGOEnsemble
        Fully initialised ensemble from ``MeGOEnsemble.from_topology``.
        ``train_matrix_tuples`` is populated in-place during this call.
    args : argparse.Namespace
        Parsed arguments; ``args.input_refs``, ``args.root_dir``, and
        ``args.system`` are used.
    custom_dict : dict
        Custom atom-name mapping dictionary forwarded to
        ``initialize_topology`` for training topologies.

    Returns
    -------
    ensemble : MeGOEnsemble
        The same ensemble object with ``train_matrix_tuples`` populated.
    matrices : dict
        Dictionary with two keys:

        - ``'reference_matrices'``: ``{flat_name: pd.DataFrame}`` enriched
          reference contact matrices with prior sigma/epsilon columns.
        - ``'train_matrices'``: ``{flat_name: pd.DataFrame}`` initialised
          training contact matrices ready for ``lj.init_LJ_datasets``.
    """
    st = time.time()
    reference_contact_matrices = {}
    train_contact_matrices = {}
    train_contact_matrices_general = {}
    computed_contact_matrices = set()  # set for O(1) membership checks
    computed_train_topologies = {}  # simulation_path → (topology_df, sbtype_dict)
    train_topology_dataframe = pd.DataFrame()

    for reference in args.input_refs:
        # Skip if this (reference, matrix) combination has already been loaded —
        # two YAML blocks can legitimately share the same reference folder.
        ref_key = f"{args.system}/{reference['reference']}/{reference['matrix']}".replace("/", "_")
        if ref_key in reference_contact_matrices:
            print("\t-", f"Reference {reference['reference']} already loaded, reusing.")
            continue
        print("\t-", f"Initializing {reference['reference']} ensemble data")
        name, contact_matrix = _load_reference_matrix(reference, ensemble, args)
        reference_contact_matrices[name] = contact_matrix
        et = time.time()
        print(f"\t- Done in: {et - st:.2f} s")
        st = et

    # TODO: enable once all test cases support intra-domain splitting
    # check_intra_domain_complementarity(reference_contact_matrices)

    reference_set = set(ensemble.topology_dataframe["name"].to_list())

    for reference in args.input_refs:
        for simulation in reference["train"]:
            print("\t-", f"Initializing {simulation} ensemble data")
            name, ref_name, contact_matrix, topology_df = _load_train_matrix(
                simulation,
                reference,
                ensemble,
                args,
                custom_dict,
                reference_contact_matrices,
                computed_contact_matrices,
                train_contact_matrices_general,
                computed_train_topologies,
            )
            train_contact_matrices[name] = contact_matrix
            train_topology_dataframe = pd.concat([train_topology_dataframe, topology_df], axis=0, ignore_index=True)
            ensemble.train_matrix_tuples.append((name, ref_name))
            et = time.time()
            print(f"\t- Done in: {et - st:.2f} s")
            st = et

    del train_contact_matrices_general
    del computed_train_topologies

    # Verify that all non-hydrogen atom names in the training topologies have
    # been mapped to multi-eGO sb_types. Unmapped names indicate missing entries
    # in from_ff_to_multiego or custom_dict.
    comparison_set = set()
    for number, molecule in enumerate(ensemble.topology.molecules, 1):
        comparison_dataframe = train_topology_dataframe.loc[
            train_topology_dataframe["molecule"] == f"{number}_{molecule}"
        ]
        if not comparison_dataframe.empty:
            comparison_set |= set(
                comparison_dataframe[
                    # TODO: replace with a cleaner filter using type_definitions
                    (~comparison_dataframe["name"].str.startswith("H"))
                    | (comparison_dataframe["name"].str == "H")
                ]["name"].to_list()
            )
        else:
            raise RuntimeError("the molecule names in the training topologies do not match those in the reference")

    difference_set = comparison_set.difference(reference_set)
    if difference_set:
        print(
            f"The following atomtypes are not converted:\n{difference_set}\n"
            f'You MUST add them in "from_ff_to_multiego" dictionary to properly merge all the contacts.'
        )
        sys.exit()

    return ensemble, {
        "reference_matrices": reference_contact_matrices,
        "train_matrices": train_contact_matrices,
    }
