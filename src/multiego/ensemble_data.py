"""
ensemble_data.py — MeGOEnsemble dataclass and topology initialisation for multi-eGO.

This module defines :class:`MeGOEnsemble`, the central data container used
throughout the multi-eGO pipeline.  It holds the parsed GROMACS topology,
per-atom DataFrames, all sb_type lookup dictionaries, and the bonded-interaction
tables generated from parmed.

Typical usage::

    ensemble = MeGOEnsemble.from_topology(args, custom_dict)

The module-level helper :func:`_initialize_topology` (private, called by
``from_topology``) converts a parmed :class:`~parmed.Structure` into the
DataFrames and dictionaries required by the rest of the pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import os
import warnings

import pandas as pd
import parmed

from . import io
from . import topology as _topology
from . import type_definitions

# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------


def _assign_molecule_type(molecule_type_dict, molecule_name, molecule_topology):
    """
    Classify a molecule as ``'protein'``, ``'nucleic_acid'``, or ``'other'``
    based on its first residue name and record the result in
    ``molecule_type_dict``.

    Parameters
    ----------
    molecule_type_dict : dict
        Mapping from molecule name to molecule type, updated in-place.
    molecule_name : str
        Key under which the classification is stored.
    molecule_topology : parmed.Structure
        Topology of the molecule; only its first residue name is inspected.

    Returns
    -------
    dict
        The updated ``molecule_type_dict``.
    """
    first_residue = molecule_topology.residues[0].name

    if first_residue in type_definitions.aminoacids_list:
        molecule_type_dict[molecule_name] = "protein"
    elif first_residue in type_definitions.nucleic_acid_list:
        molecule_type_dict[molecule_name] = "nucleic_acid"
    else:
        molecule_type_dict[molecule_name] = "other"

    return molecule_type_dict


def _build_bonded_interactions(topology):
    """
    Build the bonded-interaction tables and bond-pair lists for every molecule
    in *topology*.

    Parameters
    ----------
    topology : parmed.Structure
        Loaded GROMACS topology.

    Returns
    -------
    meGO_bonded_interactions : dict
        ``{molecule: {"bonds": …, "angles": …, "dihedrals": …,
        "impropers": …, "pairs": …}}`` — one entry per molecule.
    bond_pairs : dict
        ``{molecule: [(i, j), …]}`` — covalent bond pairs used for
        1-4 parametrization.
    user_pairs : dict
        ``{molecule: parmed adjust list}`` — non-protein 1-4 pairs.
    """
    meGO_bonded_interactions = {}
    bond_pairs = {}
    user_pairs = {}

    for molecule, topol in topology.molecules.items():
        meGO_bonded_interactions[molecule] = {
            "bonds": _topology.get_bonds(topol[0].bonds),
            "angles": _topology.get_angles(topol[0].angles),
            "dihedrals": _topology.get_dihedrals(topol[0].dihedrals),
            "impropers": _topology.get_impropers(topol[0].impropers),
            "pairs": _topology.get_pairs(topol[0].adjusts),
        }
        bond_pairs[molecule] = _topology.get_bond_pairs(topol[0].bonds)
        user_pairs[molecule] = _topology.get_pairs(topol[0].adjusts)

    return meGO_bonded_interactions, bond_pairs, user_pairs


def _initialize_topology(topology, custom_dict, args):
    """
    Build the per-atom topology DataFrame and all sb_type lookup dictionaries
    from a parmed topology object.

    Each atom is assigned a unique ``sb_type`` identifier of the form
    ``atomname_moleculename_resnum`` (e.g. ``CA_ProteinA_5``). Atom names are
    first remapped through ``type_definitions.from_ff_to_multiego`` (extended
    by ``custom_dict``) to normalise force-field-specific names to multi-eGO
    conventions. Per-atom LJ parameters (``rc_c12``, ``mg_c6``, ``mg_c12``)
    are looked up from ``type_definitions.gromos_atp`` by GROMOS atom type;
    custom c12 values from ``args.custom_c12`` override the defaults.

    Parameters
    ----------
    topology : parmed.Structure
        Loaded GROMACS topology (e.g. from ``parmed.load_file``).
    custom_dict : dict
        Additional atom-name remappings for non-standard residues, merged on
        top of ``type_definitions.from_ff_to_multiego``.
    args : argparse.Namespace
        Parsed arguments; only ``args.custom_c12`` is accessed (may be
        ``None``).

    Returns
    -------
    ensemble_topology_dataframe : pd.DataFrame
        Per-atom DataFrame with columns: ``number``, ``sb_type``, ``molecule``,
        ``molecule_number``, ``molecule_name``, ``resnum``, ``cgnr``,
        ``ptype``, ``charge``, ``rc_c6``, ``rc_c12``, ``mg_c6``, ``mg_c12``,
        ``molecule_type``, plus all columns from parmed's ``to_dataframe()``.
    ensemble_molecules_idx_sbtype_dictionary : dict
        ``{molecule_key: {atom_number_str: sb_type}}`` mapping atom indices to
        sb_types per molecule, where ``molecule_key`` is
        ``"{number}_{molecule_name}"``.
    sbtype_c12_dict : dict
        ``{sb_type: rc_c12}`` — random-coil repulsive c12 per atom.
    sbtype_mg_c12_dict : dict
        ``{sb_type: mg_c12}`` — MG attractive c12 per atom.
    sbtype_mg_c6_dict : dict
        ``{sb_type: mg_c6}`` — MG attractive c6 per atom.
    sbtype_name_dict : dict
        ``{sb_type: atom_name}`` — atom name (e.g. ``CA``) per sb_type.
    sbtype_moltype_dict : dict
        ``{sb_type: molecule_type}`` — molecule classification per sb_type.
    sbtype_number_dict : dict
        ``{sb_type: atom_number}`` — sequential atom number per sb_type.
    sbtype_type_dict : dict
        ``{sb_type: atom_type}`` — GROMOS atom type (e.g. ``NL``, ``CH2``) per sb_type.
    molecule_type_dict : dict
        ``{molecule_name: molecule_type}`` — molecule-level classification.
    """
    frames = []
    new_number, col_molecule, new_resnum = [], [], []
    ensemble_molecules_idx_sbtype_dictionary = {}
    molecule_type_dict = {}

    for molecule_number, (molecule_name, molecule_topology) in enumerate(topology.molecules.items(), 1):
        molecule_type_dict = _assign_molecule_type(molecule_type_dict, molecule_name, molecule_topology[0])
        ensemble_molecules_idx_sbtype_dictionary[f"{molecule_number}_{molecule_name}"] = {}
        frames.append(molecule_topology[0].to_dataframe())
        for atom in molecule_topology[0].atoms:
            new_number.append(str(atom.idx + 1))
            col_molecule.append(f"{molecule_number}_{molecule_name}")
            new_resnum.append(str(atom.residue.number))

    ensemble_topology_dataframe = pd.concat(frames, axis=0, ignore_index=True)
    ensemble_topology_dataframe["number"] = new_number
    ensemble_topology_dataframe["molecule"] = col_molecule
    ensemble_topology_dataframe["molecule_number"] = col_molecule
    ensemble_topology_dataframe[["molecule_number", "molecule_name"]] = ensemble_topology_dataframe.molecule.str.split(
        "_", expand=True, n=1
    )
    ensemble_topology_dataframe["resnum"] = new_resnum
    ensemble_topology_dataframe["cgnr"] = ensemble_topology_dataframe["resnum"]
    ensemble_topology_dataframe["ptype"] = "A"

    # Merge the standard name-remapping dict with any user-supplied overrides.
    from_ff_to_multiego_extended = type_definitions.from_ff_to_multiego.copy()
    from_ff_to_multiego_extended.update(custom_dict)
    ensemble_topology_dataframe = ensemble_topology_dataframe.replace({"name": from_ff_to_multiego_extended})

    ensemble_topology_dataframe["sb_type"] = (
        ensemble_topology_dataframe["name"]
        + "_"
        + ensemble_topology_dataframe["molecule_name"]
        + "_"
        + ensemble_topology_dataframe["resnum"].astype(str)
    )

    atp_c12_map = dict(zip(type_definitions.gromos_atp["name"], type_definitions.gromos_atp["rc_c12"]))
    atp_mg_c6_map = dict(zip(type_definitions.gromos_atp["name"], type_definitions.gromos_atp["mg_c6"]))
    atp_mg_c12_map = dict(zip(type_definitions.gromos_atp["name"], type_definitions.gromos_atp["mg_c12"]))

    # TODO: extend custom_c12 override to mg_c6 / mg_c12 as well
    if args.custom_c12 is not None:
        custom_c12_df = io.read_custom_c12_parameters(args.custom_c12)
        atp_c12_map.update(dict(zip(custom_c12_df.name, custom_c12_df.c12)))

    ensemble_topology_dataframe["charge"] = 0.0
    ensemble_topology_dataframe["rc_c6"] = 0.0
    ensemble_topology_dataframe["rc_c12"] = ensemble_topology_dataframe["type"].map(atp_c12_map)
    ensemble_topology_dataframe["mg_c6"] = ensemble_topology_dataframe["type"].map(atp_mg_c6_map)
    ensemble_topology_dataframe["mg_c12"] = ensemble_topology_dataframe["type"].map(atp_mg_c12_map)
    ensemble_topology_dataframe["molecule_type"] = ensemble_topology_dataframe["molecule_name"].map(molecule_type_dict)

    for molecule in ensemble_molecules_idx_sbtype_dictionary:
        tmp = ensemble_topology_dataframe.loc[ensemble_topology_dataframe["molecule"] == molecule]
        ensemble_molecules_idx_sbtype_dictionary[molecule] = (
            tmp[["number", "sb_type"]].set_index("number")["sb_type"].to_dict()
        )

    idx = ensemble_topology_dataframe.set_index("sb_type")
    sbtype_c12_dict = idx["rc_c12"].to_dict()
    sbtype_mg_c12_dict = idx["mg_c12"].to_dict()
    sbtype_mg_c6_dict = idx["mg_c6"].to_dict()
    sbtype_name_dict = idx["name"].to_dict()
    sbtype_moltype_dict = idx["molecule_type"].to_dict()
    sbtype_number_dict = idx["number"].to_dict()
    sbtype_type_dict = idx["type"].to_dict()

    return (
        ensemble_topology_dataframe,
        ensemble_molecules_idx_sbtype_dictionary,
        sbtype_c12_dict,
        sbtype_mg_c12_dict,
        sbtype_mg_c6_dict,
        sbtype_name_dict,
        sbtype_moltype_dict,
        sbtype_number_dict,
        sbtype_type_dict,
        molecule_type_dict,
    )


@dataclass
class MeGOEnsemble:
    """
    Holds all topology and parameter data for a meGO system.

    Construct via MeGOEnsemble.from_topology(args, custom_dict) rather than
    directly, so that topology loading, dataframe initialization, and bonded
    interaction generation are all performed atomically.

    Attributes
    ----------
    topology : parmed.Structure
        The base system topology loaded from topol.top.
    topology_dataframe : pd.DataFrame
        Per-atom topology table with sb_type, resnum, type, c12, mg_c6/c12, etc.
    molecules_idx_sbtype_dictionary : dict
        {molecule_key: {atom_index: sb_type}} — maps atom indices to sb_type per molecule.
    sbtype_c12_dict : dict
        {sb_type: rc_c12} — random-coil repulsive c12 per sb_type.
    sbtype_mg_c12_dict : dict
        {sb_type: mg_c12} — MG attractive c12 per sb_type.
    sbtype_mg_c6_dict : dict
        {sb_type: mg_c6} — MG attractive c6 per sb_type.
    sbtype_name_dict : dict
        {sb_type: atom_name} — atom name (e.g. N, CA) per sb_type.
    sbtype_moltype_dict : dict
        {sb_type: molecule_type} — 'protein', 'nucleic_acid', or 'other' per sb_type.
    sbtype_number_dict : dict
        {sb_type: atom_number} — sequential atom number per sb_type.
    sbtype_type_dict : dict
        {sb_type: atom_type} — GROMOS atom type (e.g. NL, CH2) per sb_type.
    molecule_type_dict : dict
        {molecule_name: molecule_type} — molecule-level type classification.
    train_matrix_tuples : list
        [(train_name, ref_name), ...] — pairs of training/reference matrix keys,
        populated during init_meGO_matrices.
    meGO_bonded_interactions : dict
        {molecule: {bonds, angles, dihedrals, impropers, pairs}} DataFrames.
    bond_pairs : dict
        {molecule: list of (i, j) bond pairs} used for 1-4 parametrization.
    user_pairs : dict
        {molecule: parmed adjust list} for non-protein 1-4 pairs.
    """

    topology: parmed.Structure
    topology_dataframe: pd.DataFrame
    molecules_idx_sbtype_dictionary: dict
    sbtype_c12_dict: dict
    sbtype_mg_c12_dict: dict
    sbtype_mg_c6_dict: dict
    sbtype_name_dict: dict
    sbtype_moltype_dict: dict
    sbtype_number_dict: dict
    sbtype_type_dict: dict
    molecule_type_dict: dict
    train_matrix_tuples: list = field(default_factory=list)
    meGO_bonded_interactions: dict = field(default_factory=dict)
    bond_pairs: dict = field(default_factory=dict)
    user_pairs: dict = field(default_factory=dict)

    @classmethod
    def from_topology(cls, args, custom_dict) -> MeGOEnsemble:
        """
        Load a topology file, initialize all per-atom parameter dictionaries,
        and populate bonded interactions in a single step.

        Parameters
        ----------
        args : argparse.Namespace
            Parsed arguments; must provide root_dir, system, and custom_c12.
        custom_dict : dict
            Custom atom-name mapping dictionary (may be empty).

        Returns
        -------
        MeGOEnsemble
            Fully initialized ensemble ready for contact matrix processing.

        Raises
        ------
        FileNotFoundError
            If the topology file does not exist.
        """
        print("\t-", "Initializing system topology")
        base_topology_path = f"{args.inputs_dir}/{args.system}/topol.top"
        if not os.path.isfile(base_topology_path):
            raise FileNotFoundError(f"{base_topology_path} not found.")

        print("\t\t-", f"Reading {base_topology_path}")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            base_reference_topology = parmed.load_file(base_topology_path, {"DISULFIDE": 1})

        (
            topology_dataframe,
            molecules_idx_sbtype_dictionary,
            sbtype_c12_dict,
            sbtype_mg_c12_dict,
            sbtype_mg_c6_dict,
            sbtype_name_dict,
            sbtype_moltype_dict,
            sbtype_number_dict,
            sbtype_type_dict,
            molecule_type_dict,
        ) = _initialize_topology(base_reference_topology, custom_dict, args)

        meGO_bonded_interactions, bond_pairs, user_pairs = _build_bonded_interactions(base_reference_topology)

        return cls(
            topology=base_reference_topology,
            topology_dataframe=topology_dataframe,
            molecules_idx_sbtype_dictionary=molecules_idx_sbtype_dictionary,
            sbtype_c12_dict=sbtype_c12_dict,
            sbtype_mg_c12_dict=sbtype_mg_c12_dict,
            sbtype_mg_c6_dict=sbtype_mg_c6_dict,
            sbtype_name_dict=sbtype_name_dict,
            sbtype_moltype_dict=sbtype_moltype_dict,
            sbtype_number_dict=sbtype_number_dict,
            sbtype_type_dict=sbtype_type_dict,
            molecule_type_dict=molecule_type_dict,
            meGO_bonded_interactions=meGO_bonded_interactions,
            bond_pairs=bond_pairs,
            user_pairs=user_pairs,
        )
