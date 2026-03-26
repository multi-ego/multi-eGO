from __future__ import annotations

from dataclasses import dataclass, field
import os
import warnings

import pandas as pd
import parmed

from . import contacts_init
from . import topology as _topology


@dataclass
class MeGOEnsemble:
    """
    Holds all topology and parameter data for a meGO system.

    Construct via MeGOEnsemble.from_topology(args, custom_dict) rather than
    directly, so that topology loading, dataframe initialization, and bonded
    interaction generation are all performed atomically.

    Fields
    ------
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

    topology: object  # parmed.Structure
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
        base_topology_path = f"{args.root_dir}/inputs/{args.system}/topol.top"
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
            molecule_type_dict,
        ) = contacts_init.initialize_topology(base_reference_topology, custom_dict, args)

        ensemble = cls(
            topology=base_reference_topology,
            topology_dataframe=topology_dataframe,
            molecules_idx_sbtype_dictionary=molecules_idx_sbtype_dictionary,
            sbtype_c12_dict=sbtype_c12_dict,
            sbtype_mg_c12_dict=sbtype_mg_c12_dict,
            sbtype_mg_c6_dict=sbtype_mg_c6_dict,
            sbtype_name_dict=sbtype_name_dict,
            sbtype_moltype_dict=sbtype_moltype_dict,
            sbtype_number_dict=topology_dataframe[["sb_type", "number"]].set_index("sb_type")["number"].to_dict(),
            sbtype_type_dict={key: name for key, name in topology_dataframe[["sb_type", "type"]].values},
            molecule_type_dict=molecule_type_dict,
        )

        ensemble._init_bonded_interactions()
        return ensemble

    def _init_bonded_interactions(self) -> None:
        """
        Populate meGO_bonded_interactions, bond_pairs, and user_pairs from
        the loaded topology. Called automatically by from_topology.
        """
        self.meGO_bonded_interactions = {}
        self.bond_pairs = {}
        self.user_pairs = {}

        for molecule, topol in self.topology.molecules.items():
            self.meGO_bonded_interactions[molecule] = {
                "bonds": _topology.get_bonds(topol[0].bonds),
                "angles": _topology.get_angles(topol[0].angles),
                "dihedrals": _topology.get_dihedrals(topol[0].dihedrals),
                "impropers": _topology.get_impropers(topol[0].impropers),
                "pairs": _topology.get_pairs(topol[0].adjusts),
            }
            self.bond_pairs[molecule] = _topology.get_bond_pairs(topol[0].bonds)
            self.user_pairs[molecule] = _topology.get_pairs(topol[0].adjusts)
