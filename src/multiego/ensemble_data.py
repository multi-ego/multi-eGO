from dataclasses import dataclass, field

import pandas as pd


@dataclass
class MeGOEnsemble:
    """
    Holds all topology and parameter data for a meGO system.

    Fields set at initialization (init_meGO_ensemble)
    --------------------------------------------------
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

    Fields set later by generate_bonded_interactions
    -------------------------------------------------
    meGO_bonded_interactions : dict
        {molecule: {bonds, angles, dihedrals, impropers, pairs}} DataFrames.
    bond_pairs : dict
        {molecule: list of (i, j) bond pairs} used for 1-4 parametrization.
    user_pairs : dict
        {molecule: parmed adjust list} for non-protein 1-4 pairs.
    """

    # --- set at init ---
    topology: object  # parmed.Structure; not typed to avoid importing parmed here
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

    # --- set later by generate_bonded_interactions ---
    meGO_bonded_interactions: dict = field(default_factory=dict)
    bond_pairs: dict = field(default_factory=dict)
    user_pairs: dict = field(default_factory=dict)
