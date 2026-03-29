"""
Unit tests for multi-eGO key functions.

Covers:
- lj.create_linearized_mask
- lj.set_sig_epsilon
- contacts.initialize_molecular_contacts
- contacts._path_to_matrix_name
- contacts.check_intra_domain_complementarity
- pairs.generate_bond_exclusions
- mg.generate_MG_LJ_pairs_rep
- mg.generate_MG_LJ_pairs_attr
- arguments.validate_args
"""

import argparse
import numpy as np
import pandas as pd
import pytest


def _make_args(**kwargs):
    """Build a minimal argparse.Namespace for testing."""
    defaults = dict(
        system="test",
        egos="production",
        p_to_learn=0.9995,
        epsilon_min=0.07,
        force_split=False,
        single_molecule=False,
        custom_dict=None,
        custom_c12=None,
        no_header=False,
        symmetry=[],
        symmetry_file="",
        relative_c12d=0.01,
        explicit_name="",
        config="",
        input_refs=[],
    )
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


def _lj_row(
    probability=0.9,
    rc_probability=0.01,
    md_threshold=0.1,
    rc_threshold=0.5,
    limit_rc_att=1.0,
    epsilon_prior=0.0,
    sigma_prior=0.35,
    epsilon_0=0.3,
    distance=0.38,
    rc_distance=0.40,
    rep=1e-6,
    bond_distance=10,
    same_chain=False,
    learned=False,
    mg_epsilon=0.085,
    mg_sigma=0.35,
):
    """Build a single-row DataFrame suitable for set_sig_epsilon."""
    return pd.DataFrame(
        [
            {
                "probability": probability,
                "rc_probability": rc_probability,
                "md_threshold": md_threshold,
                "rc_threshold": rc_threshold,
                "limit_rc_att": limit_rc_att,
                "epsilon_prior": epsilon_prior,
                "sigma_prior": sigma_prior,
                "epsilon_0": epsilon_0,
                "distance": distance,
                "rc_distance": rc_distance,
                "rep": rep,
                "bond_distance": bond_distance,
                "same_chain": same_chain,
                "learned": learned,
                "mg_epsilon": mg_epsilon,
                "mg_sigma": mg_sigma,
            }
        ]
    )


# ===========================================================================
# lj.create_linearized_mask
# ===========================================================================


class TestCreateLinearizedMask:

    @pytest.fixture(autouse=True)
    def load_module(self, stub_deps):
        from multiego.lj import create_linearized_mask

        self.create_linearized_mask = create_linearized_mask

    def test_exact_match(self):
        set1 = np.array(["O", "N", "C"])
        set2 = np.array(["O", "H", "C"])
        mask = self.create_linearized_mask(set1, set2, [("O", "O")])
        assert mask[0] is np.bool_(True)
        assert mask[1] is np.bool_(False)
        assert mask[2] is np.bool_(False)

    def test_symmetrize(self):
        set1 = np.array(["O", "N"])
        set2 = np.array(["N", "O"])
        mask_no_sym = self.create_linearized_mask(set1, set2, [("O", "N")], symmetrize=False)
        mask_sym = self.create_linearized_mask(set1, set2, [("O", "N")], symmetrize=True)
        assert mask_no_sym[0] is np.bool_(True)
        assert mask_no_sym[1] is np.bool_(False)
        assert mask_sym[0] is np.bool_(True)
        assert mask_sym[1] is np.bool_(True)

    def test_no_match_returns_all_false(self):
        set1 = np.array(["C", "C", "C"])
        set2 = np.array(["C", "C", "C"])
        mask = self.create_linearized_mask(set1, set2, [("O", "N")])
        assert not mask.any()

    def test_multiple_types(self):
        set1 = np.array(["O", "N", "H"])
        set2 = np.array(["O", "N", "H"])
        mask = self.create_linearized_mask(set1, set2, [("O", "O"), ("H", "H")])
        assert mask[0] is np.bool_(True)
        assert mask[1] is np.bool_(False)
        assert mask[2] is np.bool_(True)


# ===========================================================================
# lj.set_sig_epsilon
# ===========================================================================


class TestSetSigEpsilon:

    @pytest.fixture(autouse=True)
    def load_module(self, stub_deps):
        from multiego.lj import set_sig_epsilon

        self.set_sig_epsilon = set_sig_epsilon

    def test_attractive_contact_learned(self):
        """A contact above both md_threshold and limit_rc_att*rc_probability threshold
        should be marked as learned with positive epsilon."""
        df = _lj_row(
            probability=0.9,
            rc_probability=0.01,
            md_threshold=0.1,
            rc_threshold=0.5,
            limit_rc_att=1.0,
            epsilon_0=0.3,
            epsilon_prior=0.0,
            distance=0.38,
            bond_distance=10,
            same_chain=False,
        )
        result = self.set_sig_epsilon(df.copy(), _make_args())
        assert len(result) == 1
        assert result.iloc[0]["learned"] == 1
        assert result.iloc[0]["epsilon"] > 0
        assert abs(result.iloc[0]["sigma"] - 0.38 / 2 ** (1 / 6)) < 1e-10

    def test_repulsive_contact_below_attractive_threshold(self):
        """A contact above md_threshold but below attractive threshold and with
        rc_probability also above md_threshold should be repulsive (negative epsilon)."""
        df = _lj_row(
            probability=0.15,
            rc_probability=0.12,
            md_threshold=0.1,
            rc_threshold=0.5,
            limit_rc_att=2.0,
            epsilon_0=0.3,
            epsilon_prior=0.0,
            distance=0.38,
            rc_distance=0.40,
            rep=1e-6,
            bond_distance=10,
            same_chain=False,
        )
        result = self.set_sig_epsilon(df.copy(), _make_args())
        assert len(result) == 1
        assert result.iloc[0]["learned"] == 1
        assert result.iloc[0]["epsilon"] < 0

    def test_contact_below_md_threshold_not_learned(self):
        """A contact below md_threshold should retain prior values and learned=0."""
        df = _lj_row(
            probability=0.05,
            rc_probability=0.01,
            md_threshold=0.1,
            rc_threshold=0.5,
            limit_rc_att=1.0,
            epsilon_prior=0.05,
            sigma_prior=0.35,
            bond_distance=10,
            same_chain=False,
        )
        result = self.set_sig_epsilon(df.copy(), _make_args())
        assert len(result) == 1
        assert result.iloc[0]["learned"] == 0
        assert abs(result.iloc[0]["epsilon"] - 0.05) < 1e-10

    def test_14_interaction_forced_repulsive(self):
        """1-4 interactions (bond_distance == bond14_separation, same_chain)
        must always be repulsive regardless of probability."""
        df = _lj_row(
            probability=0.99,
            rc_probability=0.01,
            md_threshold=0.1,
            rc_threshold=0.5,
            limit_rc_att=1.0,
            epsilon_0=0.3,
            epsilon_prior=0.0,
            bond_distance=3,
            same_chain=True,
            rep=1e-6,
        )
        result = self.set_sig_epsilon(df.copy(), _make_args())
        assert len(result) == 1
        assert result.iloc[0]["learned"] == 1
        assert result.iloc[0]["epsilon"] < 0
        assert abs(result.iloc[0]["epsilon"] - (-1e-6)) < 1e-12

    def test_zero_epsilon_rows_dropped(self):
        """Rows where epsilon remains exactly 0 after processing should be dropped."""
        df = _lj_row(
            probability=0.05,
            md_threshold=0.1,
            epsilon_prior=0.0,
            sigma_prior=0.35,
            bond_distance=10,
            same_chain=False,
        )
        result = self.set_sig_epsilon(df.copy(), _make_args())
        assert len(result) == 0

    def test_repulsive_sigma_consistency(self):
        """For repulsive contacts, sigma must satisfy: epsilon = -(sigma * 2^(1/6))^12."""
        df = _lj_row(
            probability=0.15,
            rc_probability=0.12,
            md_threshold=0.1,
            rc_threshold=0.5,
            limit_rc_att=2.0,
            epsilon_0=0.3,
            epsilon_prior=0.0,
            rep=1e-6,
            distance=0.38,
            rc_distance=0.40,
            bond_distance=10,
            same_chain=False,
        )
        result = self.set_sig_epsilon(df.copy(), _make_args())
        row = result.iloc[0]
        expected_sigma = (-row["epsilon"]) ** (1 / 12) / (2 ** (1 / 6))
        assert abs(row["sigma"] - expected_sigma) < 1e-10


# ===========================================================================
# contacts.initialize_molecular_contacts
# ===========================================================================


class TestInitializeMolecularContacts:

    @pytest.fixture(autouse=True)
    def load_module(self, stub_deps):
        from multiego.contacts import initialize_molecular_contacts

        self.func = initialize_molecular_contacts

    def _make_contact_matrix(self, probabilities, learned=None):
        n = len(probabilities)
        if learned is None:
            learned = [True] * n
        return pd.DataFrame({"probability": probabilities, "learned": learned})

    def _make_prior_matrix(self, epsilon_prior_values):
        return pd.DataFrame({"epsilon_prior": epsilon_prior_values})

    def test_md_threshold_from_p_to_learn(self):
        """md_threshold should be the smallest probability that covers p_to_learn
        of the cumulative sum."""
        probs = [0.9, 0.5, 0.3, 0.1]
        result = self.func(
            self._make_contact_matrix(probs),
            self._make_prior_matrix([0.0] * 4),
            _make_args(p_to_learn=0.9, epsilon_min=0.07),
            {"reference": "ref", "epsilon": 0.3},
        )
        assert abs(result["md_threshold"].iloc[0] - 0.3) < 1e-10

    def test_zero_norm_gives_threshold_one(self):
        """When no contacts are learned, norm=0 and md_threshold should be 1."""
        result = self.func(
            self._make_contact_matrix([0.5, 0.3], learned=[False, False]),
            self._make_prior_matrix([0.0, 0.0]),
            _make_args(p_to_learn=0.9995, epsilon_min=0.07),
            {"reference": "ref", "epsilon": 0.3},
        )
        assert result["md_threshold"].iloc[0] == 1

    def test_epsilon_0_set_from_reference(self):
        result = self.func(
            self._make_contact_matrix([0.9]),
            self._make_prior_matrix([0.0]),
            _make_args(p_to_learn=0.9995, epsilon_min=0.07),
            {"reference": "ref", "epsilon": 0.25},
        )
        assert result["epsilon_0"].iloc[0] == 0.25

    def test_limit_rc_att_clamped_when_epsilon_prior_negative(self):
        """When epsilon_prior < 0 and limit_rc_att < 1, it should be set to 1."""
        result = self.func(
            self._make_contact_matrix([0.9]),
            self._make_prior_matrix([-0.5]),
            _make_args(p_to_learn=0.9995, epsilon_min=0.07),
            {"reference": "ref", "epsilon": 0.3},
        )
        assert result["limit_rc_att"].iloc[0] >= 1.0

    def test_reference_name_set(self):
        result = self.func(
            self._make_contact_matrix([0.5]),
            self._make_prior_matrix([0.0]),
            _make_args(p_to_learn=0.9995, epsilon_min=0.07),
            {"reference": "my_reference", "epsilon": 0.3},
        )
        assert result["reference"].iloc[0] == "my_reference"


# ===========================================================================
# contacts._path_to_matrix_name
# ===========================================================================


class TestPathToMatrixName:

    @pytest.fixture(autouse=True)
    def load_module(self, stub_deps):
        from multiego.contacts import _path_to_matrix_name

        self.func = _path_to_matrix_name

    def test_ndx_extension(self):
        assert self.func("/root/inputs/GB1/reference/intramat_1_1.ndx", "/root/inputs") == "GB1_reference_intramat_1_1"

    def test_ndx_gz_extension(self):
        assert (
            self.func("/root/inputs/GB1/reference/intramat_1_1.ndx.gz", "/root/inputs") == "GB1_reference_intramat_1_1"
        )

    def test_h5_extension(self):
        assert (
            self.func("/root/inputs/GB1/reference/intramat_1_1.ndx.h5", "/root/inputs") == "GB1_reference_intramat_1_1"
        )

    def test_intermat(self):
        assert self.func("/root/inputs/GB1/reference/intermat_1_2.ndx", "/root/inputs") == "GB1_reference_intermat_1_2"


# ===========================================================================
# contacts.check_intra_domain_complementarity
# ===========================================================================


class TestCheckIntraDomainComplementarity:

    @pytest.fixture(autouse=True)
    def load_module(self, stub_deps):
        from multiego.contacts import check_intra_domain_complementarity

        self.func = check_intra_domain_complementarity

    def test_complementary_flags_pass(self):
        """Two references covering non-overlapping atoms should not raise."""
        matrices = {
            "sysA_ref1_intramat_1_1": pd.DataFrame({"rc_learned": [True, False, False]}),
            "sysB_ref2_intramat_1_1": pd.DataFrame({"rc_learned": [False, True, False]}),
        }
        self.func(matrices)

    def test_overlapping_flags_raise(self):
        """Two references both claiming the same atom should raise ValueError."""
        matrices = {
            "sysA_ref1_intramat_1_1": pd.DataFrame({"rc_learned": [True, False]}),
            "sysB_ref2_intramat_1_1": pd.DataFrame({"rc_learned": [True, False]}),
        }
        with pytest.raises(ValueError, match="complementarity"):
            self.func(matrices)

    def test_single_reference_always_passes(self):
        """A single reference matrix can never violate complementarity."""
        matrices = {
            "sysA_ref1_intramat_1_1": pd.DataFrame({"rc_learned": [True, True, True]}),
        }
        self.func(matrices)


# ===========================================================================
# pairs.generate_bond_exclusions
# ===========================================================================


class TestGenerateBondExclusions:

    @pytest.fixture(autouse=True)
    def load_module(self, stub_deps):
        from multiego.pairs import generate_bond_exclusions

        self.func = generate_bond_exclusions

    def _linear_topology(self, n):
        return pd.DataFrame({"number": [str(i) for i in range(1, n + 1)]})

    def _linear_bonds(self, n):
        return [(str(i), str(i + 1)) for i in range(1, n)]

    def test_1_4_pairs_present(self):
        """Atoms exactly 3 bonds apart should appear in p14."""
        _, p14, _ = self.func(self._linear_topology(5), self._linear_bonds(5))
        assert "1_4" in p14
        assert "4_1" in p14

    def test_exclusion_bonds_within_3(self):
        """Atoms within 3 bonds should be in exclusion_bonds."""
        excl, _, _ = self.func(self._linear_topology(5), self._linear_bonds(5))
        assert "1_2" in excl
        assert "1_3" in excl
        assert "1_4" in excl

    def test_beyond_max_not_excluded(self):
        """Atoms more than max_bond_separation (5) bonds apart should not appear."""
        excl, _, nth = self.func(self._linear_topology(8), self._linear_bonds(8))
        assert "1_7" not in excl
        assert "1_7" not in nth

    def test_nth_bonds_includes_up_to_max(self):
        """nth_bonds should cover atoms up to max_bond_separation bonds away."""
        _, _, nth = self.func(self._linear_topology(7), self._linear_bonds(7))
        assert "1_6" in nth
        assert "6_1" in nth

    def test_disconnected_atom_not_in_exclusions(self):
        """An atom with no bonds should not appear in any exclusion list."""
        topo = pd.DataFrame({"number": ["1", "2", "3"]})
        excl, p14, nth = self.func(topo, [("1", "2")])
        assert "1_3" not in excl
        assert "3_1" not in excl


# ===========================================================================
# mg.generate_MG_LJ_pairs_rep
# ===========================================================================


class TestGenerateMGLJPairsRep:

    @pytest.fixture(autouse=True)
    def load_module(self, stub_deps):
        from multiego.mg import generate_MG_LJ_pairs_rep

        self.func = generate_MG_LJ_pairs_rep

    def test_self_pairs_when_equal_lists(self):
        rc_c12 = {"A_mol_1": 1e-6, "A_mol_2": 2e-6}
        result = self.func(["A_mol_1", "A_mol_2"], ["A_mol_1", "A_mol_2"], rc_c12)
        assert len(result) == 4

    def test_cross_pairs_when_different_lists(self):
        rc_c12 = {"O_mol_1": 1e-6, "N_mol_1": 2e-6}
        result = self.func(["O_mol_1"], ["N_mol_1"], rc_c12)
        pairs = set(zip(result["ai"], result["aj"]))
        assert ("O_mol_1", "N_mol_1") in pairs
        assert ("N_mol_1", "O_mol_1") in pairs

    def test_repulsive_epsilon_negative(self):
        rc_c12 = {"A_mol_1": 1e-6}
        result = self.func(["A_mol_1"], ["A_mol_1"], rc_c12)
        assert (result["epsilon"] < 0).all()

    def test_c6_is_zero(self):
        rc_c12 = {"A_mol_1": 1e-6}
        result = self.func(["A_mol_1"], ["A_mol_1"], rc_c12)
        assert (result["c6"] == 0.0).all()

    def test_fixed_c12_rep(self):
        rc_c12 = {"A_mol_1": 1e-6, "B_mol_1": 2e-6}
        result = self.func(["A_mol_1"], ["B_mol_1"], rc_c12, c12_rep=5e-7)
        assert (result["c12"] == 5e-7).all()

    def test_geometric_mean_c12_when_not_fixed(self):
        c12_a, c12_b = 1e-6, 4e-6
        rc_c12 = {"A_mol_1": c12_a, "B_mol_1": c12_b}
        result = self.func(["A_mol_1"], ["B_mol_1"], rc_c12)
        expected = np.sqrt(c12_a * c12_b)
        assert all(abs(v - expected) < 1e-15 for v in result["c12"])

    def test_sigma_consistent_with_c12(self):
        rc_c12 = {"A_mol_1": 1e-6}
        result = self.func(["A_mol_1"], ["A_mol_1"], rc_c12)
        for _, row in result.iterrows():
            expected_sigma = row["c12"] ** (1 / 12) / 2 ** (1 / 6)
            assert abs(row["sigma"] - expected_sigma) < 1e-12


# ===========================================================================
# mg.generate_MG_LJ_pairs_attr
# ===========================================================================


class TestGenerateMGLJPairsAttr:

    @pytest.fixture(autouse=True)
    def load_module(self, stub_deps):
        from multiego.mg import generate_MG_LJ_pairs_attr

        self.func = generate_MG_LJ_pairs_attr

    def test_attractive_epsilon_positive(self):
        mg_c12 = {"A_mol_1": 1e-5}
        mg_c6 = {"A_mol_1": 1e-3}
        result = self.func(["A_mol_1"], ["A_mol_1"], mg_c12, mg_c6, epsilon=0.1)
        assert (result["epsilon"] > 0).all()

    def test_c6_c12_lj_formula(self):
        mg_c12 = {"A_mol_1": 1e-5}
        mg_c6 = {"A_mol_1": 1e-3}
        eps = 0.1
        result = self.func(["A_mol_1"], ["A_mol_1"], mg_c12, mg_c6, epsilon=eps)
        for _, row in result.iterrows():
            assert abs(row["c6"] - 4 * eps * row["sigma"] ** 6) < 1e-15
            assert abs(row["c12"] - 4 * eps * row["sigma"] ** 12) < 1e-15

    def test_sigma_only_raises(self):
        with pytest.raises(ValueError):
            self.func(["A_mol_1"], ["A_mol_1"], {}, {}, epsilon=None, sigma=0.35)

    def test_sigma_from_combination_rules(self):
        c12_val, c6_val = 1e-5, 1e-3
        mg_c12 = {"A_mol_1": c12_val}
        mg_c6 = {"A_mol_1": c6_val}
        result = self.func(["A_mol_1"], ["A_mol_1"], mg_c12, mg_c6, epsilon=0.1)
        expected_sigma = (c12_val / c6_val) ** (1 / 6)
        assert abs(result.iloc[0]["sigma"] - expected_sigma) < 1e-12


# ===========================================================================
# arguments.validate_args
# ===========================================================================


class TestValidateArgs:

    @pytest.fixture(autouse=True)
    def load_module(self, stub_deps):
        # Give the io stub the specific attribute validate_args needs
        stub_deps["multiego.io"].read_custom_c12_parameters = lambda path: pd.DataFrame(
            {"name": ["CH2"], "c12": [1e-5]}
        )
        from multiego.arguments import validate_args

        self.validate = validate_args

    def _valid_ref(self):
        return {"matrix": "intramat_1_1", "epsilon": 0.3, "train": ["md_monomer"], "reference": "reference"}

    def test_valid_args_pass(self):
        self.validate(_make_args(system="GB1", egos="production", input_refs=[self._valid_ref()]))

    def test_missing_system_exits(self):
        with pytest.raises(SystemExit):
            self.validate(_make_args(system="", egos="production", input_refs=[]))

    def test_missing_egos_exits(self):
        with pytest.raises(SystemExit):
            self.validate(_make_args(system="GB1", egos=None, input_refs=[]))

    def test_epsilon_min_zero_exits(self):
        with pytest.raises(SystemExit):
            self.validate(_make_args(system="GB1", egos="production", input_refs=[], epsilon_min=0.0))

    def test_epsilon_below_epsilon_min_exits(self):
        ref = self._valid_ref()
        ref["epsilon"] = 0.01
        with pytest.raises(SystemExit):
            self.validate(_make_args(system="GB1", egos="production", input_refs=[ref], epsilon_min=0.07))

    def test_missing_required_ref_key_raises(self):
        bad_ref = {"epsilon": 0.3, "train": ["md"], "reference": "ref"}
        with pytest.raises(ValueError, match="Missing required keys"):
            self.validate(_make_args(system="GB1", egos="production", input_refs=[bad_ref]))

    def test_symmetry_conflict_exits(self):
        with pytest.raises(SystemExit):
            self.validate(
                _make_args(
                    system="GB1",
                    egos="production",
                    input_refs=[self._valid_ref()],
                    symmetry=["ARG NH1 NH2"],
                    symmetry_file="symmetry.dat",
                )
            )

    def test_p_to_learn_warning_printed(self, capsys):
        self.validate(
            _make_args(
                system="GB1",
                egos="production",
                input_refs=[self._valid_ref()],
                p_to_learn=0.5,
            )
        )
        captured = capsys.readouterr()
        assert "WARNING" in captured.out
        assert "p_to_learn" in captured.out
