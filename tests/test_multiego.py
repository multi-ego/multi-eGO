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

import sys
import types
import numpy as np
import pandas as pd
import pytest


def _make_args(**kwargs):
    """Build a minimal argparse.Namespace for testing."""
    import argparse

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

    def setup_method(self):
        # Import directly from file to avoid package init
        from multiego.lj import create_linearized_mask

        # Stub heavy deps before loading
        for dep in ["multiego.type_definitions", "multiego.mg", "multiego.model_config", "multiego.io"]:
            if dep not in sys.modules:
                stub = types.ModuleType(dep)
                # model_config needs a config object with attributes
                if dep == "multiego.model_config":
                    cfg = types.SimpleNamespace(max_bond_separation=5, bond14_separation=3)
                    stub.config = cfg
                sys.modules[dep] = stub
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
        # Without symmetrize: only (O, N) matches index 0
        mask_no_sym = self.create_linearized_mask(set1, set2, [("O", "N")], symmetrize=False)
        # With symmetrize: both (O,N) and (N,O) are in types, index 1 also matches
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
        assert mask[0] is np.bool_(True)  # O-O
        assert mask[1] is np.bool_(False)  # N-N not in types
        assert mask[2] is np.bool_(True)  # H-H


# ===========================================================================
# lj.set_sig_epsilon
# ===========================================================================


class TestSetSigEpsilon:

    @pytest.fixture(autouse=True)
    def load_module(self):
        for dep in ["multiego.type_definitions", "multiego.mg", "multiego.model_config", "multiego.io"]:
            if dep not in sys.modules:
                stub = types.ModuleType(dep)
                if dep == "multiego.model_config":
                    stub.config = types.SimpleNamespace(max_bond_separation=5, bond14_separation=3)
                sys.modules[dep] = stub
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
        args = _make_args()
        result = self.set_sig_epsilon(df.copy(), args)
        assert len(result) == 1
        assert result.iloc[0]["learned"] == 1
        assert result.iloc[0]["epsilon"] > 0
        # sigma should be derived from distance
        assert abs(result.iloc[0]["sigma"] - 0.38 / 2 ** (1 / 6)) < 1e-10

    def test_repulsive_contact_below_attractive_threshold(self):
        """A contact above md_threshold but below attractive threshold and with
        rc_probability also above md_threshold should be repulsive (negative epsilon)."""
        df = _lj_row(
            probability=0.15,
            rc_probability=0.12,
            md_threshold=0.1,
            rc_threshold=0.5,
            limit_rc_att=2.0,  # limit_rc_att * rc_probability = 0.24 > 0.15
            epsilon_0=0.3,
            epsilon_prior=0.0,
            distance=0.38,
            rc_distance=0.40,
            rep=1e-6,
            bond_distance=10,
            same_chain=False,
        )
        args = _make_args()
        result = self.set_sig_epsilon(df.copy(), args)
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
        args = _make_args()
        result = self.set_sig_epsilon(df.copy(), args)
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
            bond_distance=3,  # bond14_separation
            same_chain=True,
            rep=1e-6,
        )
        args = _make_args()
        result = self.set_sig_epsilon(df.copy(), args)
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
        args = _make_args()
        result = self.set_sig_epsilon(df.copy(), args)
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
        args = _make_args()
        result = self.set_sig_epsilon(df.copy(), args)
        row = result.iloc[0]
        eps = row["epsilon"]
        sigma = row["sigma"]
        assert eps < 0
        expected_sigma = (-eps) ** (1 / 12) / (2 ** (1 / 6))
        assert abs(sigma - expected_sigma) < 1e-10


# ===========================================================================
# contacts.initialize_molecular_contacts
# ===========================================================================


class TestInitializeMolecularContacts:

    @pytest.fixture(autouse=True)
    def load_module(self):

        # --- stub ensemble_data WITH the required symbol ---
        ensemble_stub = types.ModuleType("multiego.ensemble_data")
        ensemble_stub.MeGOEnsemble = type("MeGOEnsemble", (), {})
        sys.modules["multiego.ensemble_data"] = ensemble_stub

        # --- stub other lightweight deps ---
        for dep in ["multiego.io", "multiego.contacts_init", "parmed"]:
            sys.modules.setdefault(dep, types.ModuleType(dep))
            if dep not in sys.modules:
                sys.modules[dep] = types.ModuleType(dep)
        from multiego.contacts import initialize_molecular_contacts

        self.func = initialize_molecular_contacts

    def _make_contact_matrix(self, probabilities, learned=None):
        n = len(probabilities)
        if learned is None:
            learned = [True] * n
        return pd.DataFrame(
            {
                "probability": probabilities,
                "learned": learned,
            }
        )

    def _make_prior_matrix(self, epsilon_prior_values):
        return pd.DataFrame({"epsilon_prior": epsilon_prior_values})

    def test_md_threshold_from_p_to_learn(self):
        """md_threshold should be the smallest probability that covers p_to_learn
        of the cumulative sum."""
        probs = [0.9, 0.5, 0.3, 0.1]
        contact = self._make_contact_matrix(probs)
        prior = self._make_prior_matrix([0.0] * 4)
        args = _make_args(p_to_learn=0.9, epsilon_min=0.07)
        reference = {"reference": "ref", "epsilon": 0.3}
        result = self.func(contact.copy(), prior, args, reference)
        # p_to_learn=0.9 means threshold = smallest prob such that cumsum/total >= 0.9
        # sorted desc: [0.9, 0.5, 0.3, 0.1], total=1.8
        # cumsum: 0.9, 1.4, 1.7, 1.8  /  1.8 = 0.5, 0.78, 0.94, 1.0
        # first >= 0.9 at index 2 -> probability 0.3
        assert abs(result["md_threshold"].iloc[0] - 0.3) < 1e-10

    def test_zero_norm_gives_threshold_one(self):
        """When no contacts are learned, norm=0 and md_threshold should be 1."""
        contact = self._make_contact_matrix([0.5, 0.3], learned=[False, False])
        prior = self._make_prior_matrix([0.0, 0.0])
        args = _make_args(p_to_learn=0.9995, epsilon_min=0.07)
        reference = {"reference": "ref", "epsilon": 0.3}
        result = self.func(contact.copy(), prior, args, reference)
        assert result["md_threshold"].iloc[0] == 1

    def test_epsilon_0_set_from_reference(self):
        contact = self._make_contact_matrix([0.9])
        prior = self._make_prior_matrix([0.0])
        args = _make_args(p_to_learn=0.9995, epsilon_min=0.07)
        reference = {"reference": "ref", "epsilon": 0.25}
        result = self.func(contact.copy(), prior, args, reference)
        assert result["epsilon_0"].iloc[0] == 0.25

    def test_limit_rc_att_clamped_when_epsilon_prior_negative(self):
        """When epsilon_prior < 0 and limit_rc_att < 1, it should be set to 1."""
        contact = self._make_contact_matrix([0.9])
        prior = self._make_prior_matrix([-0.5])  # negative epsilon_prior
        args = _make_args(p_to_learn=0.9995, epsilon_min=0.07)
        reference = {"reference": "ref", "epsilon": 0.3}
        result = self.func(contact.copy(), prior, args, reference)
        # With negative epsilon_prior and the formula, limit_rc_att could go < 1;
        # the function must clamp it to 1
        assert result["limit_rc_att"].iloc[0] >= 1.0

    def test_reference_name_set(self):
        contact = self._make_contact_matrix([0.5])
        prior = self._make_prior_matrix([0.0])
        args = _make_args(p_to_learn=0.9995, epsilon_min=0.07)
        reference = {"reference": "my_reference", "epsilon": 0.3}
        result = self.func(contact.copy(), prior, args, reference)
        assert result["reference"].iloc[0] == "my_reference"


# ===========================================================================
# contacts._path_to_matrix_name
# ===========================================================================


class TestPathToMatrixName:

    @pytest.fixture(autouse=True)
    def load_module(self):

        # --- stub ensemble_data WITH the required symbol ---
        ensemble_stub = types.ModuleType("multiego.ensemble_data")
        ensemble_stub.MeGOEnsemble = type("MeGOEnsemble", (), {})
        sys.modules["multiego.ensemble_data"] = ensemble_stub

        # --- stub other lightweight deps ---
        for dep in ["multiego.io", "multiego.contacts_init", "parmed"]:
            if dep not in sys.modules:
                sys.modules[dep] = types.ModuleType(dep)
        from multiego.contacts import _path_to_matrix_name

        self.func = _path_to_matrix_name

    def test_ndx_extension(self):
        result = self.func("/root/inputs/GB1/reference/intramat_1_1.ndx", "/root")
        assert result == "GB1_reference_intramat_1_1"

    def test_ndx_gz_extension(self):
        result = self.func("/root/inputs/GB1/reference/intramat_1_1.ndx.gz", "/root")
        assert result == "GB1_reference_intramat_1_1"

    def test_h5_extension(self):
        result = self.func("/root/inputs/GB1/reference/intramat_1_1.ndx.h5", "/root")
        assert result == "GB1_reference_intramat_1_1"

    def test_intermat(self):
        result = self.func("/root/inputs/GB1/reference/intermat_1_2.ndx", "/root")
        assert result == "GB1_reference_intermat_1_2"


# ===========================================================================
# contacts.check_intra_domain_complementarity
# ===========================================================================


class TestCheckIntraDomainComplementarity:

    @pytest.fixture(autouse=True)
    def load_module(self):

        # --- stub ensemble_data WITH the required symbol ---
        ensemble_stub = types.ModuleType("multiego.ensemble_data")
        ensemble_stub.MeGOEnsemble = type("MeGOEnsemble", (), {})
        sys.modules["multiego.ensemble_data"] = ensemble_stub

        # --- stub other lightweight deps ---
        for dep in ["multiego.io", "multiego.contacts_init", "parmed"]:
            if dep not in sys.modules:
                sys.modules[dep] = types.ModuleType(dep)
        from multiego.contacts import check_intra_domain_complementarity

        self.func = check_intra_domain_complementarity

    def test_complementary_flags_pass(self):
        """Two references covering non-overlapping atoms should not raise."""
        matrices = {
            "sysA_ref1_intramat_1_1": pd.DataFrame({"rc_learned": [True, False, False]}),
            "sysB_ref2_intramat_1_1": pd.DataFrame({"rc_learned": [False, True, False]}),
        }
        self.func(matrices)  # should not raise

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
        self.func(matrices)  # should not raise


# ===========================================================================
# pairs.generate_bond_exclusions
# ===========================================================================


class TestGenerateBondExclusions:

    @pytest.fixture(autouse=True)
    def load_module(self):
        stub_td = types.ModuleType("multiego.type_definitions")
        stub_mc = types.ModuleType("multiego.model_config")
        stub_mc.config = types.SimpleNamespace(max_bond_separation=5, bond14_separation=3)
        sys.modules.setdefault("multiego.type_definitions", stub_td)
        sys.modules["multiego.model_config"] = stub_mc
        from multiego.pairs import generate_bond_exclusions

        self.func = generate_bond_exclusions

    def _linear_topology(self, n):
        """Build a reduced_topology DataFrame for a linear chain of n atoms."""
        return pd.DataFrame({"number": [str(i) for i in range(1, n + 1)]})

    def _linear_bonds(self, n):
        """Bond pairs for a linear chain 1-2-3-...-n."""
        return [(str(i), str(i + 1)) for i in range(1, n)]

    def test_1_4_pairs_present(self):
        """Atoms exactly 3 bonds apart should appear in p14."""
        topo = self._linear_topology(5)
        bonds = self._linear_bonds(5)
        _, p14, _ = self.func(topo, bonds)
        assert "1_4" in p14
        assert "4_1" in p14

    def test_exclusion_bonds_within_3(self):
        """Atoms within 3 bonds should be in exclusion_bonds."""
        topo = self._linear_topology(5)
        bonds = self._linear_bonds(5)
        excl, _, _ = self.func(topo, bonds)
        assert "1_2" in excl
        assert "1_3" in excl
        assert "1_4" in excl  # 3 bonds away, still excluded

    def test_beyond_max_not_excluded(self):
        """Atoms more than max_bond_separation (5) bonds apart should not appear."""
        topo = self._linear_topology(8)
        bonds = self._linear_bonds(8)
        excl, _, nth = self.func(topo, bonds)
        # atom 1 and atom 7 are 6 bonds apart — outside max_bond_separation=5
        assert "1_7" not in excl
        assert "1_7" not in nth

    def test_nth_bonds_includes_up_to_max(self):
        """nth_bonds should cover atoms up to max_bond_separation bonds away."""
        topo = self._linear_topology(7)
        bonds = self._linear_bonds(7)
        _, _, nth = self.func(topo, bonds)
        # atom 1 and atom 6 are exactly 5 bonds apart — should be in nth_bonds
        assert "1_6" in nth
        assert "6_1" in nth

    def test_disconnected_atom_not_in_exclusions(self):
        """An atom with no bonds should not appear in any exclusion list."""
        topo = pd.DataFrame({"number": ["1", "2", "3"]})
        bonds = [("1", "2")]  # atom 3 is disconnected
        excl, p14, nth = self.func(topo, bonds)
        assert "1_3" not in excl
        assert "3_1" not in excl


# ===========================================================================
# mg.generate_MG_LJ_pairs_rep
# ===========================================================================


class TestGenerateMGLJPairsRep:

    @pytest.fixture(autouse=True)
    def load_module(self):
        stub = types.ModuleType("multiego.type_definitions")
        stub.special_non_local = []
        sys.modules.setdefault("multiego.type_definitions", stub)
        from multiego.mg import generate_MG_LJ_pairs_rep

        self.func = generate_MG_LJ_pairs_rep

    def test_self_pairs_when_equal_lists(self):
        """When sbtype1 == sbtype2, only self-pairs should be generated."""
        rc_c12 = {"A_mol_1": 1e-6, "A_mol_2": 2e-6}
        result = self.func(["A_mol_1", "A_mol_2"], ["A_mol_1", "A_mol_2"], rc_c12)
        # itertools.product with repeat=2 gives (1,1),(1,2),(2,1),(2,2)
        assert len(result) == 4

    def test_cross_pairs_when_different_lists(self):
        """When sbtype1 != sbtype2, cross-pairs in both directions."""
        rc_c12 = {"O_mol_1": 1e-6, "N_mol_1": 2e-6}
        result = self.func(["O_mol_1"], ["N_mol_1"], rc_c12)
        pairs = set(zip(result["ai"], result["aj"]))
        assert ("O_mol_1", "N_mol_1") in pairs
        assert ("N_mol_1", "O_mol_1") in pairs

    def test_repulsive_epsilon_negative(self):
        """All epsilon values should be negative (repulsive)."""
        rc_c12 = {"A_mol_1": 1e-6}
        result = self.func(["A_mol_1"], ["A_mol_1"], rc_c12)
        assert (result["epsilon"] < 0).all()

    def test_c6_is_zero(self):
        """Repulsive pairs have no attractive component: c6 must be 0."""
        rc_c12 = {"A_mol_1": 1e-6}
        result = self.func(["A_mol_1"], ["A_mol_1"], rc_c12)
        assert (result["c6"] == 0.0).all()

    def test_fixed_c12_rep(self):
        """When c12_rep is provided, all pairs should use that value."""
        rc_c12 = {"A_mol_1": 1e-6, "B_mol_1": 2e-6}
        fixed = 5e-7
        result = self.func(["A_mol_1"], ["B_mol_1"], rc_c12, c12_rep=fixed)
        assert (result["c12"] == fixed).all()

    def test_geometric_mean_c12_when_not_fixed(self):
        """Without c12_rep, c12 should be the geometric mean of the pair."""
        c12_a, c12_b = 1e-6, 4e-6
        rc_c12 = {"A_mol_1": c12_a, "B_mol_1": c12_b}
        result = self.func(["A_mol_1"], ["B_mol_1"], rc_c12)
        # Only one direction since sbtype1 != sbtype2 but deduped via set
        expected = np.sqrt(c12_a * c12_b)
        assert all(abs(v - expected) < 1e-15 for v in result["c12"])

    def test_sigma_consistent_with_c12(self):
        """sigma should satisfy: c12 = (sigma * 2^(1/6))^12."""
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
    def load_module(self):
        stub = types.ModuleType("multiego.type_definitions")
        stub.special_non_local = []
        sys.modules.setdefault("multiego.type_definitions", stub)
        from multiego.mg import generate_MG_LJ_pairs_attr

        self.func = generate_MG_LJ_pairs_attr

    def test_attractive_epsilon_positive(self):
        """Attractive pairs should have positive epsilon."""
        mg_c12 = {"A_mol_1": 1e-5}
        mg_c6 = {"A_mol_1": 1e-3}
        result = self.func(["A_mol_1"], ["A_mol_1"], mg_c12, mg_c6, epsilon=0.1)
        assert (result["epsilon"] > 0).all()

    def test_c6_c12_lj_formula(self):
        """c6 and c12 must satisfy: c6 = 4*eps*sigma^6, c12 = 4*eps*sigma^12."""
        mg_c12 = {"A_mol_1": 1e-5}
        mg_c6 = {"A_mol_1": 1e-3}
        eps = 0.1
        result = self.func(["A_mol_1"], ["A_mol_1"], mg_c12, mg_c6, epsilon=eps)
        for _, row in result.iterrows():
            assert abs(row["c6"] - 4 * eps * row["sigma"] ** 6) < 1e-15
            assert abs(row["c12"] - 4 * eps * row["sigma"] ** 12) < 1e-15

    def test_sigma_only_raises(self):
        """Providing sigma without epsilon should raise ValueError."""
        with pytest.raises(ValueError):
            self.func(["A_mol_1"], ["A_mol_1"], {}, {}, epsilon=None, sigma=0.35)

    def test_sigma_from_combination_rules(self):
        """When epsilon is given but sigma is None, sigma should be derived from
        the geometric mean combination rules of mg_c12 and mg_c6."""
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
    def load_module(self):
        # arguments.py imports from .io — stub it
        stub_io = types.ModuleType("multiego.io")
        stub_io.read_custom_c12_parameters = lambda path: pd.DataFrame({"name": ["CH2"], "c12": [1e-5]})
        sys.modules["multiego.io"] = stub_io
        from multiego.arguments import validate_args

        self.validate = validate_args

    def _valid_ref(self):
        return {"matrix": "intramat_1_1", "epsilon": 0.3, "train": ["md_monomer"], "reference": "reference"}

    def test_valid_args_pass(self):
        args = _make_args(
            system="GB1",
            egos="production",
            input_refs=[self._valid_ref()],
            epsilon_min=0.07,
            p_to_learn=0.9995,
        )
        self.validate(args)  # should not raise or exit

    def test_missing_system_exits(self):
        args = _make_args(system="", egos="production", input_refs=[])
        with pytest.raises(SystemExit):
            self.validate(args)

    def test_missing_egos_exits(self):
        args = _make_args(system="GB1", egos=None, input_refs=[])
        with pytest.raises(SystemExit):
            self.validate(args)

    def test_epsilon_min_zero_exits(self):
        args = _make_args(system="GB1", egos="production", input_refs=[], epsilon_min=0.0)
        with pytest.raises(SystemExit):
            self.validate(args)

    def test_epsilon_below_epsilon_min_exits(self):
        ref = self._valid_ref()
        ref["epsilon"] = 0.01  # below default epsilon_min=0.07
        args = _make_args(system="GB1", egos="production", input_refs=[ref], epsilon_min=0.07)
        with pytest.raises(SystemExit):
            self.validate(args)

    def test_missing_required_ref_key_raises(self):
        bad_ref = {"epsilon": 0.3, "train": ["md"], "reference": "ref"}  # missing matrix
        args = _make_args(system="GB1", egos="production", input_refs=[bad_ref])
        with pytest.raises(ValueError, match="Missing required keys"):
            self.validate(args)

    def test_symmetry_conflict_exits(self):
        args = _make_args(
            system="GB1",
            egos="production",
            input_refs=[self._valid_ref()],
            symmetry=["ARG NH1 NH2"],
            symmetry_file="symmetry.dat",
        )
        with pytest.raises(SystemExit):
            self.validate(args)

    def test_p_to_learn_warning_printed(self, capsys):
        args = _make_args(
            system="GB1",
            egos="production",
            input_refs=[self._valid_ref()],
            p_to_learn=0.5,  # below 0.9 threshold
        )
        self.validate(args)
        captured = capsys.readouterr()
        assert "WARNING" in captured.out
        assert "p_to_learn" in captured.out
