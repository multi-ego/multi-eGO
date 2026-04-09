"""
Unit tests for the pure-computation functions in tools/make_mat/make_mat.py.

Heavy dependencies (parmed, multiego.*, h5py) are stubbed out.
scipy is used when available; a numpy-based fallback for logsumexp is
injected otherwise, so these tests can run with numpy + pandas only.
"""

import importlib.util
import os
import sys
import types

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Stubs for heavy imports
# ---------------------------------------------------------------------------


def _stub(name):
    mod = sys.modules.get(name)
    if mod is None:
        mod = types.ModuleType(name)
        sys.modules[name] = mod
    return mod


# parmed
_stub("parmed")

# h5py
_stub("h5py")


# scipy.special – use real scipy if present, otherwise a numpy fallback
def _numpy_logsumexp(a, b=None, **_kw):
    """Minimal logsumexp equivalent for environments without scipy."""
    a = np.asarray(a, dtype=float)
    a_max = np.max(a)
    if b is not None:
        b = np.asarray(b, dtype=float)
        return a_max + np.log(np.sum(b * np.exp(a - a_max)))
    return a_max + np.log(np.sum(np.exp(a - a_max)))


try:
    from scipy.special import logsumexp as _logsumexp  # noqa: F401
except ImportError:
    _scipy = _stub("scipy")
    _scipy_special = _stub("scipy.special")
    _scipy.special = _scipy_special
    _scipy_special.logsumexp = _numpy_logsumexp

# multiego package
_multiego = _stub("multiego")

# type_definitions – gromos_atp is accessed at module level
_td = _stub("multiego.type_definitions")
_td.gromos_atp = types.SimpleNamespace(
    name=["OMet", "N"],
    rc_c12=[1.0e-6, 2.0e-6],
)
_td.special_non_local = []
_td.aminoacids_list = []
_td.nucleic_acid_list = []
_td.atom_type_combinations = []
_td.from_ff_to_multiego = {}
_multiego.type_definitions = _td

_bonded = _stub("multiego.bonded")
_multiego.bonded = _bonded

_io_stub = _stub("multiego.io")
_multiego.io = _io_stub

# ---------------------------------------------------------------------------
# Import the module under test via file path
# ---------------------------------------------------------------------------

_SCRIPT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "tools", "make_mat", "make_mat.py"))
_spec = importlib.util.spec_from_file_location("make_mat_mod", _SCRIPT)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

create_matrix_mask = _mod.create_matrix_mask
get_col_params = _mod.get_col_params
calculate_probability = _mod.calculate_probability
get_cumulative_probability = _mod.get_cumulative_probability
c12_avg = _mod.c12_avg
generate_c12_values = _mod.generate_c12_values
CUTOFF_FACTOR = _mod.CUTOFF_FACTOR


# ===========================================================================
# TestConstants
# ===========================================================================


class TestConstants:
    def test_cutoff_factor_value(self):
        assert CUTOFF_FACTOR == pytest.approx(1.45)


# ===========================================================================
# TestCreateMatrixMask
# ===========================================================================


class TestCreateMatrixMask:

    def test_single_pair_match(self):
        s1 = np.array(["A", "B", "C"])
        s2 = np.array(["X", "Y"])
        mask = create_matrix_mask(s1, s2, [("A", "X")])
        assert mask[0, 0] is np.bool_(True)
        assert not mask[1, 0]
        assert not mask[0, 1]

    def test_no_match(self):
        s1 = np.array(["A", "B"])
        s2 = np.array(["X", "Y"])
        mask = create_matrix_mask(s1, s2, [("Z", "W")])
        assert not mask.any()

    def test_symmetrize_adds_reverse_pair(self):
        s1 = np.array(["A", "B"])
        s2 = np.array(["A", "B"])
        mask = create_matrix_mask(s1, s2, [("A", "B")], symmetrize=True)
        assert mask[0, 1]  # (A,B)
        assert mask[1, 0]  # (B,A) added by symmetrize

    def test_multiple_pairs(self):
        s1 = np.array(["A", "B", "C"])
        s2 = np.array(["X", "Y", "Z"])
        mask = create_matrix_mask(s1, s2, [("A", "X"), ("C", "Z")])
        assert mask[0, 0]
        assert mask[2, 2]
        assert not mask[1, 1]

    def test_output_shape(self):
        s1 = np.array(["A"] * 5)
        s2 = np.array(["B"] * 3)
        mask = create_matrix_mask(s1, s2, [])
        assert mask.shape == (5, 3)
        assert not mask.any()

    def test_empty_types_list(self):
        s1 = np.array(["A"])
        s2 = np.array(["A"])
        mask = create_matrix_mask(s1, s2, [])
        assert not mask.any()


# ===========================================================================
# TestGetColParams
# ===========================================================================


def _make_histogram(n=50, cutoff=0.5, x_min=0.1, x_max=1.0):
    """Return (values_with_sentinel, weights_with_cutoff) for a flat histogram."""
    v = np.linspace(x_min, x_max, n)
    w = np.ones(n) * 0.02
    # Append sentinel: last value of v unused; last weight = cutoff
    return np.append(v, 0.0), np.append(w, cutoff)


class TestGetColParams:

    def test_returns_five_values(self):
        v, w = _make_histogram()
        result = get_col_params(v, w)
        assert len(result) == 5

    def test_cutoff_extracted_from_last_weight(self):
        v, w = _make_histogram(cutoff=0.6)
        cutoff, i, norm, vt, wt = get_col_params(v, w)
        assert cutoff == pytest.approx(0.6)

    def test_bins_truncated_at_cutoff(self):
        v, w = _make_histogram(cutoff=0.5)
        cutoff, i, norm, vt, wt = get_col_params(v, w)
        assert np.all(vt <= cutoff + 1e-12)

    def test_norm_equals_sum_of_weights(self):
        v, w = _make_histogram(n=20, cutoff=0.6)
        cutoff, i, norm, vt, wt = get_col_params(v, w)
        assert norm == pytest.approx(np.sum(wt))

    def test_empty_when_all_bins_above_cutoff(self):
        # All bins at 0.8–1.0, cutoff = 0.3 → no valid bins
        v = np.array([0.8, 0.9, 1.0, 0.0])
        w = np.array([1.0, 1.0, 1.0, 0.3])
        result = get_col_params(v, w)
        assert result == (0, 0, 0, 0, 0)

    def test_single_valid_bin(self):
        v = np.array([0.2, 0.8, 0.0])
        w = np.array([5.0, 0.0, 0.5])  # cutoff=0.5; only bin at 0.2 qualifies
        cutoff, i, norm, vt, wt = get_col_params(v, w)
        assert len(vt) == 1
        assert vt[0] == pytest.approx(0.2)


# ===========================================================================
# TestCalculateProbability
# ===========================================================================


class TestCalculateProbability:

    def test_zero_weights_gives_zero(self):
        v, w = _make_histogram()
        w[:-1] = 0.0
        assert calculate_probability(v, w) == pytest.approx(0.0)

    def test_large_weights_capped_at_one(self):
        v = np.linspace(0.1, 0.8, 50)
        w = np.full(50, 100.0)
        v_ext = np.append(v, 0.0)
        w_ext = np.append(w, 0.8)
        assert calculate_probability(v_ext, w_ext) == pytest.approx(1.0)

    def test_result_between_zero_and_one(self):
        v, w = _make_histogram(cutoff=0.6)
        p = calculate_probability(v, w)
        assert 0.0 <= p <= 1.0

    def test_higher_weights_give_higher_probability(self):
        v_lo, w_lo = _make_histogram(cutoff=0.6)
        v_hi, w_hi = _make_histogram(cutoff=0.6)
        w_hi[:-1] *= 5.0
        p_lo = calculate_probability(v_lo, w_lo)
        p_hi = calculate_probability(v_hi, w_hi)
        assert p_hi >= p_lo


# ===========================================================================
# TestGetCumulativeProbability
# ===========================================================================


class TestGetCumulativeProbability:

    def test_returns_weight_at_cutoff_index(self):
        v = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.0])
        w = np.array([0.5, 0.7, 0.9, 0.3, 0.1, 0.35])  # cutoff = 0.35
        _, i, _, _, _ = get_col_params(v, w)
        expected = w[i]
        result = get_cumulative_probability(v, w)
        assert result == pytest.approx(expected)

    def test_empty_histogram_does_not_raise(self):
        v = np.array([0.8, 0.9, 0.0])
        w = np.array([1.0, 1.0, 0.3])  # cutoff=0.3, no valid bins
        # get_col_params returns (0,0,0,0,0) → result = w[0]
        result = get_cumulative_probability(v, w)
        assert result is not None


# ===========================================================================
# TestC12Avg
# ===========================================================================


class TestC12Avg:

    def test_zero_weights_returns_zero(self):
        v = np.linspace(0.1, 0.8, 50)
        w = np.zeros(50)
        v_ext = np.append(v, 0.0)
        w_ext = np.append(w, 0.8)
        assert c12_avg(v_ext, w_ext) == pytest.approx(0.0)

    def test_returns_positive_for_peaked_histogram(self):
        v = np.linspace(0.2, 0.6, 40)
        w = np.exp(-((v - 0.35) ** 2) / 0.01)
        v_ext = np.append(v, 0.0)
        w_ext = np.append(w, 0.6)
        assert c12_avg(v_ext, w_ext) > 0.0

    def test_result_at_most_cutoff(self):
        cutoff = 0.75
        v = np.linspace(0.1, cutoff, 60)
        w = np.ones(60)
        v_ext = np.append(v, 0.0)
        w_ext = np.append(w, cutoff)
        result = c12_avg(v_ext, w_ext)
        assert 0.0 <= result <= cutoff

    def test_peak_closer_gives_smaller_average(self):
        # A histogram peaked near the origin should give a smaller c12 average
        # (c12_avg is a c12-exp average, so closer distances get higher weight)
        v = np.linspace(0.1, 0.8, 100)
        w_near = np.exp(-((v - 0.15) ** 2) / 0.001)
        w_far = np.exp(-((v - 0.70) ** 2) / 0.001)
        v_ext = np.append(v, 0.0)
        r_near = c12_avg(v_ext, np.append(w_near, 0.8))
        r_far = c12_avg(v_ext, np.append(w_far, 0.8))
        assert r_near < r_far


# ===========================================================================
# TestGenerateC12Values
# ===========================================================================


def _df(n, c12_values=None):
    if c12_values is None:
        c12_values = np.full(n, 1.0e-6)
    return pd.DataFrame({"c12": c12_values, "resnum": np.arange(n)})


class TestGenerateC12Values:

    def test_output_shape(self):
        result = generate_c12_values(_df(5), {}, [], "other")
        assert result.shape == (5, 5)

    def test_geometric_mean_identical_values(self):
        c12 = 4.0e-6
        result = generate_c12_values(_df(3, np.full(3, c12)), {}, [], "other")
        assert result[0, 0] == pytest.approx(c12)

    def test_geometric_mean_different_values(self):
        df = _df(2, np.array([1.0e-6, 4.0e-6]))
        result = generate_c12_values(df, {}, [], "other")
        expected = np.sqrt(1.0e-6 * 4.0e-6)
        assert result[0, 1] == pytest.approx(expected)
        assert result[1, 0] == pytest.approx(expected)

    def test_no_protein_rules_for_other(self):
        """Rules in combinations are skipped when molecule_type='other'."""
        df = _df(3, np.full(3, 2.0e-6))
        types = {"BB": np.array([True, False, False])}
        combos = [("BB", "BB", 0.5, None, 0)]
        r_other = generate_c12_values(df, types, combos, "other")
        r_protein = generate_c12_values(df, types, combos, "protein")
        assert r_other[0, 0] == pytest.approx(2.0e-6)
        assert r_protein[0, 0] == pytest.approx(0.5 * 2.0e-6)

    def test_constant_rule_overrides(self):
        # shift=0 matches only pairs where resnum_i == resnum_j (diagonal).
        # Off-diagonal entries keep the geometric mean.
        df = _df(2, np.full(2, 1.0e-6))
        types = {"X": np.array([True, True])}
        combos = [("X", "X", None, 9.0e-7, 0)]
        result = generate_c12_values(df, types, combos, "protein")
        assert result[0, 0] == pytest.approx(9.0e-7)  # same residue → overridden
        assert result[1, 1] == pytest.approx(9.0e-7)  # same residue → overridden
        assert result[0, 1] == pytest.approx(1.0e-6)  # different residues → geometric mean

    def test_factor_and_constant_uses_minimum(self):
        df = _df(2, np.full(2, 1.0e-6))
        types = {"X": np.array([True, True])}
        # factor=2 → 2e-6, but constant=5e-7 caps it → result = 5e-7
        combos = [("X", "X", 2.0, 5.0e-7, 0)]
        result = generate_c12_values(df, types, combos, "protein")
        assert result[0, 0] == pytest.approx(5.0e-7)

    def test_rule_missing_both_factor_and_constant_raises(self):
        df = _df(2)
        types = {"X": np.array([True, True])}
        combos = [("X", "X", None, None, 0)]
        with pytest.raises(ValueError, match="factor or constant"):
            generate_c12_values(df, types, combos, "protein")

    def test_residue_shift_applied(self):
        """A shift of 1 means only pairs where resnum_j == resnum_i + 1 are matched."""
        df = _df(4, np.full(4, 1.0e-6))
        types = {"X": np.array([True, True, True, True])}
        combos = [("X", "X", 2.0, None, 1)]
        result = generate_c12_values(df, types, combos, "protein")
        # Only pairs (i,j) where resnum_j == resnum_i + 1 should be scaled
        assert result[0, 1] == pytest.approx(2.0e-6)  # 0→1 scaled
        assert result[0, 2] == pytest.approx(1.0e-6)  # 0→2 not scaled (shift mismatch)
