"""
Unit tests for contacts._symmetrize_reference_contacts.

The function averages rc_probability (simple mean) and rc_distance
(probability-weighted mean) for symmetry-equivalent atom pairs within
the same residue of a reference contact matrix.

sb_type format: atomname_residuename_resnum  (e.g. "OD1_ASP_5")
Contact matrix columns used: rc_ai, rc_aj, rc_probability, rc_distance,
rc_same_chain.
"""

import contextlib
import sys
import types

import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def symmetrize_fn(stub_deps):
    """
    Import _symmetrize_reference_contacts with all heavy deps stubbed.
    The shared conftest stubs most of the multiego package; we add
    _parmed_compat here since it is not in the shared fixture.
    """
    pmc = types.ModuleType("multiego._parmed_compat")
    pmc.gmxlib_in_topo_dir = contextlib.nullcontext
    sys.modules.setdefault("multiego._parmed_compat", pmc)

    from multiego.contacts import _symmetrize_reference_contacts

    return _symmetrize_reference_contacts


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _topology(sbtype_resname, molecule_name="MOL1"):
    """
    Build a minimal topology DataFrame.

    Parameters
    ----------
    sbtype_resname : dict
        Mapping {sb_type: resname}.  The residue number is parsed from the
        last ``_``-delimited segment of ``sb_type``.
    molecule_name : str
        Value written into the ``molecule_name`` column.
    """
    rows = [
        {
            "sb_type": sbt,
            "resname": rn,
            "molecule_name": molecule_name,
            "resnum": sbt.rsplit("_", 1)[-1],
        }
        for sbt, rn in sbtype_resname.items()
    ]
    return pd.DataFrame(rows)


def _matrix(rows, same_chain=True):
    """
    Build a reference contact matrix from a list of tuples.

    Parameters
    ----------
    rows : list of (rc_ai, rc_aj, rc_probability, rc_distance)
    same_chain : bool or list of bool
        If a scalar, applied to all rows.
    """
    if isinstance(same_chain, bool):
        same_chain = [same_chain] * len(rows)
    return pd.DataFrame(
        [
            {
                "rc_ai": ai,
                "rc_aj": aj,
                "rc_probability": p,
                "rc_distance": d,
                "rc_same_chain": sc,
            }
            for (ai, aj, p, d), sc in zip(rows, same_chain)
        ]
    )


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_empty_symmetry_returns_unchanged(self, symmetrize_fn):
        """No symmetry rules → matrix is returned as-is."""
        topo = _topology({"OD1_ASP_5": "ASP", "OD2_ASP_5": "ASP"})
        cm = _matrix([("X_ALA_1", "OD1_ASP_5", 0.8, 0.40), ("X_ALA_1", "OD2_ASP_5", 0.4, 0.50)])
        original_p = cm["rc_probability"].tolist()
        original_d = cm["rc_distance"].tolist()

        result = symmetrize_fn(cm, [], topo)

        assert result["rc_probability"].tolist() == original_p
        assert result["rc_distance"].tolist() == original_d

    def test_empty_matrix_returns_unchanged(self, symmetrize_fn):
        """Empty contact matrix → returns the same empty DataFrame."""
        topo = _topology({"OD1_ASP_5": "ASP"})
        cm = _matrix([])

        result = symmetrize_fn(cm, [["ASP", "OD1", "OD2"]], topo)

        assert result.empty

    def test_no_matching_atoms_returns_unchanged(self, symmetrize_fn):
        """Symmetry rule for ASP, but matrix only contains GLY → no change."""
        topo = _topology({"CA_GLY_1": "GLY", "CB_ALA_2": "ALA"})
        cm = _matrix([("CA_GLY_1", "CB_ALA_2", 0.8, 0.40)])
        symmetry = [["ASP", "OD1", "OD2"]]

        result = symmetrize_fn(cm, symmetry, topo)

        assert result["rc_probability"].iloc[0] == pytest.approx(0.8)
        assert result["rc_distance"].iloc[0] == pytest.approx(0.40)


# ---------------------------------------------------------------------------
# Probability averaging
# ---------------------------------------------------------------------------


class TestProbabilityAveraging:
    def test_simple_average_on_aj(self, symmetrize_fn):
        """Two equivalent contacts on the aj side get the mean probability."""
        topo = _topology({"X_ALA_1": "ALA", "OD1_ASP_5": "ASP", "OD2_ASP_5": "ASP"})
        cm = _matrix(
            [
                ("X_ALA_1", "OD1_ASP_5", 0.8, 0.40),
                ("X_ALA_1", "OD2_ASP_5", 0.4, 0.50),
            ]
        )
        result = symmetrize_fn(cm, [["ASP", "OD1", "OD2"]], topo)

        assert result["rc_probability"].tolist() == pytest.approx([0.6, 0.6])

    def test_simple_average_on_ai(self, symmetrize_fn):
        """Two equivalent contacts on the ai side get the mean probability."""
        topo = _topology({"OD1_ASP_5": "ASP", "OD2_ASP_5": "ASP", "X_ALA_1": "ALA"})
        cm = _matrix(
            [
                ("OD1_ASP_5", "X_ALA_1", 0.8, 0.40),
                ("OD2_ASP_5", "X_ALA_1", 0.4, 0.50),
            ]
        )
        result = symmetrize_fn(cm, [["ASP", "OD1", "OD2"]], topo)

        assert result["rc_probability"].tolist() == pytest.approx([0.6, 0.6])

    def test_average_on_both_sides(self, symmetrize_fn):
        """Symmetry on both ai and aj groups all four contacts together."""
        topo = _topology(
            {
                "OD1_ASP_5": "ASP",
                "OD2_ASP_5": "ASP",
                "NH1_ARG_10": "ARG",
                "NH2_ARG_10": "ARG",
            }
        )
        cm = _matrix(
            [
                ("OD1_ASP_5", "NH1_ARG_10", 0.8, 0.40),
                ("OD1_ASP_5", "NH2_ARG_10", 0.6, 0.50),
                ("OD2_ASP_5", "NH1_ARG_10", 0.4, 0.60),
                ("OD2_ASP_5", "NH2_ARG_10", 0.2, 0.70),
            ]
        )
        symmetry = [["ASP", "OD1", "OD2"], ["ARG", "NH1", "NH2"]]
        result = symmetrize_fn(cm, symmetry, topo)

        expected_p = (0.8 + 0.6 + 0.4 + 0.2) / 4  # 0.5
        assert result["rc_probability"].tolist() == pytest.approx([expected_p] * 4)

    def test_non_symmetric_contacts_unchanged(self, symmetrize_fn):
        """Contacts that don't involve any symmetric atom are not touched."""
        topo = _topology(
            {
                "OD1_ASP_5": "ASP",
                "OD2_ASP_5": "ASP",
                "CA_ALA_1": "ALA",
                "CB_ALA_2": "ALA",
            }
        )
        cm = _matrix(
            [
                ("CA_ALA_1", "OD1_ASP_5", 0.8, 0.40),  # symmetric
                ("CA_ALA_1", "OD2_ASP_5", 0.4, 0.50),  # symmetric
                ("CA_ALA_1", "CB_ALA_2", 0.9, 0.35),  # not symmetric → unchanged
            ]
        )
        result = symmetrize_fn(cm, [["ASP", "OD1", "OD2"]], topo)

        # Third row must not be affected
        assert result.iloc[2]["rc_probability"] == pytest.approx(0.9)
        assert result.iloc[2]["rc_distance"] == pytest.approx(0.35)


# ---------------------------------------------------------------------------
# Distance weighting
# ---------------------------------------------------------------------------


class TestDistanceWeighting:
    def test_distance_weighted_by_original_probability(self, symmetrize_fn):
        """
        rc_distance is a probability-weighted mean using the *original*
        probabilities, i.e.  Σ(p_i * d_i) / Σ(p_i).
        """
        topo = _topology({"X_ALA_1": "ALA", "OD1_ASP_5": "ASP", "OD2_ASP_5": "ASP"})
        cm = _matrix(
            [
                ("X_ALA_1", "OD1_ASP_5", 0.8, 0.40),
                ("X_ALA_1", "OD2_ASP_5", 0.6, 0.50),
            ]
        )
        result = symmetrize_fn(cm, [["ASP", "OD1", "OD2"]], topo)

        expected_d = (0.8 * 0.40 + 0.6 * 0.50) / (0.8 + 0.6)
        assert result["rc_distance"].tolist() == pytest.approx([expected_d, expected_d])

    def test_distance_not_simple_mean(self, symmetrize_fn):
        """Confirm that the result differs from a naive simple mean when
        probabilities are unequal — i.e. weighting is actually applied."""
        topo = _topology({"X_ALA_1": "ALA", "OD1_ASP_5": "ASP", "OD2_ASP_5": "ASP"})
        # Very unequal probabilities so simple vs weighted means differ clearly
        cm = _matrix(
            [
                ("X_ALA_1", "OD1_ASP_5", 0.9, 0.30),
                ("X_ALA_1", "OD2_ASP_5", 0.1, 0.70),
            ]
        )
        result = symmetrize_fn(cm, [["ASP", "OD1", "OD2"]], topo)

        simple_mean = (0.30 + 0.70) / 2  # 0.50
        weighted_mean = (0.9 * 0.30 + 0.1 * 0.70) / (0.9 + 0.1)  # 0.34

        assert result["rc_distance"].iloc[0] == pytest.approx(weighted_mean)
        assert result["rc_distance"].iloc[0] != pytest.approx(simple_mean, abs=1e-3)

    def test_zero_probability_fallback_to_simple_mean(self, symmetrize_fn):
        """When all contacts in a group have probability 0, falls back to
        the simple mean distance to avoid division by zero."""
        topo = _topology({"X_ALA_1": "ALA", "OD1_ASP_5": "ASP", "OD2_ASP_5": "ASP"})
        cm = _matrix(
            [
                ("X_ALA_1", "OD1_ASP_5", 0.0, 0.40),
                ("X_ALA_1", "OD2_ASP_5", 0.0, 0.60),
            ]
        )
        result = symmetrize_fn(cm, [["ASP", "OD1", "OD2"]], topo)

        assert result["rc_distance"].tolist() == pytest.approx([0.50, 0.50])


# ---------------------------------------------------------------------------
# Grouping behaviour
# ---------------------------------------------------------------------------


class TestGrouping:
    def test_independent_residues_averaged_separately(self, symmetrize_fn):
        """Two ASP residues with different numbers are averaged independently."""
        topo = _topology(
            {
                "X_ALA_1": "ALA",
                "OD1_ASP_5": "ASP",
                "OD2_ASP_5": "ASP",
                "OD1_ASP_9": "ASP",
                "OD2_ASP_9": "ASP",
            }
        )
        cm = _matrix(
            [
                ("X_ALA_1", "OD1_ASP_5", 0.8, 0.40),  # group A
                ("X_ALA_1", "OD2_ASP_5", 0.4, 0.50),  # group A
                ("X_ALA_1", "OD1_ASP_9", 0.6, 0.35),  # group B
                ("X_ALA_1", "OD2_ASP_9", 0.2, 0.45),  # group B
            ]
        )
        result = symmetrize_fn(cm, [["ASP", "OD1", "OD2"]], topo)

        expected_p_A = (0.8 + 0.4) / 2
        expected_d_A = (0.8 * 0.40 + 0.4 * 0.50) / (0.8 + 0.4)
        expected_p_B = (0.6 + 0.2) / 2
        expected_d_B = (0.6 * 0.35 + 0.2 * 0.45) / (0.6 + 0.2)

        assert result.iloc[0]["rc_probability"] == pytest.approx(expected_p_A)
        assert result.iloc[1]["rc_probability"] == pytest.approx(expected_p_A)
        assert result.iloc[0]["rc_distance"] == pytest.approx(expected_d_A)
        assert result.iloc[2]["rc_probability"] == pytest.approx(expected_p_B)
        assert result.iloc[3]["rc_probability"] == pytest.approx(expected_p_B)
        assert result.iloc[2]["rc_distance"] == pytest.approx(expected_d_B)
        assert result.iloc[3]["rc_distance"] == pytest.approx(expected_d_B)

    def test_same_chain_is_grouping_key(self, symmetrize_fn):
        """Contacts with the same atom pair but different same_chain values
        are averaged within their own group, not across groups."""
        topo = _topology({"X_ALA_1": "ALA", "OD1_ASP_5": "ASP", "OD2_ASP_5": "ASP"})
        cm = _matrix(
            [
                ("X_ALA_1", "OD1_ASP_5", 0.8, 0.40),  # intra
                ("X_ALA_1", "OD2_ASP_5", 0.4, 0.50),  # intra
                ("X_ALA_1", "OD1_ASP_5", 0.6, 0.30),  # inter
                ("X_ALA_1", "OD2_ASP_5", 0.2, 0.60),  # inter
            ],
            same_chain=[True, True, False, False],
        )
        result = symmetrize_fn(cm, [["ASP", "OD1", "OD2"]], topo)

        # intra group
        expected_p_intra = (0.8 + 0.4) / 2
        assert result.iloc[0]["rc_probability"] == pytest.approx(expected_p_intra)
        assert result.iloc[1]["rc_probability"] == pytest.approx(expected_p_intra)

        # inter group — different result
        expected_p_inter = (0.6 + 0.2) / 2
        assert result.iloc[2]["rc_probability"] == pytest.approx(expected_p_inter)
        assert result.iloc[3]["rc_probability"] == pytest.approx(expected_p_inter)

    def test_three_equivalent_atoms(self, symmetrize_fn):
        """A symmetry group with three atoms averages all three together."""
        topo = _topology(
            {
                "X_ALA_1": "ALA",
                "CD1_PHE_7": "PHE",
                "CD2_PHE_7": "PHE",
                "CE1_PHE_7": "PHE",
            }
        )
        cm = _matrix(
            [
                ("X_ALA_1", "CD1_PHE_7", 0.9, 0.35),
                ("X_ALA_1", "CD2_PHE_7", 0.6, 0.40),
                ("X_ALA_1", "CE1_PHE_7", 0.3, 0.50),
            ]
        )
        symmetry = [["PHE", "CD1", "CD2", "CE1"]]
        result = symmetrize_fn(cm, symmetry, topo)

        expected_p = (0.9 + 0.6 + 0.3) / 3
        expected_d = (0.9 * 0.35 + 0.6 * 0.40 + 0.3 * 0.50) / (0.9 + 0.6 + 0.3)
        assert result["rc_probability"].tolist() == pytest.approx([expected_p] * 3)
        assert result["rc_distance"].tolist() == pytest.approx([expected_d] * 3)

    def test_six_equivalent_atoms(self, symmetrize_fn):
        """A symmetry group with six atoms (e.g. lysine side-chain hydrogens)
        averages all six contacts together for both probability and distance."""
        topo = _topology(
            {
                "X_ALA_1": "ALA",
                "HZ1_LYS_12": "LYS",
                "HZ2_LYS_12": "LYS",
                "HZ3_LYS_12": "LYS",
                "HE1_LYS_12": "LYS",
                "HE2_LYS_12": "LYS",
                "HE3_LYS_12": "LYS",
            }
        )
        probs = [0.9, 0.7, 0.5, 0.4, 0.2, 0.1]
        dists = [0.30, 0.35, 0.40, 0.45, 0.55, 0.60]
        cm = _matrix(
            [
                ("X_ALA_1", "HZ1_LYS_12", probs[0], dists[0]),
                ("X_ALA_1", "HZ2_LYS_12", probs[1], dists[1]),
                ("X_ALA_1", "HZ3_LYS_12", probs[2], dists[2]),
                ("X_ALA_1", "HE1_LYS_12", probs[3], dists[3]),
                ("X_ALA_1", "HE2_LYS_12", probs[4], dists[4]),
                ("X_ALA_1", "HE3_LYS_12", probs[5], dists[5]),
            ]
        )
        symmetry = [["LYS", "HZ1", "HZ2", "HZ3", "HE1", "HE2", "HE3"]]
        result = symmetrize_fn(cm, symmetry, topo)

        expected_p = sum(probs) / 6
        sum_p = sum(probs)
        expected_d = sum(p * d for p, d in zip(probs, dists)) / sum_p
        assert result["rc_probability"].tolist() == pytest.approx([expected_p] * 6)
        assert result["rc_distance"].tolist() == pytest.approx([expected_d] * 6)


# ---------------------------------------------------------------------------
# Terminal-residue keywords
# ---------------------------------------------------------------------------


class TestTerminalKeywords:
    def _cter_topology(self):
        """Single molecule: residues 1 (non-terminal) and 2 (C-terminal)."""
        return _topology(
            {
                "X_ALA_1": "ALA",
                "O1_ALA_1": "ALA",  # not CTER — resnum 1 in a 2-residue chain
                "O2_ALA_1": "ALA",
                "O1_ALA_2": "ALA",  # CTER — last residue
                "O2_ALA_2": "ALA",
            }
        )

    def test_cter_only_affects_last_residue(self, symmetrize_fn):
        """CTER symmetry should average only the contacts at residue 2, not 1."""
        topo = self._cter_topology()
        cm = _matrix(
            [
                ("X_ALA_1", "O1_ALA_1", 0.8, 0.40),  # non-terminal → unchanged
                ("X_ALA_1", "O2_ALA_1", 0.4, 0.50),  # non-terminal → unchanged
                ("X_ALA_1", "O1_ALA_2", 0.8, 0.40),  # CTER → averaged
                ("X_ALA_1", "O2_ALA_2", 0.4, 0.50),  # CTER → averaged
            ]
        )
        result = symmetrize_fn(cm, [["CTER", "O1", "O2"]], topo)

        # Non-terminal residue: untouched
        assert result.iloc[0]["rc_probability"] == pytest.approx(0.8)
        assert result.iloc[1]["rc_probability"] == pytest.approx(0.4)

        # C-terminal residue: averaged
        expected_p = (0.8 + 0.4) / 2
        assert result.iloc[2]["rc_probability"] == pytest.approx(expected_p)
        assert result.iloc[3]["rc_probability"] == pytest.approx(expected_p)
