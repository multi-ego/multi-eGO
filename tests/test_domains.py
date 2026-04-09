"""
Unit tests for tools/domain_sectioner/domains.py.

Covers:
- dom_range: parsing, validation (decreasing / overlapping ranges)
- find_atom_start / find_atom_end: correct 0-based atom index lookup
- build_domain_mask: correct outer-product masking and inversion
"""

import importlib.util
import types
from pathlib import Path

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Import the module under test
# ---------------------------------------------------------------------------

_SCRIPT = Path(__file__).parent.parent / "tools" / "domain_sectioner" / "domains.py"

# domains.py imports parmed at the top level; stub it so we don't need it
# installed in the test environment.
_parmed_stub = types.ModuleType("parmed")
_parmed_stub.load_file = None

import sys  # noqa: E402

sys.modules.setdefault("parmed", _parmed_stub)

_spec = importlib.util.spec_from_file_location("domains", _SCRIPT)
_domains = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_domains)

dom_range = _domains.dom_range
find_atom_start = _domains.find_atom_start
find_atom_end = _domains.find_atom_end
build_domain_mask = _domains.build_domain_mask


# ---------------------------------------------------------------------------
# Mock topology helpers
# ---------------------------------------------------------------------------


def _make_residue(n_atoms):
    """Return a minimal mock residue with *n_atoms* dummy atoms."""
    res = types.SimpleNamespace()
    res.atoms = [types.SimpleNamespace()] * n_atoms
    return res


def _make_topology(*atoms_per_residue):
    """Return a mock topology whose residues have the given atom counts.

    Example: _make_topology(3, 5, 2) → three residues with 3, 5, 2 atoms.
    """
    top = types.SimpleNamespace()
    top.residues = [_make_residue(n) for n in atoms_per_residue]
    # Build a flat atoms list so build_domain_mask can access topology.atoms[i]
    atom_list = []
    for res_idx, n in enumerate(atoms_per_residue):
        for atom_idx in range(n):
            atom_list.append(types.SimpleNamespace(__str__=lambda self, r=res_idx, a=atom_idx: f"RES{r+1}:atom{a}"))
    top.atoms = atom_list
    return top


# ---------------------------------------------------------------------------
# dom_range
# ---------------------------------------------------------------------------


class TestDomRange:
    def test_single_valid_range(self):
        assert dom_range(["1-30"]) == [(1, 30)]

    def test_multiple_valid_ranges(self):
        assert dom_range(["1-30", "45-60"]) == [(1, 30), (45, 60)]

    def test_single_residue_range(self):
        assert dom_range(["4-4"]) == [(4, 4)]

    def test_decreasing_range_raises(self):
        with pytest.raises(ValueError, match="start must be"):
            dom_range(["30-20"])

    def test_overlapping_ranges_raises(self):
        with pytest.raises(ValueError, match="overlap"):
            dom_range(["1-30", "25-60"])

    def test_touching_ranges_raise(self):
        # end of first == start of second → overlap
        with pytest.raises(ValueError, match="overlap"):
            dom_range(["1-30", "30-60"])

    def test_adjacent_ranges_are_valid(self):
        # end of first < start of second → ok
        result = dom_range(["1-30", "31-60"])
        assert result == [(1, 30), (31, 60)]

    def test_out_of_order_ranges_raise(self):
        with pytest.raises(ValueError, match="overlap"):
            dom_range(["45-60", "1-30"])


# ---------------------------------------------------------------------------
# find_atom_start / find_atom_end
# ---------------------------------------------------------------------------


class TestFindAtomBoundaries:
    """Topology: residues 1-4 with 3, 5, 2, 4 atoms → cumulative 0,3,8,10,14."""

    @pytest.fixture(autouse=True)
    def topology(self):
        self.top = _make_topology(3, 5, 2, 4)

    def test_start_first_residue(self):
        assert find_atom_start(self.top, 1) == 0

    def test_start_second_residue(self):
        assert find_atom_start(self.top, 2) == 3

    def test_start_third_residue(self):
        assert find_atom_start(self.top, 3) == 8

    def test_start_fourth_residue(self):
        assert find_atom_start(self.top, 4) == 10

    def test_end_first_residue(self):
        assert find_atom_end(self.top, 1) == 2  # atoms 0,1,2

    def test_end_second_residue(self):
        assert find_atom_end(self.top, 2) == 7  # atoms 3..7

    def test_end_third_residue(self):
        assert find_atom_end(self.top, 3) == 9  # atoms 8,9

    def test_end_fourth_residue(self):
        assert find_atom_end(self.top, 4) == 13  # atoms 10..13

    def test_single_atom_residue_start_equals_end(self):
        top = _make_topology(1)
        assert find_atom_start(top, 1) == find_atom_end(top, 1) == 0

    def test_start_end_consistency(self):
        """find_atom_end - find_atom_start + 1 == atoms in that residue."""
        atoms_per_res = (3, 5, 2, 4)
        top = _make_topology(*atoms_per_res)
        for res_num, n in enumerate(atoms_per_res, start=1):
            span = find_atom_end(top, res_num) - find_atom_start(top, res_num) + 1
            assert span == n


# ---------------------------------------------------------------------------
# build_domain_mask
# ---------------------------------------------------------------------------


class TestBuildDomainMask:
    """Topology: 3 residues, 2 atoms each → 6 atoms total, 36-element mask."""

    @pytest.fixture(autouse=True)
    def setup(self):
        # residues 1,2,3 each with 2 atoms → atoms 0-1, 2-3, 4-5
        self.top = _make_topology(2, 2, 2)
        self.n = 6

    def test_full_range_all_true(self):
        mask = build_domain_mask(self.top, self.n, [(1, 3)])
        assert mask.all()

    def test_single_residue_mask_shape(self):
        mask = build_domain_mask(self.top, self.n, [(1, 1)])
        assert mask.shape == (self.n**2,)

    def test_single_first_residue(self):
        # Only atoms 0 and 1 are in the domain; only (0,0),(0,1),(1,0),(1,1)
        # should be True in the 6×6 grid.
        mask = build_domain_mask(self.top, self.n, [(1, 1)])
        mask_2d = mask.reshape(self.n, self.n)
        # domain atoms: {0, 1}
        for i in range(self.n):
            for j in range(self.n):
                expected = (i in (0, 1)) and (j in (0, 1))
                assert mask_2d[i, j] == expected, f"Unexpected value at ({i},{j})"

    def test_two_ranges(self):
        # residues 1 (atoms 0-1) and 3 (atoms 4-5)
        mask = build_domain_mask(self.top, self.n, [(1, 1), (3, 3)])
        mask_2d = mask.reshape(self.n, self.n)
        domain_atoms = {0, 1, 4, 5}
        for i in range(self.n):
            for j in range(self.n):
                expected = (i in domain_atoms) and (j in domain_atoms)
                assert mask_2d[i, j] == expected

    def test_invert(self):
        mask = build_domain_mask(self.top, self.n, [(1, 1)])
        mask_inv = build_domain_mask(self.top, self.n, [(1, 1)], invert=True)
        assert np.all(mask == ~mask_inv)

    def test_invert_full_range_all_false(self):
        mask = build_domain_mask(self.top, self.n, [(1, 3)], invert=True)
        assert not mask.any()

    def test_mask_is_symmetric(self):
        """Contact mask must be symmetric: if (i,j) is in domain, so is (j,i)."""
        mask = build_domain_mask(self.top, self.n, [(1, 2)])
        mask_2d = mask.reshape(self.n, self.n)
        assert np.array_equal(mask_2d, mask_2d.T)

    def test_number_of_true_entries(self):
        """Domain of k atoms → k² True entries in the mask."""
        # residue 2 has 2 atoms
        mask = build_domain_mask(self.top, self.n, [(2, 2)])
        assert mask.sum() == 2**2

    def test_two_residue_domain_count(self):
        # residues 1+2 = 4 atoms → 16 True entries
        mask = build_domain_mask(self.top, self.n, [(1, 2)])
        assert mask.sum() == 4**2
