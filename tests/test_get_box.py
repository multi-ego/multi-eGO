"""
Unit tests for tools/box_concentration/get_box.py.

Covers:
- Pure math functions (box_from_n_mol_and_conc, conc_from_n_mol_and_volume,
  n_mol_from_conc_and_volume, sphere_volume)
- Roundtrip consistency between the three conversion functions
- CLI behaviour: correct output, error handling, geometry-only mode
"""

import math
import subprocess
import sys
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Import the module under test (it lives outside the src/ tree)
# ---------------------------------------------------------------------------

_SCRIPT = Path(__file__).parent.parent / "tools" / "box_concentration" / "get_box.py"

import importlib.util

_spec = importlib.util.spec_from_file_location("get_box", _SCRIPT)
_get_box = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_get_box)

AVOGADRO = _get_box.AVOGADRO
NM3_PER_LITRE = _get_box.NM3_PER_LITRE
box_from_n_mol_and_conc = _get_box.box_from_n_mol_and_conc
conc_from_n_mol_and_volume = _get_box.conc_from_n_mol_and_volume
n_mol_from_conc_and_volume = _get_box.n_mol_from_conc_and_volume
sphere_volume = _get_box.sphere_volume


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _run(*args):
    """Run the script with the given CLI arguments, return CompletedProcess."""
    return subprocess.run(
        [sys.executable, str(_SCRIPT), *args],
        capture_output=True,
        text=True,
    )


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

class TestConstants:
    def test_avogadro_codata_2019(self):
        """AVOGADRO must be the 2019 CODATA exact value."""
        assert AVOGADRO == pytest.approx(6.02214076e23, rel=1e-10)

    def test_nm3_per_litre(self):
        """1 L = (10 dm)^... actually 1 L = 1e24 nm³."""
        assert NM3_PER_LITRE == 1e24


# ---------------------------------------------------------------------------
# sphere_volume
# ---------------------------------------------------------------------------

class TestSphereVolume:
    def test_unit_sphere(self):
        assert sphere_volume(1.0) == pytest.approx(4 / 3 * math.pi, rel=1e-12)

    def test_zero_radius(self):
        assert sphere_volume(0.0) == pytest.approx(0.0)

    def test_scales_with_cube_of_radius(self):
        assert sphere_volume(2.0) == pytest.approx(8 * sphere_volume(1.0), rel=1e-12)


# ---------------------------------------------------------------------------
# box_from_n_mol_and_conc
# ---------------------------------------------------------------------------

class TestBoxFromNMolAndConc:
    def test_one_molecule_one_molar(self):
        # V = 1 / (1 * NA) L → nm³ → side
        v_nm3 = NM3_PER_LITRE / AVOGADRO
        expected = v_nm3 ** (1 / 3)
        assert box_from_n_mol_and_conc(1, 1.0) == pytest.approx(expected, rel=1e-10)

    def test_box_increases_with_n_mol(self):
        assert box_from_n_mol_and_conc(100, 0.1) > box_from_n_mol_and_conc(10, 0.1)

    def test_box_decreases_with_concentration(self):
        assert box_from_n_mol_and_conc(50, 1.0) < box_from_n_mol_and_conc(50, 0.1)

    def test_known_value(self):
        # 50 molecules at 0.1 M — verified against README example
        assert box_from_n_mol_and_conc(50, 0.1) == pytest.approx(9.39881, rel=1e-5)


# ---------------------------------------------------------------------------
# conc_from_n_mol_and_volume
# ---------------------------------------------------------------------------

class TestConcFromNMolAndVolume:
    def test_one_molecule_one_litre(self):
        # 1 molecule in 1 L (= 1e24 nm³) → 1/NA mol/L
        assert conc_from_n_mol_and_volume(1, NM3_PER_LITRE) == pytest.approx(
            1 / AVOGADRO, rel=1e-10
        )

    def test_linear_in_n_mol(self):
        v = 500.0
        assert conc_from_n_mol_and_volume(20, v) == pytest.approx(
            2 * conc_from_n_mol_and_volume(10, v), rel=1e-12
        )

    def test_inversely_proportional_to_volume(self):
        n = 10
        assert conc_from_n_mol_and_volume(n, 1000.0) == pytest.approx(
            2 * conc_from_n_mol_and_volume(n, 2000.0), rel=1e-12
        )


# ---------------------------------------------------------------------------
# n_mol_from_conc_and_volume
# ---------------------------------------------------------------------------

class TestNMolFromConcAndVolume:
    def test_rounds_to_nearest_integer(self):
        result = n_mol_from_conc_and_volume(0.1, 1000.0)
        assert isinstance(result, int)

    def test_known_value(self):
        # 0.1 M in 1000 nm³: n = 0.1 * NA * (1000 / 1e24) = 0.1 * NA * 1e-21
        expected = int(round(0.1 * AVOGADRO * 1e-21))
        assert n_mol_from_conc_and_volume(0.1, 1000.0) == expected

    def test_more_volume_more_molecules(self):
        assert n_mol_from_conc_and_volume(0.1, 2000.0) > n_mol_from_conc_and_volume(
            0.1, 1000.0
        )


# ---------------------------------------------------------------------------
# Roundtrip consistency
# ---------------------------------------------------------------------------

class TestRoundtrip:
    """Each pair of functions should be exact inverses of each other."""

    def test_box_then_conc(self):
        """box_from_n_mol_and_conc → conc_from_n_mol_and_volume recovers original conc."""
        n, c = 50, 0.15
        side = box_from_n_mol_and_conc(n, c)
        volume = side**3
        assert conc_from_n_mol_and_volume(n, volume) == pytest.approx(c, rel=1e-10)

    def test_box_then_n_mol(self):
        """box_from_n_mol_and_conc → n_mol_from_conc_and_volume recovers original n_mol."""
        n, c = 100, 0.05
        side = box_from_n_mol_and_conc(n, c)
        volume = side**3
        assert n_mol_from_conc_and_volume(c, volume) == n

    def test_conc_then_box(self):
        """conc_from_n_mol_and_volume → box_from_n_mol_and_conc recovers original side."""
        n, v = 30, 800.0
        c = conc_from_n_mol_and_volume(n, v)
        side = box_from_n_mol_and_conc(n, c)
        assert side**3 == pytest.approx(v, rel=1e-10)


# ---------------------------------------------------------------------------
# CLI — correct output
# ---------------------------------------------------------------------------

class TestCLIOutput:
    def test_box_side_five_decimal_places(self):
        result = _run("--n_mol", "50", "--conc", "0.1")
        assert result.returncode == 0
        line = result.stdout.strip()
        # value should have exactly 5 decimal places
        value_str = line.split()[-2]  # e.g. "9.39881"
        assert len(value_str.split(".")[-1]) == 5

    def test_box_side_value(self):
        result = _run("--n_mol", "50", "--conc", "0.1")
        assert "9.39881" in result.stdout

    def test_concentration_output(self):
        result = _run("--n_mol", "10", "--volume", "1000.0")
        assert result.returncode == 0
        assert "concentration" in result.stdout

    def test_n_mol_output(self):
        result = _run("--conc", "0.1", "--volume", "1000.0")
        assert result.returncode == 0
        assert "n_mol" in result.stdout

    def test_sphere_r_alone_prints_volume(self):
        result = _run("--sphere_r", "5.0")
        assert result.returncode == 0
        assert "sphere volume" in result.stdout

    def test_cubic_side_alone_prints_volume(self):
        result = _run("--cubic_side", "8.0")
        assert result.returncode == 0
        assert "cubic volume" in result.stdout
        assert "512.0000" in result.stdout

    def test_sphere_r_with_conc(self):
        result = _run("--conc", "0.1", "--sphere_r", "5.0")
        assert result.returncode == 0
        assert "n_mol" in result.stdout


# ---------------------------------------------------------------------------
# CLI — error handling
# ---------------------------------------------------------------------------

class TestCLIErrors:
    def test_no_args_exits_cleanly(self):
        result = _run()
        assert result.returncode == 0  # prints help and exits 0

    def test_negative_n_mol_is_rejected(self):
        result = _run("--n_mol", "-5", "--conc", "0.1")
        assert result.returncode != 0

    def test_negative_conc_is_rejected(self):
        result = _run("--n_mol", "10", "--conc", "-0.1")
        assert result.returncode != 0

    def test_sphere_r_and_volume_are_mutually_exclusive(self):
        result = _run("--n_mol", "10", "--sphere_r", "3.0", "--volume", "100.0")
        assert result.returncode != 0

    def test_cubic_side_and_volume_are_mutually_exclusive(self):
        result = _run("--n_mol", "10", "--cubic_side", "5.0", "--volume", "100.0")
        assert result.returncode != 0

    def test_over_specified_is_rejected(self):
        result = _run("--n_mol", "10", "--conc", "0.1", "--volume", "500.0")
        assert result.returncode != 0

    def test_under_specified_is_rejected(self):
        result = _run("--n_mol", "10")
        assert result.returncode != 0
