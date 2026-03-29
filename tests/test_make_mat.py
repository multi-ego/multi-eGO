"""
Integration tests for tools/make_mat/make_mat.py and tools/make_mat/HDF52ndx.py.

Replaces run_make_mat.sh with pytest-discoverable equivalents.

Each test class covers one dataset (ttr, popc).  Within a class the fixtures
extract the compressed histo archive once per session, run the tools, and then
individual test methods assert that the output files match the references in
test_outputs/.
"""

import filecmp
import os
import subprocess
import sys
import tarfile

import pytest

TEST_ROOT = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(TEST_ROOT, ".."))
TOOLS = os.path.join(REPO_ROOT, "tools", "make_mat")
SRC = os.path.join(REPO_ROOT, "src")

MAKE_MAT = os.path.join(TOOLS, "make_mat.py")
HDF52NDX = os.path.join(TOOLS, "HDF52ndx.py")

# Ensure subprocesses find the multiego *package* (src/multiego/) and not the
# root-level multiego.py launcher, which would shadow it.
_TOOL_ENV = {
    **os.environ,
    "PYTHONPATH": SRC + (os.pathsep + os.environ["PYTHONPATH"] if "PYTHONPATH" in os.environ else ""),
}


def _run(*args, cwd=None):
    """Run a command, raising on non-zero exit."""
    result = subprocess.run(
        [sys.executable, *args],
        capture_output=True,
        text=True,
        cwd=cwd,
        env=_TOOL_ENV,
    )
    if result.returncode != 0:
        raise AssertionError(
            f"Command failed (exit {result.returncode}):\n"
            f"  {' '.join(str(a) for a in args)}\n"
            f"stdout: {result.stdout}\n"
            f"stderr: {result.stderr}"
        )
    return result


def _extract_histo(tgz_path, dest_dir):
    """Extract histo archive into dest_dir (idempotent)."""
    with tarfile.open(tgz_path) as tf:
        tf.extractall(dest_dir, filter="data")


def _assert_files_equal(actual, reference):
    """Assert two files have identical content, printing a diff on failure."""
    if not filecmp.cmp(actual, reference, shallow=False):
        with open(actual) as fh:
            actual_lines = fh.readlines()
        with open(reference) as fh:
            ref_lines = fh.readlines()
        import difflib
        diff = "".join(
            difflib.unified_diff(ref_lines, actual_lines, fromfile="reference", tofile="actual", n=3)
        )
        raise AssertionError(
            f"Files differ:\n  actual:    {actual}\n  reference: {reference}\n\n{diff}"
        )


# ===========================================================================
# TTR dataset
# ===========================================================================

TTR_IN = os.path.join(TEST_ROOT, "test_inputs", "make_mat_ttr")
TTR_OUT = os.path.join(TEST_ROOT, "test_outputs", "make_mat_ttr")


@pytest.fixture(scope="module")
def ttr_normal(tmp_path_factory):
    """
    Run make_mat.py (normal mode) + HDF52ndx.py on the TTR dataset.
    Returns the directory where output files were written.
    """
    out = tmp_path_factory.mktemp("ttr_normal")
    histo_dir = tmp_path_factory.mktemp("ttr_histo")
    _extract_histo(os.path.join(TTR_IN, "hh.tgz"), histo_dir)

    _run(
        MAKE_MAT,
        "--histo", str(histo_dir / "histo"),
        "--target_top", os.path.join(TTR_IN, "topol_mego.top"),
        "--mego_top",   os.path.join(TTR_IN, "topol.top"),
        "--cutoff", "0.75",
        "--mode", "intra+same",
        "--out", str(out) + os.sep,
    )
    _run(HDF52NDX, "--input", str(out / "intramat_1_1.ndx.h5"))
    _run(HDF52NDX, "--input", str(out / "intermat_1_1.ndx.h5"))
    return out


@pytest.fixture(scope="module")
def ttr_zero(tmp_path_factory):
    """
    Run make_mat.py (--zero mode) + HDF52ndx.py on the TTR dataset.
    Returns the directory where output files were written.
    """
    out = tmp_path_factory.mktemp("ttr_zero")

    _run(
        MAKE_MAT,
        "--zero",
        "--target_top", os.path.join(TTR_IN, "topol_mego.top"),
        "--mego_top",   os.path.join(TTR_IN, "topol.top"),
        "--cutoff", "0.75",
        "--mode", "intra+same",
        "--out", str(out) + os.sep,
    )
    _run(HDF52NDX, "--input", str(out / "intramat_1_1.ndx.h5"))
    _run(HDF52NDX, "--input", str(out / "intermat_1_1.ndx.h5"))
    return out


class TestMakeMatTTRNormal:

    def test_intramat(self, ttr_normal):
        _assert_files_equal(
            ttr_normal / "intramat_1_1.ndx",
            os.path.join(TTR_OUT, "intramat_1_1.ndx"),
        )

    def test_intermat(self, ttr_normal):
        _assert_files_equal(
            ttr_normal / "intermat_1_1.ndx",
            os.path.join(TTR_OUT, "intermat_1_1.ndx"),
        )


class TestMakeMatTTRZero:

    def test_intramat_zero(self, ttr_zero):
        _assert_files_equal(
            ttr_zero / "intramat_1_1.ndx",
            os.path.join(TTR_OUT, "intramat_1_1.ndx.zero"),
        )

    def test_intermat_zero(self, ttr_zero):
        _assert_files_equal(
            ttr_zero / "intermat_1_1.ndx",
            os.path.join(TTR_OUT, "intermat_1_1.ndx.zero"),
        )


# ===========================================================================
# POPC dataset
# ===========================================================================

POPC_IN = os.path.join(TEST_ROOT, "test_inputs", "make_mat_popc")
POPC_OUT = os.path.join(TEST_ROOT, "test_outputs", "make_mat_popc")


@pytest.fixture(scope="module")
def popc_normal(tmp_path_factory):
    """
    Run make_mat.py on the POPC dataset (--noh5 mode).
    Returns the directory where output files were written.
    """
    out = tmp_path_factory.mktemp("popc_normal")
    histo_dir = tmp_path_factory.mktemp("popc_histo")
    _extract_histo(os.path.join(POPC_IN, "hh.tgz"), histo_dir)

    _run(
        MAKE_MAT,
        "--histo", str(histo_dir / "histo"),
        "--target_top", os.path.join(POPC_IN, "topol_md.top"),
        "--mego_top",   os.path.join(POPC_IN, "topol_ref.top"),
        "--cutoff", "0.75",
        "--mode", "intra",
        "--out", str(out) + os.sep,
        "--noh5",
    )
    return out


class TestMakeMatPOPC:

    def test_runs_without_error(self, popc_normal):
        """make_mat.py should complete successfully for the POPC dataset."""
        out_file = popc_normal / "intramat_1_1.ndx.h5"
        assert out_file.exists(), f"Expected output not found: {out_file}"
        assert out_file.stat().st_size > 0, f"Output file is empty: {out_file}"
