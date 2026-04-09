"""
Regression tests for tools/domain_sectioner/domains.py.

Each fixture runs domains.py on the reference dataset in
tests/test_inputs/domain_sectioner/, then converts the HDF5 output back to
plain text with HDF52ndx.py.  The resulting text file is compared against the
reference output stored alongside the inputs.
"""

import os
import subprocess
import sys
import filecmp
import difflib
import pandas as pd

import pytest

TEST_ROOT = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(TEST_ROOT, ".."))

DOMAINS = os.path.join(REPO_ROOT, "tools", "domain_sectioner", "domains.py")
HDF52NDX = os.path.join(REPO_ROOT, "tools", "make_mat", "HDF52ndx.py")

DOM_IN = os.path.join(TEST_ROOT, "test_inputs", "domain_sectioner")
DOM_OUT = os.path.join(TEST_ROOT, "test_outputs", "domain_sectioner")

_TOOL_ENV = {
    **os.environ,
    "PYTHONPATH": os.path.join(REPO_ROOT, "src")
    + (os.pathsep + os.environ["PYTHONPATH"] if "PYTHONPATH" in os.environ else ""),
}


def _run(*args, cwd=None):
    """Run a Python script, raising AssertionError on non-zero exit."""
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


def _assert_hdf_equal(a, b):
    df_a = pd.read_hdf(a)
    df_b = pd.read_hdf(b)

    pd.testing.assert_frame_equal(df_a, df_b)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="class")
def dom_split(tmp_path_factory):
    """Run domains.py (split mode) and convert output to plain text."""
    out = tmp_path_factory.mktemp("dom_split")

    _run(
        DOMAINS,
        "--intra",
        os.path.join(DOM_IN, "intramat_1_1.ndx.h5"),
        "--top",
        os.path.join(DOM_IN, "topol.top"),
        "--dom_res",
        "2-61",
        "72-161",
        "--out",
        str(out),
    )
    return out


@pytest.fixture(scope="class")
def dom_invert(tmp_path_factory):
    """Run domains.py (invert mode) and convert output to plain text."""
    out = tmp_path_factory.mktemp("dom_invert")

    _run(
        DOMAINS,
        "--intra",
        os.path.join(DOM_IN, "intramat_1_1.ndx.h5"),
        "--top",
        os.path.join(DOM_IN, "topol.top"),
        "--dom_res",
        "2-61",
        "72-161",
        "--invert",
        "--out",
        str(out),
    )
    return out


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestDomainSplit:

    def test_split(self, dom_split):
        _assert_hdf_equal(
            dom_split / "split_2-61-72-161_intramat_1_1.h5",
            os.path.join(DOM_OUT, "ref_split_2-61-72-161_intramat_1_1.h5"),
        )


class TestDomainInvert:

    def test_invert(self, dom_invert):
        _assert_hdf_equal(
            dom_invert / "inverted_split_2-61-72-161_intramat_1_1.h5",
            os.path.join(DOM_OUT, "ref_inverted_split_2-61-72-161_intramat_1_1.h5"),
        )
