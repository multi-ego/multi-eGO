"""
Shared pytest fixtures for the multi-eGO test suite.

Provides module-scoped stubs for heavy or unavailable dependencies so that
individual test modules don't need to repeat the same boilerplate.
"""

import sys
import types
import pytest


def _install_stub(name, **attrs):
    """Register a stub module in sys.modules, optionally with extra attributes."""
    mod = sys.modules.get(name)
    if mod is None:
        mod = types.ModuleType(name)
        sys.modules[name] = mod
    for k, v in attrs.items():
        setattr(mod, k, v)
    return mod


@pytest.fixture(scope="module")
def stub_deps():
    """
    Install lightweight stub modules for every heavy dependency that the
    multiego package imports at module level.  Returns a dict mapping
    module name → stub so individual tests can add attributes if needed.

    Scope is ``module`` so the stubs are installed once per test file and
    remain for the lifetime of that file's tests.
    """
    stubs = {}

    # model_config — needs a config object with specific attributes
    stubs["multiego.model_config"] = _install_stub(
        "multiego.model_config",
        config=types.SimpleNamespace(max_bond_separation=5, bond14_separation=3),
    )

    # type_definitions — lj and mg read special_non_local at import time
    stubs["multiego.type_definitions"] = _install_stub(
        "multiego.type_definitions",
        special_non_local=[],
    )

    # ensemble_data — contacts.py references MeGOEnsemble
    stubs["multiego.ensemble_data"] = _install_stub(
        "multiego.ensemble_data",
        MeGOEnsemble=type("MeGOEnsemble", (), {}),
    )

    # multiego.io — no attributes needed at import time
    stubs["multiego.io"] = _install_stub("multiego.io")

    # parmed — third-party, may not be installed in all test environments
    stubs["parmed"] = _install_stub("parmed")

    # Note: multiego.mg and multiego.bonded are NOT stubbed here.
    # They only depend on type_definitions and model_config (both already stubbed
    # above), so they load cleanly as real modules.  Tests that import mg or
    # bonded directly (TestGenerateMGLJPairsRep, etc.) need the real
    # implementations, not empty stubs.

    return stubs
