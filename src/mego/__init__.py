"""
Entry point for the ``mego`` console command.

This package is separate from ``multiego`` so the console script can fix
``sys.path`` before importing the ``multiego`` package.  The editable install
may put the project root ahead of ``src/`` in ``sys.path``, which causes
``multiego.py`` (the root launcher script) to shadow the ``multiego`` package.
By loading this module first (it has no name conflict), we can correct the
path before any ``multiego`` import happens.
"""

import os
import sys


def main():
    """Entry point for the ``mego`` console command."""
    # This file lives at src/mego/__init__.py, so:
    #   os.path.dirname(__file__)           → src/mego/
    #   os.path.dirname(dirname(__file__))  → src/
    _src = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    _root = os.path.dirname(_src)

    # Put src/ at position 0 and strip the project root so that
    # ``import multiego`` finds src/multiego/ rather than multiego.py.
    sys.path = [_src] + [
        p for p in sys.path if os.path.realpath(p) != os.path.realpath(_root)
    ]

    from multiego._run import main as _main  # noqa: E402

    _main(root_dir=os.getcwd())
