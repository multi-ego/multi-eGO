"""
Root-level launcher script — kept for backward compatibility so that
``python multiego.py <args>`` continues to work from the repository root.

When running from the repo, ``root_dir`` is the directory that contains
this script (i.e. the repo root), which is also where ``inputs/`` and
``outputs/`` live.

Users who install the package via pip get a ``multiego`` console command
instead (defined in pyproject.toml); that command uses the current working
directory as ``root_dir`` so it works from any location.
"""

import os
import sys

# When running as a script the repo root is on sys.path, which means
# "multiego" resolves to this file itself rather than the package under src/.
# Prepend src/ so the package is found first.
_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "src")
sys.path.insert(0, _src)  # always first, so src/multiego/ beats multiego.py

from multiego._run import main  # noqa: E402

if __name__ == "__main__":
    main(root_dir=os.path.dirname(os.path.abspath(__file__)))
