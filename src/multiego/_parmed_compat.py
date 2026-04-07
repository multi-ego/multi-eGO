"""
_parmed_compat.py — ParmEd compatibility helpers for multi-eGO.

ParmEd's C preprocessor resolves ``#include`` directives by searching:
  1. The directory of the topology file being loaded.
  2. An ``includes`` kwarg added to ``GromacsTopologyFile`` in newer ParmEd versions.
  3. Some versions also read ``GMXDATA`` / ``GMXLIB``, but this is not reliable.

Because older ParmEd releases (≤ 4.2.x) do not expose the ``includes`` parameter
and do not honour ``GMXLIB``, we use a different approach: temporarily symlink every
``.ff`` directory from ``GMXLIB`` into the topology file's own directory before
calling ``parmed.load_file``, then clean up the symlinks afterwards.

This strategy works with any ParmEd version and requires no ParmEd API knowledge.
"""

import contextlib
import os


@contextlib.contextmanager
def gmxlib_in_topo_dir(topology_path: str):
    """Temporarily symlink GMXLIB force-field directories next to *topology_path*.

    ParmEd always searches the directory that contains the topology file for
    ``#include`` targets.  By symlinking ``<GMXLIB>/*.ff`` into that directory we
    make custom force fields (e.g. ``multi-ego-basic.ff``) visible without relying
    on any particular ParmEd version feature.

    Usage::

        with gmxlib_in_topo_dir(topology_path):
            struct = parmed.load_file(topology_path, ...)
    """
    gmxlib = os.environ.get("GMXLIB", "")
    created: list[str] = []

    if gmxlib and os.path.isdir(gmxlib):
        topo_dir = os.path.dirname(os.path.abspath(topology_path))
        for item in os.listdir(gmxlib):
            if item.endswith(".ff"):
                src = os.path.join(gmxlib, item)
                link = os.path.join(topo_dir, item)
                if os.path.isdir(src) and not os.path.exists(link):
                    os.symlink(src, link)
                    created.append(link)

    try:
        yield
    finally:
        for link in created:
            with contextlib.suppress(OSError):
                os.unlink(link)
