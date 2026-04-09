"""
domains.py — apply a domain mask to a multi-eGO intra-molecular contact matrix.

For each domain range supplied via ``--dom_res``, the contacts between atoms
that both fall inside the range are marked as *learned* (``True``) in the
output HDF5 file.  Contacts outside all ranges are marked as *not learned*
(``False``).  The ``--invert`` flag flips the assignment.
"""

import argparse
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import parmed as pmd

# ---------------------------------------------------------------------------
# Topology helpers
# ---------------------------------------------------------------------------


def find_atom_start(top, res_num):
    """Return the 0-based index of the first atom of residue *res_num* (1-based).

    Parameters
    ----------
    top : parmed.Structure
        Parsed topology.
    res_num : int
        1-based residue number.

    Returns
    -------
    int
    """
    return sum(len(top.residues[i].atoms) for i in range(res_num - 1))


def find_atom_end(top, res_num):
    """Return the 0-based index of the last atom of residue *res_num* (1-based).

    Parameters
    ----------
    top : parmed.Structure
        Parsed topology.
    res_num : int
        1-based residue number.

    Returns
    -------
    int
    """
    return sum(len(top.residues[i].atoms) for i in range(res_num)) - 1


def read_topology(top_path):
    """Read a GROMACS topology file and return the parsed structure and a summary DataFrame.

    Parameters
    ----------
    top_path : str
        Path to the topology file (e.g. ``topol.top``).

    Returns
    -------
    topology : parmed.Structure
        Parsed topology.
    top_df : pd.DataFrame
        Per-molecule summary with columns ``name``, ``residues``,
        ``atoms_per_res``, ``tot_atoms``, and ``atoms_name``.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        topology = pmd.load_file(top_path)

    mol_names = list(topology.molecules.keys())
    res, atoms, atoms_name, tot_num_atoms = [], [], [], []
    for name in mol_names:
        mol = topology.molecules[name][0]
        res.append([r.name for r in mol.residues])
        atoms.append([len(r.atoms) for r in mol.residues])
        atoms_name.append([[a.type for a in r.atoms] for r in mol.residues])
        tot_num_atoms.append(sum(len(r.atoms) for r in mol.residues))

    top_df = pd.DataFrame(
        {
            "name": mol_names,
            "residues": res,
            "atoms_per_res": atoms,
            "tot_atoms": tot_num_atoms,
            "atoms_name": atoms_name,
        }
    )
    return topology, top_df


# ---------------------------------------------------------------------------
# Domain range parsing
# ---------------------------------------------------------------------------


def dom_range(ranges_str):
    """Parse ``'start-end'`` strings into validated ``(start, end)`` tuples.

    Parameters
    ----------
    ranges_str : list of str
        Strings of the form ``'start-end'``, e.g. ``['1-30', '45-60']``.

    Returns
    -------
    list of (int, int)

    Raises
    ------
    ValueError
        If any range is decreasing (start > end) or if consecutive ranges
        overlap or are not in increasing order.
    """
    doms = [(int(r.split("-")[0]), int(r.split("-")[1])) for r in ranges_str]

    for start, end in doms:
        if start > end:
            raise ValueError(f"Domain range {start}-{end} is invalid: start must be ≤ end.")

    for (_, end1), (start2, _) in zip(doms[:-1], doms[1:]):
        if end1 >= start2:
            raise ValueError(
                f"Domain ranges overlap or are out of order: "
                f"end of first range ({end1}) must be < start of next ({start2})."
            )

    return doms


# ---------------------------------------------------------------------------
# Domain mask
# ---------------------------------------------------------------------------


def build_domain_mask(topology, n_atoms, ranges, invert=False):
    """Build a flat boolean mask of length *n_atoms²* marking intra-domain contacts.

    A contact (i, j) is marked ``True`` when *both* atom i and atom j fall
    inside at least one of the supplied residue ranges.

    Parameters
    ----------
    topology : parmed.Structure
        Parsed topology (used to map residue numbers to atom indices).
    n_atoms : int
        Total number of atoms in the molecule.
    ranges : list of (int, int)
        Validated residue ranges from :func:`dom_range`.
    invert : bool
        If ``True``, flip the mask so that intra-domain contacts are ``False``
        and all others are ``True``.

    Returns
    -------
    np.ndarray of bool, shape (n_atoms²,)
    """
    atom_in_domain = np.zeros(n_atoms, dtype=bool)
    for r_start, r_end in ranges:
        a_start = find_atom_start(topology, r_start)
        a_end = find_atom_end(topology, r_end)
        print(f"  Domain range: {r_start}-{r_end}")
        print(f"    Atom index range (1-based): {a_start + 1} – {a_end + 1}")
        print(f"    Number of atoms in range:   {a_end - a_start + 1}")
        print(f"    First / last atom:  " f"{topology.atoms[a_start]} – {topology.atoms[a_end]}")
        atom_in_domain[a_start : a_end + 1] = True

    # Outer product: True where both atom i and atom j are in the domain
    mask_2d = atom_in_domain[:, np.newaxis] & atom_in_domain[np.newaxis, :]
    domain_mask = mask_2d.reshape(n_atoms**2)

    if invert:
        domain_mask = ~domain_mask

    return domain_mask


# ---------------------------------------------------------------------------
# Contact matrix I/O
# ---------------------------------------------------------------------------


def readmat(intramat_path, h5=True):
    """Read an intra-molecular contact matrix and return it as a DataFrame.

    Parameters
    ----------
    intramat_path : str
        Path to the contact matrix file.
    h5 : bool
        If ``True`` (default), read as HDF5; otherwise read as whitespace-
        delimited text.

    Returns
    -------
    pd.DataFrame
        Contact matrix with columns:
        ``molecule_name_ai``, ``ai``, ``molecule_name_aj``, ``aj``,
        ``distance``, ``probability``, ``cutoff``, ``learned``.
    """
    if h5:
        return pd.read_hdf(intramat_path, key="data")

    col_names = [
        "molecule_name_ai",
        "ai",
        "molecule_name_aj",
        "aj",
        "distance",
        "probability",
        "cutoff",
        "learned",
    ]
    col_types = {
        "molecule_name_ai": int,
        "ai": int,
        "molecule_name_aj": int,
        "aj": int,
        "distance": np.float64,
        "probability": np.float64,
        "cutoff": np.float64,
        "learned": int,
    }

    raw = np.loadtxt(intramat_path, unpack=True)
    if raw.shape[0] not in (7, 8):
        raise ValueError(
            "Intramat must have 7 or 8 columns: "
            "molecule_name_ai, ai, molecule_name_aj, aj, "
            "distance, probability, cutoff, [learned]"
        )
    if raw.shape[0] == 7:
        print("  Intramat has 7 columns — domain mask column will be appended.")
        raw = np.vstack([raw, np.zeros(raw.shape[1], dtype=int)])

    df = pd.DataFrame({col: raw[i].astype(col_types[col]) for i, col in enumerate(col_names)})
    return df


def write_intramat(contact_matrix, out_path):
    """Write a contact matrix DataFrame to an HDF5 file.

    Parameters
    ----------
    contact_matrix : pd.DataFrame
        Contact matrix (must contain a ``learned`` column).
    out_path : str or Path
        Destination file path (should end in ``.h5``).
    """
    contact_matrix = contact_matrix.copy()
    contact_matrix["learned"] = contact_matrix["learned"].fillna(1).astype(bool)
    for col in ("molecule_name_ai", "ai", "molecule_name_aj", "aj"):
        contact_matrix[col] = contact_matrix[col].astype("category")

    contact_matrix.to_hdf(
        out_path,
        key="data",
        mode="w",
        format="table",
        complib="blosc:lz4",
        complevel=9,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Apply a domain mask to a multi-eGO intra-molecular contact matrix.\n\n"
            "Contacts where both atoms fall inside one of the supplied residue "
            "ranges are marked as learned; all others are marked as not learned. "
            "Use --invert to flip this assignment."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--intra",
        type=str,
        required=True,
        help="Path to the intra-molecular contact matrix (HDF5 .h5 or plain text)",
    )
    parser.add_argument(
        "--top",
        type=str,
        required=True,
        help="Path to the reference topology file (e.g. topol.top)",
    )
    parser.add_argument(
        "--dom_res",
        nargs="+",
        type=str,
        required=True,
        help=(
            "One or more residue ranges as 'start-end' pairs, e.g. --dom_res 1-30 45-60. "
            "Ranges must be non-overlapping and in increasing order."
        ),
    )
    parser.add_argument(
        "--out",
        type=str,
        default=".",
        help="Output directory (default: current directory)",
    )
    parser.add_argument(
        "--invert",
        action="store_true",
        default=False,
        help="Invert the domain mask: mark intra-domain contacts as not learned",
    )

    args = parser.parse_args()

    out_dir = Path(args.out)
    if not out_dir.is_dir():
        print(f"ERROR: output directory '{out_dir}' does not exist.")
        sys.exit(1)

    # Read topology
    topology, top_df = read_topology(args.top)
    if len(top_df) > 1:
        raise ValueError(
            "Only a single molecule species is supported; " f"topology contains {len(top_df)} molecule types."
        )

    n_atoms = top_df.tot_atoms[0]
    n_res = len(top_df.residues[0])
    print(f"Topology: {n_res} residues, {n_atoms} atoms")

    # Parse and validate domain ranges
    ranges = dom_range(args.dom_res)

    # Read contact matrix
    intra_path = Path(args.intra)
    h5 = intra_path.suffix == ".h5"
    if h5:
        print(f"Reading contact matrix: {intra_path} (HDF5)")
    intra_md = readmat(str(intra_path), h5=h5)

    dim = int(round(len(intra_md) ** 0.5))
    if dim != n_atoms:
        raise ValueError(f"Contact matrix size ({dim}²) does not match topology atom count ({n_atoms}).")

    # Build and apply domain mask
    domain_mask = build_domain_mask(topology, n_atoms, ranges, invert=args.invert)
    intra_md["learned"] = domain_mask

    # Write output
    intramat_stem = intra_path.stem
    if intramat_stem.endswith(".ndx"):
        intramat_stem = intramat_stem[: -len(".ndx")]
    range_str = "-".join(args.dom_res)
    prefix = "inverted_split" if args.invert else "split"
    out_path = out_dir / f"{prefix}_{range_str}_{intramat_stem}.h5"

    write_intramat(intra_md, out_path)
    print(f"Written: {out_path}")
