#!/usr/bin/env python3

"""
dihedral_histo.py

Normalised, periodic histograms of the dihedral time series written by
`gmx chi -all` (one xvg file per dihedral, e.g. phiALA2.xvg, psiALA2.xvg,
omegaALA2.xvg, chi1LYS12.xvg).

Every dihedral is histogrammed on the same fixed grid spanning the whole
periodic range [-180, 180) degrees, so that distributions from different
simulations can be compared bin by bin.

Sub-commands
  histo    xvg files  -> one table, bin centres + one column per dihedral
  compare  two or more `histo` tables -> per-dihedral similarity metrics
           (Jensen-Shannon divergence, Hellinger distance, overlap)
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np

# gmx chi -all file names: <type><resname><resnr>.xvg
DIH_RE = re.compile(r"^(phi|psi|omega|chi\d+)(.*?)(-?\d+)$")
FULL = 360.0
LOWER = -180.0


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def dihedral_key(name: str):
    """Sort key: residue number, then phi, psi, omega, chi1, chi2, ..."""
    m = DIH_RE.match(name)
    if not m:
        return (float("inf"), 99, name)
    kind, _, resnr = m.groups()
    if kind.startswith("chi"):
        order = 2 + int(kind[3:])
    else:
        order = {"phi": 0, "psi": 1, "omega": 2}[kind]
    return (int(resnr), order, name)


def wrap(angles_deg: np.ndarray) -> np.ndarray:
    """Map angles onto [-180, 180): +180 and -180 are the same point."""
    return np.mod(angles_deg - LOWER, FULL) + LOWER


def bin_grid(nbins: int):
    width = FULL / nbins
    edges = LOWER + width * np.arange(nbins + 1)
    centers = edges[:-1] + 0.5 * width
    return width, edges, centers


def periodic_histogram(angles_deg, nbins, weights=None):
    """Return (probability per bin, effective counts). The last bin closes on the first."""
    width = FULL / nbins
    idx = np.floor((wrap(angles_deg) - LOWER) / width).astype(np.int64)
    idx %= nbins  # guards round-off at the +180 boundary
    counts = np.bincount(idx, weights=weights, minlength=nbins)
    total = counts.sum()
    if total <= 0:
        raise ValueError("empty histogram (no frames or all weights zero)")
    return counts / total, total


def read_xvg(path: Path):
    data = np.loadtxt(path, comments=("#", "@"), ndmin=2)
    if data.shape[0] == 0 or data.shape[1] < 2:
        raise ValueError(f"{path}: expected two columns (time, angle)")
    return data[:, 0], data[:, 1]


def frame_selection(time, begin, end, stride):
    mask = time >= begin
    if end >= 0:
        mask &= time <= end
    sel = np.flatnonzero(mask)[::stride]
    return sel


def collect_files(inputs, pattern):
    files = []
    for p in inputs:
        p = Path(p)
        if p.is_dir():
            files.extend(sorted(p.glob(pattern)))
        else:
            files.append(p)
    # keep only gmx chi -all dihedral files (drops order.xvg, histo-*.xvg, ...)
    keep = [f for f in files if DIH_RE.match(f.stem)]
    skipped = sorted(set(files) - set(keep))
    return keep, skipped


def write_table(path, centers, names, columns, header_lines):
    with open(path, "w") as f:
        for line in header_lines:
            f.write(f"# {line}\n")
        f.write("# columns: center " + " ".join(names) + "\n")
        table = np.column_stack([centers] + columns)
        np.savetxt(f, table, fmt=["%10.4f"] + ["%.8e"] * len(columns))


def write_xvg(path, centers, values, name, ylabel, header_lines):
    """Single-histogram xvg file (xmgrace-readable, loads with numpy.loadtxt(comments=('#', '@')))."""
    with open(path, "w") as f:
        for line in header_lines:
            f.write(f"# {line}\n")
        f.write(f'@    title "{name}"\n')
        f.write('@    xaxis  label "angle (deg)"\n')
        f.write(f'@    yaxis  label "{ylabel}"\n')
        f.write("@TYPE xy\n")
        f.write(f'@ s0 legend "{name}"\n')
        np.savetxt(f, np.column_stack([centers, values]), fmt=["%10.4f", "%.8e"])


def read_table(path: Path):
    names = None
    norm = None
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            if line.startswith("# columns:"):
                names = line.split(":", 1)[1].split()[1:]
            if line.startswith("# normalisation:"):
                norm = line.split(":", 1)[1].split()[0]
    if names is None:
        raise ValueError(f"{path}: missing '# columns:' header, not a dihedral_histo.py table")
    data = np.loadtxt(path, comments="#", ndmin=2)
    centers = data[:, 0]
    # bring every column back to a probability per bin, whatever the normalisation
    probs = data[:, 1:]
    probs = probs / probs.sum(axis=0, keepdims=True)
    return centers, dict(zip(names, probs.T)), norm


# --------------------------------------------------------------------------- #
# histo
# --------------------------------------------------------------------------- #
def run_histo(args):
    files, skipped = collect_files(args.inputs, args.pattern)
    if not files:
        sys.exit("No gmx chi -all dihedral files found (expected names like phiALA2.xvg, chi1LYS12.xvg).")
    if skipped and args.verbose:
        print(f"skipped {len(skipped)} non-dihedral files: " + " ".join(s.name for s in skipped), file=sys.stderr)

    weights_all = None
    if args.weights is not None:
        weights_all = np.loadtxt(args.weights, comments=("#", "@"), ndmin=2)[:, -1]
        if np.any(weights_all < 0):
            sys.exit(f"{args.weights}: negative weights")

    width, _, centers = bin_grid(args.nbins)
    files = sorted(files, key=lambda f: dihedral_key(f.stem))

    names, probs, n_frames, n_eff = [], [], [], []
    ref_time = None
    for f in files:
        time, angle = read_xvg(f)
        if ref_time is None:
            ref_time = time
        elif len(time) != len(ref_time) or not np.allclose(time, ref_time):
            print(f"warning: {f.name} has a different time axis from {files[0].name}", file=sys.stderr)
        if args.rad:
            angle = np.degrees(angle)
        sel = frame_selection(time, args.begin, args.end, args.stride)
        if sel.size == 0:
            sys.exit(f"{f}: no frames left after --begin/--end/--stride")
        w = None
        if weights_all is not None:
            if len(weights_all) != len(time):
                sys.exit(f"{args.weights}: {len(weights_all)} weights but {f.name} has {len(time)} frames")
            w = weights_all[sel]
        p, tot = periodic_histogram(angle[sel], args.nbins, w)
        names.append(f.stem)
        probs.append(p)
        n_frames.append(sel.size)
        n_eff.append(tot)

    if args.norm == "density":
        columns = [p / width for p in probs]
        norm_line = "density (1/deg): sum(column) * bin_width = 1"
    else:
        columns = probs
        norm_line = "probability: sum(column) = 1"

    header = [
        "dihedral histograms from gmx chi -all, written by dihedral_histo.py",
        f"normalisation: {norm_line}",
        f"range: [-180, 180) deg, periodic (+180 == -180)",
        f"nbins: {args.nbins}",
        f"bin_width: {width:g} deg",
        f"frames: {min(n_frames)}" + ("" if min(n_frames) == max(n_frames) else f"-{max(n_frames)}"),
        f"begin: {args.begin:g} ps  end: {'last' if args.end < 0 else f'{args.end:g} ps'}  stride: {args.stride}",
        f"weights: {args.weights if args.weights else 'none'}",
        f"input: {', '.join(str(p) for p in args.inputs)}",
        f"n_dihedrals: {len(names)}",
    ]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_table(out, centers, names, columns, header)

    if args.split is not None:
        split_dir = Path(args.split) if args.split else out.parent
        split_dir.mkdir(parents=True, exist_ok=True)
        ylabel = "density (1/deg)" if args.norm == "density" else "probability"
        for name, col in zip(names, columns):
            write_xvg(split_dir / f"{args.prefix}{name}.xvg", centers, col, name, ylabel, header)
        print(f"{len(names)} files {split_dir / (args.prefix + '*.xvg')}")

    print(f"{len(names)} dihedrals, {args.nbins} bins of {width:g} deg, {min(n_frames)} frames -> {out}")


# --------------------------------------------------------------------------- #
# compare
# --------------------------------------------------------------------------- #
def jsd(p, q):
    """Jensen-Shannon divergence in bits, in [0, 1]."""
    m = 0.5 * (p + q)

    def kl(a, b):
        nz = a > 0
        return np.sum(a[nz] * np.log2(a[nz] / b[nz]))

    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def hellinger(p, q):
    return np.sqrt(max(0.0, 1.0 - np.sum(np.sqrt(p * q))))


def overlap(p, q):
    return np.sum(np.minimum(p, q))


def make_labels(paths, labels):
    if labels:
        if len(labels) != len(paths):
            sys.exit("--labels needs one label per table")
        return labels
    stems = [p.stem for p in paths]
    if len(set(stems)) == len(stems):
        return stems
    return [f"{p.parent.name}/{p.stem}" if p.parent.name else p.stem for p in paths]


def run_compare(args):
    paths = [Path(p) for p in args.tables]
    labels = make_labels(paths, args.labels)
    tables = [read_table(p) for p in paths]

    ref_centers, ref_hist, _ = tables[0]
    for p, (c, _, _) in zip(paths[1:], tables[1:]):
        if len(c) != len(ref_centers) or not np.allclose(c, ref_centers):
            sys.exit(f"{p} and {paths[0]} use different bins: rerun histo with the same --nbins")

    common = set(ref_hist)
    for p, (_, h, _) in zip(paths[1:], tables[1:]):
        missing = sorted(set(ref_hist) ^ set(h), key=dihedral_key)
        if missing:
            print(f"warning: {len(missing)} dihedrals not in both {paths[0].name} and {p.name}: {' '.join(missing)}", file=sys.stderr)
        common &= set(h)
    common = sorted(common, key=dihedral_key)
    if not common:
        sys.exit("no dihedrals in common")

    metrics = {"jsd": jsd, "hellinger": hellinger, "overlap": overlap}
    cols, col_names = [], []
    for lab, (_, h, _) in zip(labels[1:], tables[1:]):
        for mname, fn in metrics.items():
            cols.append([fn(ref_hist[d], h[d]) for d in common])
            col_names.append(f"{mname}:{lab}")

    with open(args.out, "w") as f:
        f.write(f"# reference: {labels[0]} ({paths[0]})\n")
        for lab, p in zip(labels[1:], paths[1:]):
            f.write(f"# compared: {lab} ({p})\n")
        f.write("# jsd: Jensen-Shannon divergence, bits, 0 = identical, 1 = disjoint\n")
        f.write("# hellinger: Hellinger distance, 0 = identical, 1 = disjoint\n")
        f.write("# overlap: sum_i min(p_i, q_i), 1 = identical, 0 = disjoint\n")
        f.write("# rows sorted by residue number, then phi, psi, omega, chi1, chi2, ...\n")
        f.write("resnr dihedral " + " ".join(col_names) + "\n")
        for i, d in enumerate(common):
            resnr = dihedral_key(d)[0]
            f.write(f"{resnr:>5d} {d:<14s} " + " ".join(f"{c[i]:.6f}" for c in cols) + "\n")

    print(f"{len(common)} dihedrals compared against {labels[0]} -> {args.out}")
    for lab, (_, h, _) in zip(labels[1:], tables[1:]):
        j = np.array([jsd(ref_hist[d], h[d]) for d in common])
        worst = np.argsort(j)[::-1][: args.top]
        print(f"  {lab}: mean JSD {j.mean():.4f}, median {np.median(j):.4f}; largest: "
              + ", ".join(f"{common[k]} {j[k]:.3f}" for k in worst))


# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Periodic, normalised histograms of gmx chi -all dihedrals.")
    sub = parser.add_subparsers(dest="command", required=True)

    h = sub.add_parser("histo", help="histogram all dihedrals of one simulation")
    h.add_argument("inputs", nargs="+", help="directory with the gmx chi -all xvg files, or the files themselves")
    h.add_argument("--out", default="dihedral_histograms.dat", help="output table (default: dihedral_histograms.dat)")
    h.add_argument("--nbins", type=int, default=72, help="number of bins over 360 deg (default: 72, i.e. 5 deg)")
    h.add_argument("--norm", choices=["density", "probability"], default="density",
                   help="density: integral = 1 (1/deg); probability: sum over bins = 1 (default: density)")
    h.add_argument("--begin", type=float, default=0.0, help="first time to use, ps (default: 0)")
    h.add_argument("--end", type=float, default=-1.0, help="last time to use, ps (default: all)")
    h.add_argument("--stride", type=int, default=1, help="use every n-th frame (default: 1)")
    h.add_argument("--weights", type=Path, default=None,
                   help="per-frame weights, one per frame of the xvg files (last column is used)")
    h.add_argument("--rad", action="store_true", help="xvg angles are in radians (gmx chi -rad)")
    h.add_argument("--pattern", default="*.xvg", help="glob used inside directories (default: *.xvg)")
    h.add_argument("--split", nargs="?", const="", default=None, metavar="DIR",
                   help="also write one <prefix><dihedral>.xvg file per dihedral, in DIR (default: next to --out)")
    h.add_argument("--prefix", default="histo-", help="file name prefix for --split (default: histo-)")
    h.add_argument("-v", "--verbose", action="store_true")

    c = sub.add_parser("compare", help="compare histogram tables, first one is the reference")
    c.add_argument("tables", nargs="+", help="two or more tables written by `histo` with the same --nbins")
    c.add_argument("--labels", nargs="+", default=None, help="one label per table")
    c.add_argument("--out", default="dihedral_compare.dat", help="output table (default: dihedral_compare.dat)")
    c.add_argument("--top", type=int, default=5, help="how many most different dihedrals to print (default: 5)")

    args = parser.parse_args(argv)
    if args.command == "histo":
        if args.nbins < 1:
            parser.error("--nbins must be >= 1")
        if args.stride < 1:
            parser.error("--stride must be >= 1")
    if args.command == "compare" and len(args.tables) < 2:
        parser.error("compare needs at least two tables")
    return args


def main(argv=None):
    args = parse_args(argv)
    if args.command == "histo":
        run_histo(args)
    else:
        run_compare(args)


if __name__ == "__main__":
    main()
