#!/usr/bin/env python3

"""
dssp_content.py

Compute the per-residue secondary structure content from the dssp.dat file
written by `gmx dssp -o dssp.dat` (GROMACS >= 2023).

Each line of dssp.dat is one frame: a string with one character per residue,
with '=' marking chain breaks. The output is a text table with, for each
residue, the fraction of frames spent in each DSSP state and in the grouped
helix / strand / coil classes.
"""

import argparse
import sys
from pathlib import Path

# DSSP states as written by gmx dssp
DSSP_STATES = {
    "H": "alpha-helix",
    "G": "3_10-helix",
    "I": "pi-helix",
    "P": "PPII-helix",
    "E": "beta-strand",
    "B": "beta-bridge",
    "T": "turn",
    "S": "bend",
    "~": "loop",
}
GROUPS = {
    "helix": "HGI",
    "strand": "EB",
    "coil": "PTS~",
}
CHAIN_BREAK = "="


def read_dssp(path: Path, begin: int = 0, end: int = -1, stride: int = 1):
    """Read dssp.dat and return the list of frames (strings) to analyse."""
    with open(path) as f:
        frames = [line.strip() for line in f if line.strip() and not line.startswith(("#", "@"))]
    if not frames:
        raise ValueError(f"No frames found in {path}")
    stop = None if end < 0 else end
    frames = frames[slice(begin, stop, stride)]
    if not frames:
        raise ValueError("No frames left after applying --begin/--end/--stride")
    lengths = {len(fr) for fr in frames}
    if len(lengths) != 1:
        raise ValueError(f"Frames have different lengths {sorted(lengths)}: is {path} a gmx dssp output?")
    return frames


def residue_map(frame: str):
    """Return (string position, chain index) for each residue, skipping chain-break characters."""
    mapping = []
    chain = 1
    for pos, ch in enumerate(frame):
        if ch == CHAIN_BREAK:
            chain += 1
            continue
        mapping.append((pos, chain))
    return mapping


def compute_content(frames, first_res: int = 1, restart: bool = False):
    """Return a list of per-residue dicts with the fraction of each state."""
    nframes = len(frames)
    mapping = residue_map(frames[0])
    unknown = set()
    rows = []
    resid = first_res
    prev_chain = 1
    for pos, chain in mapping:
        if chain != prev_chain:
            prev_chain = chain
            if restart:
                resid = first_res
        counts = dict.fromkeys(DSSP_STATES, 0)
        for fr in frames:
            ch = fr[pos]
            if ch in counts:
                counts[ch] += 1
            else:
                unknown.add(ch)
        row = {"chain": chain, "resid": resid}
        row.update({s: counts[s] / nframes for s in DSSP_STATES})
        row.update({g: sum(counts[s] for s in states) / nframes for g, states in GROUPS.items()})
        rows.append(row)
        resid += 1
    if unknown:
        print(f"WARNING: unrecognised DSSP characters ignored: {''.join(sorted(unknown))}", file=sys.stderr)
    return rows


def write_content(rows, out: Path, nframes: int, source: Path):
    cols = ["chain", "resid"] + list(GROUPS) + list(DSSP_STATES)
    with open(out, "w") as f:
        f.write(f"# per-residue secondary structure content from {source} ({nframes} frames)\n")
        f.write("# helix = H+G+I ; strand = E+B ; coil = P+T+S+~\n")
        f.write("# " + "  ".join(f"{s}={n}" for s, n in DSSP_STATES.items()) + "\n")
        f.write("\t".join("loop" if c == "~" else c for c in cols) + "\n")
        for r in rows:
            f.write(f"{r['chain']}\t{r['resid']}\t" + "\t".join(f"{r[c]:.4f}" for c in cols[2:]) + "\n")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Per-residue secondary structure content from gmx dssp output.")
    parser.add_argument("--dssp", type=Path, required=True, help="dssp.dat file from gmx dssp -o")
    parser.add_argument(
        "--out", type=Path, default=Path("ss_content.dat"), help="Output text file (default: ss_content.dat)"
    )
    parser.add_argument("--begin", type=int, default=0, help="First frame to analyse (0-based, default: 0)")
    parser.add_argument("--end", type=int, default=-1, help="Last frame (exclusive, default: all)")
    parser.add_argument("--stride", type=int, default=1, help="Analyse every n-th frame (default: 1)")
    parser.add_argument("--first-res", type=int, default=1, help="Number of the first residue (default: 1)")
    parser.add_argument(
        "--restart-numbering", action="store_true", help="Restart residue numbering at --first-res in each chain"
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.stride < 1:
        sys.exit("--stride must be >= 1")
    frames = read_dssp(args.dssp, args.begin, args.end, args.stride)
    rows = compute_content(frames, args.first_res, args.restart_numbering)
    write_content(rows, args.out, len(frames), args.dssp)
    nchains = rows[-1]["chain"] if rows else 0
    print(f"{len(frames)} frames, {len(rows)} residues, {nchains} chain(s) -> {args.out}")


if __name__ == "__main__":
    main()
