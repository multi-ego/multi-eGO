#!/usr/bin/env python3
"""
Numerical comparator for cmdata histogram output files.

Usage:
    compare_histo.py <file_a> <file_b> [--rtol RTOL] [--atol ATOL]

Exit codes:
    0  all values agree within tolerance
    1  at least one value differs beyond tolerance
    2  file shape mismatch (different number of lines or columns)

Tolerance check (numpy-compatible):
    |a - b| <= atol + rtol * max(|a|, |b|)

Defaults:  rtol=1e-4  atol=1e-7
"""

import sys
import argparse


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("file_a")
    p.add_argument("file_b")
    p.add_argument("--rtol", type=float, default=1e-4, help="relative tolerance (default: 1e-4)")
    p.add_argument("--atol", type=float, default=1e-7, help="absolute tolerance (default: 1e-7)")
    return p.parse_args()


def read_rows(path):
    rows = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                rows.append((lineno, [float(v) for v in line.split()]))
            except ValueError as exc:
                print(f"ERROR: {path}:{lineno}: cannot parse float — {exc}", file=sys.stderr)
                sys.exit(2)
    return rows


def compare(file_a, file_b, rtol, atol):
    rows_a = read_rows(file_a)
    rows_b = read_rows(file_b)

    if len(rows_a) != len(rows_b):
        print(f"FAIL  {file_a}: line count mismatch " f"({len(rows_a)} vs {len(rows_b)})", file=sys.stderr)
        return False

    failures = []
    for (lineno, vals_a), (_, vals_b) in zip(rows_a, rows_b):
        if len(vals_a) != len(vals_b):
            print(
                f"FAIL  {file_a}:{lineno}: column count mismatch " f"({len(vals_a)} vs {len(vals_b)})", file=sys.stderr
            )
            return False
        for col, (a, b) in enumerate(zip(vals_a, vals_b), 1):
            tol = atol + rtol * max(abs(a), abs(b))
            if abs(a - b) > tol:
                failures.append((lineno, col, a, b, abs(a - b), tol))

    if failures:
        for lineno, col, a, b, diff, tol in failures[:10]:  # show at most 10
            print(
                f"FAIL  {file_a}:{lineno} col {col}: " f"{a:.6g} vs {b:.6g}  diff={diff:.3e}  tol={tol:.3e}",
                file=sys.stderr,
            )
        if len(failures) > 10:
            print(f"      ... and {len(failures) - 10} more failures", file=sys.stderr)
        return False

    return True


def main():
    args = parse_args()
    ok = compare(args.file_a, args.file_b, args.rtol, args.atol)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
