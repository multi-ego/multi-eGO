# domains.py

Apply a domain mask to a multi-eGO intra-molecular contact matrix.

For each residue range supplied via `--dom_res`, contacts where **both** atoms
fall inside the range are marked as *learned* (`True`) in the output HDF5 file.
Contacts outside all ranges are marked as *not learned* (`False`).
Use `--invert` to flip this assignment.

## Usage

```
python domains.py --intra <matrix> --top <topology> --dom_res <ranges> [--out <dir>] [--invert]
```

## Arguments

| Argument | Required | Description |
|---|---|---|
| `--intra` | yes | Path to the intra-molecular contact matrix (`.h5` HDF5 or plain text) |
| `--top` | yes | Path to the reference topology file (e.g. `topol.top`) |
| `--dom_res` | yes | One or more residue ranges as `start-end` pairs (see below) |
| `--out` | no | Output directory (default: current directory) |
| `--invert` | no | Invert the mask: intra-domain contacts become *not learned* |

### `--dom_res` rules

Ranges are given as `start-end` pairs and must satisfy:
- Each range is non-decreasing (`start ≤ end`).
- Ranges are non-overlapping and in increasing order.

```
--dom_res 1-30 45-60    # valid: two non-overlapping ranges
--dom_res 4-4           # valid: single residue
--dom_res 1-30 25-60    # invalid: overlapping
--dom_res 30-20         # invalid: decreasing
```

## Output filename

The output file is written to `--out` with the name:

```
split_<ranges>_<input_stem>.h5          # normal mode
inverted_split_<ranges>_<input_stem>.h5  # with --invert
```

For example, running with `--dom_res 1-30 45-60` on `intramat_1_1.ndx.h5` produces:

```
split_1-30-45-60_intramat_1_1.h5
```

## Examples

**Mark contacts within two domains as learned:**
```
python domains.py \
  --intra reference/intramat_1_1.ndx.h5 \
  --top   reference/topol.top \
  --dom_res 1-30 45-60 \
  --out   output/
```
```
Topology: 80 residues, 623 atoms
  Domain range: 1-30
    Atom index range (1-based): 1 – 234
    Number of atoms in range:   234
    First / last atom:  MET1 – ALA30
  Domain range: 45-60
    Atom index range (1-based): 352 – 469
    Number of atoms in range:   118
    First / last atom:  GLY45 – VAL60
Written: output/split_1-30-45-60_intramat_1_1.h5
```

**Invert — mark everything *outside* the domain as learned:**
```
python domains.py \
  --intra reference/intramat_1_1.ndx.h5 \
  --top   reference/topol.top \
  --dom_res 1-30 45-60 \
  --out   output/ \
  --invert
```
Output file: `output/inverted_split_1-30-45-60_intramat_1_1.h5`
