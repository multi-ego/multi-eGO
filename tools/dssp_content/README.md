# dssp_content.py

Per-residue secondary structure content from the `dssp.dat` file written by `gmx dssp` (GROMACS >= 2023).

`dssp.dat` has one line per frame and one character per residue (`=` marks chain breaks).
For every residue the tool reports the fraction of analysed frames spent in each DSSP state,
and in three grouped classes:

| group  | DSSP states                                   |
|--------|-----------------------------------------------|
| helix  | H (alpha), G (3_10), I (pi)                   |
| strand | E (beta-strand), B (beta-bridge)              |
| coil   | P (PPII), T (turn), S (bend), ~ (loop)        |

Only Python standard library is required.

## Usage

```
gmx dssp -f traj.xtc -s topol.tpr -o dssp.dat
python dssp_content.py --dssp dssp.dat --out ss_content.dat
```

## Arguments

| argument              | default          | description                                          |
|-----------------------|------------------|------------------------------------------------------|
| `--dssp`              | (required)       | `dssp.dat` from `gmx dssp -o`                        |
| `--out`               | `ss_content.dat` | output text file                                     |
| `--begin`             | 0                | first frame to analyse (0-based)                     |
| `--end`               | all              | last frame (exclusive)                               |
| `--stride`            | 1                | analyse every n-th frame                             |
| `--first-res`         | 1                | number assigned to the first residue                 |
| `--restart-numbering` | off              | restart numbering at `--first-res` in each chain     |

## Output

Tab-separated, `#` comment header, one row per residue:

```
chain  resid  helix  strand  coil  H  G  I  P  E  B  T  S  loop
1      1      0.0000 0.0000  1.0000 ...
```

Unrecognised characters are ignored (with a warning) and are not counted in any column.
