# dihedral_histo.py

Normalised, periodic histograms of all dihedrals written by `gmx chi -all`, on a fixed grid,
so distributions from different simulations can be compared bin by bin.

- every dihedral uses the same grid over the full range [-180, 180) deg
- periodic: angles are wrapped, +180 and -180 fall in the same bin, the last bin closes on the first
- normalised: `density` (integral = 1, units 1/deg, default) or `probability` (sum over bins = 1)
- optional per-frame weights, time window and stride

Requires numpy.

## Usage

```
gmx chi -f traj.xtc -s topol.tpr -all -phi -psi -omega -maxchi 6    # writes phiALA2.xvg, chi1LYS12.xvg, ...

# one table per simulation: bin centres + one column per dihedral
python dihedral_histo.py histo simA/chi --out simA.dat
python dihedral_histo.py histo simB/chi --out simB.dat

# compare, first table is the reference
python dihedral_histo.py compare simA.dat simB.dat --labels ref test --out compare.dat
```

`histo` takes a directory (or a list of xvg files). Only files named like the `gmx chi -all`
output (`phi|psi|omega|chiN` + residue name + residue number) are used, so `order.xvg`,
`histo-*.xvg` and other files in the same directory are ignored.
Columns are ordered by residue number, then phi, psi, omega, chi1, chi2, ...

## `histo` arguments

| argument      | default                   | description                                                    |
|---------------|---------------------------|----------------------------------------------------------------|
| `inputs`      | (required)                | directory with the xvg files, or the files themselves          |
| `--out`       | `dihedral_histograms.dat` | output table                                                   |
| `--nbins`     | 72                        | bins over 360 deg (72 = 5 deg)                                 |
| `--norm`      | `density`                 | `density` (integral = 1) or `probability` (sum = 1)            |
| `--begin`     | 0                         | first time, ps                                                 |
| `--end`       | -1                        | last time, ps (-1 = all)                                       |
| `--stride`    | 1                         | use every n-th frame                                           |
| `--weights`   | none                      | per-frame weights, one line per frame (last column is used)    |
| `--rad`       | off                       | xvg angles are in radians (`gmx chi -rad`); output is in deg   |
| `--pattern`   | `*.xvg`                   | glob used inside directories                                   |
| `--split DIR` | off                       | also write `DIR/<dihedral>.dat` (centre, density, probability) |

Output header (`# ...` lines) records normalisation, bins, frames, time window and weights;
the `# columns:` line names the columns. Tables load with `numpy.loadtxt`.

## `compare`

Takes two or more `histo` tables made with the same `--nbins` (either normalisation) and,
for each dihedral present in all of them, compares each table with the first one:

| metric      | definition                    | identical | disjoint |
|-------------|-------------------------------|-----------|----------|
| `jsd`       | Jensen-Shannon divergence, bits | 0       | 1        |
| `hellinger` | Hellinger distance            | 0         | 1        |
| `overlap`   | sum_i min(p_i, q_i)           | 1         | 0        |

It also prints mean and median JSD and the `--top` most different dihedrals.
Dihedrals missing from one table are reported and skipped.
