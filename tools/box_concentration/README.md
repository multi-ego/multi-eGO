# get_box.py

Interconvert number of molecules, molar concentration, and box geometry.
Provide any two of the three quantities and the third is computed.

```
python get_box.py --n_mol <N> --conc <C>
python get_box.py --n_mol <N> --volume <V>
python get_box.py --conc <C>  --volume <V>
python get_box.py --conc <C>  --sphere_r <R>
python get_box.py --conc <C>  --cubic_side <L>
```

## Arguments

| Argument | Unit | Description |
|---|---|---|
| `--n_mol` | — | Number of molecules |
| `--conc` | M | Molar concentration |
| `--volume` | nm³ | Box volume |
| `--sphere_r` | nm | Radius of a sphere (converted to volume) |
| `--cubic_side` | nm | Side length of a cube (converted to volume) |

`--volume`, `--sphere_r`, and `--cubic_side` are mutually exclusive.
Passing a geometry argument alone prints the corresponding volume without
requiring `--n_mol` or `--conc`.

## Examples

**Cubic box side from molecule count and concentration**
```
$ python get_box.py --n_mol 50 --conc 0.1
cubic box side: 9.39881 nm
```

**Concentration from molecule count and box volume**
```
$ python get_box.py --n_mol 10 --volume 1000.0
concentration:  0.016605 M
```

**Number of molecules from concentration and sphere radius**
```
$ python get_box.py --conc 0.1 --sphere_r 5.0
sphere volume: 523.5988 nm³  (cubic-equivalent side: 8.05996 nm)
n_mol:          32
```

**Volume of a cube (geometry only)**
```
$ python get_box.py --cubic_side 8.0
cubic volume:  512.0000 nm³
```
