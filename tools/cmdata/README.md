# cmdata

`cmdata` calculates interatomic distance histograms from GROMACS trajectories and outputs probability density and CDF files used by multi-eGO to build its input contact matrices.

## Installation

### Prerequisites

GROMACS 2023, 2024, or 2026+ must be compiled and installed with:

```bash
cmake .. -DGMX_INSTALL_LEGACY_API=ON \
         -DBUILD_SHARED_LIBS=ON \
         -DGMX_INSTALL_NBLIB_API=ON
make -j$(nproc) && make install
```

HDF5 is optional. If found by CMake, output files can be written in HDF5 format (`.h5`) in addition to plain text (`.dat`).

### Building cmdata

```bash
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install \
         -DGROMACS_DIR=/path/to/gromacs
make -j$(nproc) && make install
```

> **macOS note:** the build system automatically strips the stale `libm.tbd` reference that GROMACS exports in its CMake config — this stub was removed from the macOS 26+ SDK (math is now part of `libSystem`). No manual workaround is needed.

Set `GROMACS_DIR` to the root of your GROMACS installation if it is not in a standard location.

To build with tests:

```bash
cmake .. -DCMDATA_BUILD_TESTS=ON
```

## Usage

```
cmdata -f TRAJ -s TOP [OPTIONS]
```

### Required options

| Option | Description |
|--------|-------------|
| `-f`, `--traj=FILE` | Input trajectory file (`.xtc`, `.trr`, `.gro`, or `.pdb`) |
| `-s`, `--top=FILE`  | Input topology file (`.tpr`) |
| `--mode=STRING`     | Calculation mode (see below) |

### Optional options

| Option | Default | Description |
|--------|---------|-------------|
| `-o`, `--out=STRING`         | `""` (current dir) | Output file prefix |
| `-b`, `--t_begin=FLOAT`      | `0.0` | Start time (ps) |
| `-e`, `--t_end=FLOAT`        | `-1` (end of traj) | End time (ps); `-1` means read to the end |
| `--dt=INT`                   | `0` (every frame) | Only process frames at multiples of this time (ps) |
| `--cutoff=FLOAT`             | `0.75` | Distance cutoff for atom pairs (nm) |
| `--mol_cutoff=FLOAT`         | `6.0` | Centre-of-mass cutoff for molecule pairs (nm) |
| `--nskip=INT`                | `0` | Skip every N frames (0 = no skipping) |
| `--num_threads=INT`          | `1` | Number of threads for the mindist accumulation step |
| `--mol_threads=INT`          | `1` | Number of concurrent molecule threads |
| `--weights=FILE`             | `""` (none) | Per-frame weight file (one float per line) |
| `--bkbn_H=STRING`            | `""` (none) | Backbone hydrogen atom name to include (all other H atoms are skipped) |
| `--no_pbc`                   | off | Disable periodic boundary corrections |
| `--noh5`                     | off | Write plain-text `.dat` files instead of HDF5 `.h5` files (ignored if HDF5 was not found at build time) |

### Mode

`--mode` controls which molecule pairs are analysed. It is a `+`-separated combination of:

| Token | Description |
|-------|-------------|
| `intra` | Intra-molecular distances (atom pairs within the same molecule) |
| `same`  | Inter-molecular distances between molecules of the **same** species |
| `cross` | Inter-molecular distances between molecules of **different** species |

Examples: `--mode intra`, `--mode same+cross`, `--mode intra+same+cross`.

## Output files

Each output file is named `<prefix><type>_<i>_<j>_aa_<ii>.h5` (or `.dat`), where:

- `<type>` is `intra_mol`, `inter_mol`, or `inter_mol_c` (cross)
- `<i>`, `<j>` are zero-based molecule-type indices
- `<ii>` is the atom index within molecule type `i`

Each file contains one column per atom pair `(ii, jj)`. The rows are histogram bins from `dx/2` to `cutoff` in steps of `dx = cutoff / n_bins`.

For `same` and `cross` modes two files are written per atom: one for the full density histogram and one for the per-molecule maximum-CDF (`maxcdf`).

### Trajectory format notes

| Format | Progress bar | Notes |
|--------|-------------|-------|
| `.xtc` | Percentage | Seek table read at startup; full progress display |
| `.trr` | Frame count | No seek table; TRR files include velocities and forces (only positions used) |
| `.gro` | Frame count | Single-frame or multi-frame structure files |
| `.pdb` | Frame count | Multi-model PDB files supported |

For TRR, GRO, and PDB files the progress bar shows the running frame count rather than a percentage because the total number of frames cannot be determined cheaply without reading the whole file.

If you have a GRO or PDB single-structure file you want to treat as a one-frame trajectory, pass it directly — cmdata will read the single frame and process it normally.

## Examples

Compute intra-molecular histograms only, using 8 threads:

```bash
cmdata -f traj.xtc -s topol.tpr --mode intra --num_threads 8 -o output/
```

Compute all three modes with a 1.0 nm cutoff, skipping the first 10 ns:

```bash
cmdata -f traj.xtc -s topol.tpr \
       --mode intra+same+cross \
       --cutoff 1.0 \
       --t_begin 10000 \
       --num_threads 4 --mol_threads 4 \
       -o run1/
```

Use a per-frame weight file (e.g., from replica-exchange or metadynamics):

```bash
cmdata -f traj.xtc -s topol.tpr --mode intra+same \
       --weights weights.dat -o weighted/
```
