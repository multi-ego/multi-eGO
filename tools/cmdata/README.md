# cmdata

`cmdata` calculates interatomic distance histograms from GROMACS trajectories and outputs probability density and cumulative distribution function (CDF) files used by multi-eGO to build its contact matrices.

## Installation

### Prerequisites

**GROMACS** (2023+) must be built with the legacy API and shared libraries enabled:

```bash
cmake .. -DGMX_INSTALL_LEGACY_API=ON \
         -DBUILD_SHARED_LIBS=ON \
         -DGMX_INSTALL_NBLIB_API=ON
make -j$(nproc) && make install
```

**OpenMP** is used for parallelism and is detected automatically at build time. On Linux it is provided by GCC or Clang out of the box. On macOS, Apple Clang does not include OpenMP; install it via Homebrew:

```bash
brew install libomp
```

then point CMake at it:

```bash
cmake .. -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp" \
         -DOpenMP_CXX_LIB_NAMES="omp" \
         -DOpenMP_omp_LIBRARY=$(brew --prefix libomp)/lib/libomp.dylib
```

If OpenMP is not found, cmdata builds and runs correctly but uses only one core.

**HDF5** is optional but recommend. If found by CMake, output files can be written in HDF5 format (`.h5`) in addition to plain text (`.dat`).

### Building cmdata

```bash
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install \
         -DGROMACS_DIR=/path/to/gromacs
make -j$(nproc) && make install
```

Set `GROMACS_DIR` to the root of your GROMACS installation if it is not in a standard location.

To build with unit tests (no GROMACS required for the tests themselves):

```bash
cmake .. -DCMDATA_BUILD_TESTS=ON
make && ctest
```

## Usage

```
cmdata -f TRAJ -s TOP [OPTIONS]
```

### Required options

| Option | Description |
|--------|-------------|
| `-f`, `--traj=FILE` | Input trajectory file (`.xtc`) |
| `-s`, `--top=FILE`  | Input topology file (`.tpr`) |
| `--mode=STRING`     | Calculation mode (see below) |

### Optional options

| Option | Default | Description |
|--------|---------|-------------|
| `-o`, `--out=STRING`         | `""` (current dir)    | Output file prefix |
| `-b`, `--t_begin=FLOAT`      | `0.0`                 | Start time (ps) |
| `-e`, `--t_end=FLOAT`        | `-1` (end of traj)    | End time (ps); `-1` means read to the end |
| `--dt=INT`                   | `0` (every frame)     | Only process frames at multiples of this time (ps) |
| `--cutoff=FLOAT`             | `0.75`                | Distance cutoff for atom pairs (nm) |
| `--mol_cutoff=FLOAT`         | `6.0`                 | Centre-of-mass cutoff for molecule pairs (nm) |
| `--nskip=INT`                | `0`                   | Skip every N frames (0 = no skipping) |
| `--weights=FILE`             | `""` (none)           | Per-frame weight file (one float per line) |
| `--bkbn_H=STRING`            | `""` (none)           | Backbone hydrogen atom name to include (all other H are skipped) |
| `--no_pbc`                   | off                   | Disable periodic boundary corrections |
| `--noh5`                     | off                   | Write plain-text `.dat` files instead of HDF5 `.h5` (ignored if HDF5 was not found at build time) |

### Mode

`--mode` controls which molecule pairs are analysed. It is a `+`-separated combination of:

| Token | Description |
|-------|-------------|
| `intra` | Intra-molecular distances (atom pairs within the same molecule) |
| `same`  | Inter-molecular distances between molecules of the **same** species |
| `cross` | Inter-molecular distances between molecules of **different** species |

Examples: `--mode intra`, `--mode same+cross`, `--mode intra+same+cross`.

### Parallelism

cmdata uses **OpenMP** to parallelise the per-molecule pairwise distance loop. The number of threads is controlled at runtime via the standard `OMP_NUM_THREADS` environment variable — no build-time flags are needed:

```bash
OMP_NUM_THREADS=8 cmdata -f traj.xtc -s topol.tpr --mode intra+same
```

If `OMP_NUM_THREADS` is not set, OpenMP defaults to using all available cores. Set it to `1` to force single-threaded execution.

## Output files

Each output file is named `<prefix><type>_<i>_<j>_aa_<ii>.dat` (or `.h5`), where:

- `<type>` is `intra_mol`, `inter_mol`, or `inter_mol_c`
- `<i>`, `<j>` are 1-based molecule-type indices
- `<ii>` is the 1-based atom index within molecule type `i`

Each file contains one column per atom in molecule type `j`. The rows correspond to histogram bins uniformly spaced from `dx/2` to `cutoff`, with `dx = cutoff / n_bins`.

For `same` and `cross` modes, two files are written per atom:

- `inter_mol_*` — probability density histogram
- `inter_mol_c_*` — per-molecule maximum CDF (used directly by multi-eGO)

## Examples

Compute intra-molecular histograms only, using 8 threads:

```bash
OMP_NUM_THREADS=8 cmdata -f traj.xtc -s topol.tpr --mode intra -o output/
```

Compute all three modes with a 1.0 nm cutoff, skipping the first 10 ns:

```bash
OMP_NUM_THREADS=8 cmdata -f traj.xtc -s topol.tpr \
    --mode intra+same+cross \
    --cutoff 1.0 \
    --t_begin 10000 \
    -o run1/
```

Use a per-frame weight file (e.g. from replica-exchange or metadynamics):

```bash
OMP_NUM_THREADS=8 cmdata -f traj.xtc -s topol.tpr \
    --mode intra+same \
    --weights weights.dat \
    -o weighted/
```
