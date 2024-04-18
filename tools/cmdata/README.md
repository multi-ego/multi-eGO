# cmdata

`cmdata` is a tool that calculates interatomic distances and generates probability distribution and CDF used in multi-eGO to generate the input matrices. It works using the GROMACS API.

## Installation

To install it you need to compile and install GROMACS (2023/2024) using `-DGMX_INSTALL_LEGACY_API=ON -DBUILD_SHARED_LIBS=ON 
-DGMX_INSTALL_NBLIB_API=ON`. Then in a local subfolder (e.g., `build`, run `cmake ..` then `make` and `make install`).

## Usage
cmdata calculates a histogram for interatomic distances for multi-eGO.

```
Usage: cmdata [OPTION...]
  -f, --traj=FILE               Trajectory file
  -s, --top=FILE                Topology file
  -b, --t_begin[=FLOAT]         Start time
  -e, --t_end[=FLOAT]           End time
  -o, --out[=STRING]            Output prefix
      --dt[=INT]                Time step
      --cutoff[=DOUBLE]         Cutoff distance
      --mol_cutoff[=DOUBLE]     Molecule cutoff distance
      --nskip[=INT]             Number of frames to skip
      --num_threads[=INT]       Number of threads
      --mode[=STRING]           Mode of operation
      --weights[=FILE]          Weights file
      --no_pbc                  Ignore pbcs

Help options:
  -?, --help                    Show this help message
      --usage                   Display brief usage message
```
