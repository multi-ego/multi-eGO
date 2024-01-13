# cmdata

`cmdata` is a GROMACS tool that implements the contact analysis that is required to get the multi-eGO imput files
it works with GROMACS 2021-2023.

## Installation

To install it you need to patch your GROMACS source code and recompile it.
```
./patch_gromacs.sh $GROMACS_ROOT_DIR
```

## Usage
cmdata calculates a histogram for interatomic distances for multi-eGO.

```
gmx cmdata [-f [<.xtc/.trr/...>]] [-s [<.tpr/.gro/...>]] [-n [<.ndx>]]
           [-b <time>] [-e <time>] [-dt <time>] [-tu <enum>]
           [-fgroup <selection>] [-xvg <enum>] [-[no]rmpbc] [-[no]pbc]
           [-sf <file>] [-selrpos <enum>] [-seltype <enum>] [-cutoff <real>]
           [-mol_cutoff <real>] [-ointra <string>] [-ointer <string>]
           [-nskip <int>] [-reference <selection>] [-sym <file>]
           [-[no]write_sym] [-num_threads <int>] [-mode <string>]
           [-weights <string>]
```

OPTIONS

Options to specify input files:
```
 -f      [<.xtc/.trr/...>]  (traj.xtc)       (Opt.)
           Input trajectory or single configuration: xtc trr cpt gro g96 pdb
           tng
 -s      [<.tpr/.gro/...>]  (topol.tpr)      (Opt.)
           Input structure: tpr gro g96 pdb brk ent
 -n      [<.ndx>]           (index.ndx)      (Opt.)
           Extra index groups

Other options:

 -b      <time>             (0)
           First frame (ps) to read from trajectory
 -e      <time>             (0)
           Last frame (ps) to read from trajectory
 -dt     <time>             (0)
           Only use frame if t MOD dt == first time (ps)
 -tu     <enum>             (ps)
           Unit for time values: fs, ps, ns, us, ms, s
 -fgroup <selection>
           Atoms stored in the trajectory file (if not set, assume first N
           atoms)
 -xvg    <enum>             (xmgrace)
           Plot formatting: xmgrace, xmgr, none
 -[no]rmpbc                 (yes)
           Make molecules whole for each frame
 -[no]pbc                   (yes)
           Use periodic boundary conditions for distance calculation
 -sf     <file>
           Provide selections from files
 -selrpos <enum>            (atom)
           Selection reference positions: atom, res_com, res_cog, mol_com,
           mol_cog, whole_res_com, whole_res_cog, whole_mol_com,
           whole_mol_cog, part_res_com, part_res_cog, part_mol_com,
           part_mol_cog, dyn_res_com, dyn_res_cog, dyn_mol_com, dyn_mol_cog
 -seltype <enum>            (atom)
           Default selection output positions: atom, res_com, res_cog,
           mol_com, mol_cog, whole_res_com, whole_res_cog, whole_mol_com,
           whole_mol_cog, part_res_com, part_res_cog, part_mol_com,
           part_mol_cog, dyn_res_com, dyn_res_cog, dyn_mol_com, dyn_mol_cog
 -cutoff <real>             (0.75)
           Cutoff in which to consider contacts
 -mol_cutoff <real>         (6)
           Molecular cutoff in which to consider contacts intermolecularly
 -ointra <string>           (intramat.ndx)
           Output of the intra-molecular contacts
 -ointer <string>           (intramat.ndx)
           Output of the intra-molecular contacts
 -nskip  <int>              (0)
           Use every nskip-th frame
 -reference <selection>
           Groups to calculate distances to
 -sym    <file>                              (Opt.)
           Atoms symmetry file path
 -[no]write_sym             (no)
           Write symmetry list
 -num_threads <int>         (1)
           Number of threads to run task on
 -mode   <string>           (intra+same+cross)
           Select the mode of analysis from intra, same or cross (combine with
           +)
 -weights <string>
           Test
```
