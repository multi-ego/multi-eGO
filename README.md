# Multi-*e*GO: data-driven force fields for molecular dynamics simulations

[![Version](https://img.shields.io/badge/Version-1.0-blue)](https://github.com/multi-ego/multi-eGO/releases)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Multi-eGO test](https://github.com/multi-ego/multi-eGO/actions/workflows/test.yml/badge.svg)](https://github.com/multi-ego/multi-eGO/actions/workflows/test.yml)
[![cmdata](https://github.com/multi-ego/multi-eGO/actions/workflows/cmdata.yml/badge.svg)](https://github.com/multi-ego/multi-eGO/actions/workflows/cmdata.yml)
[![CodeQL](https://github.com/multi-ego/multi-eGO/actions/workflows/github-code-scanning/codeql/badge.svg)](https://github.com/multi-ego/multi-eGO/actions/workflows/github-code-scanning/codeql)

Multi-*e*GO is a framework for building data-driven, atomic-resolution force fields for molecular dynamics simulations. It learns pairwise Lennard-Jones interactions directly from contact probability distributions, producing force fields that faithfully reproduce structural ensembles, thermodynamics and kinetics — including protein folding, aggregation, and ligand binding — at a fraction of the computational cost of all-atom simulations.

## Table of Contents

- [How it works](#how-it-works)
- [Requirements](#requirements)
- [Installation](#installation)
- [Workflow overview](#workflow-overview)
- [Step-by-step guide](#step-by-step-guide)
  - [1. Prepare the system topology](#1-prepare-the-system-topology)
  - [2. Analyse a training simulation](#2-analyse-a-training-simulation)
  - [3. Generate the molten-globule prior (`mg`)](#3-generate-the-molten-globule-prior-mg)
  - [4. Run and analyse the reference simulation](#4-run-and-analyse-the-reference-simulation)
  - [5. Generate the production force field](#5-generate-the-production-force-field)
- [Advanced usage](#advanced-usage)
  - [Multiple training simulations via config file](#multiple-training-simulations-via-config-file)
  - [Key parameters](#key-parameters)
- [Developers](#developers)
- [Cite us](#cite-us)

---

## How it works

Multi-*e*GO constructs an atomistic Lennard-Jones force field in two stages:

1. **Molten-globule (`mg`) prior** — a baseline transferable force field derived from a combination of physico-chemical and statistical properties of proteins.
2. **Production force field** — the prior is refined following Bayes, learning contact probabilities from one or more training dataset (obtained from structures, MD simulations, conformational ensemble resulting from integrative/AI approches).

The result is a force field encoded in standard GROMACS `ffnonbonded.itp` and `topol_mego.top` files, compatible with any GROMACS workflow.

---

## Requirements

- Python ≥ 3.10
- [GROMACS](https://www.gromacs.org) > 2022 (required for running simulations)
- The `cmdata` tool for extracting contact histograms from trajectories — this requires recompiling GROMACS from source; see [its README](tools/cmdata/README.md) for instructions

---

## Installation

Create and activate the conda environment:

```bash
conda env create -f conda/environment.yml
conda activate meGO
```

Alternatively, with pip:

```bash
pip install -r requirements.txt
```

For the `cmdata` trajectory analysis tool, follow the separate [installation instructions](tools/cmdata/README.md).

---

## Workflow overview

![Multi-eGO workflow](img/mego_workflow_black.png)

The full workflow from a PDB structure to a production simulation involves five steps:

1. Prepare a GROMACS topology using the `multi-ego-basic.ff` force field
2. Run a training simulation with your all-atom force field and extract contact data
3. Generate the `mg` prior force field and run a reference (molten-globule) simulation
4. Analyse the reference simulation to obtain contact probabilities
5. Generate the production force field by learning from the training data against the reference

---

## Step-by-step guide

### 1. Prepare the system topology

Copy your PDB file and the `multi-ego-basic.ff/` directory into a working folder, then generate a GROMACS topology:

```bash
gmx pdb2gmx -f your_structure.pdb -ignh
```

Select the `multi-ego-basic` force field when prompted. This produces a `.gro` coordinate file and a `topol.top` topology.

Create the input directory structure:

```
inputs/
└── SYSTEM_NAME/
    └── topol.top
```

### 2. Analyse a training simulation

Given an existing all-atom MD trajectory, first extract contact histograms using `cmdata`:

```bash
cmdata -f trajectory.xtc -s topology.tpr
```

Then convert the histograms into contact matrices using `make_mat.py`:

```bash
python tools/make_mat/make_mat.py \
    --histo     MD_DIRECTORY/histo/ \
    --target_top MD_DIRECTORY/topol.top \
    --mego_top   inputs/SYSTEM_NAME/topol.top \
    --out        inputs/SYSTEM_NAME/md_ensemble
```

Copy the topology and force field into the training folder:

```
inputs/
└── SYSTEM_NAME/
    ├── topol.top
    └── md_ensemble/
        ├── topol.top
        ├── intramat_1_1.ndx
        └── all-atom.ff/
```

> **Note:**  The `moleculetype` name in all topology files must be consistent. The program will exit with an error if names do not match.

### 3. Generate the molten-globule prior (`mg`)

The `mg` prior provides generic local interactions and is used as the starting point for the reference simulation. Generate it with:

```bash
python multiego.py --system SYSTEM_NAME --egos mg
```
> **Note:** before running, ensure that the `moleculetype` name in all topology files is consistent. The program will exit with an error if names do not match.

Output files are written to `outputs/SYSTEM_NAME/mg_#/` and contain:

- `ffnonbonded.itp` — copy into your simulation directory
- `topol_mego.top` — copy into your simulation directory
- `meGO.log` — a log including the detail of the generated force field 

Run the reference simulation using the provided `mdp` files in order:

```
1. ff_em.mdp       (energy minimisation)
2. ff_cg.mdp       (coarse-grained equilibration)
3. ff_aa-posre.mdp (restrained all-atom equilibration)
4. ff_aa.mdp       (production reference run)
```

### 4. Run and analyse the reference simulation

Once the reference simulation is complete, extract its contact data exactly as in step 2, but write the output into the `reference/` folder:

```bash
cmdata -f reference_trajectory.xtc -s reference_topology.tpr

python tools/make_mat/make_mat.py \
    --histo      RC_DIRECTORY/histo/ \
    --target_top RC_DIRECTORY/topol.top \
    --mego_top   inputs/SYSTEM_NAME/topol.top \
    --out        inputs/SYSTEM_NAME/reference
```

The final input directory should look like:

```
inputs/
└── SYSTEM_NAME/
    ├── topol.top
    ├── reference/
    │   ├── topol.top
    │   ├── intramat_1_1.ndx
    │   ├── ffnonbonded.itp
    │   └── topol_mego.top 
    └── md_ensemble/
        ├── topol.top
        ├── intramat_1_1.ndx
        └── all-atom.ff/
```

### 5. Generate the production force field

```bash
python multiego.py \
    --system SYSTEM_NAME \
    --egos production \
    --epsilon 0.3 \
    --train md_ensemble
```

`--epsilon` sets the maximum interaction energy (in kJ/mol) for any learned contact pair; 0.3 kJ/mol is a reasonable starting value for most proteins.

Output is written to `outputs/SYSTEM_NAME/production_#/`. Copy `ffnonbonded.itp` and `topol_mego.top` to your simulation directory as before, and run using `ff_aa.mdp` as the production mdp file.

> **Note:** before running, ensure that the `moleculetype` name in all topology files is consistent. The program will exit with an error if names do not match.

---

## Advanced usage

### Multiple training simulations via config file

When learning from multiple training simulations or multiple reference conditions simultaneously, use a YAML configuration file instead of command-line arguments:

```bash
python multiego.py --config config.yml
```

A typical config file looks like:

```yaml
- system: SYSTEM_NAME
- egos: production
- symmetry: 
  - ARG NH1 NH2
  - ASP OD1 OD2
  - GLU OE1 OE2
  - PHE CD1 CD2
  - PHE CE1 CE2
  - TYR CD1 CD2
  - TYR CE1 CE2
  - ALA O1 O2
- input_refs:
  - reference: reference_a
    train: native_MD
    matrix: intramat_1_1
    epsilon: 0.25
  - reference: reference_b
    train: fibril
    matrix: intramat_1_1
    epsilon: 0.25
  - reference: reference_a
    train: fibril
    matrix: intermat_1_1
    epsilon: 0.25
```

Each entry under `input_refs` defines one training/reference pair. Intra- and inter-molecular contacts can be learned simultaneously by listing both matrix types.

### Key parameters

| Parameter | Default | Description |
|---|---|---|
| `--epsilon` | — | Maximum interaction energy per contact pair (kJ/mol). Required for `production`. |
| `--epsilon_min` | 0.07 | Minimum meaningful epsilon; contacts below this are discarded. |
| `--p_to_learn` | 0.9995 | Fraction of the cumulative contact probability used to set the learning threshold. Should be close to 1. |
| `--relative_c12d` | 0.01 | Tolerance for removing contacts that are effectively identical to the MG prior. |
| `--force_split` | False | Write intra- and inter-molecular interactions to separate files. |
| `--single_molecule` | False | Enable optimisations valid when simulating a single molecule in isolation. |
| `--symmetry` | — | List of equivalent atoms, e.g. `ARG NH1 NH2`. Contacts are averaged over equivalent atoms. |
| `--custom_dict` | — | JSON file for custom atom-name mappings (non-standard residues). |
| `--custom_c12` | — | CSV file for custom repulsive C12 parameters. |

---

## Developers

- Fran Bacic Toplek
- Carlo Camilloni
- Riccardo Capelli
- Bruno Stegani

---

## Cite us

If you use multi-*e*GO in your work, please cite the relevant papers:

### Method

1. Scalone, E., et al. [Multi-eGO: An in silico lens to look into protein aggregation kinetics at atomic resolution.](https://www.pnas.org/doi/10.1073/pnas.2203181119) *Proc. Natl. Acad. Sci. USA* 119, e2203181119 (2022). ([preprint](https://www.biorxiv.org/content/10.1101/2022.02.18.481033v2))

2. Bacic Toplek, F., Scalone, E., et al. [Multi-eGO: model improvements towards the study of complex self-assembly processes.](https://doi.org/10.1021/acs.jctc.3c01182) *J. Chem. Theory Comput.* 20, 459–468 (2024). ([preprint](https://doi.org/10.26434/chemrxiv-2023-67255))

### Applications

3. Pennacchietti, V., et al. [A PDZ tandem repeat folds and unfolds via different pathways.](https://dx.doi.org/10.1002/pro.5203) *Protein Sci.* 33, e5203 (2024).

4. Stegani, B., et al. [Estimation of Ligand Binding Free Energy Using Multi-eGO.](https://doi.org/10.1021/acs.jcim.4c01545) *J. Chem. Inf. Model.* 65, 351–362 (2025). ([preprint](https://doi.org/10.26434/chemrxiv-2024-jcmgc-v2))

5. Stegani, B. and Capelli, R. [Kinetic rates calculation via non-equilibrium dynamics.](https://doi.org/10.1063/5.0277524) *J. Chem. Phys.* 163, 104103 (2025).

6. Kulshrestha, A., et al. [iMultiscale Simulations Elucidate the Mechanism of Polyglutamine Aggregation and the Role of Flanking Domains in Fibril Polymorphism](https://pubs.acs.org/doi/10.1021/acs.jpcb.5c06627) *J. Phys. Chem. B* 129, 43, 11205 (2025)
