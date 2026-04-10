# fasta2pdb.py 

A simple command-line tool to generate an extended peptide structure from an amino acid sequence and export it as a PDB file.

## Features

- Generates peptide structures using `PeptideBuilder`
- Outputs standard PDB files via `Bio.PDB`
- Lightweight and easy to integrate into workflows

## Requirements

- Python ≥ 3.8
- PeptideBuilder
- Biopython

Install dependencies:

pip install PeptideBuilder biopython

## Usage

python peptide_from_fasta.py --sequence ACDEFGHIK --output peptide.pdb

## Arguments

Argument	Required	Description
--sequence	Yes	Amino acid sequence (single-letter code)
--output	No	Output PDB file (default: output.pdb)

## Example

python peptide_from_fasta.py --sequence MKTFFVLLL

