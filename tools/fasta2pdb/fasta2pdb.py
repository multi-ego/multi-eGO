#!/usr/bin/env python3

"""
fasta2pdb.py

Generate an extended peptide structure from a FASTA sequence
and save it as a PDB file.
"""

import argparse
from pathlib import Path

import PeptideBuilder
from Bio.PDB import PDBIO


def build_structure(sequence: str):
    """Generate an extended peptide structure from an amino acid sequence."""
    normalized_sequence = "".join(sequence.split()).upper() if sequence is not None else ""
    if not normalized_sequence:
        raise ValueError("Input sequence must not be empty.")
    if not normalized_sequence.isalpha():
        raise ValueError("Input sequence must contain only amino acid letters.")

    valid_residues = set("ACDEFGHIKLMNPQRSTVWY")
    invalid_residues = sorted(set(normalized_sequence) - valid_residues)
    if invalid_residues:
        raise ValueError("Input sequence contains invalid amino acid residue(s): " + ", ".join(invalid_residues))

    return PeptideBuilder.make_extended_structure(normalized_sequence)


def save_pdb(structure, output_path: Path):
    """Save a structure to a PDB file."""
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path))


def parse_args():
    parser = argparse.ArgumentParser(description="Generate a PDB structure from a FASTA sequence.")
    parser.add_argument(
        "--sequence",
        type=str,
        required=True,
        help="Amino acid sequence (single-letter code).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output.pdb"),
        help="Output PDB file path (default: output.pdb)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    structure = build_structure(args.sequence)
    save_pdb(structure, args.output)

    print(f"PDB file successfully written to: {args.output}")


if __name__ == "__main__":
    main()
