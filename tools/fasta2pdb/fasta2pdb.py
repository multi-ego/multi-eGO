import PeptideBuilder
from Bio.PDB import PDBIO
import argparse

parser = argparse.ArgumentParser(description="Generate PDB files from a fasta file")
parser.add_argument("--fasta", type=str, help="Fasta sequence as input")
parser.add_argument("--output", type=str, default="./output.pdb", help="Path to the output PDB file")
args = parser.parse_args()

structure = PeptideBuilder.make_extended_structure(args.fasta)
io = PDBIO()
io.set_structure(structure)
io.save(args.output)
