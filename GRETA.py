from pymol import cmd


cmd.load("https://files.rcsb.org/download/2KB8.pdb")

cmd.get_distance(atom1='pk1', atom2='pk2', state=0)