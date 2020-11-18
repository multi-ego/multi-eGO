from pymol import cmd


cmd.load("https://files.rcsb.org/download/2KB8.pdb")


prova = cmd.select(name = 'prova', selection="name CA")

for at1 in cmd.index("resi 10"):
   for at2 in cmd.index("resi 11"):
       cmd.distance(None, "%s`%d"%at1, "%s`%d"%at2)