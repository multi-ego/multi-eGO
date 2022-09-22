import string

printable = list(string.printable)
#z = list(map(chr, range(0, 64)))
for h in printable[:172]: # number of chains
    for i in range(0, 144):# 0 to atoms in chain
        print(h)

#print(len(string.printable))



#5440 linee
#64 chains
#85 atoms
