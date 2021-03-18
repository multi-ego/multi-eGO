import string

#z = range(1, 65)
printable = list(string.printable)
#print(printable)
#z = list(map(chr, range(0, 64)))
for h in printable[:22]: # number of chains
    for i in range(0, 185):# 0 to atoms in chain
        print(h)

#print(len(string.printable))



#5440 linee
#64 chains
#85 atoms