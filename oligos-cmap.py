# Python script to calculate oligomerization state for each monomer chain and inter-contact map of protein chains residues

import numpy as np
import copy as cp
import subprocess as sp
import os as os
import shutil as sh
import MDAnalysis as mdana
import sys
from MDAnalysis.analysis.distances import distance_array
import networkx as nx
import pandas as pd
import mdtraj as md
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt

#input parameters

ref_structure=sys.argv[1]
traj=sys.argv[2]
Min_Distance=int(sys.argv[3])

#structure parameters

topology = md.load(ref_structure).topology
trajectory = md.load(traj, top=ref_structure)
frames=trajectory.n_frames				#Number of frames
chains=topology.n_chains				#Number of chains

print(f'Number of chains', chains)

atoms=int(topology.n_atoms/chains)			#Number of atoms in each monomer 
AminoAcids = int(topology.n_residues/chains)#-2		#Number of residues per chain ('-2' avoid the N- and C- cap residues as individual residues)

isum=1
atoms_list=[]
atomsperAminoAcid=[]
residue_list=[]

for residue in topology.chain(0).residues:
    print(f'Residues in topology', residue)
    atoms_list.append(residue.n_atoms)
    residue_list.append(residue)
    ', '.join(map(lambda x: "'" + x + "'", str(residue_list)))
#The N- and C- cap residues are part of the 1st and last residue index. If no N- and C- cap residues for the protein, comment the line below using "#"
#del residue_list[0]; del residue_list[-1] 

for ii in range(len(atoms_list)):
    isum = isum + atoms_list[ii]
    atomsperAminoAcid.append(isum)
atomsperAminoAcid.insert(0, 1)

#The N- and C- cap residues are part of the 1st and last residue index. If no N- and C- cap residues for the protein, comment the line below using "#"
#del atomsperAminoAcid[1]; del atomsperAminoAcid[-2]      

# Create Universe

uni = mdana.Universe(ref_structure,traj)
n,t = list(enumerate(uni.trajectory))[0]
box = t.dimensions[:6]



atom_Groups = [[] for x in range(chains)]
m_start=0
for m in range(0,chains):
    m_end = atoms * (m+1)
    atom_Groups[m].extend([uni.select_atoms('bynum '+ str(m_start) + ':' + str(m_end))])
    m_start = m_end + 1
fileout1 =  open('oligomer-groups.dat','w')
fileout2 =  open('oligomer-states.dat','w')

for tt in uni.trajectory[1:]:
    print(tt,f'out of', len(uni.trajectory[1:]))
    fileout1.write (str(tt.frame) + '\t')
    fileout2.write (str(tt.frame) + '\t')
    
    mySet = set([])
    graph = []
    atom_Groups = [[] for x in range(chains)]
    m_start=0
    for m in range(0,chains):
        m_end = atoms * (m+1)
        atom_Groups[m].extend([uni.select_atoms('bynum '+ str(m_start) + ':' + str(m_end))])
        m_start = m_end + 1
    count = 0
    for i in range(chains-1):
        for j in range(i+1,chains):
            dist = distance_array(atom_Groups[i][0].positions,atom_Groups[j][0].positions,box).min()
            if(dist <= Min_Distance):
                my_tuple = (i,j)
                mySet.add(my_tuple)
                
    graph = nx.from_edgelist(mySet)   
    for i in range(chains):
        if i not in list(graph):
            fileout1.write ('['+ str(i)+']' + '\t')
            fileout2.write ('1' + '\t')
        else:
            fileout1.write (str(list(nx.node_connected_component(graph, i))) + '\t')
            fileout2.write (str(len(list(nx.node_connected_component(graph, i)))) + '\t')
    fileout1.write ('\n')
    fileout2.write ('\n')
   
    
fileout1.close()
fileout2.close()


# Get oligomerization data

OligoStates = [[0 for z in range(chains)] for x in range(frames+1)]
print('Oligomerization part')
file = open("oligomer-groups.dat",'r')
line = file.readline()
j = 0
while line:
    temp = line.split('\t')
    for k in range(chains):
        OligoStates[j][k] = temp[k + 1][1:-1].split(',')
    j += 1
    line = file.readline()
file.close


# Create contact matrix

ContactMap = [[0 for x in range(AminoAcids)] for y in range(AminoAcids)]
print('Matrix Contact')
# Create atom groups for each amino acid of each monomer

AtomGroups = [[] for x in range(chains)]

for m in range(0,chains):
    for aa in range(0,AminoAcids):
        AtomGroups[m].extend([uni.select_atoms('bynum '+str(atoms*m + atomsperAminoAcid[aa])+':'+str(atoms*m + atomsperAminoAcid[aa + 1] - 1 ))])



count = 0
for n,t in enumerate(uni.trajectory[1:]):
    on = 0
    Groups = []
    for i in OligoStates[n]:
        if len(i) > 1:
            on = 1
            Groups.extend([i])
    Set = set(tuple(x) for x in Groups)
    Groups = [ list(x) for x in Set ]
    
    if on == 1:
    # Calculate dimension of the box to considered PBC
        box = t.dimensions[:6]
        for g in Groups:
        # Calculate contacts
            for n1,i in enumerate(g):
                for j in g[(n1 + 1)::]:
                    count += 1
                    for n2,atoms1 in enumerate(AtomGroups[int(float(i))]):
                        for n3,atoms2 in enumerate(AtomGroups[int(float(j))]):
                            if ((distance_array(atoms1.positions,atoms2.positions,box).min()) <= Min_Distance):
                                ContactMap[n2][n3] +=1
                                ContactMap[n3][n2] +=1 

#print(count)
Norm_ContactMap = np.true_divide(ContactMap,float(count)*2.0)

# Save contact map in a file

fileout = open ('contact-map.dat','w')
for i in Norm_ContactMap:
	for j in i:
		fileout.write (str(j) + '\t')
	fileout.write ('\n')
fileout.close()



#Highest Oligomer size in each frame

states=open('oligomer-states.dat', 'r')
ter=states.readlines()[0:frames+1]

result=[]
for freq in (ter):
    	result.append([int(hist) for hist in freq.strip().split('\t')[1:]])

fileout3 = open ('oligo-highest-size.dat', 'w')
for oli_count in range(len(ter)):
	fileout3.write("{} {} {}\n".format(oli_count, '\t', np.max(result[oli_count])))
fileout3.close()



# Block Average


size_data = np.loadtxt('oligomer-states.dat')
window = 25                                         	# specify over how many frames the running average is to be calculated
weights = np.repeat(1.0,window)/window
size_data_m = np.convolve(size_data[:,1],weights,'valid')

fileout4 = open('oligo-block-average.dat', 'w')
for t,b in enumerate(size_data_m):
    	fileout4.write("{} {} {}\n".format(t, '\t', b))
fileout4.close()



