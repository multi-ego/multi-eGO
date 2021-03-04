# Python script to calculate oligomerization state for each monomer chain and inter-contact map of protein chains residues

from networkx.readwrite.graph6 import data_to_n
import numpy as np
import copy as cp
import subprocess as sp
import os as os
import shutil as sh
import MDAnalysis as mdana
from MDAnalysis.analysis import distances
import sys
from MDAnalysis.analysis.distances import distance_array
import networkx as nx
import pandas as pd
import mdtraj as md
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import itertools
import pandas as pd


import time

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

print('atomsperAminoAcid', atomsperAminoAcid) # Non mi e' chiaro



#The N- and C- cap residues are part of the 1st and last residue index. If no N- and C- cap residues for the protein, comment the line below using "#"
#del atomsperAminoAcid[1]; del atomsperAminoAcid[-2]      

# Create Universe
print('Create Universe')

uni = mdana.Universe(ref_structure,traj)
n,t = list(enumerate(uni.trajectory))[0]
box = t.dimensions[:6]
print('Box dimensions', box)


# Why is this function repeated below??
atom_Groups = [[] for x in range(chains)]
m_start=0
for m in range(0,chains):
    # With this functions every chain is defined as a group of atoms
    # m is the chain number, atom is the number of atoms in every chain
    m_end = atoms * (m+1) # potevano evitare di fare questa cacata
    atom_Groups[m].extend([uni.select_atoms('bynum '+ str(m_start) + ':' + str(m_end))])
    m_start = m_end + 1


#fileout1 =  open('oligomer-groups.dat','w')
#fileout2 =  open('oligomer-states.dat','w')



cluster_for_traj = {}

for tt in uni.trajectory:
    print(tt.frame,f'out of', len(uni.trajectory))

    chains_interactions = pd.DataFrame(columns=['chain_ai', 'chain_aj', 'min_distance'])

    #fileout1.write (str(tt.frame) + '\t')
    #fileout2.write (str(tt.frame) + '\t')
    
    mySet = set([])
    graph = []

    chain_list = []
    for i in range(chains):
        chain_list.append(i)
    chain_comb = list(itertools.combinations(chain_list, 2))

    start_time = time.time()
    df_chain_ai, df_chain_aj, df_distance = [], [], []
    for n in range(0, len(chain_comb)): # Circa 16 secondi per iterazione
        i = chain_comb[n][0]
        j = chain_comb[n][1]
        dist = distance_array(atom_Groups[i][0].positions,atom_Groups[j][0].positions,box).min()
        if dist <= Min_Distance:
            my_tuple = (i, j)
            mySet.add(my_tuple)
            df_chain_ai.append(i)
            df_chain_aj.append(j)
            df_distance.append(dist)

    chains_interactions['chain_ai'] = df_chain_ai
    chains_interactions['chain_aj'] = df_chain_aj
    chains_interactions['min_distance'] = df_distance
    
    cluster_for_traj[tt.frame] = chains_interactions

    print("--- %s seconds ---" % (time.time() - start_time))



    # Le lunghezze sono uguali ma uno e' un set e dunque sempre in disordine
    #print(len(chains_interactions))
    #print(len(mySet))



    #graph = nx.from_edgelist(mySet) 

    #print(list(graph)) # Sono 584 chains nel sistema dove alla fine ci sono 4 molecole fuori 
                        # Graph di fatto e' l'elenco delle catene che fanno parte della fibrilla
    #print(len(graph))
    #for g in graph:
    #    print('g', g)

    #coso = nx.petersen_graph()
    #plt.subplot(121)
    #nx.draw(coso)
    #plt.show()


    # Di fatto quello che mi serve e' la lista dei clusters dunque basterebbe ciclare tutto e prendere tutto quello che e' piu grande di 1 (monomero)
    # Ogni linea e' un frame e per ogni linea ci sono tutti i clusters.
    # Se c'e' solo il monomero viene salvato, per ogni catena viene salvato tutto il cluster ogni volta
    #for i in range(chains):
    #    if i not in list(graph): # Graph sono le chain in fibrilla
    #        #print('i', i)
    #        fileout1.write ('['+ str(i)+']' + '\t') # Questi sono i monomeri
    #        fileout2.write ('1' + '\t')
    #    else: 
    #        fileout1.write (str(list(nx.node_connected_component(graph, i))) + '\t') # Per ogni catena in una fibrilla
                                                                                    # mi deve printare tutti i suoi amici vicini
                                                                                    # dunque ho 580 liste tutte uguali.
    #        #print(list(nx.node_connected_component(graph, i)))
    #        fileout2.write (str(len(list(nx.node_connected_component(graph, i)))) + '\t')
    #fileout1.write ('\n')
    #fileout2.write ('\n')
   
    
#fileout1.close()
#fileout2.close()


print(cluster_for_traj)


exit()


# Get oligomerization data

# OligoStates e' una lista 
print('FRAAAAAAAAAAMEEEEEES',frames+1)
OligoStates = [[0 for z in range(chains)] for x in range(frames+1)]

# OligoStates 9 liste contententi uno zero per ogni catena
# Dentro la lista frame va ad inserire
print('Oligomerization part')
file = open("oligomer-groups.dat",'r')
line = file.readline()
j = 0
while line:
    temp = line.split('\t') # Una linea equivale al frame
    for k in range(chains):
        #print(temp[k+1]) # il +1 e' perche' il primo elemento della riga e' il numero di frame
        OligoStates[j][k] = temp[k + 1][1:-1].split(',') # OligoStates[frame][chain#] = la lista corrispondente alla chain con tutto il contenuto e separato dalla virgola
    j += 1
    line = file.readline()
file.close

print(OligoStates[0][0])

# Create contact matrix

ContactMap = [[0 for x in range(AminoAcids)] for y in range(AminoAcids)]
print(ContactMap)
print('Matrix Contact')
# Create atom groups for each amino acid of each monomer

AtomGroups = [[] for x in range(chains)]

for m in range(0,chains):
    for aa in range(0,AminoAcids):
        #print(uni.select_atoms('bynum '+str(atoms*m + atomsperAminoAcid[aa])+':'+str(atoms*m + atomsperAminoAcid[aa + 1] - 1 )))
        AtomGroups[m].extend([uni.select_atoms('bynum '+str(atoms*m + atomsperAminoAcid[aa])+':'+str(atoms*m + atomsperAminoAcid[aa + 1] - 1 ))])



fibril_clust_prova = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574, 575, 576, 577, 578, 579, 580, 581, 582, 583]
fibril_clust = [int(x) for x in fibril_clust_prova]
atom_start = 0
for m in range(0, fibril_clust_prova):
    print(m)
    atom_end = atoms*m
    print(uni.select_atoms(f'bynum {atoms*m}'))
    atom_start = atom_end + 1




print('AtomGroups created')

count = 0
for n,t in enumerate(uni.trajectory[1:]):
    #print(n) # 0 che boh... al prossimo ciclo si vede
    print('frame', t) # timestep n
    on = 0
    Groups = []
    for i in OligoStates[n]: # Seleziona gli OligoStates in base al numero di frame
            # Se ho capito bene, gli oligostates sono le chains in ogni singolo cluster
        if len(i) > 1: # Dunque qui mi esclude i monomeri
            on = 1
            Groups.extend([i]) # Qui mi aggiunge le chain degli oligomeri in Groups
    Set = set(tuple(x) for x in Groups)
    Groups = [ list(x) for x in Set ]
    if on == 1:
    # Calculate dimension of the box to considered PBC
        box = t.dimensions[:6]
        for g in Groups: # seleziono il cluster
            print('len g',len(g))
        # Calculate contacts
            for n1,i in enumerate(g):
                #print('n1', n1)
                #print('i', i)
                for j in g[(n1 + 1)::]:
                    #print('g', g)
                    #print('j',j)
                    #print(g[(n1 + 1)::])
                    count += 1
                    for n2,atoms1 in enumerate(AtomGroups[int(float(i))]):
                        #print(i)
                        #print(n2)
                        #print(atoms1)
                        for n3,atoms2 in enumerate(AtomGroups[int(float(j))]):
                            #print('i', i, 'j',j, atoms1, atoms2, (distance_array(atoms1.positions,atoms2.positions,box).min()))
                            if ((distance_array(atoms1.positions,atoms2.positions,box).min()) <= Min_Distance):
                                ContactMap[n2][n3] +=1
                                ContactMap[n3][n2] +=1


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



