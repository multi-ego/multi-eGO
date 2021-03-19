import MDAnalysis
from MDAnalysis.analysis import distances
import itertools
from numpy.lib.function_base import average
import pandas as pd
import time
import sys

distance_cutoff = 5.5
reference_structure = sys.argv[1]
reference_trajectory = sys.argv[2]
#u = MDAnalysis.Universe('box_peptide_greta.gro', 'prod_peptide_greta.xtc')
u = MDAnalysis.Universe(reference_structure, reference_trajectory)
#peptides = u.select_atoms('resid 1:11 and not name H*') # ho dovuto selezionare il rage di aa perche senno non mi toglieva gli H
peptides = u.select_atoms('all')
print(u.residues)
print(len(peptides))

atomtypes = []
for atom in peptides:
    atp = str(atom.name) + '_' + str(atom.resnum)
    atomtypes.append(atp)

pairs_list = list(itertools.combinations(atomtypes, 2))
pairs_ai, pairs_aj = [], []
for n in range(0, len(pairs_list)):
    i = pairs_list[n][0]
    pairs_ai.append(i)
    j = pairs_list[n][1]
    pairs_aj.append(j)

print(len(pairs_list))
print(len(u.trajectory))

extended_ai, extended_aj, monomer_distances = [], [], []
for frame in u.trajectory:
    self_distance = distances.self_distance_array(peptides.positions)
    if len(self_distance) == len(pairs_list):
        # Which is a list of arrays
        monomer_distances.extend(self_distance)
        extended_ai.extend(pairs_ai)
        extended_aj.extend(pairs_aj)
    else:
        print('Pairs list and distance array have different length')
        break

print(len(pairs_list))
print(len(u.trajectory))
print(len(monomer_distances))

if ((len(pairs_list) * len(u.trajectory)) == len(monomer_distances)) and (len(extended_ai) == len(extended_aj) == len(monomer_distances)):
    print('BRAVOH')

# exclusion list degli aminoacidi

monomer_pairs_df = pd.DataFrame(columns=['ai', 'aj','ai_name', 'aj_name', 'ai_resnum', 'aj_resnum', 'distances'])
monomer_pairs_df['ai'] = extended_ai
monomer_pairs_df['aj'] = extended_aj
monomer_pairs_df['distances'] = monomer_distances

print(len(monomer_pairs_df))
monomer_pairs_df = monomer_pairs_df[monomer_pairs_df['distances'] < distance_cutoff]
print(len(monomer_pairs_df))
monomer_pairs_df[['ai_name','ai_resnum']] = monomer_pairs_df.ai.str.split("_", expand=True)
monomer_pairs_df[['aj_name','aj_resnum']] = monomer_pairs_df.aj.str.split("_", expand=True)
monomer_pairs_df = monomer_pairs_df.astype({"ai_resnum": int, "aj_resnum": int})
monomer_pairs_df.drop(monomer_pairs_df[abs(monomer_pairs_df['aj_resnum'] - monomer_pairs_df['ai_resnum']) < 3].index, inplace=True)
print(len(monomer_pairs_df))
print(monomer_pairs_df)


count_ai, count_aj, count_distance, count_ratio, average_distance = [], [], [], [], []
start_time, c_num =time.time(), 1
for pair in pairs_list:
    # filtering the data frame based on the pairs values
    #print(c_num, len(pairs_list))
    count_ai.append(pair[0])
    count_aj.append(pair[1])
    # salvati il df che serve per la media delle distanze e del sigma
    counts_df = monomer_pairs_df[(monomer_pairs_df['ai'] == pair[0]) & (monomer_pairs_df['aj'] == pair[1])]
    average_distance.append(counts_df['distances'].mean())
    count_distance.append(len(counts_df))
    count_ratio.append(len(counts_df)/len(u.trajectory))
    c_num += 1
    #print(f'{time.time() - start_time}')

pairs_count = pd.DataFrame(columns=['ai', 'aj', 'count', 'ratio', 'average_distance'])
pairs_count['ai'] = count_ai
pairs_count['aj'] = count_aj
pairs_count['count'] = count_distance
pairs_count['ratio'] = count_ratio
pairs_count['average_distance'] = average_distance
pairs_count.sort_values(by = ['ratio'], inplace = True, ascending=False)

file = open('monomer_pairs.txt', 'w')
file.write(pairs_count.to_string(index=False, header=False))
file.close()