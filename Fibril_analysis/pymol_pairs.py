import pandas as pd
import numpy as np

# Import the pairs section of FFnonbonded.itp


pairs = pd.read_csv('analysis/pairs', sep = "\s+", header = None)
pairs.columns = [";ai", "aj", "type", "A", "B"]
pairs = pairs.drop(['type', 'A', 'B'], axis = 1)




print(pairs)


























# This function is a basic pairs analyzing tool which can come handy and to develop in the future

def pairs_separator(pairs):
    # Poi commenta un po' quello che serve che neanche ora ti ricordi come funziona questo script

    ai_in = pairs['ai'].str.split('_', n = 1, expand = True)
    aj_jn = pairs['aj'].str.split('_', n = 1, expand = True)
    pairs['ain'] = ai_in[1]
    pairs['air'] = ai_in[0]
    pairs['ajn'] = aj_jn[1]
    pairs['ajr'] = aj_jn[0]
    pairs.astype(str)

    # Here we start to select the backbone and the sidechains
    # The mixed pairs between the two are not included atm
    back_atomtype = ['CA', 'N', 'O', 'C', 'OXT']
    # Backbone
    drop_sidechains = pairs[pairs['air'].isin(back_atomtype)]
    backbone = drop_sidechains[drop_sidechains['ajr'].isin(back_atomtype)]
    # Sidechains
    drop_backbone = pairs[~ pairs['air'].isin(back_atomtype)]
    sidechains = drop_backbone[~ drop_backbone['ajr'].isin(back_atomtype)]
    # mixed #con questo metodo non funziona
    mix_air = pairs[pairs['air'].isin(back_atomtype)]
    mix_ajr = pairs[pairs['ajr'].isin(back_atomtype)]
    mix_full = mix_air.append(mix_ajr, sort = False, ignore_index = True)
    # mixed poi ci penso
    pymol_backbone = backbone.copy()
    pymol_sidechains = sidechains.copy()
    pymol_mixed = mix_full.copy()

    # cleaning the dataframes
    backbone = backbone.drop(['air', 'ajr', 'ain', 'ajn'], axis = 1)
    sidechains = sidechains.drop(['air', 'ajr', 'ain', 'ajn'], axis = 1)
    mix_full = mix_full.drop(['air', 'ajr', 'ain', 'ajn'], axis = 1)

    file = open("pymol/backbone", "w")
    file.write(str(backbone.to_string(index = False, header = False)))
    file.close()

    file = open("pymol/sidechains", "w")
    file.write(str(sidechains.to_string(index = False, header = False)))
    file.close()

    file = open("pymol/mixed", "w")
    file.write(str(mix_full.to_string(index = False, header = False)))
    file.close()

    # backbone
    pymol_backbone = pymol_backbone.drop(['type', 'A', 'B'], axis = 1)
    pymol_backbone['interaction'] = np.where((pymol_backbone['ai'] == pymol_backbone['aj']), 'v', 'h')
    # PROVA AD UNTILIZZARE .JOIN
    name = pymol_backbone['interaction'].apply(str) + '_' + pymol_backbone['ai'].apply(str) + ':' + pymol_backbone[
        'aj'].apply(str)
    pymol_backbone = pymol_backbone.drop(['ai', 'aj'], axis = 1)
    pymol_backbone.insert(0, 'distance', 'distance')
    pymol_backbone.insert(1, 'name', name)
    pymol_backbone.insert(2, 'resi', ', resi')
    pymol_backbone.insert(4, 'and', 'and name')
    pymol_backbone.insert(6, 'resi2', ', resi')
    pymol_backbone.insert(8, 'and2', 'and name')
    pymol_backbone.insert(10, 'cutoff', ', 6')
    pymol_backbone = pymol_backbone.sort_values(by = 'interaction')
    pymol_backbone = pymol_backbone.drop(['interaction'], axis = 1)

    # sidechains
    pymol_sidechains = pymol_sidechains.drop(['type', 'A', 'B'], axis = 1)
    pymol_sidechains['interaction'] = np.where((pymol_sidechains['ai'] == pymol_sidechains['aj']), 'v', 'h')
    name = pymol_sidechains['interaction'].apply(str) + '_' + pymol_sidechains['ai'].apply(str) + ':' + \
           pymol_sidechains['aj'].apply(str)
    pymol_sidechains = pymol_sidechains.drop(['ai', 'aj'], axis = 1)
    pymol_sidechains.insert(0, 'distance', 'distance')
    pymol_sidechains.insert(1, 'name', name)
    pymol_sidechains.insert(2, 'resi', ', resi')
    pymol_sidechains.insert(4, 'and', 'and name')
    pymol_sidechains.insert(6, 'resi2', ', resi')
    pymol_sidechains.insert(8, 'and2', 'and name')
    pymol_sidechains.insert(10, 'cutoff', ', 6')
    pymol_sidechains = pymol_sidechains.sort_values(by = 'interaction')
    pymol_sidechains = pymol_sidechains.drop(['interaction'], axis = 1)

    # mixed
    pymol_mixed = pymol_mixed.drop(['type', 'A', 'B'], axis = 1)
    pymol_mixed['interaction'] = np.where((pymol_mixed['ai'] == pymol_mixed['aj']), 'v', 'h')
    name = pymol_mixed['interaction'].apply(str) + '_' + pymol_mixed['ai'].apply(str) + ':' + pymol_mixed['aj'].apply(
        str)
    pymol_mixed = pymol_mixed.drop(['ai', 'aj'], axis = 1)
    pymol_mixed.insert(0, 'distance', 'distance')
    pymol_mixed.insert(1, 'name', name)
    pymol_mixed.insert(2, 'resi', ', resi')
    pymol_mixed.insert(4, 'and', 'and name')
    pymol_mixed.insert(6, 'resi2', ', resi')
    pymol_mixed.insert(8, 'and2', 'and name')
    pymol_mixed.insert(10, 'cutoff', ', 6')
    pymol_mixed = pymol_mixed.sort_values(by = 'interaction')
    pymol_mixed = pymol_mixed.drop(['interaction'], axis = 1)

    file = open("pymol/pymol_backbone", "w")
    file.write(str(pymol_backbone.to_string(index = False, header = False)))
    file.close()

    file = open("pymol/pymol_sidechains", "w")
    file.write(str(pymol_sidechains.to_string(index = False, header = False)))
    file.close()

    file = open("pymol/pymol_mixed", "w")
    file.write(str(pymol_mixed.to_string(index = False, header = False)))
    file.close()