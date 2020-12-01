import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals, distances
import numpy as np
from MDAnalysis.analysis.dihedrals import Dihedral
from read_input import read_gro_bonds, read_gro_angles, read_gro_dihedrals
import pandas as pd
import itertools
from itertools import product, combinations

#native_pdb = mda.Universe('GRETA/native/pep.pdb', guess_bonds = True) # Da spostare su read pdb
native_bonds =  read_gro_bonds()
native_angles = read_gro_angles()
native_dihedrals = read_gro_dihedrals()

def make_exclusion_list (structure_pdb, native_bonds, native_angles, native_dihedrals, native_impropers):
    
    exclusion_list = []
    for index, row in native_bonds.iterrows():
        # For every bonds two atoms are defined and for every atom it is retrieved the atomtype
        exclusion_list.append(str(structure_pdb.atoms[row['ai'] - 1].name) + '_' + str(structure_pdb.atoms[row['ai'] - 1].resnum) + '_' + str(structure_pdb.atoms[row['aj'] - 1].name) + '_' + str(structure_pdb.atoms[row['aj'] - 1].resnum))   
    
    print('Exclusion List from bonds:               ', len(exclusion_list)) # 87

    for index, row in native_angles.iterrows():
        exclusion_list.append(str(structure_pdb.atoms[row['ai'] - 1].name) + '_' + str(structure_pdb.atoms[row['ai'] - 1].resnum) + '_' + str(structure_pdb.atoms[row['aj'] - 1].name) + '_' + str(structure_pdb.atoms[row['aj'] - 1].resnum))
        exclusion_list.append(str(structure_pdb.atoms[row['ai'] - 1].name) + '_' + str(structure_pdb.atoms[row['ai'] - 1].resnum) + '_' + str(structure_pdb.atoms[row['ak'] - 1].name) + '_' + str(structure_pdb.atoms[row['ak'] - 1].resnum))
        exclusion_list.append(str(structure_pdb.atoms[row['aj'] - 1].name) + '_' + str(structure_pdb.atoms[row['aj'] - 1].resnum) + '_' + str(structure_pdb.atoms[row['ak'] - 1].name) + '_' + str(structure_pdb.atoms[row['ak'] - 1].resnum))
    
    print('Addition of angles to exclusion list:    ', len(exclusion_list)) # 447


    for index, row in native_dihedrals.iterrows():
        exclusion_list.append(str(structure_pdb.atoms[row['ai'] - 1].name) + '_' + str(structure_pdb.atoms[row['ai'] - 1].resnum) + '_' + str(structure_pdb.atoms[row['aj'] - 1].name) + '_' + str(structure_pdb.atoms[row['aj'] - 1].resnum))
        exclusion_list.append(str(structure_pdb.atoms[row['ai'] - 1].name) + '_' + str(structure_pdb.atoms[row['ai'] - 1].resnum) + '_' + str(structure_pdb.atoms[row['ak'] - 1].name) + '_' + str(structure_pdb.atoms[row['ak'] - 1].resnum))
        exclusion_list.append(str(structure_pdb.atoms[row['ai'] - 1].name) + '_' + str(structure_pdb.atoms[row['ai'] - 1].resnum) + '_' + str(structure_pdb.atoms[row['al'] - 1].name) + '_' + str(structure_pdb.atoms[row['al'] - 1].resnum))
        exclusion_list.append(str(structure_pdb.atoms[row['aj'] - 1].name) + '_' + str(structure_pdb.atoms[row['aj'] - 1].resnum) + '_' + str(structure_pdb.atoms[row['ak'] - 1].name) + '_' + str(structure_pdb.atoms[row['ak'] - 1].resnum))
        exclusion_list.append(str(structure_pdb.atoms[row['aj'] - 1].name) + '_' + str(structure_pdb.atoms[row['aj'] - 1].resnum) + '_' + str(structure_pdb.atoms[row['al'] - 1].name) + '_' + str(structure_pdb.atoms[row['al'] - 1].resnum))
        exclusion_list.append(str(structure_pdb.atoms[row['ak'] - 1].name) + '_' + str(structure_pdb.atoms[row['ak'] - 1].resnum) + '_' + str(structure_pdb.atoms[row['al'] - 1].name) + '_' + str(structure_pdb.atoms[row['al'] - 1].resnum))

    print('Addition of dihedrals to exclusion list: ', len(exclusion_list)) # 861

    for index, row in native_impropers.iterrows():
        exclusion_list.append(str(structure_pdb.atoms[row['ai'] - 1].name) + '_' + str(structure_pdb.atoms[row['ai'] - 1].resnum) + '_' + str(structure_pdb.atoms[row['aj'] - 1].name) + '_' + str(structure_pdb.atoms[row['aj'] - 1].resnum))
        exclusion_list.append(str(structure_pdb.atoms[row['ai'] - 1].name) + '_' + str(structure_pdb.atoms[row['ai'] - 1].resnum) + '_' + str(structure_pdb.atoms[row['ak'] - 1].name) + '_' + str(structure_pdb.atoms[row['ak'] - 1].resnum))
        exclusion_list.append(str(structure_pdb.atoms[row['ai'] - 1].name) + '_' + str(structure_pdb.atoms[row['ai'] - 1].resnum) + '_' + str(structure_pdb.atoms[row['al'] - 1].name) + '_' + str(structure_pdb.atoms[row['al'] - 1].resnum))
        exclusion_list.append(str(structure_pdb.atoms[row['aj'] - 1].name) + '_' + str(structure_pdb.atoms[row['aj'] - 1].resnum) + '_' + str(structure_pdb.atoms[row['ak'] - 1].name) + '_' + str(structure_pdb.atoms[row['ak'] - 1].resnum))
        exclusion_list.append(str(structure_pdb.atoms[row['aj'] - 1].name) + '_' + str(structure_pdb.atoms[row['aj'] - 1].resnum) + '_' + str(structure_pdb.atoms[row['al'] - 1].name) + '_' + str(structure_pdb.atoms[row['al'] - 1].resnum))
        exclusion_list.append(str(structure_pdb.atoms[row['ak'] - 1].name) + '_' + str(structure_pdb.atoms[row['ak'] - 1].resnum) + '_' + str(structure_pdb.atoms[row['al'] - 1].name) + '_' + str(structure_pdb.atoms[row['al'] - 1].resnum))
    
    print('Addition of impropers to exclusion list: ', len(exclusion_list)) # 1119

    # Keep only unique values
    exclusion_list = list(set(exclusion_list))
    print('Drop duplicates in the exclusion list:   ', len(exclusion_list)) # 350
    return exclusion_list


########################## DIHEDRALS

def sb_dihedrals (structure_pdb):
    phi_dihedrals = []
    for index, row in native_dihedrals.iterrows():
        # Here are selected only the atom numbers for every dihedral from the pdb structure
        # ANCORA DA FINIRE PERCHE' I DIEDRI DA PDB NON STANNO CORRISPONDENDO
        atom_selection = structure_pdb.atoms[row['ai'] - 1] + structure_pdb.atoms[row['aj'] - 1] + structure_pdb.atoms[row['ak'] - 1] + structure_pdb.atoms[row['al'] - 1]
        phi = atom_selection.dihedral.value()
        phi_dihedrals.append(phi)

    native_dihedrals['func'] = 9
    native_dihedrals['phi'] = phi_dihedrals
    native_dihedrals['kd'] = ''
    native_dihedrals['mult'] = ''
    #print(native_dihedrals)
    return native_dihedrals

################################ PAIRS

def make_pairs (structure_pdb, exclusion_list):

    # Selection of all the atoms required to compute LJ
    atom_sel = structure_pdb.select_atoms('all')
    # Calculating all the distances between atoms
    # The output is a combination array 
    self_distances = distances.self_distance_array(atom_sel.positions)

    # Making a list of the atomtypes of the protein
    # Da spostare su atomtypes_definitions
    atomtype = []
    for atom in atom_sel:
        atp = str(atom.name) + '_' + str(atom.resnum)
        atomtype.append(atp)

    # Combining all the atomtypes in the list to create a pair list corresponding to the distance array
    pairs_list = list(itertools.combinations(atomtype, 2))

    pairs_ai = []
    pairs_aj = []
    # But the combinations are list of list and we need to separate them
    for n in range(0, len(pairs_list)):
        i = pairs_list[n][0]
        pairs_ai.append(i)
        j = pairs_list[n][1]
        pairs_aj.append(j)

    # Creation of the dataframe containing the ffnonbonded.itp
    structural_LJ = pd.DataFrame(columns = ['ai', 'aj', 'distance', 'sigma', 'epsilon', 'check'])
    structural_LJ['ai'] = pairs_ai
    structural_LJ['aj'] = pairs_aj
    structural_LJ['distance'] = self_distances
    structural_LJ['check'] = structural_LJ['ai'] + '_' + structural_LJ['aj']
    # Here we keep only the one without the exclusions
    structural_LJ = structural_LJ[~structural_LJ['check'].isin(exclusion_list)]
    # Keep only the atoms within 6 A
    structural_LJ = structural_LJ[structural_LJ.distance < 6]
    # Drop of the pairs which are included in bonds, angles and dihedrals.
    structural_LJ['sigma'] = (structural_LJ['distance']/10) / (2**(1/6))
    structural_LJ['epsilon'] = 1

    # This part is to filter more the LJ like in smog: if two pairs are made by aminoacids closer than
    # 3 they'll be deleted. Therefore aminoacids 1, 2, 3 and 4 does not make any contacts.
    # Therefore I copy the LJ dataframe and apply some more filters
    structural_LJ[['type_ai', 'resnum_ai']] = structural_LJ.ai.str.split("_", expand = True)
    structural_LJ[['type_aj', 'resnum_aj']] = structural_LJ.aj.str.split("_", expand = True)
    # And to do that it is necessary to convert the two columns into integer
    structural_LJ = structural_LJ.astype({"resnum_ai": int, "resnum_aj": int})
    structural_LJ['diff'] = ''
    structural_LJ.drop(structural_LJ[abs((structural_LJ['resnum_aj'] - structural_LJ['resnum_ai'])) < 4].index, inplace = True)
    structural_LJ['diff'] = abs(structural_LJ['resnum_aj'] - structural_LJ['resnum_ai'])
    
    #print(structural_LJ.to_string()) 
    print(len(structural_LJ))
    #print(exclusion_list)

    return structural_LJ



# EXCLUSION SMOG togli i contatti a < di 3 a di distanza. 




























##### SEPARA LA TOPOLOGIA A MANO PER FARE LE EXCLUSIONS


#native_bonds = pd.DataFrame(columns = [';ai', 'aj'])
#b_ai = []
#b_aj = []
#for b in native_pdb.bonds:
#    b_ai.append(str((b.atoms[0].name) + '_' + str(b.atoms[0].resnum))) 
#    b_aj.append(str((b.atoms[1].name) + '_' + str(b.atoms[1].resnum)))


#bonds_combinations = list(product(b_ai, b_aj))
#print(len(bonds_combinations))

#native_bonds[';ai'] = b_ai
#native_bonds['aj'] = b_aj


#print(native_bonds) # OK


#native_angles = pd.DataFrame(columns = [';ai', 'aj', 'ak'])
#a_ai = []
#a_aj = []
#a_ak = []
#for a in native_pdb.angles:
#    a_ai.append(str((a.atoms[0].name) + '_' + str(a.atoms[0].resnum))) 
#    a_aj.append(str((a.atoms[1].name) + '_' + str(a.atoms[1].resnum)))
#    a_ak.append(str((a.atoms[2].name) + '_' + str(a.atoms[2].resnum)))

#native_angles[';ai'] = a_ai
#native_angles['aj'] = a_aj
#native_angles['ak'] = a_ak

#native_dihe = pd.DataFrame(columns = [';ai', 'aj', 'ak', 'al'])
#d_ai = []
#d_aj = []
#d_ak = []
#d_al = []
#for d in native_pdb.dihedrals:
    #d_ai.append(str((d.atoms[0].name) + '_' + str(d.atoms[0].resnum))) 
    #d_aj.append(str((d.atoms[1].name) + '_' + str(d.atoms[1].resnum)))
    #d_ak.append(str((d.atoms[2].name) + '_' + str(d.atoms[2].resnum)))
    #d_al.append(str((d.atoms[3].name) + '_' + str(d.atoms[3].resnum)))
 #   d_ai.append(d.atoms[0].index + 1)
 #   d_aj.append(d.atoms[1].index + 1)
 #   d_ak.append(d.atoms[2].index + 1)
 #   d_al.append(d.atoms[3].index + 1)

#print(len(native_pdb.dihedrals))
#
##native_dihe[';ai'] = d_ai
#native_dihe['aj'] = d_aj
#native_dihe['ak'] = d_ak
#native_dihe['al'] = d_ak

#print(native_dihe)