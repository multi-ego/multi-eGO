from read_input import read_native_pairs
from MDAnalysis.analysis import distances
import numpy as np
from pandas.core.frame import DataFrame
import pandas as pd
import itertools
from protein_configuration import distance_cutoff, distance_residue, epsilon_input, idp, ratio_treshold, protein, N_terminal, sigma_method
from topology_definitions import topology_atoms, gromos_atp, gro_to_amb_dict, topology_bonds, atom_topology_num, first_resid


def make_pdb_atomtypes (native_pdb, fibril_pdb):


    print('\t Native atomtypes')
    native_sel = native_pdb.select_atoms('all')
    native_atomtypes, ffnb_sb_type = [], []

    for atom in native_sel:
        # Also the chain ID is printed along because it might happen an interaction between two atoms of different 
        # chains but of the same residue which would be deleted based only by the eclusion list as example
        # CA_1 N_1 is in the exclusion list but in fibril might interact the CA_1 chain 1 with N_1 chain 2
        
        # Native atomtypes will be used to create the pairs list
        atp = str(atom.name) + '_' + str(atom.resnum) + ':' + str(atom.segid)
        native_atomtypes.append(atp)

        # This part is for attaching to FFnonbonded.itp
        # ffnb_sb_type is the first column
        # check gromologist
        sb_type = str(atom.name) + '_' + str(atom.resnum)
        ffnb_sb_type.append(sb_type)

    check_topology = DataFrame(ffnb_sb_type, columns=['sb_type'])
    check_topology['check'] = np.where(topology_atoms.sb_type == check_topology.sb_type, 'same', 'different')
    
    if len(np.unique(check_topology.check)) == 1:
        print('\n Atoms of topol.top and pdb have the same number')

    else:
        print('\n Check PDB and topology numeration')
        
    print('\t Fibril atomtypes')
    fibril_sel = fibril_pdb.select_atoms('all')
    fibril_atomtypes = []
    for atom in fibril_sel:
        atp = str(atom.name) + '_' + str(atom.resnum) + ':' + str(atom.segid)
        fibril_atomtypes.append(atp)

    # ffnonbonded making
    # Making a dictionary with atom number and type
    print('\t FFnonbonded atomtypes section creation')
    ffnb_atomtype = pd.DataFrame(columns = ['; type', 'chem', 'at.num', 'mass', 'charge', 'ptype', 'c6', 'c12'])
    ffnb_atomtype['; type'] = topology_atoms['sb_type']
    ffnb_atomtype['chem'] = topology_atoms['type']
    ffnb_atomtype['at.num'] = ffnb_atomtype['chem'].map(gromos_atp['at.num'])
    ffnb_atomtype['mass'] = topology_atoms['mass']
    ffnb_atomtype['charge'] = '0.000000'
    ffnb_atomtype['ptype'] = 'A'
    ffnb_atomtype['c6'] = '0.00000e+00'
    ffnb_atomtype['c12'] = ffnb_atomtype['chem'].map(gromos_atp['c12'])
    
    residue_list = topology_atoms['residue'].to_list()

    if 'PRO' in residue_list:
        print('There are prolines in the structure. The c12 of N should be the half')
        proline_n = topology_atoms.loc[(topology_atoms['residue'] == 'PRO') & (topology_atoms['atom'] == 'N'), 'sb_type'].to_list()
        ffnb_atomtype.loc[(ffnb_atomtype['; type'].isin(proline_n)), 'c12'] = ffnb_atomtype['c12']/2
    else:
        print('There not are prolines in the structure. The c12 of N should be the half')
    

    if N_terminal == True:
        print('Changing the c12 value of N-terminal')
        #first_resid = 'N_'+str(atom_topology_resid[0])
        ffnb_atomtype.loc[(ffnb_atomtype['; type'] == first_resid), 'c12'] = ffnb_atomtype['c12']*5 # Harp 2
        

    # This will be needed for exclusion and pairs to paste in topology
    type_c12_dict = ffnb_atomtype.set_index('; type')['c12'].to_dict()
    
    ffnb_atomtype['c12'] = ffnb_atomtype["c12"].map(lambda x:'{:.6e}'.format(x))
    ffnb_atomtype.drop(columns = ['chem'], inplace = True)

    print('\t Topology atomtypes section creation')
    topology_atoms['type'] = topology_atoms['sb_type']
    topology_atoms.insert(5, 'cgnr', 1)
    topology_atoms.insert(6, 'charge', '')
    topology_atoms['mass'] = ''
    topology_atoms['typeB'] = ''
    topology_atoms['chargeB'] = ''
    topology_atoms['massB'] = ''
    topology_atoms.rename(columns={'nr':'; nr'}, inplace=True)
    topology_atoms.drop(columns=['sb_type'], inplace=True)

    print('\t Atomtypes.atp file creation')
    atomtypes_atp = ffnb_atomtype[['; type', 'mass']].copy()

    return native_atomtypes, fibril_atomtypes, ffnb_atomtype, atomtypes_atp, topology_atoms, type_c12_dict

################################ PAIRS

def make_pairs (structure_pdb, atomtypes):

    print('\n\t Measuring distances between all atom in the pdb')
    # Selection of all the atoms required to compute LJ
    atom_sel = structure_pdb.select_atoms('all')
    # Calculating all the distances between atoms
    # The output is a combination array 
    self_distances = distances.self_distance_array(atom_sel.positions)
    print('\t Number of distances measured :', len(self_distances))
    
    print('\n\t Preparing the atomtype array')
    # Combining all the atomtypes in the list to create a pair list corresponding to the distance array
    pairs_list = list(itertools.combinations(atomtypes, 2))
    pairs_ai, pairs_aj = [], []

    # But the combinations are list of list and we need to separate them
    for n in range(0, len(pairs_list)):
        i = pairs_list[n][0]
        pairs_ai.append(i)
        j = pairs_list[n][1]
        pairs_aj.append(j)
    print('\t Atomtype array ready')

    print('\n\t Creating the pairs dataframes')
    # Creation of the dataframe containing the ffnonbonded.itp
    structural_LJ = pd.DataFrame(columns = ['ai', 'aj', 'distance', 'sigma', 'epsilon', 'check'])
    structural_LJ['ai'] = pairs_ai
    structural_LJ['aj'] = pairs_aj
    structural_LJ['distance'] = self_distances
    print('\t Raw pairs list ', len(structural_LJ))
    
    print('\n\t Extracting chain label for every atom pair')
    structural_LJ[['ai', 'chain_ai']] = structural_LJ.ai.str.split(":", expand = True)
    structural_LJ[['aj', 'chain_aj']] = structural_LJ.aj.str.split(":", expand = True)

    # Create the check column for the exclusion list
    structural_LJ['check'] = structural_LJ['ai'] + '_' + structural_LJ['aj']
    same_chain = np.where(structural_LJ['chain_ai'] == structural_LJ['chain_aj'], 'Yes', 'No')
    structural_LJ['same_chain'] = same_chain

    print('\t Tagging the chain association')
    are_same = structural_LJ.loc[structural_LJ['same_chain'] == 'Yes']
    not_same = structural_LJ.loc[structural_LJ['same_chain'] == 'No']
    print('\t Pairs within the same chain: ', len(are_same))
    print('\t Pairs not in the same chain: ', len(not_same))
    print('\t Raw pairs list ', len(structural_LJ))

    print('\t Tagging pairs included in bonded exclusion list')
    # Here we keep only the one without the exclusions
    structural_LJ['exclude'] = ''
    not_same = structural_LJ.loc[structural_LJ['same_chain'] == 'No']

    print(f'\n\t Applying distance cutoff of {distance_cutoff} A')
    # Keep only the atoms within 6 A
    structural_LJ = structural_LJ[structural_LJ.distance < distance_cutoff] # PROTEIN CONFIGURATION
    print(f'\t Pairs below cutoff {distance_cutoff}: ', len(structural_LJ))

    print(f'\t Applying residue number cutoff of {distance_residue}')
    # This part is to filter more the LJ like in smog: if two pairs are made by aminoacids closer than
    # 3 they'll be deleted. Therefore aminoacids 1, 2, 3 and 4 does not make any contacts.
    # Therefore I copy the LJ dataframe and apply some more filters
    structural_LJ[['type_ai', 'resnum_ai']] = structural_LJ.ai.str.split("_", expand = True)
    structural_LJ[['type_aj', 'resnum_aj']] = structural_LJ.aj.str.split("_", expand = True)
    # And to do that it is necessary to convert the two columns into integer
    structural_LJ = structural_LJ.astype({"resnum_ai": int, "resnum_aj": int})
    structural_LJ['diff'] = ''
    # Da riattivare successivamente, test tenendo solo exlcusion list bonded e non SB
    structural_LJ.drop(structural_LJ[(abs(structural_LJ['resnum_aj'] - structural_LJ['resnum_ai']) < distance_residue) & (structural_LJ['same_chain'] == 'Yes')].index, inplace = True)    
    structural_LJ['diff'] = abs(structural_LJ['resnum_aj'] - structural_LJ['resnum_ai'])
    print(f'\t All the pairs further than {distance_residue} aminoacids and not in the same chain: ', len(structural_LJ))

    print('\n\t Making the reverse duplicate dataframe')
    # Inverse pairs calvario
    inv_LJ = structural_LJ[['aj', 'ai', 'distance', 'sigma', 'epsilon', 'check', 'chain_ai', 'chain_aj', 'same_chain', 'exclude', 'type_ai', 'resnum_ai', 'type_aj', 'resnum_aj', 'diff']].copy()
    inv_LJ.columns = ['ai', 'aj', 'distance', 'sigma', 'epsilon', 'check', 'chain_ai', 'chain_aj', 'same_chain', 'exclude', 'type_ai', 'resnum_ai', 'type_aj', 'resnum_aj', 'diff']
    structural_LJ = structural_LJ.append(inv_LJ, sort = False, ignore_index = True)
    print('\t Doubled pairs list: ', len(structural_LJ))


    if sigma_method == 'minimum':

        print('\t Sorting and dropping all the duplicates')
        # Sorting the pairs
        structural_LJ.sort_values(by = ['ai', 'aj', 'distance'], inplace = True)
        # Cleaning the duplicates
        structural_LJ = structural_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
        # Removing the reverse duplicates
        cols = ['ai', 'aj']
        structural_LJ[cols] = np.sort(structural_LJ[cols].values, axis=1)
        structural_LJ = structural_LJ.drop_duplicates()
        print('\t Cleaning Complete ', len(structural_LJ))
    
    print('\n\t Calculating sigma and epsilon')
    structural_LJ['sigma'] = (structural_LJ['distance']/10) / (2**(1/6))
    structural_LJ['epsilon'] = epsilon_input #2.49 # PROTEIN CONFIGURATION # 0.41 epsilon MAGROS

    print('\n\n\t Sigma and epsilon completed ', len(structural_LJ))
    return structural_LJ


def merge_GRETA(native_pdb_pairs, fibril_pdb_pairs):
    # Merging native and fibril LJ pairs and cleaning all the duplicates among them
    if idp == True:
        # PROVA CON SOLO LA FIBRILLA forse il copy non e' necessario
        greta_LJ = fibril_pdb_pairs.copy()
    else:
        greta_LJ = native_pdb_pairs.append(fibril_pdb_pairs, sort = False, ignore_index = True)

    # Harp test, we don't have the fibril structure
    if protein == 'harp':
        greta_LJ = native_pdb_pairs.copy()

    # Sorting the pairs
    greta_LJ.sort_values(by = ['ai', 'aj', 'distance'], inplace = True)


    if sigma_method == 'NMR':

        # Alternative sigma distances choice
        check_pairs = greta_LJ['check'].to_list()
        check_pairs = list(dict.fromkeys(check_pairs))
        #print(len(greta_LJ['check']))
        #print(len(check_pairs))
        new_greta_LJ = []
        for c in check_pairs:
            filter = (greta_LJ['check'] == c)
            contact_subset = greta_LJ[filter]  
            contact_subset.sort_values(by = ['ai', 'aj', 'distance'], inplace = True)
            # Drop only the duplicates with the same distance due to the reverse copy of greta_LJ
            # Which used to be deleted before merge
            contact_subset = contact_subset.drop_duplicates(subset = ['distance'], keep = 'first')
            contact_subset.insert(3, 'new_sigma', 1)
            contact_subset.insert(3, 'new_distance', 1)
            contact_subset['new_distance'] = 1/(((sum((1/contact_subset['distance'])**6))/len(contact_subset))**(1/6))
            contact_subset['new_sigma'] = (contact_subset['new_distance']/10) / (2**(1/6))
            new_greta_LJ.append(contact_subset)
        cols = ['ai', 'aj']
        new_greta_LJ = pd.concat(new_greta_LJ, ignore_index=True)
        new_greta_LJ = new_greta_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
        new_greta_LJ[cols] = np.sort(new_greta_LJ[cols].values, axis=1)
        new_greta_LJ = new_greta_LJ.drop_duplicates()
    
    # Cleaning the duplicates
    greta_LJ = greta_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')

    # Removing the reverse duplicates
    cols = ['ai', 'aj']
    greta_LJ[cols] = np.sort(greta_LJ[cols].values, axis=1)
    greta_LJ = greta_LJ.drop_duplicates()
    
    if sigma_method == 'NMR':
        if len(new_greta_LJ) == len(greta_LJ):
            #\n\n\n\n\n\n NEW_GRETA_LJ == GRETA_LJ \n\n\n\n\n\n\n\n\n\n\n\n\n') 
            new_greta_LJ['sigma'] = new_greta_LJ['new_sigma']
            new_greta_LJ = new_greta_LJ.drop(columns = ['new_distance', 'new_sigma'])
            greta_LJ = new_greta_LJ.copy()
        else:
            print(len(new_greta_LJ))
            print(len(greta_LJ))
            print('New sigmas dont match with old ones')
            exit()

    #check_GRETA = greta_LJ[['ai', 'aj']].copy()
    greta_LJ.insert(2, 'type', 1)
    greta_LJ.insert(3, 'c12', '')
    greta_LJ['c12'] = 4 * greta_LJ['epsilon'] * (greta_LJ['sigma'] ** 12)
    greta_LJ.insert(3, 'c6', '')
    greta_LJ['c6'] = 4 * greta_LJ['epsilon'] * (greta_LJ['sigma'] ** 6)
    greta_LJ.insert(5, '', ';')
    greta_LJ.drop(columns = ['distance', 'check', 'chain_ai', 'chain_aj', 'same_chain', 'exclude', 'type_ai', 'resnum_ai', 'type_aj', 'resnum_aj', 'diff'], inplace = True)

    # SELF INTERACTIONS
    print('\n GRETA - Self interactions')
    atomtypes = set(greta_LJ['ai'])
    greta_LJ['double'] = ''

    print('\t Checking how many atoms do not self interact')
    for i in atomtypes:
        # Questo funziona e riesco a fargli dire quello che voglio.
        # Cioe' flaggare solo i valori che hanno un loro corrispettivo: N_1 N_1, CA_1 CA_1 ...
        greta_LJ.loc[(greta_LJ['ai'] == i) & (greta_LJ['aj'] == i), 'double'] = 'True'

    # Create a subset of the main dataframe of the self interactions.
    doubles = greta_LJ.loc[(greta_LJ['double'] == 'True')]
    atp_doubles = list(doubles['ai'])
    # The list is used to obtain all the atomtypes which does not make a self interaction
    atp_notdoubles = list(set(set(atomtypes) - set(atp_doubles)))
    atp_notdoubles.sort()

    if len(atp_notdoubles) == 0:
        print('\n\t All atoms interacts with themself')
        
    else:
        print('\n\t There are', len(atp_notdoubles), 'self interactions to add:\n\n\t')
        #print('\n\t There are', len(atp_notdoubles), 'self interactions to add:\n\n\t', atp_notdoubles, '\n')

        #print(doubles)
        #print(atp_doubles) # 84
        #print(atp_notdoubles) # 1 che in totale fanno 85, nel caso di TTR e' giusto
                                # perche' la prima riga l'ho tolta

        # From the list of atomtypes to add, a new dataframe is created to append to the main one
        pairs_toadd = pd.DataFrame(columns = ['ai', 'aj', 'type', 'c6', 'c12', 'sigma', 'epsilon'])
        pairs_toadd['ai'] = atp_notdoubles
        pairs_toadd['aj'] = atp_notdoubles
        pairs_toadd['type'] = '1'

        # Here i want to check every value for all the atom type and if they're similar
        # make an average and paste into the main dataframe
        # I am checking every doubles based on the atomtype (except the information of the residue number) and make an average of the sigma
        # since all the epsilon are equal
        atomtypes_toadd = pairs_toadd['ai'].str.split('_', n = 1, expand = True)
        atomtypes_toadd = atomtypes_toadd[0].drop_duplicates()
        atomtypes_toadd = atomtypes_toadd.to_list()
        atomtypes_toadd = [x + '_' for x in atomtypes_toadd]

        for a in atomtypes_toadd:
            doubles_a = doubles.loc[(doubles['ai'].str.contains(a)) & (doubles['aj'].str.contains(a))]
        
            # All the epsilon are the same, therefore the average sigma will be added on the self interaction
            sigma = doubles_a['sigma']
                        #for s in sigma:
            #    print('\t\t','sigma of ', a, '= ', s)
            
            if len(sigma) == 1:
                print('\n\t Only one self interacting pair has been found for', a, '==> Skip')

            elif len(sigma) == 0:
                print('\n\t There are not self interactions for', a, '==> Skip')

            else:
                print('\n\t There are', len(sigma), 'of', a, 'contributing to the average of self interacting sigma')
                media_sigma = sigma.mean()
                print('\n\t\t', 'Average Sigma for', a, ':', '\t', media_sigma)
                
                
                # Nuovi c6 e c12
                new_c6 = 4 * epsilon_input * (media_sigma ** 6)
                new_c12 = 4 * epsilon_input * (media_sigma ** 12)

                print('\t\t New c6 for ', a, '=\t', new_c6)
                print('\t\t New c12 for ', a, '=\t', new_c12)
                pairs_toadd.loc[(pairs_toadd['ai'].str.contains(a)) & (pairs_toadd['aj'].str.contains(a)), 'c6'] = new_c6#.mean()
                pairs_toadd.loc[(pairs_toadd['ai'].str.contains(a)) & (pairs_toadd['aj'].str.contains(a)), 'c12'] = new_c12#.mean()
                pairs_toadd.loc[(pairs_toadd['ai'].str.contains(a)) & (pairs_toadd['aj'].str.contains(a)), 'sigma'] = media_sigma
                pairs_toadd.loc[(pairs_toadd['ai'].str.contains(a)) & (pairs_toadd['aj'].str.contains(a)), 'epsilon'] = epsilon_input
        
        pairs_toadd.insert(5, '', ';')

        # Drop NaN: SD_1 SD_100 and OXT_100 -> in case of B2m
        pairs_toadd.dropna(inplace = True)
        greta_LJ = greta_LJ.append(pairs_toadd, sort = False, ignore_index = True)
        print('\n\t Self interactions added to greta_LJ\n')

    # Drop columns
    greta_LJ.drop(columns = ['double'], inplace = True)

    if idp == True:
        print('Addition of reweighted native pairs')
        # Here i join ai and aj of greta_LJ to compare with the monomer pairs
        pairs_check = (greta_LJ['ai'] + '_' + greta_LJ['aj']).to_list()
        native_pairs = read_native_pairs()
        native_pairs = native_pairs[native_pairs.ratio > ratio_treshold]
        native_pairs = native_pairs.replace({'ai':gro_to_amb_dict})
        native_pairs = native_pairs.replace({'aj':gro_to_amb_dict})
        native_pairs['pairs_check'] = native_pairs['ai'] + '_' + native_pairs['aj']
        # I keep only the one which are NOT included in pairs_check
        native_pairs = native_pairs[~native_pairs['pairs_check'].isin(pairs_check)]
        native_pairs['pairs_check'] = native_pairs['aj'] + '_' + native_pairs['ai']
        native_pairs = native_pairs[~native_pairs['pairs_check'].isin(pairs_check)]
        
        # Seems that all the native contacts are not included in the fibril, which is perfect!

        # Add sigma, add epsilon reweighted, add c6 and c12
        native_pairs['sigma'] = (native_pairs['distance']/10) / (2**(1/6))
        
        # Reweight 1
        #native_pairs['epsilon'] = native_pairs['ratio'] * epsilon_input

        # Reweight 2
        native_pairs['epsilon'] = epsilon_input*(1-((np.log(native_pairs['ratio']))/(np.log(ratio_treshold))))
        
        native_pairs['c12'] = 4 * native_pairs['epsilon'] * (native_pairs['sigma'] ** 12)   
        native_pairs['c6'] = 4 * native_pairs['epsilon'] * (native_pairs['sigma'] ** 6)
        native_pairs['type'] = 1
        native_pairs.drop(columns = ['counts', 'ratio', 'distance', 'pairs_check'], inplace = True)
        
        native_pairs = native_pairs[['ai', 'aj', 'type', 'c6', 'c12', 'sigma', 'epsilon']]
        native_pairs.insert(5, '', ';')
        greta_LJ = greta_LJ.append(native_pairs, ignore_index = True)

        # This sort and drop step is repeated since it might happen that native contacts might be duplicated with the fibril ones
        # This step has been added after the harp2 FF but ofc for harp2 FF we checked that all contacts were unique (comment above!)
        # 886 contacts before the addition of this second drop, 886 contacts (sorted differently) after this addition

        # Sorting the pairs
        greta_LJ.sort_values(by = ['ai', 'aj', 'sigma'], inplace = True)
        # Cleaning the duplicates
        greta_LJ = greta_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
        # Removing the reverse duplicates
        cols = ['ai', 'aj']
        greta_LJ[cols] = np.sort(greta_LJ[cols].values, axis=1)
        greta_LJ = greta_LJ.drop_duplicates()


    if N_terminal == True:
        print('Removing N_1 N_1 pair')
        greta_LJ.loc[(greta_LJ['ai'] == first_resid) & (greta_LJ['aj'] == first_resid), 'ai'] = ';'+greta_LJ['ai'] 

    greta_LJ = greta_LJ.rename(columns = {'ai':'; ai'})
    greta_LJ['sigma'] = greta_LJ["sigma"].map(lambda x:'{:.6e}'.format(x))
    greta_LJ['c6'] = greta_LJ["c6"].map(lambda x:'{:.6e}'.format(x))
    greta_LJ['c12'] = greta_LJ["c12"].map(lambda x:'{:.6e}'.format(x))
    print('\t GRETA FF COMPLETE: ', len(greta_LJ))

    if idp == False:
        return greta_LJ
    else:
        return greta_LJ, native_pairs


def make_pairs_exclusion_topology(greta_merge, type_c12_dict):
    '''
    This function prepares the [ exclusion ] and [ pairs ] section to paste in topology.top
    '''

    greta_merge = greta_merge.rename(columns = {'; ai': 'ai'})
    atnum_type_top = topology_atoms[['; nr', 'type']]
    atnum_type_top = atnum_type_top.rename(columns = {'; nr': 'nr'})
    atnum_type_dict = atnum_type_top.set_index('type')['nr'].to_dict()

    atnum_topology_bonds = topology_bonds.copy()
    atnum_topology_bonds['ai'] = atnum_topology_bonds['ai'].map(atnum_type_dict)
    atnum_topology_bonds['aj'] = atnum_topology_bonds['aj'].map(atnum_type_dict)
    atnum_topology_bonds['ai'] = atnum_topology_bonds['ai'].astype(int)
    atnum_topology_bonds['aj'] = atnum_topology_bonds['aj'].astype(int)

    bond_tuple = list(map(tuple, atnum_topology_bonds.to_numpy()))
    #print(bond_tuple)

    ex, exclusion_bonds = [], []
    for atom in atom_topology_num:
        for t in bond_tuple:
            if t[0] == atom:
                first = t[1]
                ex.append(t[1])
            elif t[1] == atom:
                first = t[0]
                ex.append(t[0])
            else: continue
            for tt in bond_tuple:
                if (tt[0] == first) & (tt[1] != atom):
                    second = tt[1]
                    ex.append(tt[1])
                elif (tt[1] == first) & (tt[0] != atom):
                    second = tt[0]
                    ex.append(tt[0])
                else: continue
                for ttt in bond_tuple:
                    if (ttt[0] == second) & (ttt[1] != first):
                        ex.append(ttt[1])
                    elif (ttt[1] == second) & (ttt[0] != first):
                        ex.append(ttt[0])
        for e in ex:
            exclusion_bonds.append((str(str(atom) + '_' + str(e))))
            exclusion_bonds.append((str(str(e) + '_' + str(atom))))
        ex = []
    

# Questa si puo prendere direttamente durante il merge per evitare di fare calcoli ridondanti
    pairs = greta_merge[['ai', 'aj']].copy()
    pairs['c12_ai'] = pairs['ai']
    pairs['c12_aj'] = pairs['aj']

    # Keeping based on resnum
    pairs[['type_ai', 'resnum_ai']] = pairs.ai.str.split("_", expand = True)
    pairs[['type_aj', 'resnum_aj']] = pairs.aj.str.split("_", expand = True)
    pairs['resnum_ai'] = pairs['resnum_ai'].astype(int)
    pairs['resnum_aj'] = pairs['resnum_aj'].astype(int)
    pairs = pairs.loc[(abs(pairs['resnum_aj'] - pairs['resnum_ai']) < distance_residue)]
    pairs = pairs[pairs['ai'] != pairs['aj']]

    pairs['ai'] = pairs['ai'].map(atnum_type_dict)
    pairs['aj'] = pairs['aj'].map(atnum_type_dict)

    pairs['check'] = pairs['ai'] + '_' + pairs['aj']

    # Here we keep only the one without the exclusions
    pairs['exclude'] = ''
    pairs.loc[(pairs['check'].isin(exclusion_bonds)), 'exclude'] = 'Yes' 
    mask = pairs.exclude == 'Yes'
    pairs = pairs[~mask]
    


    pairs['c12_ai'] = pairs['c12_ai'].map(type_c12_dict)
    pairs['c12_aj'] = pairs['c12_aj'].map(type_c12_dict)
    pairs['func'] = 1
    pairs['c6'] = 0.00000e+00
    pairs['c6'] = pairs["c6"].map(lambda x:'{:.6e}'.format(x))
    pairs['c12'] = np.sqrt(pairs['c12_ai'] * pairs['c12_aj'])
    pairs['c12'] = pairs["c12"].map(lambda x:'{:.6e}'.format(x))

    pairs.drop(columns = ['type_ai', 'resnum_ai', 'type_aj', 'resnum_aj', 'c12_ai', 'c12_aj', 'check', 'exclude'], inplace = True)
    
    exclusion = pairs[['ai', 'aj']].copy()
    
    pairs = pairs.rename(columns = {'ai': '; ai'})
    exclusion = exclusion.rename(columns = {'ai': '; ai'})

    return pairs, exclusion
    
########################## DIHEDRALS

#def sb_dihedrals (structure_pdb):
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
    native_dihedrals['mult'] = ''  # manca questa
    #print(native_dihedrals)
    return native_dihedrals