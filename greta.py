from operator import concat

from MDAnalysis.lib.util import parse_residue
from read_input import read_native_pairs
from MDAnalysis.analysis import distances
import numpy as np
from pandas.core.frame import DataFrame
import pandas as pd
import itertools
from protein_configuration import distance_cutoff, distance_residue, epsilon_input, ratio_treshold, lj_reduction, multiply_c6
from topology_definitions import raw_topology_atoms, gromos_atp, gromos_atp_c6, gro_to_amb_dict, topology_bonds, atom_topology_num

# TODO
#def make_pdb_atomtypes (native_pdb, *fibril_pdb):
def make_pdb_atomtypes (native_pdb, fibril_pdb):
    '''
    This function defines the SB based atomtypes to add in topology.top, atomtypes.atp and ffnonbonded.itp.
    '''
    topology_atoms = raw_topology_atoms.copy()
    print('\tBuilding native atomtypes')
    native_sel = native_pdb.select_atoms('all')
    native_atomtypes, ffnb_sb_type = [], []

    for atom in native_sel:
        '''
        This for loop reads the coordinates file and starts building the SB atomtypes.
        The native atomtypes will be used for topology.top, atomtypes.atp and ffnonbonded.itp.
        We define atps as atomname_resnumber:chainnumber, as we need this three informations to create atom pairs.
        '''
        # TODO add the exclusion
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
    
    # Just checking that the pdb and the topology have the same number of atoms
    if len(np.unique(check_topology.check)) == 1:
        print('\n\tAtoms of topol.top and pdb have the same number')
    else:
        print('\n\tCheck PDB and topology numeration')
        exit()
        
    print('\tBuilding fibril atomtypes')
    fibril_sel = fibril_pdb.select_atoms('all')
    fibril_atomtypes = []
    for atom in fibril_sel:
        atp = str(atom.name) + '_' + str(atom.resnum) + ':' + str(atom.segid)
        fibril_atomtypes.append(atp)

    # ffnonbonded making
    # Making a dictionary with atom number and type
    print('\tFFnonbonded atomtypes section creation')
    ffnb_atomtype = pd.DataFrame(columns = ['; type', 'chem', 'at.num', 'mass', 'charge', 'ptype', 'c6', 'c12'])
    ffnb_atomtype['; type'] = topology_atoms['sb_type']
    ffnb_atomtype['chem'] = topology_atoms['type']
    ffnb_atomtype['at.num'] = ffnb_atomtype['chem'].map(gromos_atp['at.num'])
    ffnb_atomtype['mass'] = topology_atoms['mass']
    ffnb_atomtype['charge'] = '0.000000'
    ffnb_atomtype['ptype'] = 'A'
    ffnb_atomtype['c6'] = '0.00000e+00'
    ffnb_atomtype['c12'] = ffnb_atomtype['chem'].map(gromos_atp['c12'])
    
    # This will be used to check if there are prolines in the structure and half their N c12
    residue_list = topology_atoms['residue'].to_list()

    # This is used to remove the proline N interactions in Pairs and Exclusions
    proline_n = []
    if 'PRO' in residue_list:
        print('\tThere are prolines in the structure. The c12 of N should be the half')
        proline_n = topology_atoms.loc[(topology_atoms['residue'] == 'PRO') & (topology_atoms['atom'] == 'N'), 'nr'].to_list()
    #    ffnb_atomtype.loc[(ffnb_atomtype['; type'].isin(proline_n)), 'c12'] = ffnb_atomtype['c12']/20
    #else:
    #    print('\tThere not are prolines in the structure. The c12 of N should be the half')
    

    # The N terminal of the structure should be bigger than the others since it has an H more and charged
    #if N_terminal == True:
    #    print('\nChanging the c12 value of N-terminal')
        # In this case we multiply by 5 since the c12 are already doubled
    #    ffnb_atomtype.loc[(ffnb_atomtype['; type'] == first_resid), 'c12'] = ffnb_atomtype['c12']*5 # Harp 2
        

    # This will be needed for exclusion and pairs to paste in topology
    # A dictionary with the c12 of each atom in the system
    type_c12_dict = ffnb_atomtype.set_index('; type')['c12'].to_dict()
    
    ffnb_atomtype['c12'] = ffnb_atomtype["c12"].map(lambda x:'{:.6e}'.format(x))
    ffnb_atomtype.drop(columns = ['chem'], inplace = True)

    print('\tTopology atomtypes section creation')
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

    return native_atomtypes, fibril_atomtypes, ffnb_atomtype, atomtypes_atp, topology_atoms, type_c12_dict, proline_n


def make_pairs(structure_pdb, atomtypes):
    '''
    This function measures all the distances between all atoms using MDAnalysis.
    It works on both native and fibril in the same manner.
    '''

    print('\n\tMeasuring distances between all atom in the pdb')
    # Selecting all atoms in the system
    atom_sel = structure_pdb.select_atoms('all')

    # Calculating all the distances between atoms.
    # The output is a combination array.
    self_distances = distances.self_distance_array(atom_sel.positions)
    print('\tNumber of distances measured :', len(self_distances))
    
    print('\n\tPreparing the atomtype array')
    # The MDAnalysis contains only distances, so we rebuilt atom pairs in the same manner
    # using the atomtypes list of native and fibril which will match with the distance array.

    # TODO create directly the two separated lists
    pairs_list = list(itertools.combinations(atomtypes, 2))

    # But the combinations are list of list and we need to separate them.
    pairs_ai, pairs_aj = [], []
    for n in range(0, len(pairs_list)):
        i = pairs_list[n][0]
        pairs_ai.append(i)
        j = pairs_list[n][1]
        pairs_aj.append(j)
    print('\tAtomtype array ready')

    print('\n\tCreating the pairs dataframes')
    # Creation of the dataframe containing the atom pairs and the distances.
    # Also, it will be prepared for sigma and epsilon.
    structural_LJ = pd.DataFrame(columns = ['ai', 'aj', 'distance', 'sigma', 'epsilon'])
    structural_LJ['ai'] = pairs_ai
    structural_LJ['aj'] = pairs_aj
    structural_LJ['distance'] = self_distances
    raw_structural_LJ = len(structural_LJ)
    print('\tRaw pairs list ', raw_structural_LJ)
    
    print(f'\n\tApplying distance cutoff of {distance_cutoff} A')
    # Keep only the atoms within cutoff
    structural_LJ = structural_LJ[structural_LJ.distance < distance_cutoff] # PROTEIN CONFIGURATION
    print(f'\tPairs below cutoff {distance_cutoff}: ', len(structural_LJ))
    print(f'\tDeleted {raw_structural_LJ - len(structural_LJ)} pairs')

    print('\n\tExtracting chain label for every atom pair')
    # That is name_resname:resid made from the previous function.
    # Extracting the resid information to check if the atom pair is on the same chain.
    structural_LJ[['ai', 'chain_ai']] = structural_LJ.ai.str.split(":", expand = True)
    structural_LJ[['aj', 'chain_aj']] = structural_LJ.aj.str.split(":", expand = True)
    
    structural_LJ['same_chain'] = np.where(structural_LJ['chain_ai'] == structural_LJ['chain_aj'], 'Yes', 'No')
    
    print('\tPairs within the same chain: ', len(structural_LJ.loc[structural_LJ['same_chain'] == 'Yes']))
    print('\tPairs not in the same chain: ', len(structural_LJ.loc[structural_LJ['same_chain'] == 'No']))
    print('\tRaw pairs list ', len(structural_LJ))
    print(f'\tApplying residue number cutoff of {distance_residue}')

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
    print(f'\tAll the pairs further than {distance_residue} aminoacids and not in the same chain: ', len(structural_LJ))

    #print('\n\tMaking the reverse duplicate dataframe')
    ## Inverse pairs calvario
    #inv_LJ = structural_LJ[['aj', 'ai', 'distance', 'sigma', 'epsilon', 'chain_ai', 'chain_aj', 'same_chain', 'type_ai', 'resnum_ai', 'type_aj', 'resnum_aj', 'diff']].copy()
    #inv_LJ.columns = ['ai', 'aj', 'distance', 'sigma', 'epsilon', 'chain_ai', 'chain_aj', 'same_chain', 'type_ai', 'resnum_ai', 'type_aj', 'resnum_aj', 'diff']
    #structural_LJ = structural_LJ.append(inv_LJ, sort = False, ignore_index = True)
    #print('\tDoubled pairs list: ', len(structural_LJ))

    ## Here we sort all the atom pairs based on the distance and we keep the closer ones.
    #print('\tSorting and dropping all the duplicates')
    ## Sorting the pairs
    #structural_LJ.sort_values(by = ['ai', 'aj', 'distance'], inplace = True)
    ## Cleaning the duplicates
    #structural_LJ = structural_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    ## Removing the reverse duplicates
    #cols = ['ai', 'aj']
    #structural_LJ[cols] = np.sort(structural_LJ[cols].values, axis=1)
    #structural_LJ = structural_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    #print('\tCleaning Complete ', len(structural_LJ))
    
    print('\n\tCalculating sigma and epsilon')
    structural_LJ['sigma'] = (structural_LJ['distance']/10) / (2**(1/6))
    # As declared in protein_configuration.py
    structural_LJ['epsilon'] = epsilon_input

    print('\n\n\tSigma and epsilon completed ', len(structural_LJ))
    structural_LJ.drop(columns = ['chain_ai', 'chain_aj', 'same_chain', 'type_ai', 'resnum_ai', 'type_aj', 'resnum_aj', 'diff'], inplace = True)

    return structural_LJ


def make_idp_epsilon():
    #TODO Vedere se mettere questa parte all'inizio del merge
    # Controlla con Carlo

    # In the case of an IDP, it is possible to add dynamical informations based on a simulation
    print('Addition of reweighted native pairs')
    # Here i join ai and aj of greta_LJ to compare with the monomer pairs
    #pairs_check = (greta_LJ['ai'] + '_' + greta_LJ['aj']).to_list()
    
    native_pairs = read_native_pairs()
    # The ratio treshold considers only pairs occurring at a certain probability
    native_pairs = native_pairs[native_pairs.ratio > ratio_treshold]
    # This dictionary was made to link amber and greta atomtypes
    native_pairs = native_pairs.replace({'ai':gro_to_amb_dict})
    native_pairs = native_pairs.replace({'aj':gro_to_amb_dict})
    
    # Add sigma, add epsilon reweighted, add c6 and c12
    native_pairs['sigma'] = (native_pairs['distance']/10) / (2**(1/6))
    
    # Epsilon reweight based on probability
    native_pairs['epsilon'] = epsilon_input#*(1-((np.log(native_pairs['ratio']))/(np.log(ratio_treshold))))
    native_pairs.drop(columns = ['counts', 'ratio'], inplace = True)

    return native_pairs


def merge_GRETA(greta_LJ):
    '''
    This function merges the atom contacts from native and fibril.
    '''
    
    
    # Inverse pairs calvario
    inv_LJ = greta_LJ[['aj', 'ai', 'distance', 'sigma', 'epsilon']].copy()
    inv_LJ.columns = ['ai', 'aj', 'distance', 'sigma', 'epsilon']
    greta_LJ = greta_LJ.append(inv_LJ, sort = False, ignore_index = True)
    print('\tDoubled pairs list: ', len(greta_LJ))

    # Here we sort all the atom pairs based on the distance and we keep the closer ones.
    print('\tSorting and dropping all the duplicates')
    # Sorting the pairs
    greta_LJ.sort_values(by = ['ai', 'aj', 'distance'], inplace = True)
    # Cleaning the duplicates
    greta_LJ = greta_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    # Removing the reverse duplicates
    cols = ['ai', 'aj']
    greta_LJ[cols] = np.sort(greta_LJ[cols].values, axis=1)
    greta_LJ = greta_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    print('\tCleaning Complete ', len(greta_LJ))



    ## Sorting the pairs
    #greta_LJ.sort_values(by = ['ai', 'aj', 'distance'], inplace = True)
    ## Cleaning the duplicates
    #greta_LJ = greta_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    ## Removing the reverse duplicates
    #cols = ['ai', 'aj']
    #greta_LJ[cols] = np.sort(greta_LJ[cols].values, axis=1)
    #greta_LJ = greta_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    
    
    
    
    
    
    
    greta_LJ.insert(2, 'type', 1)
    greta_LJ.insert(3, 'c12', '')
    greta_LJ['c12'] = 4 * greta_LJ['epsilon'] * (greta_LJ['sigma'] ** 12)
    greta_LJ.insert(3, 'c6', '')
    greta_LJ['c6'] = 4 * greta_LJ['epsilon'] * (greta_LJ['sigma'] ** 6)
    greta_LJ.insert(5, '', ';')
    greta_LJ.drop(columns = ['distance'], inplace = True)

    # SELF INTERACTIONS
    # In the case of fibrils which are not fully modelled we add self interactions which is a feature of amyloids
    # So that the balance between native and fibril is less steep.
    print('\n GRETA - Self interactions')
    atomtypes = set(greta_LJ['ai'])
    greta_LJ['double'] = ''

    print('\tChecking how many atoms do not self interact')
    for i in atomtypes:
        # Selection of already known atoms which contacts with themself
        greta_LJ.loc[(greta_LJ['ai'] == i) & (greta_LJ['aj'] == i), 'double'] = 'True'

    # Create a subset of the main dataframe of the self interactions.
    doubles = greta_LJ.loc[(greta_LJ['double'] == 'True')]
    atp_doubles = list(doubles['ai'])
    # The list is used to obtain all the atomtypes which does not make a self interaction
    atp_notdoubles = list(set(set(atomtypes) - set(atp_doubles)))
    atp_notdoubles.sort()

    if len(atp_notdoubles) == 0:
        print('\n\tAll atoms interacts with themself')
        
    else:
        print('\n\tThere are', len(atp_notdoubles), 'self interactions to add')
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
            # Selects the atom pairs from the double pairs 
            doubles_a = doubles.loc[(doubles['ai'].str.contains(a)) & (doubles['aj'].str.contains(a))]
            # All the epsilon are the same, therefore the average sigma will be added on the self interaction
            sigma = doubles_a['sigma']
            
            if len(sigma) == 1:
                # If there is only onw sigma for the averages it will be skipped
                print('\n\t Only one self interacting pair has been found for', a, '==> Skip')
            elif len(sigma) == 0:
                # If the missing atom pairs is not represented in the strcture there are not
                # sigmas to average
                print('\n\t There are not self interactions for', a, '==> Skip')
            else:
                # If there are enough sigmas to make an average then it creates the missing atom pairs
                print('\n\t There are', len(sigma), 'of', a, 'contributing to the average of self interacting sigma')
                media_sigma = sigma.mean()
                print('\n\t\t', 'Average Sigma for', a, ':', '\t', media_sigma)
                
                # Creation of new c6 and c12
                new_c6 = 4 * epsilon_input * (media_sigma ** 6)
                new_c12 = 4 * epsilon_input * (media_sigma ** 12)

                print('\t\t New c6 for ', a, '=\t', new_c6)
                print('\t\t New c12 for ', a, '=\t', new_c12)

                # In the pairs to add dataframe all those new information are inserted

                pairs_toadd.loc[(pairs_toadd['ai'].str.contains(a)) & (pairs_toadd['aj'].str.contains(a)), 'c6'] = new_c6
                pairs_toadd.loc[(pairs_toadd['ai'].str.contains(a)) & (pairs_toadd['aj'].str.contains(a)), 'c12'] = new_c12
                pairs_toadd.loc[(pairs_toadd['ai'].str.contains(a)) & (pairs_toadd['aj'].str.contains(a)), 'sigma'] = media_sigma
                pairs_toadd.loc[(pairs_toadd['ai'].str.contains(a)) & (pairs_toadd['aj'].str.contains(a)), 'epsilon'] = epsilon_input

        pairs_toadd.insert(5, '', ';')
        pairs_toadd.dropna(inplace = True)
        # Appending the missing atom pairs to the main dataframe
        greta_LJ = greta_LJ.append(pairs_toadd, sort = False, ignore_index = True)
        print('\n\t Self interactions added to greta_LJ\n')

    # Drop double, we don't need it anymore
    greta_LJ.drop(columns = ['double'], inplace = True)
    greta_LJ = greta_LJ.rename(columns = {'ai':'; ai'})
    greta_LJ['epsilon'] = greta_LJ["epsilon"].map(lambda x:'{:.6f}'.format(x))
    greta_LJ['sigma'] = greta_LJ["sigma"].map(lambda x:'{:.6e}'.format(x))
    greta_LJ['c6'] = greta_LJ["c6"].map(lambda x:'{:.6e}'.format(x))
    greta_LJ['c12'] = greta_LJ["c12"].map(lambda x:'{:.6e}'.format(x))
    print('\t GRETA FF COMPLETE: ', len(greta_LJ))

    return greta_LJ


def make_pairs_exclusion_topology(greta_merge, type_c12_dict, proline_n):
    '''
    This function prepares the [ exclusion ] and [ pairs ] section to paste in topology.top
    Here we define the GROMACS exclusion list and drop from the LJ list made using GRETA so that all the remaining
    contacts will be defined in pairs and exclusions as particular cases.
    Since we are not defining explicit H, the 1-4 list is defined by 2 bonds and not 3 bonds.
    '''


    topology_atoms = raw_topology_atoms.copy()

    greta_merge = greta_merge.rename(columns = {'; ai': 'ai'})
    #atnum_type_top = topology_atoms[['; nr', 'type']]
    atnum_type_top = topology_atoms[['nr', 'sb_type', 'resnr', 'atom', 'type', 'residue']]
    #atnum_type_top = atnum_type_top.rename(columns = {'; nr': 'nr'})
    atnum_type_top['resnr'] = atnum_type_top['resnr'].astype(int)

    # Dictionaries definitions to map values
    atnum_type_dict = atnum_type_top.set_index('sb_type')['nr'].to_dict()
    type_atnum_dict = atnum_type_top.set_index('nr')['sb_type'].to_dict()

    # Bonds from topology
    atnum_topology_bonds = topology_bonds.copy()
    atnum_topology_bonds['ai'] = atnum_topology_bonds['ai'].map(atnum_type_dict)
    atnum_topology_bonds['aj'] = atnum_topology_bonds['aj'].map(atnum_type_dict)
    atnum_topology_bonds['ai'] = atnum_topology_bonds['ai'].astype(int)
    atnum_topology_bonds['aj'] = atnum_topology_bonds['aj'].astype(int)
    bond_tuple = list(map(tuple, atnum_topology_bonds.to_numpy()))


    #TODO this should be in topology_definitions.py

    # Building the exclusion bonded list
    ex, ex14, p14, exclusion_bonds = [], [], [], []
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
                        ex14.append(ttt[1])

                    elif (ttt[1] == second) & (ttt[0] != first):
                        ex.append(ttt[0])
                        ex14.append(ttt[0])
        for e in ex:
            exclusion_bonds.append((str(str(atom) + '_' + str(e))))
            exclusion_bonds.append((str(str(e) + '_' + str(atom))))
        ex = []
        for e in ex14:
            p14.append((str(str(atom) + '_' + str(e))))
            p14.append((str(str(e) + '_' + str(atom))))
        ex14 = []
    

    # TODO Questa si puo prendere direttamente durante il merge per evitare di fare calcoli ridondanti
    pairs = greta_merge[['ai', 'aj']].copy()
    pairs['c12_ai'] = pairs['ai']
    pairs['c12_aj'] = pairs['aj']
    pairs[['type_ai', 'resnum_ai']] = pairs.ai.str.split("_", expand = True)
    pairs[['type_aj', 'resnum_aj']] = pairs.aj.str.split("_", expand = True)
    pairs['resnum_ai'] = pairs['resnum_ai'].astype(int)
    pairs['resnum_aj'] = pairs['resnum_aj'].astype(int)
    
    # We keep the pairs we dropped from the make_pairs: those are from the fibril interactions exclusively
    # The incoming pairs are obtained from atoms close in space. Some are from different chains but from residues close in sequence
    # Here we keep pairs between residues close in sequence but from different chains, thus close in space
    # If we keep such pairs we cause severe frustration to the system and artifacts
    pairs = pairs.loc[(abs(pairs['resnum_aj'] - pairs['resnum_ai']) < distance_residue)]
    
    # We remove the contact with itself
    pairs = pairs[pairs['ai'] != pairs['aj']]

    # The exclusion list was made based on the atom number
    pairs['ai'] = pairs['ai'].map(atnum_type_dict)
    pairs['aj'] = pairs['aj'].map(atnum_type_dict)
    pairs['check'] = pairs['ai'] + '_' + pairs['aj']

    # Here the drop the contacts which are already defined by GROMACS, including the eventual 1-4 exclusion defined in the LJ_pairs
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
    pairs.drop(columns = ['type_ai', 'resnum_ai', 'type_aj', 'resnum_aj', 'c12_ai', 'c12_aj', 'check', 'exclude'], inplace = True)    

    # Only 1-4 exclusions are fully reintroduced
    pairs_14 = pd.DataFrame(columns=['ai', 'aj', 'exclusions'])
    pairs_14['exclusions'] = p14
    pairs_14[['ai', 'aj']] = pairs_14.exclusions.str.split("_", expand = True)

    # TODO migliora questo passaggio mettendo direttamente il map
    # TODO sistemare direttamente c12_ax con type e poi chiamare la colonna quando serve con quel nome
    pairs_14['c12_ai'] = pairs_14['ai']
    pairs_14['c12_aj'] = pairs_14['aj']
    pairs_14['c12_ai'] = pairs_14['c12_ai'].map(type_atnum_dict)
    pairs_14['c12_aj'] = pairs_14['c12_aj'].map(type_atnum_dict)

    # Adding an atom column because we want to flag NOT N N interactions, this because the N should have an explicit H we removed
    # TODO chiedere a Carlo perche stiamo considerando solo le interazioni con N e non con O ad esempio, o tutti gli altri
    pairs_14[['ai_type', 'ai_resid']] = pairs_14.c12_ai.str.split("_", expand = True)
    pairs_14[['aj_type', 'aj_resid']] = pairs_14.c12_aj.str.split("_", expand = True)

    # NOT 1_4 N X interactions will be dropped
    # Only the pairs with an N involved are retained
    pairs_14.loc[(pairs_14['ai_type'] == 'N') | (pairs_14['aj_type'] == 'N'), 'c12_tozero'] = False
    pairs_14.drop(pairs_14[pairs_14.c12_tozero != False].index, inplace=True)

    # TODO da mettere sotto
    # Here we take a particular interest of the interaction between two N, because both should have an explicit H
    pairs_14.loc[(pairs_14['ai_type'] == 'N') & (pairs_14['aj_type'] == 'N'), 'c12_tozero'] = True

    # Thus, only N with X LJ 1_4 interactions will be kept
    # All the other 1_4 interactions will NOT interact with each others
    pairs_14['c12_ai'] = pairs_14['c12_ai'].map(type_c12_dict)
    pairs_14['c12_aj'] = pairs_14['c12_aj'].map(type_c12_dict)
    pairs_14['func'] = 1
    pairs_14['c6'] = 0.00000e+00
    pairs_14['c6'] = pairs_14["c6"].map(lambda x:'{:.6e}'.format(x))
    pairs_14['c12'] = (np.sqrt(pairs_14['c12_ai'] * pairs_14['c12_aj']))*lj_reduction
    
    # The N N interactions are still too close, so we double the c12
    pairs_14.loc[(pairs_14['c12_tozero'] == True), 'c12'] *= 2



    pairs_14.drop(columns = ['exclusions', 'c12_ai', 'c12_aj', 'ai_type', 'ai_resid','aj_type', 'aj_resid', 'c12_tozero'], inplace = True)    

    # Exclusions 1-4
    pairs = pairs.append(pairs_14)

    # Drop duplicates
    pairs.sort_values(by = ['ai', 'aj', 'c12'], inplace = True)
    # Cleaning the duplicates
    pairs = pairs.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')

    # Removing the reverse duplicates
    cols = ['ai', 'aj']
    pairs[cols] = np.sort(pairs[cols].values, axis=1)
    pairs = pairs.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    pairs['c12'] = pairs["c12"].map(lambda x:'{:.6e}'.format(x))
    pairs.sort_values(by = ['ai', 'aj'], inplace = True)
    exclusion = pairs[['ai', 'aj']].copy()
    pairs = pairs.rename(columns = {'ai': '; ai'})
    exclusion = exclusion.rename(columns = {'ai': '; ai'})
    

    # Left alpha helix
    # Here we introduce the c6 and the c12 between the O and the CB of the next residue
    # to optimize the left alpha region in the Ramachandran

    # Adding the c6 and c12 (I do it now because later is sbatti)
    atnum_type_top['c6'] = atnum_type_top['type'].map(gromos_atp_c6['c6'])
    atnum_type_top['c12'] = atnum_type_top['sb_type'].map(type_c12_dict)

    # Here we make a dictionary of the backbone oxygen as atom number
    backbone_oxygen = atnum_type_top.loc[atnum_type_top['atom'] == 'O']
    sidechain_cb = atnum_type_top.loc[atnum_type_top['atom'] == 'CB']
    # Left alphas do not occur in GLY and PRO
    sidechain_cb = sidechain_cb[sidechain_cb.residue != 'PRO']
    sidechain_cb = sidechain_cb[sidechain_cb.residue != 'GLY']

    # For each backbone oxygen take the CB of the next residue and save in a pairs tuple
    left_alpha_ai, left_alpha_aj, left_alpha_c6, left_alpha_c12 = [], [], [], []
    for index, line_backbone_oxygen in backbone_oxygen.iterrows():
        line_sidechain_cb = sidechain_cb.loc[sidechain_cb['resnr'] == (line_backbone_oxygen['resnr'])+1].squeeze(axis=None)
        if not line_sidechain_cb.empty:
            left_alpha_ai.append(line_backbone_oxygen['nr'])
            left_alpha_aj.append(line_sidechain_cb['nr'])
            left_alpha_c6.append(np.sqrt(line_backbone_oxygen['c6']*line_sidechain_cb['c6'])*multiply_c6)
            left_alpha_c12.append(np.sqrt(line_backbone_oxygen['c12']*line_sidechain_cb['c12']))


    left_alpha_pairs = pd.DataFrame(columns=['ai', 'aj', 'c6', 'c12'])
    left_alpha_pairs['ai'] = left_alpha_ai
    left_alpha_pairs['aj'] = left_alpha_aj
    left_alpha_pairs['c6'] = left_alpha_c6
    left_alpha_pairs['c12'] = left_alpha_c12
    left_alpha_pairs['c6'] = left_alpha_pairs["c6"].map(lambda x:'{:.6e}'.format(x))
    left_alpha_pairs['c12'] = left_alpha_pairs["c12"].map(lambda x:'{:.6e}'.format(x))
    left_alpha_pairs['func'] = 1

    left_alpha_exclusions = left_alpha_pairs[['ai', 'aj']].copy()
    left_alpha_pairs = left_alpha_pairs.rename(columns = {'ai': '; ai'})
    left_alpha_exclusions = left_alpha_exclusions.rename(columns = {'ai': '; ai'})


    pairs = pairs.append(left_alpha_pairs)
    exclusion = exclusion.append(left_alpha_exclusions)

    # Removing the interactions with the prolin N

    pairs = pairs[~pairs['; ai'].isin(proline_n)]
    pairs = pairs[~pairs['aj'].isin(proline_n)]

    exclusion = exclusion[~exclusion['; ai'].isin(proline_n)] 
    exclusion = exclusion[~exclusion['aj'].isin(proline_n)]


    return pairs, exclusion


    # 3 bonds versione
    ##TODO this should be in topology_definitions.py
    #ex, exclusion_bonds = [], []
    #for atom in atom_topology_num:
    #    for t in bond_tuple:
    #        if t[0] == atom:
    #            first = t[1]
    #        elif t[1] == atom:
    #            first = t[0]
    #            ex.append(t[0])
    #        else: continue
    #        for tt in bond_tuple:
    #            if (tt[0] == first) & (tt[1] != atom):
    #                second = tt[1]
    #                ex.append(tt[1])
    #            elif (tt[1] == first) & (tt[0] != atom):
    #                second = tt[0]
    #                ex.append(tt[0])
    #            else: continue
    #            for ttt in bond_tuple:
    #                if (ttt[0] == second) & (ttt[1] != first):
    #                    ex.append(ttt[1])
    #                elif (ttt[1] == second) & (ttt[0] != first):
    #                    ex.append(ttt[0])
    #    for e in ex:
    #        exclusion_bonds.append((str(str(atom) + '_' + str(e))))
    #        exclusion_bonds.append((str(str(e) + '_' + str(atom))))
    #    ex = []
    #