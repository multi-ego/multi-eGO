from operator import concat
from MDAnalysis.lib.util import parse_residue
from MDAnalysis.analysis import distances
import numpy as np
from pandas.core.frame import DataFrame
import pandas as pd
import itertools
pd.options.mode.chained_assignment = None  # default='warn'


gromos_atp = pd.DataFrame(
    {'name': ['O', 'OA', 'N', 'C', 'CH1', 
            'CH2', 'CH3', 'CH2r', 'NT', 'S',
            'NR', 'OM', 'NE', 'NL', 'NZ'],
     'at.num': [8, 8, 7, 6, 6, 6, 6, 6, 7, 16, 7, 8, 7, 7, 7],
     'c12': [1e-06, 1.505529e-06, 2.319529e-06, 4.937284e-06, 9.70225e-05, # CH1
            3.3965584e-05, 2.6646244e-05, 2.8058209e-05, 5.0625e-06, 1.3075456e-05,
            3.389281e-06, 7.4149321e-07, 2.319529e-06, 2.319529e-06, 2.319529e-06],
     # here the 0 should be corrected with the correct c6 (anyway they are not used now)
     'c6': [0.0022619536, 0, 0, 0, 0.00606841, 0.0074684164, 0.0096138025, 0, 0, 0, 0, 0, 0, 0, 0]
     }
)
gromos_atp.to_dict()
gromos_atp.set_index('name', inplace=True)

# TODO make it more general (without the residue number)
# native reweight for TTR and ABeta. This dictionary will rename the amber topology to gromos topology
#gro_to_amb_dict = {'OC1_11' : 'O1_11', 'OC2_11':'O2_11'}
gro_to_amb_dict = {'OT1_42' : 'O1_42', 'OT2_42':'O2_42'}

def make_pdb_atomtypes(native_pdb, topology_atoms, parameters):
    '''
    This function prepares the ffnonbonded.itp section of Multi-ego.ff.
    The topology of the plainMD is read, and custom c12 are added.
    The native .pdb is read and a list of the atoms is prepared to use in PDB_LJ_pairs.
    It also provides a dictionary based on atom_type and c12 used in make_pairs_exclusion_topology.
    '''
    native_sel = native_pdb.select_atoms('all')
    native_atomtypes, ffnb_sb_type = [], []

    for atom in native_sel:
        '''
        This loop is required for the atom list to be used in PDB_LJ_pairs.
        '''
        # Native atomtypes will be used to create the pairs list
        #TODO print another for check
        atp = str(atom.name) + '_' + str(atom.resnum) + ':' + str(atom.segid)
        native_atomtypes.append(atp)

        # This part is for attaching to FFnonbonded.itp
        # ffnb_sb_type is the first column
        # check gromologist
        sb_type = str(atom.name) + '_' + str(atom.resnum)
        ffnb_sb_type.append(sb_type)

    check_topology = DataFrame(ffnb_sb_type, columns=['sb_type'])
    check_topology = check_topology.drop_duplicates(subset=['sb_type'])
    check_topology['check'] = np.where(topology_atoms.sb_type == check_topology.sb_type, 'same', 'different')
    
    # Just checking that the pdb and the topology have the same number of atoms
    if len(np.unique(check_topology.check)) != 1:
        print('\n\tCheck PDB and topology because they have different numbers of atoms')
        exit()
        
    # ffnonbonded making
    # Making a dictionary with atom number and type
    ffnb_atomtype = pd.DataFrame(columns = ['; type', 'chem', 'at.num', 'mass', 'charge', 'ptype', 'c6', 'c12'])
    ffnb_atomtype['; type'] = topology_atoms['sb_type']
    ffnb_atomtype['chem'] = topology_atoms['atom_type']
    ffnb_atomtype['at.num'] = ffnb_atomtype['chem'].map(gromos_atp['at.num'])
    ffnb_atomtype['mass'] = topology_atoms['mass']
    ffnb_atomtype['charge'] = '0.000000'
    ffnb_atomtype['ptype'] = 'A'
    ffnb_atomtype['c6'] = '0.00000e+00'
    ffnb_atomtype['c12'] = ffnb_atomtype['chem'].map(gromos_atp['c12'])
    
    
    # This will be needed for exclusion and pairs to paste in topology
    # A dictionary with the c12 of each atom in the system
    type_c12_dict = ffnb_atomtype.set_index('; type')['c12'].to_dict()
    
    ffnb_atomtype['c12'] = ffnb_atomtype["c12"].map(lambda x:'{:.6e}'.format(x))
    ffnb_atomtype.drop(columns = ['chem'], inplace = True)

    atomtypes_atp = ffnb_atomtype[['; type', 'mass']].copy()

    return native_atomtypes, ffnb_atomtype, atomtypes_atp, type_c12_dict

def make_more_atomtypes(fibril_pdb):
    '''
    Like the previous one, this part is needed in PDB_LJ_pairs when computing the pairs.
    '''
    fibril_sel = fibril_pdb.select_atoms('all')
    fibril_atomtypes = []
    for atom in fibril_sel:
        atp = str(atom.name) + '_' + str(atom.resnum) + ':' + str(atom.segid)
        fibril_atomtypes.append(atp)

    return fibril_atomtypes

def PDB_LJ_pairs(structure_pdb, atomic_mat_random_coil, atomtypes, parameters):
    '''
    This function measures all the distances between all atoms using MDAnalysis.
    It works on both native and fibril in the same manner.
    Pairs are filtered based on the distance_cutoff which is fixed at 5.5A (in main parameters).
    A second filter is based on the distance_residue. Pairs between atoms closer than two residues are removed 
    if the contact is in the same chain.
    The function also provides the sigma and the epsilon of each pair.
    In case of intramolecular contacts, pairs are reweighted based on the random_coil probability
    '''
    print('\tAddition of PDB derived native LJ-pairs')

    print('\t\tMeasuring distances between all atom in the structure')
    # Selecting all atoms in the system
    atom_sel = structure_pdb.select_atoms('all')

    # Calculating all the distances between atoms.
    # The output is a combination array.
    self_distances = distances.self_distance_array(atom_sel.positions)
    print('\t\tNumber of distances measured :', len(self_distances))
    
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

    # Creation of the dataframe containing the atom pairs and the distances.
    # Also, it will be prepared for sigma and epsilon.
    structural_LJ = pd.DataFrame(columns = ['ai', 'aj', 'distance', 'sigma', 'epsilon'])
    structural_LJ['ai'] = pairs_ai
    structural_LJ['aj'] = pairs_aj
    structural_LJ['distance'] = self_distances
    raw_structural_LJ = len(structural_LJ)
    print('\t\tRaw pairs list ', raw_structural_LJ)
    
    # Keep only the atoms within cutoff
    structural_LJ = structural_LJ[structural_LJ.distance < parameters["distance_cutoff"]] # PROTEIN CONFIGURATION
    print(f'\t\tPairs below cutoff {parameters["distance_cutoff"]}: ', len(structural_LJ))

    # That is name_resname:resid made from the previous function.
    # Extracting the resid information to check if the atom pair is on the same chain.
    structural_LJ[['ai', 'chain_ai']] = structural_LJ.ai.str.split(":", expand = True)
    structural_LJ[['aj', 'chain_aj']] = structural_LJ.aj.str.split(":", expand = True)
   
    structural_LJ['same_chain'] = np.where(structural_LJ['chain_ai'] == structural_LJ['chain_aj'], 'Yes', 'No')
    
    print('\t\tPairs within the same chain: ', len(structural_LJ.loc[structural_LJ['same_chain'] == 'Yes']))
    print('\t\tPairs not in the same chain: ', len(structural_LJ.loc[structural_LJ['same_chain'] == 'No']))

    # if two pairs are made by aminoacids closer than X they'll be deleted. 
    structural_LJ[['type_ai', 'resnum_ai']] = structural_LJ.ai.str.split("_", expand = True)
    structural_LJ[['type_aj', 'resnum_aj']] = structural_LJ.aj.str.split("_", expand = True)
    # And to do that it is necessary to convert the two columns into integer
    structural_LJ = structural_LJ.astype({"resnum_ai": int, "resnum_aj": int})
    structural_LJ['diff'] = ''
    structural_LJ.drop(structural_LJ[(abs(structural_LJ['resnum_aj'] - structural_LJ['resnum_ai']) < parameters['distance_residue']) & (structural_LJ['same_chain'] == 'Yes')].index, inplace = True)    
    structural_LJ['diff'] = abs(structural_LJ['resnum_aj'] - structural_LJ['resnum_ai'])
    print(f'\t\tAll the pairs further than {parameters["distance_residue"]} aminoacids or not in the same chain: ', len(structural_LJ))
    
    # sigma copied below
    structural_LJ['sigma'] = (structural_LJ['distance']/10) / (2**(1/6))
    structural_LJ['epsilon'] = parameters['epsilon_structure']

    # Inverse pairs calvario
    # this must list ALL COLUMNS!
    inv_LJ = structural_LJ[['aj', 'ai', 'distance', 'sigma', 'epsilon', 'chain_ai', 'chain_aj', 'same_chain', 'type_ai', 'resnum_ai', 'type_aj', 'resnum_aj', 'diff']].copy()
    inv_LJ.columns = ['ai', 'aj', 'distance', 'sigma', 'epsilon', 'chain_ai', 'chain_aj', 'same_chain', 'type_ai', 'resnum_ai', 'type_aj', 'resnum_aj', 'diff']
    structural_LJ = pd.concat([structural_LJ, inv_LJ], axis=0, sort = False, ignore_index = True)
    

    # Here we sort all the atom pairs based on the distance and we keep the closer ones.
    # Sorting the pairs
    structural_LJ.sort_values(by = ['ai', 'aj', 'distance'], inplace = True)
    # Cleaning the duplicates
    structural_LJ = structural_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    # Removing the reverse duplicates
    cols = ['ai', 'aj']
    structural_LJ[cols] = np.sort(structural_LJ[cols].values, axis=1)
    structural_LJ = structural_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')

    structural_LJ[['idx_ai', 'idx_aj']] = structural_LJ[['ai', 'aj']]
    structural_LJ.set_index(['idx_ai', 'idx_aj'], inplace=True)

    print(f'\t\tAll the pairs after removing duplicates: ', len(structural_LJ))

    # Take the contact from different chains 
    inter_mask = structural_LJ.same_chain == 'No'
    structural_LJ_inter = structural_LJ[inter_mask]

    # Intramolecular contacts will be reweighted based on the Random coil simulation
    structural_LJ_intra = structural_LJ[~inter_mask]

    atomic_mat_random_coil[['idx_ai', 'idx_aj']] = atomic_mat_random_coil[['rc_ai', 'rc_aj']]
    atomic_mat_random_coil.set_index(['idx_ai', 'idx_aj'], inplace = True)

    # Using inner ho un dataframe vuoto, dunque vuol dire che i contatti tra nativa e fibrilla sono completamente diversi
    # E' un caso generico da prevedere
    structural_LJ_intra = pd.concat([structural_LJ_intra, atomic_mat_random_coil], axis=1, join='inner')

    # OLD prima di usare la formula
    #structural_LJ_intra['epsilon'].loc[structural_LJ_intra['rc_probability'] == 1] = 0
    pd.options.mode.chained_assignment = None
    
    u_threshold = 1-parameters['ratio_threshold']

    # Paissoni Equation 2.0
    # Attractive pairs
    #structural_LJ_intra['epsilon'].loc[u_threshold >=  structural_LJ_intra['rc_probability'])] = epsilon_structure*(1-((np.log(u_threshold))/(np.log(structural_LJ_intra['rc_probability']))))
    #structural_LJ_intra['epsilon'].loc[(u_threshold <  structural_LJ_intra['rc_probability'])] = 0
    # Paissoni Equation 2.1
    # Attractive pairs
    structural_LJ_intra['epsilon'].loc[(u_threshold >=  structural_LJ_intra['rc_probability'])] = -(parameters['epsilon_structure']/np.log(0.1*parameters['ratio_threshold']))*(np.log(u_threshold/structural_LJ_intra['rc_probability']))
    structural_LJ_intra['epsilon'].loc[(u_threshold <  structural_LJ_intra['rc_probability'])] = 0
    # Too small epsilon will be removed
    structural_LJ_intra['epsilon'].loc[abs(structural_LJ_intra['epsilon']) < 0.01*parameters['epsilon_structure']] = 0
    structural_LJ_intra.dropna(inplace=True)

    # This is included in the old before using the formula
    structural_LJ_intra = structural_LJ_intra[structural_LJ_intra.epsilon != 0]
    structural_LJ = pd.concat([structural_LJ_intra,structural_LJ_inter], axis=0, sort = False, ignore_index = True)

    pd.options.mode.chained_assignment = 'warn' 
    print('\t\tSigma and epsilon completed ', len(structural_LJ))
    structural_LJ.drop(columns = ['distance', 'chain_ai', 'chain_aj', 'same_chain', 'type_ai', 'resnum_ai', 'type_aj', 'resnum_aj', 'diff', 'rc_ai',  'rc_aj',  'rc_distance', 'rc_probability', 'rc_residue_ai', 'rc_residue_aj'], inplace = True)

    return structural_LJ


def MD_LJ_pairs(atomic_mat_plainMD, atomic_mat_random_coil, parameters):
    '''
    This function reads the probabilities obtained using mdmat on the plainMD and the random coil simulations.
    For each atom contact the sigma and epsilon are obtained.
    '''
    # In the case of an IDP, it is possible to add dynamical informations based on a simulation
    print('\tAddition of MD derived LJ-pairs')
    # The ratio threshold considers only pairs occurring at a certain probability
    # This dictionary was made to link amber and greta atomtypes
    atomic_mat_plainMD = atomic_mat_plainMD.replace({'ai':gro_to_amb_dict})
    atomic_mat_plainMD = atomic_mat_plainMD.replace({'aj':gro_to_amb_dict})

    atomic_mat_random_coil = atomic_mat_random_coil.replace({'rc_ai':gro_to_amb_dict})
    atomic_mat_random_coil = atomic_mat_random_coil.replace({'rc_aj':gro_to_amb_dict})

    # Add sigma, add epsilon reweighted, add c6 and c12
    atomic_mat_plainMD['sigma'] = (atomic_mat_plainMD['distance']) / (2**(1/6))
    # Merge the two dataframes by ai and aj which are also indexes now
    atomic_mat_plainMD[['idx_ai', 'idx_aj']] = atomic_mat_plainMD[['ai', 'aj']]
    atomic_mat_plainMD.set_index(['idx_ai', 'idx_aj'], inplace = True)

    atomic_mat_random_coil[['idx_ai', 'idx_aj']] = atomic_mat_random_coil[['rc_ai', 'rc_aj']]
    atomic_mat_random_coil.set_index(['idx_ai', 'idx_aj'], inplace = True)

    atomic_mat_merged = pd.concat([atomic_mat_plainMD, atomic_mat_random_coil], axis=1)

    # Epsilon reweight based on probability
    atomic_mat_merged['epsilon'] = ''    

    pd.options.mode.chained_assignment = None

    # Paissoni Equation 2.0
    # Attractive pairs
    #atomic_mat_merged['epsilon'].loc[(atomic_mat_merged['probability'] >=  atomic_mat_merged['rc_probability'])] = epsilon_md*(1-((np.log(atomic_mat_merged['probability']))/(np.log(atomic_mat_merged['rc_probability']))))
    # Repulsive pairs test
    #atomic_mat_merged['epsilon'].loc[(atomic_mat_merged['probability'] <  atomic_mat_merged['rc_probability'])] = -(epsilon_md*(1-((np.log(atomic_mat_merged['rc_probability']))/(np.log(atomic_mat_merged['probability'])))))
    #atomic_mat_merged['sigma'].loc[(atomic_mat_merged['probability'] <  atomic_mat_merged['rc_probability'])] = atomic_mat_merged['rc_distance']/(2**(1/6))

    # Paissoni Equation 2.1
    # Attractive
    atomic_mat_merged['epsilon'].loc[(atomic_mat_merged['probability'] >=  atomic_mat_merged['rc_probability'])] = -(parameters['epsilon_md']/np.log(0.1*parameters['ratio_threshold']))*(np.log(atomic_mat_merged['probability']/atomic_mat_merged['rc_probability']))
    # Repulsive
    atomic_mat_merged['epsilon'].loc[(atomic_mat_merged['probability'] <  atomic_mat_merged['rc_probability'])] = (parameters['epsilon_md']/np.log(parameters['ratio_threshold']))*(np.log(atomic_mat_merged['rc_probability']/atomic_mat_merged['probability']))
    atomic_mat_merged['sigma'].loc[(atomic_mat_merged['probability'] <  atomic_mat_merged['rc_probability'])] = atomic_mat_merged['rc_distance']/(2**(1/6))

    # Treshold vari ed eventuali
    atomic_mat_merged['epsilon'].loc[(atomic_mat_merged['probability'] <  parameters['ratio_threshold'])] = 0
    atomic_mat_merged['epsilon'].loc[abs(atomic_mat_merged['epsilon']) < 0.01*parameters['epsilon_md']] = 0
    pd.options.mode.chained_assignment = 'warn' 


    atomic_mat_merged.drop(columns = ['distance', 'rc_residue_ai', 'rc_residue_aj', 'residue_ai', 'residue_aj', 'probability', 'rc_ai', 'rc_aj', 'rc_probability', 'rc_distance'], inplace = True)
    atomic_mat_merged.dropna(inplace=True)
    atomic_mat_merged = atomic_mat_merged[atomic_mat_merged.epsilon != 0]

    print("\t\t",len(atomic_mat_merged), " interactions added")
    print("\t\t average epsilon is ", atomic_mat_merged['epsilon'].mean())

    return atomic_mat_merged


def merge_and_clean_LJ(greta_LJ, parameters):
    '''
    This function merges the atom contacts from native and fibril and removed eventual duplicates.
    '''
    # Inverse pairs calvario
    inv_LJ = greta_LJ[['aj', 'ai', 'sigma', 'epsilon']].copy()
    inv_LJ.columns = ['ai', 'aj', 'sigma', 'epsilon']
    greta_LJ = pd.concat([greta_LJ,inv_LJ], axis=0, sort = False, ignore_index = True)
    print('\tDoubled pairs list: ', len(greta_LJ))

    # Here we sort all the atom pairs based on the distance and we keep the closer ones.
    print('\tSorting and dropping all the duplicates')
    # Sorting the pairs
    greta_LJ.sort_values(by = ['ai', 'aj', 'sigma'], inplace = True)

    # Cleaning the duplicates
    greta_LJ = greta_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    # Removing the reverse duplicates
    cols = ['ai', 'aj']
    greta_LJ[cols] = np.sort(greta_LJ[cols].values, axis=1)
    greta_LJ = greta_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    print('\tCleaning Complete ', len(greta_LJ))
    
    greta_LJ.insert(2, 'type', 1)
    greta_LJ.insert(3, 'c6', '')
    greta_LJ['c6'] = 4 * greta_LJ['epsilon'] * (greta_LJ['sigma'] ** 6)
    greta_LJ.insert(4, 'c12', '')
    greta_LJ['c12'] = abs(4 * greta_LJ['epsilon'] * (greta_LJ['sigma'] ** 12))

    # SELF INTERACTIONS
    # In the case of fibrils which are not fully modelled we add self interactions which is a feature of amyloids
    # So that the balance between native and fibril is less steep.
    print('- Self interactions')
    atomtypes = set(greta_LJ['ai'])
    greta_LJ['double'] = ''

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
        print('\tAll atoms interacts with themself')
        
    else:
        print('\tThere are', len(atp_notdoubles), 'self interactions to add')
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
                print('\t\tOnly one self interacting pair available for {} ==> {}'.format((str(a)[:-1]), 'Skip'))
            elif len(sigma) == 0:
                # If the missing atom pairs is not represented in the strcture there are not
                # sigmas to average
                print('\t\tThere are not self interactions for {:<12} ==> {}'.format((str(a)[:-1]), 'Skip'))
            else:
                # If there are enough sigmas to make an average then it creates the missing atom pairs
                media_sigma = sigma.mean()
                sd_sigma = sigma.std()
                print('\t\tThere are {:<3} {:<3} with an average Sigma of: {:>17.10f} +/- {}'.format((len(sigma)), (str(a)[:-1]), media_sigma, sd_sigma))
                
                # Creation of new c6 and c12
                # Epsilon structure because those are self
                new_c6 = 4 * parameters['epsilon_structure'] * (media_sigma ** 6)
                new_c12 = 4 * parameters['epsilon_structure'] * (media_sigma ** 12)

                # In the pairs to add dataframe all those new information are inserted

                pairs_toadd.loc[(pairs_toadd['ai'].str.contains(a)) & (pairs_toadd['aj'].str.contains(a)), 'c6'] = new_c6
                pairs_toadd.loc[(pairs_toadd['ai'].str.contains(a)) & (pairs_toadd['aj'].str.contains(a)), 'c12'] = new_c12
                pairs_toadd.loc[(pairs_toadd['ai'].str.contains(a)) & (pairs_toadd['aj'].str.contains(a)), 'sigma'] = media_sigma
                pairs_toadd.loc[(pairs_toadd['ai'].str.contains(a)) & (pairs_toadd['aj'].str.contains(a)), 'epsilon'] = parameters['epsilon_structure']

        pairs_toadd.dropna(inplace = True)
        # Appending the missing atom pairs to the main dataframe
        greta_LJ = pd.concat([greta_LJ,pairs_toadd], axis=0, sort = False, ignore_index = True)
        print('\tSelf interactions added to greta_LJ ', len(pairs_toadd))

    # Drop double, we don't need it anymore
    greta_LJ.drop(columns = ['double'], inplace = True)
    greta_LJ.insert(5, '', ';')
    greta_LJ = greta_LJ.rename(columns = {'ai':'; ai'})
    greta_LJ['epsilon'] = greta_LJ["epsilon"].map(lambda x:'{:.6f}'.format(x))
    greta_LJ['sigma'] = greta_LJ["sigma"].map(lambda x:'{:.6e}'.format(x))
    greta_LJ['c6'] = greta_LJ["c6"].map(lambda x:'{:.6e}'.format(x))
    greta_LJ['c12'] = greta_LJ["c12"].map(lambda x:'{:.6e}'.format(x))

    print('\tLJ Merging completed: ', len(greta_LJ))

    return greta_LJ


def make_pairs_exclusion_topology(type_c12_dict, topology_atoms, topology_bonds, parameters, greta_merge=pd.DataFrame()):
    '''
    This function prepares the [ exclusion ] and [ pairs ] section to paste in topology.top
    Here we define the GROMACS exclusion list and drop from the LJ list made using GRETA so that all the remaining
    contacts will be defined in pairs and exclusions as particular cases.
    Since we are not defining explicit H, the 1-4 list is defined by 2 bonds and not 3 bonds.
    '''
    if not greta_merge.empty:
        greta_merge = greta_merge.rename(columns = {'; ai': 'ai'})
    
    atnum_type_top = topology_atoms[['atom_number', 'sb_type', 'residue_number', 'atom', 'atom_type', 'residue']].copy()
    atnum_type_top['residue_number'] = atnum_type_top['residue_number'].astype(int)
    atnum_type_top['atom_number'] = atnum_type_top['atom_number'].astype(str)

    # Dictionaries definitions to map values
    atnum_type_dict = atnum_type_top.set_index('sb_type')['atom_number'].to_dict()
    type_atnum_dict = atnum_type_top.set_index('atom_number')['sb_type'].to_dict()
    # Bonds from topology
    bond_tuple = topology_bonds.bond_pairs

    #TODO this should be in topology_definitions.py
    # Building the exclusion bonded list
    # exclusion_bonds are all the interactions within 3 bonds
    # p14 are specifically the interactions at exactly 3 bonds
    ex, ex14, p14, exclusion_bonds = [], [], [], []
    for atom in topology_atoms['atom_number'].to_list():
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
    

    if not greta_merge.empty:
        # pairs from greta does not have duplicates because these have been cleaned before
        pairs = greta_merge[['ai', 'aj']].copy()
        pairs['c12_ai'] = pairs['ai']
        pairs['c12_aj'] = pairs['aj']
        pairs[['type_ai', 'resnum_ai']] = pairs.ai.str.split("_", expand = True)
        pairs[['type_aj', 'resnum_aj']] = pairs.aj.str.split("_", expand = True)
        pairs['resnum_ai'] = pairs['resnum_ai'].astype(int)
        pairs['resnum_aj'] = pairs['resnum_aj'].astype(int)

	# When generating LJ interactions we kept intermolecular interactions between atoms belonging to residues closer than distance residues
        # Now we neeed to be sure that these are excluded intramolecularly
        # If we keep such LJ they cause severe frustration to the system and artifacts
        pairs = pairs.loc[(abs(pairs['resnum_aj'] - pairs['resnum_ai']) < parameters['distance_residue'])]

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
        pairs['c12'] = np.sqrt(pairs['c12_ai'] * pairs['c12_aj'])
        pairs.drop(columns = ['type_ai', 'resnum_ai', 'type_aj', 'resnum_aj', 'c12_ai', 'c12_aj', 'check', 'exclude'], inplace = True)    

    else:
        pairs = pd.DataFrame()

    # Only 1-4 exclusions are fully reintroduced
    pairs_14 = pd.DataFrame(columns=['ai', 'aj', 'exclusions'])
    pairs_14['exclusions'] = p14
    pairs_14[['ai', 'aj']] = pairs_14.exclusions.str.split("_", expand = True)

    pairs_14['c12_ai'] = pairs_14['ai']
    pairs_14['c12_aj'] = pairs_14['aj']
    pairs_14['c12_ai'] = pairs_14['c12_ai'].map(type_atnum_dict)
    pairs_14['c12_aj'] = pairs_14['c12_aj'].map(type_atnum_dict)

    # Adding an atom column because we want to flag NOT N N interactions, this because the N should have an explicit H we removed
    pairs_14[['ai_type', 'ai_resid']] = pairs_14.c12_ai.str.split("_", expand = True)
    pairs_14[['aj_type', 'aj_resid']] = pairs_14.c12_aj.str.split("_", expand = True)

    # NOT 1-4 N-X interactions will be dropped
    pairs_14.loc[(pairs_14['ai_type'] == 'N') | (pairs_14['aj_type'] == 'N'), 'c12_tozero'] = False
    # Here we take a particular interest of the interaction between two N, because both should have an explicit H
    pairs_14.loc[(pairs_14['ai_type'] == 'N') & (pairs_14['aj_type'] == 'N'), 'c12_tozero'] = True
    # Only the pairs with an N involved are retained
    pairs_14.drop(pairs_14[pairs_14.c12_tozero != False].index, inplace=True)

    # Thus, only N with X LJ 1-4 interactions will be kept
    # All the other 1-4 interactions will NOT interact with each others
    pairs_14['c12_ai'] = pairs_14['c12_ai'].map(type_c12_dict)
    pairs_14['c12_aj'] = pairs_14['c12_aj'].map(type_c12_dict)
    pairs_14['func'] = 1
    pairs_14['c6'] = 0.00000e+00
    # in general 1-4 interactions are excluded, N-X 1-4 interactions are retained but scaled down
    pairs_14['c12'] = (np.sqrt(pairs_14['c12_ai'] * pairs_14['c12_aj']))*parameters['lj_reduction']
    
    # The N-N interactions are less scaled down, double the c12
    pairs_14.loc[(pairs_14['c12_tozero'] == True), 'c12'] *= 2

    # Removing the interactions with the proline N becasue this does not have the H
    residue_list = topology_atoms['residue'].to_list()
    proline_n = []
    if 'PRO' in residue_list:
        proline_n = topology_atoms.loc[(topology_atoms['residue'] == 'PRO') & (topology_atoms['atom'] == 'N'), 'nr'].to_list()

    pairs_14 = pairs_14[~pairs_14['ai'].isin(proline_n)]
    pairs_14 = pairs_14[~pairs_14['aj'].isin(proline_n)]

    pairs_14.drop(columns = ['exclusions', 'c12_ai', 'c12_aj', 'ai_type', 'ai_resid','aj_type', 'aj_resid', 'c12_tozero'], inplace = True)    

    # Exclusions 1-4
    pairs = pd.concat([pairs,pairs_14], axis=0, sort=False, ignore_index=True)

    # Drop duplicates
    pairs.sort_values(by = ['ai', 'aj', 'c12'], inplace = True)
    # Cleaning the duplicates (in case of doubt keep the smallest c12)
    pairs = pairs.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')

    # Removing the reverse duplicates
    cols = ['ai', 'aj']
    pairs[cols] = np.sort(pairs[cols].values, axis=1)
    pairs = pairs.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')

    # Left alpha helix
    # Here we introduce the c6 and the c12 between the O and the CB of the next residue
    # to optimize the left alpha region in the Ramachandran

    # Adding the c6 and c12 (I do it now because later is sbatti)
    atnum_type_top['c6'] = atnum_type_top['atom_type'].map(gromos_atp['c6'])
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
        line_sidechain_cb = sidechain_cb.loc[sidechain_cb['residue_number'] == (line_backbone_oxygen['residue_number'])+1].squeeze(axis=None)
        if not line_sidechain_cb.empty:
            left_alpha_ai.append(line_backbone_oxygen['atom_number'])
            left_alpha_aj.append(line_sidechain_cb['atom_number'])
            #left_alpha_c6.append(np.sqrt(line_backbone_oxygen['c6']*line_sidechain_cb['c6'])*parameters['multiply_c6'])
            #left_alpha_c12.append(np.sqrt(line_backbone_oxygen['c12']*line_sidechain_cb['c12']))
            #we use the parameters of Alanine because dihedrals have been optimised with these
            left_alpha_c6.append(0.007115485)
            left_alpha_c12.append(0.000005162090)

    left_alpha_pairs = pd.DataFrame(columns=['ai', 'aj', 'c6', 'c12'])
    left_alpha_pairs['ai'] = left_alpha_ai
    left_alpha_pairs['aj'] = left_alpha_aj
    left_alpha_pairs['c6'] = left_alpha_c6
    left_alpha_pairs['c12'] = left_alpha_c12
    left_alpha_pairs['func'] = 1

    pairs = pd.concat([pairs,left_alpha_pairs], axis=0, sort=False, ignore_index=True)
    # Cleaning the duplicates (the left alpha pairs win on pairs that may be previously defined)
    pairs.sort_values(by = ['ai', 'aj', 'c6'], inplace = True)
    pairs = pairs.drop_duplicates(subset = ['ai', 'aj'], keep = 'last')

    pairs['ai'] = pairs['ai'].astype(int)
    pairs['aj'] = pairs['aj'].astype(int)

    # Here we want to sort so that ai is smaller than aj
    inv_pairs = pairs[['aj', 'ai', 'func', 'c6', 'c12']].copy()
    inv_pairs.columns = ['ai', 'aj', 'func', 'c6', 'c12']
    pairs = pd.concat([pairs,inv_pairs], axis=0, sort = False, ignore_index = True)
    pairs = pairs[pairs['ai']<pairs['aj']]
    
    pairs.sort_values(by = ['ai', 'aj'], inplace = True)

    pairs = pairs.rename(columns = {'ai': '; ai'})
    pairs['c6'] = pairs["c6"].map(lambda x:'{:.6e}'.format(x))
    pairs['c12'] = pairs["c12"].map(lambda x:'{:.6e}'.format(x))
    exclusion = pairs[['; ai', 'aj']].copy()

    return pairs, exclusion
