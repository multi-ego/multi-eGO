from operator import concat
import MDAnalysis as mda
from MDAnalysis.lib.util import parse_residue
from MDAnalysis.analysis import distances
import numpy as np
from pandas.core.frame import DataFrame
import pandas as pd
import itertools
import parmed as pmd
from read_input import random_coil_mdmat, plainMD_mdmat
import mdtraj as md

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

# native reweight for TTR and ABeta. This dictionary will rename the amber topology to gromos topology
#gro_to_amb_dict = {'OC1_11' : 'O1_11', 'OC2_11':'O2_11'}
#gro_to_amb_dict = {'OT1_42' : 'O1_42', 'OT2_42':'O2_42'}

class ensemble:
    '''
    Ensemble class: aggregates all the parameters used in the script.
    '''
    #def __init__(self, topology_file, structure_file, parameters, get_structure_pairs=False, is_MD=False, not_matching_native=False, use_RC=False, get_pairs_exclusions=False):
    def __init__(self, parameters, ensemble_parameters):
        
        # RC is needed for native pairs
        if ensemble_parameters['get_structure_pairs'] == True:
            ensemble_parameters['use_RC'] = True
        
        # Topology Section
        # Atoms
        print('\t- Generating ensemble Atomtypes')
        print('\t- Reading topology and structure')
        topology = pmd.load_file(ensemble_parameters['topology_file'], parametrize=False)
        structure = pmd.load_file(ensemble_parameters['structure_file'])
        ensemble_top, renamed_structure = prepare_ensemble_topology(topology, structure, ensemble_parameters, parameters)
        print('\t- Ensemble topology generated')

        if ensemble_parameters['not_matching_native']:
            sbtype_idx_dict = ensemble_parameters['not_matching_native'].sbtype_idx_dict           
            # Fibril might not be modelled, therefore the numbering does not match the native structure
            # Here a dictionary is supplied to renumber the atoms
            print('\t- Renumbering the fibril atom numbers')
            ensemble_top['atom_number'] = ensemble_top['sb_type'].map(sbtype_idx_dict)
        
        self.multiego_topology = ensemble_top

        native_atomtypes = (ensemble_top['sb_type'] +':'+ ensemble_top['chain']).tolist()
        self.native_atomtypes = native_atomtypes

        atomtypes_atp = ensemble_top[['sb_type', 'mass']].copy()
        atomtypes_atp.rename(columns={'sb_type':'; type'}, inplace=True)
        self.atomtypes_atp = atomtypes_atp

        type_c12_dict = ensemble_top[['sb_type', 'c12']].copy()
        type_c12_dict.rename(columns={'sb_type':'; type'}, inplace=True)
        type_c12_dict = type_c12_dict.set_index('; type')['c12'].to_dict()
        self.type_c12_dict = type_c12_dict

        idx_sbtype_dict = ensemble_top[['atom_number', 'sb_type']].copy()
        idx_sbtype_dict = idx_sbtype_dict.set_index('atom_number')['sb_type'].to_dict()
        self.idx_sbtype_dict = idx_sbtype_dict

        sbtype_idx_dict = ensemble_top[['atom_number', 'sb_type']].copy()
        sbtype_idx_dict = sbtype_idx_dict.set_index('sb_type')['atom_number'].to_dict()
        self.sbtype_idx_dict = sbtype_idx_dict

        ffnonbonded_atp = ensemble_top[['sb_type', 'atomic_number', 'mass', 'charge', 'ptype', 'c6', 'c12']].copy()
        ffnb_colnames = ['; type', 'at.num', 'mass', 'charge', 'ptype', 'c6', 'c12']
        ffnonbonded_atp.columns = ffnb_colnames
        ffnonbonded_atp['c12'] = ffnonbonded_atp['c12'].map(lambda x:'{:.6e}'.format(x))
        # TODO there's a problem with the atom number, but the topol.top ones seems correct
        self.ffnonbonded_atp = ffnonbonded_atp

        # ACID pH
        # Selection of the aminoacids and the charged atoms (used for B2m)
        # TODO add some options for precise pH setting
        acid_ASP = ensemble_top[(ensemble_top['residue'] == "ASP") & ((ensemble_top['atom'] == "OD1") | (ensemble_top['atom'] == "OD2") | (ensemble_top['atom'] == "CG"))]
        acid_GLU = ensemble_top[(ensemble_top['residue'] == "GLU") & ((ensemble_top['atom'] == "OE1") | (ensemble_top['atom'] == "OE2") | (ensemble_top['atom'] == "CD"))]
        acid_HIS = ensemble_top[(ensemble_top['residue'] == "HIS") & ((ensemble_top['atom'] == "ND1") | (ensemble_top['atom'] == "CE1") | (ensemble_top['atom'] == "NE2") | (ensemble_top['atom'] == "CD2") | (ensemble_top['atom'] == "CG"))]
        frames = [acid_ASP, acid_GLU, acid_HIS]
        acid_atp = pd.concat(frames, ignore_index = True)
        #this is used
        self.acid_atp = acid_atp['sb_type'].tolist()

        if ensemble_parameters['get_pairs_exclusions'] == True:
            print('\t- Generating bonds, pairs and exclusions')
            ensemble_bonds = get_topology_bonds(topology)
            bonds_pairs = list([(str(ai), str(aj)) for ai, aj in zip(ensemble_bonds['ai'].to_list(), ensemble_bonds['aj'].to_list())])
            self.multiego_bonds = ensemble_bonds
            self.bonds_pairs = bonds_pairs

            # Pairs and Exclusions
            topology_pairs, topology_exclusion = make_pairs_exclusion_topology(ensemble_top, bonds_pairs, type_c12_dict, parameters)
            self.topology_pairs = topology_pairs
            self.topology_exclusion = topology_exclusion

        # The random coil should come from the same native ensemble
        if ensemble_parameters['use_RC']==True:
            if not ensemble_parameters['not_matching_native']:
                atomic_mat_random_coil = random_coil_mdmat(parameters, idx_sbtype_dict)
                self.atomic_mat_random_coil = atomic_mat_random_coil
        
        # MD_contacts
        if ensemble_parameters['is_MD']==True:
            atomic_mat_plainMD = plainMD_mdmat(parameters, idx_sbtype_dict)
            self.atomic_mat_plainMD = atomic_mat_plainMD
        
        # Structure based contacts
        if ensemble_parameters['get_structure_pairs']==True:
            if ensemble_parameters['not_matching_native']:
                native_ensemble = ensemble_parameters['not_matching_native']
                #mda_structure = mda.Universe(ensemble_parameters['structure_file'], guess_bonds = True)
                mda_structure = mda.Universe(renamed_structure, guess_bonds = True, topology_format='PDB')
                structure_pairs = PDB_LJ_pairs(mda_structure, native_ensemble.atomic_mat_random_coil, native_atomtypes, parameters)
            else:
                #mda_structure = mda.Universe(ensemble_parameters['structure_file'], guess_bonds = True)
                mda_structure = mda.Universe(renamed_structure, guess_bonds = True, topology_format='PDB')
                structure_pairs = PDB_LJ_pairs(mda_structure, atomic_mat_random_coil, native_atomtypes, parameters)
            self.structure_pairs = structure_pairs


def prepare_ensemble_topology(topology, structure, ensemble_parameters, parameters):
    print('\t\t- Checking the atoms in both Topology and Structure')
    for atom_top, atom_struct in zip(topology.atoms, structure.atoms):
        if str(atom_top) != str(atom_struct):
            print('- Topology and structure have different atom definitions\n\n')
            print(atom_top, atom_struct)
            exit()
    print('\t\t- Atoms between Topology and Structure are corresponding')

    topology_df = topology.to_dataframe()
    structure_df = structure.to_dataframe()
    hydrogen_number = 0
    if ensemble_parameters['is_MD'] == True:
        print('\t\t- Removing Hydrogens')
        # MD has hydrogen which we don't use
        hydrogen = topology['@H=']
        hydrogen_number = len(hydrogen.atoms)
        hydrogen_namelist = []
        for h in hydrogen.atoms:
            hydrogen_namelist.append(h.name)
        hydrogen_namelist = list(set(hydrogen_namelist))
        mask = (topology_df['name']).isin(hydrogen_namelist)
        topology_df = topology_df[~mask]
        topology_df.reset_index(inplace=True)
        mask = (structure_df['name']).isin(hydrogen_namelist)
        structure_df = structure_df[~mask]
        structure_df['number'] = list(range(1, len(structure_df)+1))
        structure_df.reset_index(inplace=True)

    print('\t\t- Chain renaming')
    chain_labels, chain_ids = [], []
    if ensemble_parameters['is_MD'] == True:
        n_molecules = ['1']
    else:
        n_molecules = topology.molecules
    for chain in n_molecules:
        chain_labels.append(chain[-1]) 
    atoms_per_molecule = int((len(topology.atoms)-hydrogen_number)/len(n_molecules))
    for label in chain_labels:
        for atom in range(atoms_per_molecule):
            chain_ids.append(label)

    print('\t\t- Generating multi-eGO topology')
    multiego_top = pd.DataFrame()
    multiego_top['atom_number'] = structure_df['number']
    multiego_top['atom_type'] = topology_df['type']
    multiego_top['residue_number'] = structure_df['resnum']
    multiego_top['residue'] = structure_df['resname']
    multiego_top['atom'] = topology_df['name']
    multiego_top['cgnr'] = structure_df['resnum']
    multiego_top['mass'] = topology_df['mass']
    multiego_top['atomic_number'] = topology_df['atomic_number']
    multiego_top['chain'] = chain_ids
    multiego_top['charge'] = '0.000000'
    multiego_top['ptype'] = 'A'
    multiego_top['c6'] = '0.00000e+00'
    multiego_top['c12'] = multiego_top['atom_type'].map(gromos_atp['c12'])

    print('\t\t- Topology fixes')
    # Removing an extra H to PRO
    mask = ((multiego_top['residue'] == "PRO") & (multiego_top['atom_type'] == 'N'))
    multiego_top['mass'][mask] = multiego_top['mass'][mask].astype(float).sub(1)
    # Adding an extra H to the N terminal
    mask = ((multiego_top['residue_number'] == multiego_top['residue_number'].min()) & (multiego_top['atom_type'] == 'N'))
    multiego_top['mass'][mask] = multiego_top['mass'][mask].astype(float).add(2)

    # Aromatic carbons dictionary
    aromatic_carbons_dict = {
        'PHE': ['CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TYR': ['CD1', 'CD2', 'CE1', 'CE2'],
        'HIS': ['CE1', 'CD2'],
        'TRP': ['CD1', 'CE3', 'CZ2', 'CZ3', 'CH2']
    }

    for resname, atomnames in aromatic_carbons_dict.items():
        for atom in atomnames:
            mask = ((multiego_top['residue'] == resname) & (multiego_top['atom'] == atom))
            multiego_top['mass'][mask] = multiego_top['mass'][mask].astype(float).add(1)

    print('\t\t- Defining multi-eGO atomtypes')
    multiego_top['sb_type'] = multiego_top['atom'] + '_' + multiego_top['residue_number'].astype(str)

    print('\t\t- Writing a new pdb with renamed chains')
    structure_dict = make_file_dictionary(ensemble_parameters['structure_file'])
    renamed_structure = rename_chains(structure_dict, chain_ids, ensemble_parameters['structure_file'])

    return multiego_top, renamed_structure


def make_file_dictionary(filename):
    file_dict = {}
    with open(filename) as f:
        for line_number, line in enumerate(f):
            if line.startswith('ATOM'):
                file_dict[line_number+1] = line.strip()
    return file_dict


def rename_chains(structure_dict, chain_ids, file_name):
    # Change the column number 22 which is the chain id in the pdb
    column = 22 - 1
    renamed_pdb = []
    file_name = file_name.split('/')
    file_path = file_name[:-1]
    file_name = file_name[-1]
    file_name = file_name.split('.')
    file_name = '/'.join(file_path)+'/'+file_name[0]+'_renamed.'+file_name[1]
    with open(file_name, 'w') as scrivi:
        for (line, string), chain_id in zip(structure_dict.items(), chain_ids):
            new_string = string[0:column]+str(chain_id)+string[column+1:]
            scrivi.write(f"{new_string}\n")
    
    return file_name


def get_topology_bonds(topology):

    bond_atom1, bond_atom2, bond_funct = [], [], []
    for bond in topology.bonds:
        bond_atom1.append(bond.atom1.idx+1)
        bond_atom2.append(bond.atom2.idx+1)
        bond_funct.append(bond.funct)
        # Here is missing the gd_definition, maybe in writing the input will solve this
        #print(dir(bond))
        #print(bond.type)

    bonds_df = pd.DataFrame()
    bonds_df['ai'] = bond_atom1
    bonds_df['aj'] = bond_atom2
    bonds_df['funct'] = bond_funct

    return bonds_df

def sb_type_conversion(multiego_ensemble, md_ensemble):
    '''
    This functions is needed to convert the structure based atomtypes from a force field to gromos.
    It is tested using charmm were the only different atomtypes are OT1 and OT2 which has to be renamed to O1 and O2.
    Other differences are on the atom index which is solved using the structure based atomtype.
    Here a dictionary is made by charmm key: gromos value.
    '''

    multiego_topology = multiego_ensemble.multiego_topology
    md_topology = md_ensemble.multiego_topology

    convert_dict = {}
    multiego_atoms = set(multiego_topology['atom'].to_list())
    md_atoms = set(md_topology['atom'].to_list())
    diff_atoms = list(multiego_atoms - md_atoms)
    merged_atoms = pd.DataFrame()
    merged_atoms['atoms_multiego'] = multiego_topology['atom']
    merged_atoms['multiego_resnum'] = multiego_topology['residue_number']
    merged_atoms['atoms_md'] = md_topology['atom']
    merged_atoms['md_resnum'] = md_topology['residue_number']
    merged_atoms = merged_atoms.loc[merged_atoms['atoms_multiego'].isin(diff_atoms)]
    merged_atoms['sb_multiego'] = merged_atoms['atoms_multiego']+'_'+merged_atoms['multiego_resnum'].astype(str)
    merged_atoms['sb_md'] = merged_atoms['atoms_md']+'_'+merged_atoms['md_resnum'].astype(str)
    merged_atoms_dict = merged_atoms.set_index('sb_md')['sb_multiego'].to_dict()
    
    updated_mat_plainMD = md_ensemble.atomic_mat_plainMD
    updated_mat_plainMD = updated_mat_plainMD.replace({'ai':merged_atoms_dict})
    updated_mat_plainMD = updated_mat_plainMD.replace({'aj':merged_atoms_dict})

    return updated_mat_plainMD, merged_atoms_dict













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

def topology_check(top1, top2):
    if top1 == top2:
        print('- Same topology definitions')
    else:
        difference = set(top1).symmetric_difference(set(top2))
        atom_difference = list(difference)
        #print(atom_difference)
    

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
    #atomic_mat_plainMD = atomic_mat_plainMD.replace({'ai':gro_to_amb_dict})
    #atomic_mat_plainMD = atomic_mat_plainMD.replace({'aj':gro_to_amb_dict})

    #atomic_mat_random_coil = atomic_mat_random_coil.replace({'rc_ai':gro_to_amb_dict})
    #atomic_mat_random_coil = atomic_mat_random_coil.replace({'rc_aj':gro_to_amb_dict})

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
    Also, in case of missing residues in the structure, predicts the self contacts based on the contacts available.
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


def make_pairs_exclusion_topology(ego_topology, bond_tuple, type_c12_dict, parameters, greta_merge=pd.DataFrame()):
    '''
    This function prepares the [ exclusion ] and [ pairs ] section to paste in topology.top
    Here we define the GROMACS exclusion list and drop from the LJ list made using GRETA so that all the remaining
    contacts will be defined in pairs and exclusions as particular cases.
    Since we are not defining explicit H, the 1-4 list is defined by 2 bonds and not 3 bonds.
    This function also fixes the dihedral issue of the left alpha to be explored.
    '''
    if not greta_merge.empty:
        greta_merge = greta_merge.rename(columns = {'; ai': 'ai'})

    ego_topology['atom_number'] = ego_topology['atom_number'].astype(str)
    atnum_type_top = ego_topology[['atom_number', 'sb_type', 'residue_number', 'atom', 'atom_type', 'residue']].copy()
    atnum_type_top['residue_number'] = atnum_type_top['residue_number'].astype(int)

    # Dictionaries definitions to map values
    atnum_type_dict = atnum_type_top.set_index('sb_type')['atom_number'].to_dict()
    type_atnum_dict = atnum_type_top.set_index('atom_number')['sb_type'].to_dict()
    # Bonds from topology

    #TODO this should be in topology_definitions.py
    # Building the exclusion bonded list
    # exclusion_bonds are all the interactions within 3 bonds
    # p14 are specifically the interactions at exactly 3 bonds
    ex, ex14, p14, exclusion_bonds = [], [], [], []
    for atom in ego_topology['atom_number'].to_list():
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
    pairs_14.dropna(inplace=True)#(pairs_14[pairs_14.c12_tozero != False].index, inplace=True)

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
    residue_list = ego_topology['residue'].to_list()
    proline_n = []
    if 'PRO' in residue_list:
        proline_n = ego_topology.loc[(ego_topology['residue'] == 'PRO') & (ego_topology['atom'] == 'N'), 'atom_number'].to_list()

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

    # For each backbone oxygen take the CB of the next residue and save in a pairs tuple
    alpha_beta_rift_ai, alpha_beta_rift_aj, alpha_beta_rift_c6, alpha_beta_rift_c12 = [], [], [], []
    for index, line_backbone_oxygen in backbone_oxygen.iterrows():
        line_sidechain_cb = sidechain_cb.loc[sidechain_cb['residue_number'] == (line_backbone_oxygen['residue_number'])].squeeze(axis=None)
        if not line_sidechain_cb.empty:
            alpha_beta_rift_ai.append(line_backbone_oxygen['atom_number'])
            alpha_beta_rift_aj.append(line_sidechain_cb['atom_number'])
            #alpha_beta_rift_c6.append(np.sqrt(line_backbone_oxygen['c6']*line_sidechain_cb['c6'])*parameters['multiply_c6'])
            #alpha_beta_rift_c12.append(np.sqrt(line_backbone_oxygen['c12']*line_sidechain_cb['c12']))
            #we use the parameters of Alanine because dihedrals have been optimised with these
            alpha_beta_rift_c6.append(0.0)
            alpha_beta_rift_c12.append((0.000005162090)*0.1)

    alpha_beta_rift_pairs = pd.DataFrame(columns=['ai', 'aj', 'c6', 'c12'])
    alpha_beta_rift_pairs['ai'] = alpha_beta_rift_ai
    alpha_beta_rift_pairs['aj'] = alpha_beta_rift_aj
    alpha_beta_rift_pairs['c6'] = alpha_beta_rift_c6
    alpha_beta_rift_pairs['c12'] = alpha_beta_rift_c12
    alpha_beta_rift_pairs['func'] = 1

    pairs = pd.concat([pairs,alpha_beta_rift_pairs], axis=0, sort=False, ignore_index=True)

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
