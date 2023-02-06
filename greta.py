from multiprocessing.dummy import Pool
from turtle import distance
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import pandas as pd
import itertools
import parmed as pmd
from read_input import read_ensemble_mdmat_contacs
from topology_parser import read_topology, topology_parser, extra_topology_ligands
import plotly.express as px
from tqdm import tqdm
from multiprocessing import pool

pd.options.mode.chained_assignment = None  # default='warn'
pd.options.mode.chained_assignment = 'warn' 


gromos_atp = pd.DataFrame(
    {'name': ['O', 'OA', 'N', 'C', 'CH1', 
            'CH2', 'CH3', 'CH2r', 'NT', 'S',
            'NR', 'OM', 'NE', 'NL', 'NZ',
            'CH3p', 'P', 'OE', 'CR1'],            
     'at.num': [8, 8, 7, 6, 6, 6, 6, 6, 7, 16, 7, 8, 7, 7, 7, 6, 15, 8, 6],
     #'c12': [1e-06/3.8, 1.505529e-06/3, 2.319529e-06/2.65, 4.937284e-06/1.9, 9.70225e-05/1.48), # CH1
     #       3.3965584e-05/2.2, 2.6646244e-05/3.1, 2.8058209e-05/2.35, 5.0625e-06/1.95, 1.3075456e-05/4.8,
     #       3.389281e-06/2.25, 7.4149321e-07/4.3, 2.319529e-06/2.65, 2.319529e-06/2.65, 2.319529e-06/2.65],
     #       2.6646244e-05/3.05, 2.2193521e-05/5.7, 1.21e-06/3.4, 1.5116544e-05/2.4],     
     'c12': [2.631580e-07, 5.018430e-07, 8.752940e-07, 2.598570e-06, 6.555574e-05, # CH1
             1.543890e-05, 8.595562e-06, 1.193966e-05, 2.596154e-06, 2.724050e-06, 
             1.506347e-06, 1.724403e-07, 8.752940e-07, 8.752940e-07, 8.752940e-07,
             8.736473e-06, 3.893600e-06, 3.558824e-07, 6.29856e-06]
     }
)
gromos_atp.to_dict()
gromos_atp.set_index('name', inplace=True)

from_ff_to_multiego = {
    'OC1' : 'O1',
    'OC2' : 'O2',
    'OT1' : 'O1',
    'OT2' : 'O2',
    'C13' :'CN1',
    'C14' :'CN2',
    'C15' :'CN3',
    'N'   :'N',
    'C12' :'CA',
    'C11' :'CB',
    'O12' :'OA',
    'P'   :'P',
    'O13' :'OB',
    'O14' :'OC',
    'O11' :'OD',
    'C1'  :'CC',
    'C2'  :'CD',
    'O21' :'OE',
    'C21' :'C1A',
    'O22' :'OF',
    'C22' :'C1B',
    'C23' :'C1C',
    'C24' :'C1D',
    'C25' :'C1E',
    'C26' :'C1F',
    'C27' :'C1G',
    'C28' :'C1H',
    'C29' :'C1I',
    'C210':'C1J',
    'C211':'C1K',
    'C212':'C1L',
    'C213':'C1M',
    'C214':'C1N',
    'C215':'C1O',
    'C216':'C1P',
    'C217':'C1Q',
    'C218':'C1R',
    'C3'  :'CE',
    'O31' :'OG',
    'C31' :'C2A',
    'O32' :'OH',
    'C32' :'C2B',
    'C33' :'C2C',
    'C34' :'C2D',
    'C35' :'C2E',
    'C36' :'C2F',
    'C37' :'C2G',
    'C38' :'C2H',
    'C39' :'C2I',
    'C310':'C2J',
    'C311':'C2K',
    'C312':'C2L',
    'C313':'C2M',
    'C314':'C2N',
    'C315':'C2O',
    'C316':'C2P'
}

class multiego_ensemble:
    #TODO insert the ego_native part since it is redundant
    '''
    This ensemble type gathers different topologies to make a single one.
    '''
    # Checklist of topology sections we need
    multiego_ensemble_top = pd.DataFrame()

    moleculetype = ''
    bonds = pd.DataFrame()
    bond_pairs = pd.DataFrame()
    angles = pd.DataFrame()
    dihedrals = pd.DataFrame()
    impropers = pd.DataFrame()
        
    ligand_bonds = pd.DataFrame()
    ligand_bond_pairs = []
    ligand_angles = pd.DataFrame()
    ligand_dihedrals = pd.DataFrame()

    pairs = pd.DataFrame()
    exclusions = pd.DataFrame()

    structure_based_contacts_dict = {
    }

    greta_LJ = pd.DataFrame(columns=['; ai', 'aj', 'type', 'c6', 'c12', 'sigma', 'epsilon'])
    
    def __init__(self, parameters):
        self.parameters = parameters


    def add_ensemble_top(self, ensemble_toadd):
        '''
        This method allow the addition of atoms into the multi-eGO ensemble
        '''
        # ATOMTYPES
        ensemble_top = ensemble_toadd.ensemble_top.copy()
        ensemble_top['idx_sbtype'] = ensemble_top['sb_type']
        ensemble_top.set_index(['idx_sbtype'], inplace = True)

        multiego_idx = self.multiego_ensemble_top.index
        ensemble_idx = ensemble_top.index
        diff_index = ensemble_idx.difference(multiego_idx)
        if not diff_index.empty:
            #print(f'\t- The following atoms are being inserted in multiego topology: {list(diff_index)}')
            print(f'\t- Inserting atoms in multiego ensemble')
        ensemble_top = ensemble_top.loc[diff_index]
        ensemble_top.sort_values(by=['atom_number'], inplace=True)
        self.multiego_ensemble_top = pd.concat([self.multiego_ensemble_top, ensemble_top], axis=0, sort=False)
    
        # Those will be recreated after the addition of all atoms, including the ligand ones
        type_c12_dict = self.multiego_ensemble_top[['sb_type', 'c12']].copy()
        type_c12_dict.rename(columns={'sb_type':'; type'}, inplace=True)
        type_c12_dict = type_c12_dict.set_index('; type')['c12'].to_dict()
        self.type_c12_dict = type_c12_dict

        type_q_dict = self.multiego_ensemble_top[['sb_type', 'charge']].copy()
        type_q_dict.rename(columns={'sb_type':'; type'}, inplace=True)
        type_q_dict = type_q_dict.set_index('; type')['charge'].to_dict()
        self.type_q_dict = type_q_dict

        return self
    
    
    def add_parsed_topology(self, ensemble_toadd):
        self.moleculetype = ensemble_toadd.parsed_topology.moleculetype.copy()
        self.bonds = ensemble_toadd.bonds.copy()
        self.bond_pairs = ensemble_toadd.bond_pairs.copy()
        self.angles = ensemble_toadd.angles.copy()
        self.dihedrals = ensemble_toadd.dihedrals.copy()
        self.impropers = ensemble_toadd.impropers.copy()
        self.system = ensemble_toadd.parsed_topology.system.copy()
        self.molecules = ensemble_toadd.parsed_topology.molecules.copy()
        
        return self

        
    def add_structure_based_contacts(self, name, contacts):
        '''
        Names to use: random_coil, atomic_mat_plainMD, native_pairs, fibril_pairs, ligand_MD_pairs
        '''
        self.structure_based_contacts_dict[name] = contacts
        # TODO questa funzione sotto si potrebbe anche togliere, mi sa che e' ridondante a quella sopra
        self.structure_based_contacts_dict = self.structure_based_contacts_dict
        return self

        
    def generate_multiego_LJ(self): # TODO nel self bisogna mettere i nomi degli MD ensembles
        
        # TODO ACID DEMMERDA
        #if parameters['acid_ff'] == True and top.acid_atp !=0:
        #        greta_LJ = greta_LJ[~greta_LJ.ai.isin(top.acid_atp)]
        #        greta_LJ = greta_LJ[~greta_LJ.aj.isin(top.acid_atp)]

        greta_LJ = pd.DataFrame()
        greta_MD_LJ = pd.DataFrame()
        ligand_MD_LJ = pd.DataFrame()
        
        # Get the md ensembles names
        for param, value in self.parameters.items():
            if 'md_ensemble' in param:
                greta_MD_LJ = MD_LJ_pairs(self.structure_based_contacts_dict[value], self.structure_based_contacts_dict['random_coil'], self.parameters, self.parameters[param])
                greta_LJ = pd.concat([greta_LJ, greta_MD_LJ], axis=0, sort=False, ignore_index=True)

        #rc_oxyg_LJ = make_oxygens_LJ(self.multiego_ensemble_top, self.type_c12_dict)
        #greta_LJ = pd.concat([greta_LJ, rc_oxyg_LJ], axis=0, sort=False, ignore_index=True)

        if greta_LJ.empty:
            greta_ffnb = greta_LJ 
            greta_lj14 = greta_LJ
        else:
            greta_ffnb, greta_lj14 = merge_and_clean_LJ(self.multiego_ensemble_top, greta_LJ, self.type_c12_dict, self.parameters)

        self.greta_ffnb = greta_ffnb
        self.greta_lj14 = greta_lj14

        return self
        
    def generate_pairs_exclusions(self):
        # Here different bonds_pairs should be added:
        # from native and MD they should be the same, the ligand will be added.
        # Then here the pairs and exclusions will be made.

        bond_pairs = self.bond_pairs# + self.ligand_bond_pairs
        topology_pairs, topology_exclusions = make_pairs_exclusion_topology(self.multiego_ensemble_top, bond_pairs, self.type_c12_dict, self.parameters, 
                                                                            self.greta_lj14) 
        self.pairs = topology_pairs
        self.exclusions = topology_exclusions

        return self

    def add_parsed_ligand_topology(self, ensemble_toadd):
        '''
        This one will be kept separated by the protein parsed topology since the definitions are different
        '''
        self.ligand_moleculetype = ensemble_toadd.ligand_moleculetype.copy()
        self.ligand_topology = ensemble_toadd.ensemble_top.copy()
        self.ligand_bonds = ensemble_toadd.ligand_bonds.copy()
        self.ligand_bond_pairs = ensemble_toadd.ligand_pair_bonds.copy()
        self.ligand_angles = ensemble_toadd.ligand_angles.copy()
        self.ligand_dihedrals = ensemble_toadd.ligand_dihedrals.copy()
        
        # This is used when when want to read ligand pairs from the original topology
        # We might want to remove this part
        self.ligand_pairs = ensemble_toadd.ligand_pairs

        return self

    def list_acid_pH(self):
        # ACID pH
        # Selection of the aminoacids and the charged atoms (used for B2m)
        # TODO add some options for precise pH setting
        acid_ASP = self.ensemble_top[(self.ensemble_top['residue'] == "ASP") & ((self.ensemble_top['atom'] == "OD1") | (self.ensemble_top['atom'] == "OD2") | (self.ensemble_top['atom'] == "CG"))]
        acid_GLU = self.ensemble_top[(self.ensemble_top['residue'] == "GLU") & ((self.ensemble_top['atom'] == "OE1") | (self.ensemble_top['atom'] == "OE2") | (self.ensemble_top['atom'] == "CD"))]
        acid_HIS = self.ensemble_top[(self.ensemble_top['residue'] == "HIS") & ((self.ensemble_top['atom'] == "ND1") | (self.ensemble_top['atom'] == "CE1") | (self.ensemble_top['atom'] == "NE2") | (self.ensemble_top['atom'] == "CD2") | (self.ensemble_top['atom'] == "CG"))]
        frames = [acid_ASP, acid_GLU, acid_HIS]
        acid_atp = pd.concat(frames, ignore_index = True)
        #this is used
        self.acid_atp = acid_atp['sb_type'].tolist()
        return self


    def generate_outputs_toWrite(self):
        # Single and merge are right
        # Topol.top is left
        #pd.set_option('display.colheader_justify', 'left')
        pd.set_option('display.colheader_justify', 'right')

        self.moleculetype_toWrite = self.moleculetype.to_string(index=False)

        ffnonbonded_atp = self.multiego_ensemble_top[['sb_type', 'atomic_number', 'mass', 'charge', 'ptype', 'c6', 'c12']].copy()
        ffnb_colnames = ['; type', 'at.num', 'mass', 'charge', 'ptype', 'c6', 'c12']
        ffnonbonded_atp.columns = ffnb_colnames
        # set charges to zero because they are only used internally
        ffnonbonded_atp['charge'] = 0.0
        ffnonbonded_atp['c12'] = ffnonbonded_atp['c12'].map(lambda x:'{:.6e}'.format(x))
        self.ffnonbonded_atp_toWrite = ffnonbonded_atp.to_string(index = False)
        
        atomtypes_top = self.multiego_ensemble_top[['atom_number', 'sb_type', 'residue_number', 'residue', 'atom', 'cgnr']].copy()
        atomtypes_top.rename(columns = {'atom_number':'; nr', 'sb_type':'type', 'residue_number':'resnr'}, inplace=True)
        self.atomtypes_top_toWrite = atomtypes_top.to_string(index=False)
        
        atomtypes_atp = self.multiego_ensemble_top[['sb_type', 'mass']].copy()
        atomtypes_atp.rename(columns={'sb_type':'; type'}, inplace=True)
        self.atomtypes_atp_toWrite = atomtypes_atp.to_string(index = False, header = False)

        bonds = self.bonds
        bonds.rename(columns = {'ai':'; ai'}, inplace=True)
        self.bonds_toWrite = bonds.to_string(index=False)
        angles = self.angles
        angles.rename(columns = {'ai':'; ai'}, inplace=True)
        self.angles_toWrite = angles.to_string(index=False)
        dihedrals = self.dihedrals
        dihedrals.rename(columns = {'ai':'; ai'}, inplace=True)
        self.dihedrals_toWrite = dihedrals.to_string(index=False)
        impropers = self.impropers
        impropers.rename(columns = {'ai':'; ai'}, inplace=True)
        self.impropers_toWrite = impropers.to_string(index=False)
        pairs = self.pairs
        pairs.insert(5, '', ';')
        pairs.rename(columns = {'ai':'; ai'}, inplace=True)
        if(not pairs.empty):
            self.pairs_toWrite = pairs.to_string(index=False)
        else:
            self.pairs_toWrite = "" 
        exclusions = self.exclusions
        exclusions.rename(columns = {'ai':'; ai'}, inplace=True)
        if(not exclusions.empty):
            self.exclusions_toWrite = exclusions.to_string(index=False)
        else:
            self.exclusions_toWrite = "" 
        #self.system_toWrite = self.system.to_string(index=False)
        self.system_toWrite = self.parameters['protein']
        self.molecules_toWrite = self.molecules.to_string(index=False)

        #if not self.parameters['egos'] == 'rc':
        atomtypes_top['; nr'] = atomtypes_top['; nr'].astype(str)
        atomtypes_top['resnr'] = atomtypes_top['resnr'].astype(int)
        # Dictionaries definitions to map values
        atnum_type_dict = atomtypes_top.set_index('type')['; nr'].to_dict()
        type_atnum_dict = atomtypes_top.set_index('; nr')['type'].to_dict()

        if self.parameters['ligand'] == True:
        
            ligand_atnum_type_df = self.ligand_topology[['sb_type']].copy()
            ligand_atnum_type_df['new_index'] = self.ligand_topology.index+1
            ligand_atnum_type_dict = ligand_atnum_type_df.set_index('sb_type')['new_index'].to_dict()
            atnum_type_dict = {**atnum_type_dict, **ligand_atnum_type_dict}

        if not self.greta_ffnb.empty:
            self.greta_ffnb['ai_n'] = self.greta_ffnb['ai'].map(atnum_type_dict)
            self.greta_ffnb['aj_n'] = self.greta_ffnb['aj'].map(atnum_type_dict)
            self.greta_ffnb['ai_n'] = self.greta_ffnb['ai_n'].astype(int)
            self.greta_ffnb['aj_n'] = self.greta_ffnb['aj_n'].astype(int)
            # Here we want to sort so that ai is smaller than aj
            inv_greta_ffnb = self.greta_ffnb[['aj', 'ai', 'type', 'c6', 'c12', 'sigma', 'epsilon', 'same_chain', 'probability', 'rc_probability', 'source', 'aj_n', 'ai_n']].copy()
            inv_greta_ffnb.columns = ['ai', 'aj', 'type', 'c6', 'c12', 'sigma', 'epsilon', 'same_chain', 'probability', 'rc_probability', 'source', 'ai_n', 'aj_n']
            self.greta_ffnb = pd.concat([self.greta_ffnb,inv_greta_ffnb], axis=0, sort = False, ignore_index = True)
            self.greta_ffnb = self.greta_ffnb[self.greta_ffnb['ai_n']<=self.greta_ffnb['aj_n']]
            self.greta_ffnb.sort_values(by = ['ai_n', 'aj_n'], inplace = True)
            self.greta_ffnb = self.greta_ffnb.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
            
            self.greta_ffnb.insert(5, '', ';')
            self.greta_ffnb = self.greta_ffnb.rename(columns = {'ai':'; ai'})
            self.greta_ffnb['epsilon'] = self.greta_ffnb["epsilon"].map(lambda x:'{:.6f}'.format(x))
            self.greta_ffnb['sigma'] = self.greta_ffnb["sigma"].map(lambda x:'{:.6e}'.format(x))
            self.greta_ffnb['c6'] = self.greta_ffnb["c6"].map(lambda x:'{:.6e}'.format(x))
            self.greta_ffnb['c12'] = self.greta_ffnb["c12"].map(lambda x:'{:.6e}'.format(x))
            self.greta_ffnb = self.greta_ffnb[['; ai', 'aj', 'type', 'c6', 'c12', '', 'sigma', 'epsilon', 'same_chain', 'probability', 'rc_probability', 'source', 'ai_n', 'aj_n']]
            self.greta_ffnb_toWrite = self.greta_ffnb.to_string(index = False)

        if self.parameters['ligand'] == True:
            self.ligand_moleculetype_toWrite = self.ligand_moleculetype.to_string(index=False)
            
            ligand_ffnonbonded_atp = self.ligand_topology[['sb_type', 'atomic_number', 'mass', 'charge', 'ptype', 'c6', 'c12']].copy()
            ffnb_colnames = ['; type', 'at.num', 'mass', 'charge', 'ptype', 'c6', 'c12']
            ligand_ffnonbonded_atp.columns = ffnb_colnames
            ligand_ffnonbonded_atp['c12'] = ligand_ffnonbonded_atp['c12'].map(lambda x:'{:.6e}'.format(x))
            ffnonbonded_atp = pd.concat([ffnonbonded_atp, ligand_ffnonbonded_atp], axis=0, sort=False, ignore_index=True)
            self.ffnonbonded_atp_toWrite = ffnonbonded_atp.to_string(index = False)

            ligand_atomtypes_top = self.ligand_topology[['atom_number', 'sb_type', 'residue_number', 'residue', 'atom', 'cgnr']].copy()
            #ligand_atomtypes_top['atom_number'] = list(range(1, len(ligand_atomtypes_top['atom_number'])+1))
            ligand_atomtypes_top['residue'] = 1
            ligand_atomtypes_top['cgnr'] = 1
            ligand_atomtypes_top.rename(columns = {'atom_number':'; nr', 'sb_type':'type', 'residue_number':'resnr'}, inplace=True)
            self.ligand_atomtypes_top_toWrite = ligand_atomtypes_top.to_string(index=False)
            ligand_bonds = self.ligand_bonds
            ligand_bonds.rename(columns = {'ai':'; ai'}, inplace=True)
            self.ligand_bonds_toWrite = ligand_bonds.to_string(index=False)
            ligand_angles = self.ligand_angles
            ligand_angles.rename(columns = {'ai':'; ai'}, inplace=True)
            self.ligand_angles_toWrite = ligand_angles.to_string(index=False)
            ligand_dihedrals = self.ligand_dihedrals
            ligand_dihedrals.rename(columns = {'ai':'; ai'}, inplace=True)
            self.ligand_dihedrals_toWrite = ligand_dihedrals.to_string(index=False)
            ligand_pairs = self.ligand_pairs
            ligand_pairs.rename(columns = {'ai':'; ai'}, inplace=True)
            ligand_pairs['c6'] = ligand_pairs['c6'].map(lambda x:'{:.6e}'.format(x))
            ligand_pairs['c12'] = ligand_pairs['c12'].map(lambda x:'{:.6e}'.format(x))
            self.ligand_pairs_toWrite = ligand_pairs.to_string(index=False)
            ligand_exclusions = self.ligand_pairs[['; ai', 'aj']].copy()
            self.ligand_exclusions_toWrite = ligand_exclusions.to_string(index=False)

        return self


class ensemble:
    '''
    Ensemble class: aggregates all the parameters used in the script.
    '''


    def __init__(self, parameters, ensemble_parameters, name):        
        # Topology Section
        # Atoms
        print('\t- Generating ensemble Atomtypes')
        print('\t- Reading topology and structure')
        self.name = name
        self.parameters = parameters
        self.ensemble_parameters = ensemble_parameters
        self.topology = pmd.load_file(ensemble_parameters[f'{name}_topology'], parametrize=False)
        #self.structure = pmd.load_file(ensemble_parameters[f'{name}_structure'])

    def prepare_ensemble(self, add_native_ensemble = False):
        '''
        This function creates a dataframe from the topol.top file read using ParmEd.
        The ensemble dataframe contains all the topology details we need in one place and it will be used to
        create dictionaries. Here we add some details such as our custom c12 values.
        Additionally, we also increase the mass of certain atoms to include the H we removed.
        Finally we create the following dictionaries:
        ...
        ...
        ...
        
        '''

        print('\t\t- Generating multi-eGO topology')

        # Removing solvent from the dataframe
        atom_selection = self.topology["!((:TIP3)|(:SOL)|(:WAT)|(:NA))"]
        topology_df = atom_selection.to_dataframe()

        topology_df['number'] = list(range(1, len(atom_selection.atoms)+1))
        topology_df['resid'] = topology_df['resid'] + 1
        topology_df['resnum'] = topology_df['resid']
        topology_df['cgnr'] = topology_df['resid']
        topology_df['ptype'] = 'A'
        topology_df['c6'] = '0.00000e+00'
        topology_df['c12'] = topology_df['type'].map(gromos_atp['c12'])

        # TODO metti qui il dizionario degli atomtype basati sui diversi force-fields
        topology_df.rename(columns={
            'number':'atom_number',
            'type':'atom_type',
            'resnum':'residue_number',
            'resname':'residue',
            'name':'atom',
        }, inplace=True)

        print('\t\t- Applying topology fixes')
        # Removing an extra H to PRO
        pd.options.mode.chained_assignment = None 
        mask = ((topology_df['residue'] == "PRO") & (topology_df['atom_type'] == 'N'))
        topology_df['mass'][mask] = topology_df['mass'][mask].astype(float).sub(1)
        # Adding an extra H to the N terminal
        mask = ((topology_df['residue_number'] == topology_df['residue_number'].min()) & (topology_df['atom_type'] == 'N'))
        topology_df['mass'][mask] = topology_df['mass'][mask].astype(float).add(2)

        # Aromatic carbons dictionary
        aromatic_carbons_dict = {
            'PHE': ['CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TYR': ['CD1', 'CD2', 'CE1', 'CE2'],
            'HIS': ['CE1', 'CD2'],
            'TRP': ['CD1', 'CE3', 'CZ2', 'CZ3', 'CH2']
        }

        for resname, atomnames in aromatic_carbons_dict.items():
            for atom in atomnames:
                mask = ((topology_df['residue'] == resname) & (topology_df['atom'] == atom))
                topology_df['mass'][mask] = topology_df['mass'][mask].astype(float).add(1)

        # Converting the different atomtypes from different forcefields
        topology_df = topology_df.replace({'atom':from_ff_to_multiego})

        print('\t\t- Defining multi-eGO atomtypes')
        topology_df['sb_type'] = topology_df['atom'] + '_' + topology_df['residue_number'].astype(str)        
                
        # Definition of the different dictionaries
        self.ensemble_top = topology_df

        self.atoms_size = len(topology_df)

        sbtype_idx_dict = topology_df[['atom_number', 'sb_type']].copy()
        sbtype_idx_dict = sbtype_idx_dict.set_index('sb_type')['atom_number'].to_dict()
        self.sbtype_idx_dict = sbtype_idx_dict

        type_c12_dict = topology_df[['sb_type', 'c12']].copy()
        type_c12_dict.rename(columns={'sb_type':'; type'}, inplace=True)
        type_c12_dict = type_c12_dict.set_index('; type')['c12'].to_dict()
        self.type_c12_dict = type_c12_dict
        
        type_q_dict = topology_df[['sb_type', 'charge']].copy()
        type_q_dict.rename(columns={'sb_type':'; type'}, inplace=True)
        type_q_dict = type_q_dict.set_index('; type')['charge'].to_dict()
        self.type_q_dict = type_q_dict

        idx_sbtype_dict = topology_df[['atom_number', 'sb_type']].copy()
        idx_sbtype_dict = idx_sbtype_dict.set_index('atom_number')['sb_type'].to_dict()
        self.idx_sbtype_dict = idx_sbtype_dict
    
        print('\t- Ensemble topology generated')


    def get_parsed_topology(self):
        '''
        Topol.top sort of things except atoms. Namely bonds, angles, dihedrals, impropers, pairs and exclusions.
        This method uses the parser i wrote and not ParmEd.
        '''
        parsed_topology = topology_parser(read_topology(self.ensemble_parameters[f'{self.name}_topology']))

        # This one is self so i wrote less things in multiego and here
        self.parsed_topology = parsed_topology
        
        self.bonds = parsed_topology.bonds
        bond_pairs = parsed_topology.bond_pairs
        self.bond_pairs = bond_pairs
        self.angles = parsed_topology.angles
        self.dihedrals = parsed_topology.dihedrals
        self.impropers = parsed_topology.impropers


    def convert_topology(self, ego_native):
        '''
        This functions is needed to convert the structure based atomtypes from a force field to gromos.
        It is tested using charmm were the only different atomtypes are OT1 and OT2 which has to be renamed to O1 and O2.
        Other differences are on the atom index which is solved using the structure based atomtype.
        Here a dictionary is made by charmm key: gromos value.
        '''
        multiego_topology = ego_native.ensemble_top
        md_topology = self.ensemble_top

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

        md_topology = md_topology.replace({'sb_type':merged_atoms_dict})
        self.ensemble_top = md_topology
        updated_mat_plainMD = self.atomic_mat_MD
        updated_mat_plainMD = updated_mat_plainMD.replace({'ai':merged_atoms_dict})
        updated_mat_plainMD = updated_mat_plainMD.replace({'aj':merged_atoms_dict})
        self.atomic_mat_MD = updated_mat_plainMD
        self.conversion_dict = merged_atoms_dict


    # TODO queste due sotto possono essere unite se ci piace avere anche un inter rc
    def add_random_coil(self):
        # The random coil should come from the same native ensemble
        #if not ensemble_parameters['not_matching_native']:
        atomic_mat_random_coil = read_ensemble_mdmat_contacs(self.ensemble_parameters[f'{self.name}_contacts'], self.idx_sbtype_dict)
        self.atomic_mat_random_coil = atomic_mat_random_coil
        
        
    def add_MD_contacts(self):
        # MD_contacts
        #atomic_mat_MD = plainMD_mdmat(self.parameters, self.ensemble_parameters[f'{self.name}_contacts'], self.idx_sbtype_dict, self.idx_chain_dict)
        atomic_mat_MD = read_ensemble_mdmat_contacs(self.ensemble_parameters[f'{self.name}_contacts'], self.idx_sbtype_dict)
        self.atomic_mat_MD = atomic_mat_MD


    def get_ligand_ensemble(self): # TODO change name
        '''
        Reading ligand .itp and .prm to get topology parameters to add in topol_ligand and ffnonbonded.itp.
        Parameters are ligand C6 and C12, bonds, angles, dihedrals and pairs.
        '''
        # ATOMS
        # Here is just filtering by the ligand 
        ligand_residue_number = self.ensemble_top['residue_number'].max()
        self.ligand_residue_number = ligand_residue_number

        # extra atomtypes for c12 definitions
        print('\t - Retrieving ligand extra atom definitions')
        itp = read_topology(self.ensemble_parameters['itp_file'])
        prm = read_topology(self.ensemble_parameters['prm_file'])
        extra_ligand_top = extra_topology_ligands(itp, prm, ligand_residue_number)

        self.ligand_moleculetype = extra_ligand_top.ligand_moleculetype

        # Inserting the new c12 in ffnonbonded.itp
        ligand_ensemble_top = self.ensemble_top.loc[self.ensemble_top['residue_number'] == ligand_residue_number]
        ligand_ensemble_top['c12'] = ligand_ensemble_top['sb_type'].map(extra_ligand_top.ligand_sbtype_c12_dict)
        ligand_ensemble_top['atom_number'] = ligand_ensemble_top['sb_type'].map(extra_ligand_top.sbtype_ligand_number_dict)
        #ligand_ensemble_top['c12'] = ligand_ensemble_top['c12'].map(lambda x:'{:.6e}'.format(x))
        self.ensemble_top = ligand_ensemble_top
        
        ligand_sbtype_number_dict = ligand_ensemble_top[['sb_type', 'atom_number']].copy()
        ligand_sbtype_number_dict = ligand_sbtype_number_dict.set_index('sb_type')['atom_number'].to_dict()   
        self.ligand_sbtype_number_dict = ligand_sbtype_number_dict

        update_ligand_bonds = extra_ligand_top.ligand_bonds
        update_ligand_bonds['ai'] = update_ligand_bonds['ai'].map(ligand_sbtype_number_dict)
        update_ligand_bonds['aj'] = update_ligand_bonds['aj'].map(ligand_sbtype_number_dict)
        update_ligand_bonds.dropna(inplace=True)
        update_ligand_bonds['ai'] = update_ligand_bonds['ai'].astype(int)
        update_ligand_bonds['aj'] = update_ligand_bonds['aj'].astype(int)
        # This is for gromacs which is bitchy
        update_ligand_bonds['c1'] = update_ligand_bonds['c1']/self.parameters['ligand_reduction']
        self.ligand_bonds = update_ligand_bonds
        bond_pairs = list([(str(ai), str(aj)) for ai, aj in zip(update_ligand_bonds['ai'].to_list(), update_ligand_bonds['aj'].to_list())])
        self.ligand_pair_bonds = bond_pairs

        update_ligand_angles = extra_ligand_top.ligand_angles
        update_ligand_angles['ai'] = update_ligand_angles['ai'].map(ligand_sbtype_number_dict)
        update_ligand_angles['aj'] = update_ligand_angles['aj'].map(ligand_sbtype_number_dict)
        update_ligand_angles['ak'] = update_ligand_angles['ak'].map(ligand_sbtype_number_dict)
        update_ligand_angles.dropna(inplace=True)
        update_ligand_angles['ai'] = update_ligand_angles['ai'].astype(int)
        update_ligand_angles['aj'] = update_ligand_angles['aj'].astype(int)
        update_ligand_angles['ak'] = update_ligand_angles['ak'].astype(int)
        self.ligand_angles = update_ligand_angles

        update_ligand_dihedrals = extra_ligand_top.ligand_dihedrals
        update_ligand_dihedrals['ai'] = update_ligand_dihedrals['ai'].map(ligand_sbtype_number_dict)
        update_ligand_dihedrals['aj'] = update_ligand_dihedrals['aj'].map(ligand_sbtype_number_dict)
        update_ligand_dihedrals['ak'] = update_ligand_dihedrals['ak'].map(ligand_sbtype_number_dict)
        update_ligand_dihedrals['al'] = update_ligand_dihedrals['al'].map(ligand_sbtype_number_dict)
        update_ligand_dihedrals.dropna(inplace=True)
        update_ligand_dihedrals['ai'] = update_ligand_dihedrals['ai'].astype(int)
        update_ligand_dihedrals['aj'] = update_ligand_dihedrals['aj'].astype(int)
        update_ligand_dihedrals['ak'] = update_ligand_dihedrals['ak'].astype(int)
        update_ligand_dihedrals['al'] = update_ligand_dihedrals['al'].astype(int)
        self.ligand_dihedrals = update_ligand_dihedrals

        # This is used when when want to read ligand pairs from the original topology
        # We might want to remove this part
        update_ligand_pairs = extra_ligand_top.ligand_pairs
        update_ligand_pairs['ai'] = update_ligand_pairs['ai'].map(ligand_sbtype_number_dict)
        update_ligand_pairs['aj'] = update_ligand_pairs['aj'].map(ligand_sbtype_number_dict)
        update_ligand_pairs.dropna(inplace=True)
        update_ligand_pairs['ai'] = update_ligand_pairs['ai'].astype(int)
        update_ligand_pairs['aj'] = update_ligand_pairs['aj'].astype(int)
        self.ligand_pairs = update_ligand_pairs
     

    def ligand_MD_LJ_pairs(self):
        # TODO this one will be moved in multiego_ensemble
        print('\t- Adding the ligand LJ potential')

        # Following the idea that the ligand is considered as a residue in the same chain,
        # we select ligand pairs by selecting the spare residue compared with the ligand
        atomic_ligand_mat = self.atomic_mat_MD
        # Drop values below md_treshold, otherwise the epsilon will be negative
        atomic_ligand_mat.drop(atomic_ligand_mat[atomic_ligand_mat['probability'] < self.parameters['md_threshold']].index, inplace=True)
        #atomic_ligand_mat.drop(atomic_ligand_mat[atomic_ligand_mat['probability'] <= 0].index, inplace=True)
        # Filtering for the ligand
        #atomic_ligand_mat.dropna(inplace=True)
        atomic_ligand_mat = atomic_ligand_mat.loc[(atomic_ligand_mat['residue_ai'] == self.ligand_residue_number) | (atomic_ligand_mat['residue_aj'] == self.ligand_residue_number)]
        self_interactions = atomic_ligand_mat.loc[(atomic_ligand_mat['residue_ai'] == self.ligand_residue_number) & (atomic_ligand_mat['residue_aj'] == self.ligand_residue_number)]
        
        # TODO to be fixed. Here I am removing the self contacts of the ligand which with the latest updates causes the molecules to overlap.
        # Namely, with the new update more self ligand contacts are retained with a very high probability.
        # A quick and dirty fix is to remove such self contacts, but those should be added in pairs and exclusions properly.
        mask = ((atomic_ligand_mat['residue_ai'] == self.ligand_residue_number) & (atomic_ligand_mat['residue_aj'] == self.ligand_residue_number))
        #print(atomic_ligand_mat[mask].to_string())
        atomic_ligand_mat = atomic_ligand_mat[~mask]

        atomic_ligand_mat['sigma'] = (atomic_ligand_mat['distance']) / (2**(1/6))
        atomic_ligand_mat[['idx_ai', 'idx_aj']] = atomic_ligand_mat[['ai', 'aj']]
        atomic_ligand_mat.set_index(['idx_ai', 'idx_aj'], inplace = True)
        
        # We are using the old equation (prior to the RC version)
        atomic_ligand_mat['epsilon'] = self.parameters['epsilon_ligand']*(1-((np.log(atomic_ligand_mat['probability']))/(np.log(self.parameters['md_threshold']))))
        atomic_ligand_mat.drop(columns = ['distance', 'residue_ai', 'residue_aj', 'probability'], inplace = True)
        atomic_ligand_mat.dropna(inplace=True)
        atomic_ligand_mat = atomic_ligand_mat[atomic_ligand_mat.epsilon != 0]
        self.ligand_atomic_mat_MD = atomic_ligand_mat

# END CLASS


def MD_LJ_pairs(atomic_mat_plainMD, atomic_mat_random_coil, parameters, name):
    '''
    This function reads the probabilities obtained using mdmat on the plainMD and the random coil simulations.
    For each atom contact the sigma and epsilon are obtained.
    '''
    print('\t- Addition of MD derived LJ-pairs')
    atomic_mat_plainMD[['idx_ai', 'idx_aj']] = atomic_mat_plainMD[['ai', 'aj']]
    atomic_mat_plainMD.set_index(['idx_ai', 'idx_aj'], inplace = True)

    atomic_mat_random_coil[['idx_ai', 'idx_aj']] = atomic_mat_random_coil[['rc_ai', 'rc_aj']]
    atomic_mat_random_coil.set_index(['idx_ai', 'idx_aj'], inplace = True)

    intra_mat_probability = pd.DataFrame()
    inter_mat_probability = pd.DataFrame()

    # Retrieving the intramolecular contacts
    intra_atomic_mat_plainMD = atomic_mat_plainMD.loc[atomic_mat_plainMD['same_chain'] == 'Yes']
    intra_atomic_mat_rc = atomic_mat_random_coil.loc[atomic_mat_random_coil['rc_same_chain'] == 'Yes']
    inter_atomic_mat_plainMD = atomic_mat_plainMD.loc[atomic_mat_plainMD['same_chain'] == 'No']
    inter_atomic_mat_rc = atomic_mat_random_coil.loc[atomic_mat_random_coil['rc_same_chain'] == 'No']

    # Le sto tenendo separate perché non sto mettendo regole specifiche intra-inter all'interno della funzione
    # Adesso per comodità sto usando intramolecular per entrambi
    # TODO controllare con Carlo dove mettere bene la separazione di same_chain all'interno di reweight intramolecular 
    # E rinominarla di conseguenza


    if not intra_atomic_mat_plainMD.empty:
        intra_mat_probability = reweight_contacts(intra_atomic_mat_plainMD, intra_atomic_mat_rc, parameters, name)
        print(f'\t\t- Intramolecular contact pairs added after removing duplicates: ', len(intra_mat_probability))
        print("\t\t\t- Epsilon input:", parameters['epsilon_md']) 
        print("\t\t\t- Average epsilon is", intra_mat_probability['epsilon'].loc[(intra_mat_probability['epsilon']>0.)].mean())
        print("\t\t\t- Maximum epsilon is", intra_mat_probability['epsilon'].max())

    if not inter_atomic_mat_plainMD.empty:
        inter_mat_probability = reweight_contacts(inter_atomic_mat_plainMD, inter_atomic_mat_rc, parameters, name)
        print(f'\t\t- Intermolecular contact pairs added after removing duplicates: ', len(inter_mat_probability))
        print("\t\t\t- Epsilon input:", parameters['epsilon_md']) 
        print("\t\t\t- Average epsilon is", inter_mat_probability['epsilon'].loc[(inter_mat_probability['epsilon']>0.)].mean())
        print("\t\t\t- Maximum epsilon is", inter_mat_probability['epsilon'].max())

    atomic_mat_merged = pd.concat([intra_mat_probability, inter_mat_probability], axis=0, sort = False, ignore_index = True)

    return atomic_mat_merged


def reweight_contacts(atomic_mat_plainMD, atomic_mat_random_coil, parameters, name):
    
    rew_mat = atomic_mat_plainMD.copy()
    rew_mat.to_csv(f'analysis/rew_mat_{name}')

    rew_mat = pd.concat([rew_mat, atomic_mat_random_coil], axis=1)
    rew_mat.drop(columns = ['rc_ai', 'rc_aj'], inplace=True)
    #rew_mat = rew_mat.loc[(rew_mat['probability']>parameters['md_threshold'])]
    #rew_mat = rew_mat.loc[(rew_mat['probability']>parameters['md_threshold'])|((rew_mat['probability']<parameters['md_threshold'])&(rew_mat['rc_probability']>parameters['md_threshold']))]
    rew_mat = rew_mat.loc[(rew_mat['probability']>parameters['md_threshold'])|((rew_mat['probability']<parameters['md_threshold'])&(rew_mat['rc_probability']>parameters['md_threshold']))|((rew_mat['probability']<parameters['rc_threshold'])&(rew_mat['rc_probability']>parameters['rc_threshold']))]

    # Add sigma, add epsilon reweighted, add c6 and c12
    rew_mat['sigma'] = (rew_mat['distance']) / (2**(1/6))
    rew_mat['epsilon'] = np.nan 

    # Epsilon reweight based on probability
    # Paissoni Equation 2.1
    # Attractive intramolecular
    rew_mat['epsilon'].loc[(rew_mat['probability']>1.1*np.maximum(rew_mat['rc_probability'],parameters['rc_threshold']))&(rew_mat['same_chain']=='Yes')] = -(parameters['epsilon_md']/np.log(parameters['rc_threshold']))*(np.log(rew_mat['probability']/np.maximum(rew_mat['rc_probability'],parameters['rc_threshold'])))
    # Attractive intermolecular
    rew_mat['epsilon'].loc[(rew_mat['probability']>1.1*np.maximum(rew_mat['rc_probability'],parameters['rc_threshold']))&(rew_mat['same_chain']=='No')] = -(parameters['epsilon_amyl']/np.log(parameters['rc_threshold']))*(np.log(rew_mat['probability']/np.maximum(rew_mat['rc_probability'],parameters['rc_threshold'])))
    # Repulsive
    # case 1: rc > md > md_t
    rew_mat['epsilon'].loc[(rew_mat['probability']>parameters['md_threshold'])&(rew_mat['probability']<0.9*rew_mat['rc_probability'])] = np.log(rew_mat['probability']/rew_mat['rc_probability'])*(np.minimum(rew_mat['distance'],rew_mat['rc_distance'])**12)
    # case 2: rc > md_t > md
    rew_mat['epsilon'].loc[((rew_mat['probability']<parameters['md_threshold'])&(rew_mat['rc_probability']>parameters['rc_threshold']))&(np.maximum(rew_mat['probability'],parameters['md_threshold'])<0.9*rew_mat['rc_probability'])] = np.log(np.maximum(rew_mat['probability'],parameters['md_threshold'])/rew_mat['rc_probability'])*(np.minimum(rew_mat['distance'],rew_mat['rc_distance'])**12)
    # case 3: rc > rc_t > md
    rew_mat['epsilon'].loc[((rew_mat['probability']<parameters['rc_threshold'])&(rew_mat['rc_probability']>parameters['rc_threshold']))&(np.maximum(rew_mat['probability'],parameters['rc_threshold'])<0.9*rew_mat['rc_probability'])] = np.log(np.maximum(rew_mat['probability'],parameters['rc_threshold'])/rew_mat['rc_probability'])*(np.minimum(rew_mat['distance'],rew_mat['rc_distance'])**12)

    # clean NaN and zeros 
    rew_mat.dropna(inplace=True)
    rew_mat = rew_mat[rew_mat.epsilon != 0]

    # Retrieving the ai and aj information from the index and reindexing back again
    # TODO maybe there's a better way to do the same thing
    rew_mat = rew_mat.reset_index()
    rew_mat[['ai', 'aj']] = rew_mat[['idx_ai', 'idx_aj']]
    rew_mat.set_index(['idx_ai', 'idx_aj'], inplace = True)

    # Changing the columns order
    rew_mat = rew_mat[['ai', 'aj', 'sigma', 'epsilon', 'probability', 'rc_probability', 'same_chain']]
    
    # Inverse pairs calvario
    # this must list ALL COLUMNS!
    inv_mat = rew_mat[['aj', 'ai', 'sigma', 'epsilon', 'probability', 'rc_probability', 'same_chain']].copy()
    inv_mat.columns = ['ai', 'aj', 'sigma', 'epsilon', 'probability', 'rc_probability', 'same_chain']
    rew_mat = pd.concat([rew_mat, inv_mat], axis=0, sort = False, ignore_index = True)
    # Here we sort all the atom pairs based on the distance and we keep the closer ones.
    # Sorting the pairs
    rew_mat.sort_values(by = ['ai', 'aj', 'sigma'], inplace = True)
    # Cleaning the duplicates
    rew_mat = rew_mat.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    # Removing the reverse duplicates
    cols = ['ai', 'aj']
    rew_mat[cols] = np.sort(rew_mat[cols].values, axis=1)
    rew_mat = rew_mat.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    rew_mat[['idx_ai', 'idx_aj']] = rew_mat[['ai', 'aj']]
    rew_mat.set_index(['idx_ai', 'idx_aj'], inplace=True)
    rew_mat['source'] = name

    return rew_mat

def make_oxygens_LJ(ego_topology, type_c12_dict):
    ego_topology['atom_number'] = ego_topology['atom_number'].astype(str)
    atnum_type_top = ego_topology[['atom_number', 'sb_type', 'residue_number', 'atom', 'atom_type', 'residue']].copy()
    atnum_type_top['residue_number'] = atnum_type_top['residue_number'].astype(int)

    # Dictionaries definitions to map values
    atnum_type_dict = atnum_type_top.set_index('sb_type')['atom_number'].to_dict()
    type_atnum_dict = atnum_type_top.set_index('atom_number')['sb_type'].to_dict()

    # Adding the c12 and charge
    atnum_type_top['c12'] = atnum_type_top['sb_type'].map(type_c12_dict)

    # Here we add a special c12 only for oxygen-oxygen interactions
    neg_oxygen = atnum_type_top.loc[(atnum_type_top['atom'].astype(str).str[0]=='O')]
    neg_oxygen['key'] = 1
    rc_oxyg_LJ = pd.merge(neg_oxygen,neg_oxygen,on='key').drop('key',axis=1)
    rc_oxyg_LJ = rc_oxyg_LJ[['sb_type_x', 'sb_type_y', 'c12_x', 'c12_y']]
    rc_oxyg_LJ['sigma'] = 0.55 
    rc_oxyg_LJ['epsilon'] = -11.4*np.sqrt(rc_oxyg_LJ['c12_x']*rc_oxyg_LJ['c12_y']) 
    rc_oxyg_LJ.drop(columns=['c12_y', 'c12_x'], inplace=True)
    rc_oxyg_LJ.columns=['ai','aj','sigma','epsilon']
    rc_oxyg_LJ['probability'] = np.nan
    rc_oxyg_LJ['rc_probability'] = np.nan
    rc_oxyg_LJ['same_chain'] = 'Yes'
    rc_oxyg_LJ['source'] = 'rc'
    # drop inverse duplicates
    cols = ['ai', 'aj']
    rc_oxyg_LJ[cols] = np.sort(rc_oxyg_LJ[cols].values, axis=1)
    rc_oxyg_LJ = rc_oxyg_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
 
    return rc_oxyg_LJ 

def merge_and_clean_LJ(ego_topology, greta_LJ, type_c12_dict, parameters):
    '''
    This function merges the atom contacts from native and fibril and removed eventual duplicates.
    Also, in case of missing residues in the structure, predicts the self contacts based on the contacts available.
    '''
    ego_topology['atom_number'] = ego_topology['atom_number'].astype(str)
    atnum_type_top = ego_topology[['atom_number', 'sb_type', 'residue_number', 'atom', 'atom_type', 'residue']].copy()
    atnum_type_top['residue_number'] = atnum_type_top['residue_number'].astype(int)

    # Dictionaries definitions to map values
    atnum_type_dict = atnum_type_top.set_index('sb_type')['atom_number'].to_dict()
    type_atnum_dict = atnum_type_top.set_index('atom_number')['sb_type'].to_dict()

    # Adding the c12 and charge
    atnum_type_top['c12'] = atnum_type_top['sb_type'].map(type_c12_dict)

    print('- Merging Inter and Intra molecular interactions')
    print('\t- Merged pairs list: ', len(greta_LJ))
    print('\t- Sorting and dropping all the duplicates')
    # Inverse pairs calvario
    inv_LJ = greta_LJ[['aj', 'ai', 'sigma', 'epsilon', 'same_chain', 'probability', 'rc_probability', 'source']].copy()
    inv_LJ.columns = ['ai', 'aj', 'sigma', 'epsilon', 'same_chain', 'probability', 'rc_probability', 'source']

    greta_LJ = pd.concat([greta_LJ, inv_LJ], axis=0, sort = False, ignore_index = True)

    if parameters['egos'] == 'split':
        # in this case we use intra and inter molecular contacts from specific simulations
        # yet we check the compatibility of the distances
        # we evaluate the minimum sigma for each contact
        greta_LJ['new_sigma'] = greta_LJ.groupby(by=['ai', 'aj', 'same_chain'])['sigma'].transform('min')
        greta_LJ['energy_at_new_sigma'] = 4.*greta_LJ['epsilon']*((greta_LJ['sigma']/(greta_LJ['new_sigma']*(2.**(1./6.))))**12-(greta_LJ['sigma']/(greta_LJ['new_sigma']*(2.**(1./6.))))**6)
        greta_LJ['energy_at_new_sigma'].loc[(greta_LJ['epsilon']<0.)] = -greta_LJ['epsilon']/(greta_LJ['new_sigma']*(2.**(1./6.)))**12
        # not use  interaction if at new_sigma the repulsion would be too strong 
        greta_LJ = greta_LJ.loc[~((greta_LJ['energy_at_new_sigma']>2.49)&(greta_LJ['source']!='rc'))]
        greta_LJ.drop('new_sigma', axis=1, inplace=True)
        greta_LJ.drop('energy_at_new_sigma', axis=1, inplace=True)
        # split inter and intra depending from the source
        greta_LJ = greta_LJ.loc[(((greta_LJ['same_chain']=='Yes')&((greta_LJ['source']==parameters['intra'])|(greta_LJ['source']=='rc')))|((greta_LJ['same_chain']=='No')&(greta_LJ['source']==parameters['inter'])))]


    pairs_LJ = greta_LJ.copy()

    # Greta prioritise intermolecular interactions and shorter length ones
    greta_LJ.sort_values(by = ['ai', 'aj', 'same_chain', 'sigma'], ascending = [True, True, True, True], inplace = True)
    greta_LJ = greta_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    # Removing the reverse duplicates
    cols = ['ai', 'aj']
    greta_LJ[cols] = np.sort(greta_LJ[cols].values, axis=1)
    greta_LJ = greta_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')

    # Pairs prioritise intramolecular interactions
    pairs_LJ.sort_values(by = ['ai', 'aj', 'same_chain', 'sigma'], ascending = [True, True, False, True], inplace = True)
    pairs_LJ = pairs_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    pairs_LJ[cols] = np.sort(pairs_LJ[cols].values, axis=1)
    pairs_LJ = pairs_LJ.drop_duplicates(subset = ['ai', 'aj'], keep = 'first')
    # where pairs_LJ is the same of greta_LJ and same_chain is yes that the line can be dropped
    # that is I want to keep lines with same_chain no or lines with same chain yes that have same_chain no in greta_LJ
    test = pd.merge(pairs_LJ, greta_LJ, how="right", on=["ai", "aj"])
    pairs_LJ = test.loc[(test['same_chain_x']=='No')|((test['same_chain_x']=='Yes')&(test['same_chain_y']=='No'))]
    pairs_LJ.drop(columns = ['sigma_y', 'epsilon_y', 'same_chain_y', 'probability_y', 'rc_probability_y', 'source_y'], inplace = True)
    pairs_LJ.rename(columns = {'sigma_x': 'sigma', 'probability_x': 'probability', 'rc_probability_x': 'rc_probability', 'epsilon_x': 'epsilon', 'same_chain_x': 'same_chain', 'source_x': 'source'}, inplace = True)

    greta_LJ.insert(2, 'type', 1)
    greta_LJ.insert(3, 'c6', '')
    greta_LJ['c6'] = 4 * greta_LJ['epsilon'] * (greta_LJ['sigma'] ** 6)
    greta_LJ.insert(4, 'c12', '')
    greta_LJ['c12'] = abs(4 * greta_LJ['epsilon'] * (greta_LJ['sigma'] ** 12))
    # repulsive interactions have just a very large C12
    greta_LJ['c6'].loc[(greta_LJ['epsilon']<0.)] = 0.
    greta_LJ['c12'].loc[(greta_LJ['epsilon']<0.)] = np.maximum(-greta_LJ['epsilon'],np.sqrt(greta_LJ['ai'].map(type_c12_dict)*greta_LJ['aj'].map(type_c12_dict)))

    pairs_LJ.insert(2, 'type', 1)
    pairs_LJ.insert(3, 'c6', '')
    pairs_LJ['c6'] = 4 * pairs_LJ['epsilon'] * (pairs_LJ['sigma'] ** 6)
    pairs_LJ.insert(4, 'c12', '')
    pairs_LJ['c12'] = abs(4 * pairs_LJ['epsilon'] * (pairs_LJ['sigma'] ** 12))
    # repulsive interactions have just a very large C12
    pairs_LJ['c6'].loc[(pairs_LJ['epsilon']<0.)] = 0.
    pairs_LJ['c12'].loc[(pairs_LJ['epsilon']<0.)] = np.maximum(-pairs_LJ['epsilon'],np.sqrt(pairs_LJ['ai'].map(type_c12_dict)*pairs_LJ['aj'].map(type_c12_dict)))  

    print('\t- Cleaning and Merging Complete, pairs count:', len(greta_LJ))
    print('\t- LJ Merging completed: ', len(greta_LJ))
    print("\t\t- Average epsilon is ", greta_LJ['epsilon'].mean())
    print("\t\t- Maximum epsilon is ", greta_LJ['epsilon'].max())

    return greta_LJ, pairs_LJ


def make_pairs_exclusion_topology(ego_topology, bond_tuple, type_c12_dict, parameters, greta_merge):
    '''
    This function prepares the [ exclusion ] and [ pairs ] section to paste in topology.top
    '''

    ego_topology['atom_number'] = ego_topology['atom_number'].astype(str)
    atnum_type_top = ego_topology[['atom_number', 'sb_type', 'residue_number', 'atom', 'atom_type', 'residue']].copy()
    atnum_type_top['residue_number'] = atnum_type_top['residue_number'].astype(int)
    # Dictionaries definitions to map values
    atnum_type_dict = atnum_type_top.set_index('sb_type')['atom_number'].to_dict()
    type_atnum_dict = atnum_type_top.set_index('atom_number')['sb_type'].to_dict()

    # Adding the c12
    atnum_type_top['c12'] = atnum_type_top['sb_type'].map(type_c12_dict)

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
    print(exclusion_bonds)
    exit()
    if not greta_merge.empty:
        greta_merge = greta_merge.rename(columns = {'; ai': 'ai'})
        # pairs from greta does not have duplicates because these have been cleaned before
        pairs = greta_merge[['ai', 'aj', 'c6', 'c12', 'epsilon', 'same_chain', 'probability', 'rc_probability', 'source']].copy()
        pairs['c12_ai'] = pairs['ai']
        pairs['c12_aj'] = pairs['aj']
        
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
        # Intermolecular interactions are excluded 
        pairs['c6'].loc[(pairs['same_chain']=='No')] = 0.
        pairs['c12'].loc[(pairs['same_chain']=='No')] = np.sqrt(pairs['c12_ai'] * pairs['c12_aj'])  
        # this is a safety check 
        pairs = pairs[pairs['c12']>0.]
 
        pairs.drop(columns = ['same_chain', 'c12_ai', 'c12_aj', 'check', 'exclude', 'epsilon'], inplace = True)
        pairs = pairs[['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source']]
    else:
        pairs = pd.DataFrame()
    
    # Drop NaNs. This is an issue when adding the ligand ensemble.
    # pairs.dropna(inplace=True)

    # Here we make a dictionary of the atoms used for local geometry 
    backbone_nitrogen = atnum_type_top.loc[atnum_type_top['atom'] == 'N']
    backbone_carbonyl = atnum_type_top.loc[atnum_type_top['atom'] == 'C']
    backbone_oxygen = atnum_type_top.loc[atnum_type_top['atom']=='O']
    ct_oxygen = atnum_type_top.loc[(atnum_type_top['atom']=='O1')|(atnum_type_top['atom']=='O2')]
    sidechain_cb = atnum_type_top.loc[atnum_type_top['atom'] == 'CB']
    pro_cd = atnum_type_top.loc[(atnum_type_top['atom'] == 'CD')&(atnum_type_top['residue'] == 'PRO')]
    sidechain_cgs = atnum_type_top.loc[(atnum_type_top['atom'] == 'CG')|(atnum_type_top['atom'] == 'CG1')|(atnum_type_top['atom'] == 'CG2')|(atnum_type_top['atom'] == 'SG')|(atnum_type_top['atom'] == 'OG')|(atnum_type_top['atom'] == 'OG1')&(atnum_type_top['residue'] != 'PRO')]

    # For proline CD take the CB of the previous residue and save in a pairs tuple
    # CB-1-CD is related to the extended region of the ramachandran
    pairs_14_ai, pairs_14_aj, pairs_14_c6, pairs_14_c12 = [], [], [], []
    for index, line_pro_cd in pro_cd.iterrows():
        line_sidechain_cb = sidechain_cb.loc[(sidechain_cb['residue_number'] == line_pro_cd['residue_number']-1)].squeeze(axis=None)
        if not line_sidechain_cb.empty:
            pairs_14_ai.append(line_pro_cd['atom_number'])
            pairs_14_aj.append(line_sidechain_cb['atom_number'])
            pairs_14_c6.append(0.0)
            pairs_14_c12.append(2.715402e-06)

    pairs_14 = pd.DataFrame(columns=['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source'])
    pairs_14['ai'] = pairs_14_ai
    pairs_14['aj'] = pairs_14_aj
    pairs_14['func'] = 1
    pairs_14['c6'] = pairs_14_c6
    pairs_14['c12'] = pairs_14_c12
    pairs_14['source'] = '1-4' 
    pairs = pd.concat([pairs,pairs_14], axis=0, sort=False, ignore_index=True)

    # For backbone carbonyl take the CB of the next residue and save in a pairs tuple
    pairs_14_ai, pairs_14_aj, pairs_14_c6, pairs_14_c12 = [], [], [], []
    for index, line_backbone_carbonyl in backbone_carbonyl.iterrows():
        line_sidechain_cb = sidechain_cb.loc[(sidechain_cb['residue_number'] == (line_backbone_carbonyl['residue_number']+1))].squeeze(axis=None)
        if not line_sidechain_cb.empty:
            pairs_14_ai.append(line_backbone_carbonyl['atom_number'])
            pairs_14_aj.append(line_sidechain_cb['atom_number'])
            pairs_14_c6.append(0.0)
            pairs_14_c12.append(0.275*np.sqrt(line_sidechain_cb['c12']*line_backbone_carbonyl['c12']))

    pairs_14 = pd.DataFrame(columns=['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source'])
    pairs_14['ai'] = pairs_14_ai
    pairs_14['aj'] = pairs_14_aj
    pairs_14['func'] = 1
    pairs_14['c6'] = pairs_14_c6
    pairs_14['c12'] = pairs_14_c12
    pairs_14['source'] = '1-4' 
    pairs = pd.concat([pairs,pairs_14], axis=0, sort=False, ignore_index=True)

    # For backbone oxygen take the CB of the same residue and save in a pairs tuple
    pairs_14_ai, pairs_14_aj, pairs_14_c6, pairs_14_c12 = [], [], [], []
    for index, line_backbone_oxygen in backbone_oxygen.iterrows():
        line_sidechain_cb = sidechain_cb.loc[(sidechain_cb['residue_number'] == line_backbone_oxygen['residue_number'])].squeeze(axis=None)
        if not line_sidechain_cb.empty:
            pairs_14_ai.append(line_backbone_oxygen['atom_number'])
            pairs_14_aj.append(line_sidechain_cb['atom_number'])
            pairs_14_c6.append(0.0)
            pairs_14_c12.append(0.1*np.sqrt(line_sidechain_cb['c12']*line_backbone_oxygen['c12']))

    pairs_14 = pd.DataFrame(columns=['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source'])
    pairs_14['ai'] = pairs_14_ai
    pairs_14['aj'] = pairs_14_aj
    pairs_14['func'] = 1
    pairs_14['c6'] = pairs_14_c6
    pairs_14['c12'] = pairs_14_c12
    pairs_14['source'] = '1-4' 
    pairs = pd.concat([pairs,pairs_14], axis=0, sort=False, ignore_index=True)

    # now we add the pair between the last CB and the two OCT ones
    pairs_14_ai, pairs_14_aj, pairs_14_c6, pairs_14_c12 = [], [], [], []
    for index, line_ct_oxygen in ct_oxygen.iterrows():
        line_last_cb = sidechain_cb.loc[sidechain_cb['residue_number'] == line_ct_oxygen['residue_number']].squeeze(axis=None)
        if not line_last_cb.empty:
            pairs_14_ai.append(line_ct_oxygen['atom_number'])
            pairs_14_aj.append(line_last_cb['atom_number'])
            pairs_14_c6.append(0.0)
            pairs_14_c12.append(0.1*np.sqrt(line_last_cb['c12']*line_ct_oxygen['c12']))

    pairs_14 = pd.DataFrame(columns=['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source'])
    pairs_14['ai'] = pairs_14_ai
    pairs_14['aj'] = pairs_14_aj
    pairs_14['func'] = 1
    pairs_14['c6'] = pairs_14_c6
    pairs_14['c12'] = pairs_14_c12
    pairs_14['source'] = '1-4' 
    pairs = pd.concat([pairs,pairs_14], axis=0, sort=False, ignore_index=True)

    # For each backbone nitrogen take the CB of the previuos residue and save in a pairs tuple
    pairs_14_ai, pairs_14_aj, pairs_14_c6, pairs_14_c12 = [], [], [], []
    for index, line_backbone_nitrogen in backbone_nitrogen.iterrows():
        line_sidechain_cb = sidechain_cb.loc[sidechain_cb['residue_number'] == (line_backbone_nitrogen['residue_number']-1)].squeeze(axis=None)
        if not line_sidechain_cb.empty:
            pairs_14_ai.append(line_backbone_nitrogen['atom_number'])
            pairs_14_aj.append(line_sidechain_cb['atom_number'])
            pairs_14_c6.append(0.0)
            pairs_14_c12.append(0.65*np.sqrt(line_sidechain_cb['c12']*line_backbone_nitrogen['c12']))

    pairs_14 = pd.DataFrame(columns=['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source'])
    pairs_14['ai'] = pairs_14_ai
    pairs_14['aj'] = pairs_14_aj
    pairs_14['func'] = 1
    pairs_14['c6'] = pairs_14_c6
    pairs_14['c12'] = pairs_14_c12
    pairs_14['source'] = '1-4' 
    pairs = pd.concat([pairs,pairs_14], axis=0, sort=False, ignore_index=True)

    # For each backbone nitrogen take the N of the next residue and save in a pairs tuple
    pairs_14_ai, pairs_14_aj, pairs_14_c6, pairs_14_c12 = [], [], [], []
    for index, line_backbone_nitrogen in backbone_nitrogen.iterrows():
        line_next_n = backbone_nitrogen.loc[backbone_nitrogen['residue_number'] == (line_backbone_nitrogen['residue_number']+1)].squeeze(axis=None)
        if not line_next_n.empty:
            pairs_14_ai.append(line_backbone_nitrogen['atom_number'])
            pairs_14_aj.append(line_next_n['atom_number'])
            pairs_14_c6.append(0.0)
            pairs_14_c12.append(0.343*np.sqrt(line_next_n['c12']*line_backbone_nitrogen['c12']))

    pairs_14 = pd.DataFrame(columns=['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source'])
    pairs_14['ai'] = pairs_14_ai
    pairs_14['aj'] = pairs_14_aj
    pairs_14['func'] = 1
    pairs_14['c6'] = pairs_14_c6
    pairs_14['c12'] = pairs_14_c12
    pairs_14['source'] = '1-4' 
    pairs = pd.concat([pairs,pairs_14], axis=0, sort=False, ignore_index=True)

    # For each backbone oxygen take the O of the next residue and save in a pairs tuple
    pairs_14_ai, pairs_14_aj, pairs_14_c6, pairs_14_c12 = [], [], [], []
    for index, line_backbone_oxygen in backbone_oxygen.iterrows():
        line_next_o = backbone_oxygen.loc[backbone_oxygen['residue_number'] == (line_backbone_oxygen['residue_number']+1)].squeeze(axis=None)
        if not line_next_o.empty:
            pairs_14_ai.append(line_backbone_oxygen['atom_number'])
            pairs_14_aj.append(line_next_o['atom_number'])
            pairs_14_c6.append(0.0)
            pairs_14_c12.append(11.4*np.sqrt(line_backbone_oxygen['c12']*line_next_o['c12']))

    pairs_14 = pd.DataFrame(columns=['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source'])
    pairs_14['ai'] = pairs_14_ai
    pairs_14['aj'] = pairs_14_aj
    pairs_14['func'] = 1
    pairs_14['c6'] = pairs_14_c6
    pairs_14['c12'] = pairs_14_c12
    pairs_14['source'] = '1-4' 
    pairs = pd.concat([pairs,pairs_14], axis=0, sort=False, ignore_index=True)

    # now we add the pair between the penultimate oxygen and the two CT ones
    pairs_14_ai, pairs_14_aj, pairs_14_c6, pairs_14_c12 = [], [], [], []
    for index, line_ct_oxygen in ct_oxygen.iterrows():
        line_prev_o = backbone_oxygen.loc[backbone_oxygen['residue_number'] == (line_ct_oxygen['residue_number']-1)].squeeze(axis=None)
        if not line_prev_o.empty:
            pairs_14_ai.append(line_ct_oxygen['atom_number'])
            pairs_14_aj.append(line_prev_o['atom_number'])
            pairs_14_c6.append(0.0)
            pairs_14_c12.append(11.4*np.sqrt(line_ct_oxygen['c12']*line_prev_o['c12']))

    pairs_14 = pd.DataFrame(columns=['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source'])
    pairs_14['ai'] = pairs_14_ai
    pairs_14['aj'] = pairs_14_aj
    pairs_14['func'] = 1
    pairs_14['c6'] = pairs_14_c6
    pairs_14['c12'] = pairs_14_c12
    pairs_14['source'] = '1-4' 
    pairs = pd.concat([pairs,pairs_14], axis=0, sort=False, ignore_index=True)

    # For each backbone carbonyl take the carbonyl of the next residue and save in a pairs tuple
    pairs_14_ai, pairs_14_aj, pairs_14_c6, pairs_14_c12 = [], [], [], []
    for index, line_backbone_carbonyl in backbone_carbonyl.iterrows():
        line_prev_c = backbone_carbonyl.loc[backbone_carbonyl['residue_number'] == (line_backbone_carbonyl['residue_number']-1)].squeeze(axis=None)
        if not line_prev_c.empty:
            pairs_14_ai.append(line_backbone_carbonyl['atom_number'])
            pairs_14_aj.append(line_prev_c['atom_number'])
            pairs_14_c6.append(0.0)
            pairs_14_c12.append(0.5*np.sqrt(line_backbone_carbonyl['c12']*line_prev_c['c12']))

    pairs_14 = pd.DataFrame(columns=['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source'])
    pairs_14['ai'] = pairs_14_ai
    pairs_14['aj'] = pairs_14_aj
    pairs_14['func'] = 1
    pairs_14['c6'] = pairs_14_c6
    pairs_14['c12'] = pairs_14_c12
    pairs_14['source'] = '1-4' 
    pairs = pd.concat([pairs,pairs_14], axis=0, sort=False, ignore_index=True)

    # For each backbone carbonyl take the CGs of the same residue and save in a pairs tuple
    pairs_14_ai, pairs_14_aj, pairs_14_c6, pairs_14_c12 = [], [], [], []
    for index, line_sidechain_cgs in sidechain_cgs.iterrows():
        line_c = backbone_carbonyl.loc[backbone_carbonyl['residue_number'] == (line_sidechain_cgs['residue_number'])].squeeze(axis=None)
        if not line_c.empty:
            pairs_14_ai.append(line_sidechain_cgs['atom_number'])
            pairs_14_aj.append(line_c['atom_number'])
            pairs_14_c6.append(0.0)
            pairs_14_c12.append(0.078*np.sqrt(line_c['c12']*line_sidechain_cgs['c12']))

    pairs_14 = pd.DataFrame(columns=['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source'])
    pairs_14['ai'] = pairs_14_ai
    pairs_14['aj'] = pairs_14_aj
    pairs_14['func'] = 1
    pairs_14['c6'] = pairs_14_c6
    pairs_14['c12'] = pairs_14_c12
    pairs_14['source'] = '1-4' 
    pairs = pd.concat([pairs,pairs_14], axis=0, sort=False, ignore_index=True)

    # For each backbone nitrogen take the CGs of the same residue and save in a pairs tuple
    pairs_14_ai, pairs_14_aj, pairs_14_c6, pairs_14_c12 = [], [], [], []
    for index, line_sidechain_cgs in sidechain_cgs.iterrows():
        line_c = backbone_nitrogen.loc[backbone_nitrogen['residue_number'] == (line_sidechain_cgs['residue_number'])].squeeze(axis=None)
        if not line_c.empty:
            pairs_14_ai.append(line_sidechain_cgs['atom_number'])
            pairs_14_aj.append(line_c['atom_number'])
            pairs_14_c6.append(0.0)
            pairs_14_c12.append(0.087*np.sqrt(line_c['c12']*line_sidechain_cgs['c12']))

    pairs_14 = pd.DataFrame(columns=['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source'])
    pairs_14['ai'] = pairs_14_ai
    pairs_14['aj'] = pairs_14_aj
    pairs_14['func'] = 1
    pairs_14['c6'] = pairs_14_c6
    pairs_14['c12'] = pairs_14_c12
    pairs_14['source'] = '1-4' 
    pairs = pd.concat([pairs,pairs_14], axis=0, sort=False, ignore_index=True)

    # remove duplicates
    inv_LJ = pairs[['aj', 'ai', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source']].copy()
    inv_LJ.columns = ['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source']
    # in case of duplicates we keep the last occurrence this because in the input pairs there are not duplicates,
    # duplicates can results due to the addition of local geometry pairs that then should superseed the former
    pairs = pd.concat([pairs, inv_LJ], axis=0, sort = False, ignore_index = True)
    pairs = pairs.drop_duplicates(subset = ['ai', 'aj'], keep = 'last')
    # drop inverse duplicates
    cols = ['ai', 'aj']
    pairs[cols] = np.sort(pairs[cols].values, axis=1)
    pairs = pairs.drop_duplicates(subset = ['ai', 'aj'], keep = 'last')
    pairs['ai'] = pairs['ai'].astype(int)
    pairs['aj'] = pairs['aj'].astype(int)

    # Here we want to sort so that ai is smaller than aj
    inv_pairs = pairs[['aj', 'ai', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source']].copy()
    inv_pairs.columns = ['ai', 'aj', 'func', 'c6', 'c12', 'probability', 'rc_probability', 'source']
    pairs = pd.concat([pairs,inv_pairs], axis=0, sort = False, ignore_index = True)
    pairs = pairs[pairs['ai']<pairs['aj']]
    
    pairs.sort_values(by = ['ai', 'aj'], inplace = True)

    pairs = pairs.rename(columns = {'ai': '; ai'})
    pairs['c6'] = pairs["c6"].map(lambda x:'{:.6e}'.format(x))
    pairs['c12'] = pairs["c12"].map(lambda x:'{:.6e}'.format(x))
    exclusion = pairs[['; ai', 'aj']].copy()


    return pairs, exclusion


def convert_topology(from_ff_to_multiego, *ensembles):
    '''
    This function is required to check the different atomtypes between different force fields.
    The atom types MUST match otherwise a proper ffnobonded cannot be created.
    '''

    dict_set = set(from_ff_to_multiego.keys())


    diff_atoms_set = set([])
    for topology in ensembles:
        topology = topology[~topology['atom_type'].astype(str).str.startswith('H')]
        atoms_set = set(topology['atom'].to_list())
        diff_atoms_set = atoms_set - diff_atoms_set
        
    check_dictionary = diff_atoms_set - dict_set

    if check_dictionary:
        print(f'The following atomtypes are not converted. \nYou MUST add them in "from_ff_to_multiego" dictionary to properly merge all the contacts \n{check_dictionary}')
        exit()
