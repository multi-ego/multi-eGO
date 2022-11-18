import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'

def read_topology(topology_path):
    sections_dict = get_sections_dict(make_file_dictionary(topology_path))
    return sections_dict

def make_file_dictionary(topology_path):
    file_dict = {}
    with open(topology_path) as f:
        for line_number, line in enumerate(f):
            if not (line.startswith(';') or line.startswith('#')):
                if not line.strip():
                    continue
                else:
                    file_dict[line_number+1] = line.strip()
    return file_dict

def get_sections_dict(file_dict):
    section_numbers = []
    sections_dict = {}
    for line_number, line in file_dict.items():
        if '[' in line:
            section_numbers.append(line_number)

    end_line = (max(file_dict.keys())+1)
    section_numbers.append(end_line)
    counter = list(range(0, len(section_numbers)))

    for number in counter:
        if number == counter[-1]:
            break
        else:
            lines = list(range(section_numbers[number], section_numbers[number+1]))
            if (file_dict[lines[0]] =='[ dihedrals ]' and '[ dihedrals ]' in list(sections_dict.keys())):
                section_name = '[ impropers ]'
            else:
                section_name = file_dict[lines[0]]
            subset_file_dict = {key: value for key, value in file_dict.items() if key in lines[1:]}
            for k, v in subset_file_dict.items():
                subset_file_dict[k] = v.split()
            sections_dict[section_name] = subset_file_dict
            subset_file_dict = {}
    
    return sections_dict

class topology_parser:
    def __init__(self, sections_dict):
        self.sections_dict = sections_dict

        colnames = ['; Name', 'nrexcl']
        section_dict = sections_dict['[ moleculetype ]']
        moleculetype_df = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        self.moleculetype = moleculetype_df
        #self.string_moleculetype = moleculetype_df.to_string(index=False)

        pd.options.mode.chained_assignment = None  # default='warn'
        colnames = ['atom_number', 'atom_type', 'residue_number', 'residue', 'atom', 'cgnr', 'charge', 'mass', 'typeB', 'chargeB', 'massB']

        section_dict = self.sections_dict['[ atoms ]']
        df_topology_atoms = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        df_topology_atoms.drop(columns=['charge', 'typeB', 'chargeB', 'massB'], inplace=True)
        df_topology_atoms.reset_index(inplace=True, drop=True)   
        
        # Changing the mass of the atoms section by adding the H
        df_topology_atoms['atom_number'] = df_topology_atoms['atom_number'].astype(int)
        df_topology_atoms['residue_number'] = df_topology_atoms['residue_number'].astype(int)
        df_topology_atoms['cgnr'] = df_topology_atoms['cgnr'].astype(int)
        df_topology_atoms['mass'] = df_topology_atoms['mass'].astype(float)

        # Removing an extra H to PRO 
        mask = ((df_topology_atoms['residue'] == "PRO") & (df_topology_atoms['atom_type'] == 'N'))
        df_topology_atoms['mass'][mask] = df_topology_atoms['mass'][mask].astype(float).sub(1)
        # Adding an extra H to the N terminal
        mask = ((df_topology_atoms['residue_number'] == df_topology_atoms['residue_number'].min()) & (df_topology_atoms['atom_type'] == 'N'))
        df_topology_atoms['mass'][mask] = df_topology_atoms['mass'][mask].astype(float).add(2)

        # Aromatic carbons dictionary
        aromatic_carbons_dict = {
            'PHE': ['CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TYR': ['CD1', 'CD2', 'CE1', 'CE2'],
            'HIS': ['CE1', 'CD2'],
            'TRP': ['CD1', 'CE3', 'CZ2', 'CZ3', 'CH2']
        }

        for resname, atomnames in aromatic_carbons_dict.items():
            for atom in atomnames:
               mask = ((df_topology_atoms['residue'] == resname) & (df_topology_atoms['atom'] == atom))
               df_topology_atoms['mass'][mask] = df_topology_atoms['mass'][mask].astype(float).add(1)

        ## Structure based atomtype definition
        df_topology_atoms['sb_type'] = df_topology_atoms['atom'] + '_' + df_topology_atoms['residue_number'].astype(str)
        self.df_topology_atoms = df_topology_atoms
        # Considering that we use only the native topology, only one chain per file is considered
        self.topology_atoms_check = (df_topology_atoms['sb_type']+':1').to_list()
        # ACID pH
        # Selection of the aminoacids and the charged atoms (used for B2m)
        # TODO add some options for precise pH setting
        acid_ASP = df_topology_atoms[(df_topology_atoms['residue'] == "ASP") & ((df_topology_atoms['atom'] == "OD1") | (df_topology_atoms['atom'] == "OD2") | (df_topology_atoms['atom'] == "CG"))]
        acid_GLU = df_topology_atoms[(df_topology_atoms['residue'] == "GLU") & ((df_topology_atoms['atom'] == "OE1") | (df_topology_atoms['atom'] == "OE2") | (df_topology_atoms['atom'] == "CD"))]
        acid_HIS = df_topology_atoms[(df_topology_atoms['residue'] == "HIS") & ((df_topology_atoms['atom'] == "ND1") | (df_topology_atoms['atom'] == "CE1") | (df_topology_atoms['atom'] == "NE2") | (df_topology_atoms['atom'] == "CD2") | (df_topology_atoms['atom'] == "CG"))]
        frames = [acid_ASP, acid_GLU, acid_HIS]
        acid_atp = pd.concat(frames, ignore_index = True)
        #this is used
        self.acid_atp = acid_atp['sb_type'].tolist()

        # BONDS
        colnames = ['ai', 'aj', 'func', 'func_type']
        section_dict = self.sections_dict['[ bonds ]']
        bonds_df = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        bonds_df.reset_index(inplace=True, drop=True)
        self.bonds = bonds_df

        bonds_pairs = list([(ai, aj) for ai, aj in zip(bonds_df['ai'].to_list(), bonds_df['aj'].to_list())])
        self.bond_pairs = bonds_pairs

        #bonds_df.rename(columns={'ai':'; ai'}, inplace=True)
        #self.string_bonds = bonds_df.to_string(index=False)

        colnames = ['ai', 'aj', 'ak', 'funct', 'c0']
        section_dict = self.sections_dict['[ angles ]']
        angles_df = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        angles_df.reset_index(inplace=True, drop=True) 
        self.angles = angles_df

        #angles_df.rename(columns={'ai':'; ai'}, inplace=True)
        #self.string_angles = angles_df.to_string(index=False)

        colnames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0']
        section_dict = self.sections_dict['[ dihedrals ]']
        dihedrals_df = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        dihedrals_df.reset_index(inplace=True, drop=True)
        self.dihedrals = dihedrals_df
        
        #dihedrals_df.rename(columns={'ai':'; ai'}, inplace=True)
        #self.string_dihedrals = dihedrals_df.to_string(index=False)

        colnames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0']
        section_dict = self.sections_dict['[ impropers ]']
        impropers_df = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        impropers_df.reset_index(inplace=True, drop=True)
        self.impropers = impropers_df
        
        #impropers_df.rename(columns={'ai':'; ai'}, inplace=True)
        #self.string_impropers = impropers_df.to_string(index=False)

        section_dict = self.sections_dict['[ system ]']
        system_df = pd.DataFrame.from_dict(section_dict, orient='index')
        system_df.reset_index(inplace=True, drop=True)
        system_df['; Name'] = 'Protein'
        system_df = system_df['; Name']
        self.system = system_df
        #self.string_system = system_df.to_string(index=False)

        colnames = ['; Compound', '#mols']
        section_dict = self.sections_dict['[ molecules ]']
        molecules_df = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        molecules_df.reset_index(inplace=True, drop=True)
        self.molecules = molecules_df
        #self.string_molecules = molecules_df.to_string(index=False)


class extra_topology_ligands:
    '''
    Cose
    '''
    def __init__(self, sections_dict_itp, sections_dict_prm, ligand_residue_number):
        pd.options.mode.chained_assignment = None  # default='warn'
        
        # PRM to get the c12
        section_dict_prm = sections_dict_prm['[ atomtypes ]']
        colnames = ['atom_type', 'atom_number', 'mass', 'charge', 'ptype', 'sigma', 'epsilon']
        df_ligand_atomtypes = pd.DataFrame.from_dict(section_dict_prm, orient='index', columns=colnames)
        df_ligand_atomtypes['sigma'] = df_ligand_atomtypes['sigma'].astype(float)
        df_ligand_atomtypes['epsilon'] = df_ligand_atomtypes['epsilon'].astype(float)
        df_ligand_atomtypes['c12'] = abs(4 * df_ligand_atomtypes['epsilon'] * (df_ligand_atomtypes['sigma'] ** 12))
        # Dictionary for ligand c12
        ligand_c12_dict = df_ligand_atomtypes[['atom_type', 'c12']].copy()
        ligand_c12_dict = ligand_c12_dict.set_index('atom_type')['c12'].to_dict()   
        self.ligand_c12_dict = ligand_c12_dict
    
        # ITP
        
        colnames = ['; Name', 'nrexcl']
        section_dict = sections_dict_itp['[ moleculetype ]']
        ligand_moleculetype_df = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        self.ligand_moleculetype = ligand_moleculetype_df

        # ATOMS
        section_dict_itp = sections_dict_itp['[ atoms ]']
        colnames = ['atom_number', 'atom_type', 'residue_number', 'residue', 'atom', 'cgnr', 'charge', 'mass', 'typeB', 'chargeB', 'massB']
        df_ligand_atoms = pd.DataFrame.from_dict(section_dict_itp, orient='index', columns=colnames)        
        # Changing the mass of the atoms section by adding the H
        df_ligand_atoms['atom_number'] = df_ligand_atoms['atom_number'].astype(int)
        df_ligand_atoms['residue_number'] = ligand_residue_number
        df_ligand_atoms['cgnr'] = df_ligand_atoms['cgnr'].astype(int)
        df_ligand_atoms['mass'] = df_ligand_atoms['mass'].astype(float)
        ## Structure based atomtype definition
        df_ligand_atoms['sb_type'] = df_ligand_atoms['atom'] + '_' + df_ligand_atoms['residue_number'].astype(str)
        df_ligand_atoms['c12'] = df_ligand_atoms['atom_type'].map(ligand_c12_dict)
        self.df_topology_ligand = df_ligand_atoms

        ligand_sbtype_c12_dict = df_ligand_atoms[['sb_type', 'c12']].copy()
        ligand_sbtype_c12_dict = ligand_sbtype_c12_dict.set_index('sb_type')['c12'].to_dict()   
        self.ligand_sbtype_c12_dict = ligand_sbtype_c12_dict
        
        #ligand_number_c12_dict = df_ligand_atoms[['atom_number', 'c12']].copy()
        #ligand_number_c12_dict = ligand_number_c12_dict.set_index('atom_number')['c12'].to_dict()   
        #self.ligand_number_c12_dict = ligand_number_c12_dict
        
        sbtype_ligand_number_dict = df_ligand_atoms[['sb_type', 'atom_number']]
        sbtype_ligand_number_dict = sbtype_ligand_number_dict.set_index('sb_type')['atom_number'].to_dict()   
        self.sbtype_ligand_number_dict = sbtype_ligand_number_dict

        ligand_sbtype_number_dict = df_ligand_atoms[['atom_number', 'sb_type']]
        ligand_sbtype_number_dict = ligand_sbtype_number_dict.set_index('atom_number')['sb_type'].to_dict()   
        self.ligand_sbtype_number_dict = ligand_sbtype_number_dict

        # BONDS
        section_dict_itp =  sections_dict_itp['[ bonds ]']
        colnames = ['ai', 'aj', 'funct', 'c0', 'c1']
        df_ligand_bonds = pd.DataFrame.from_dict(section_dict_itp, orient='index', columns=colnames)
        df_ligand_bonds['ai'] = df_ligand_bonds['ai'].astype(int)
        df_ligand_bonds['aj'] = df_ligand_bonds['aj'].astype(int)
        df_ligand_bonds['c1'] = df_ligand_bonds['c1'].astype(float)
        df_ligand_bonds['ai'] = df_ligand_bonds['ai'].map(ligand_sbtype_number_dict)
        df_ligand_bonds['aj'] = df_ligand_bonds['aj'].map(ligand_sbtype_number_dict)
        self.ligand_bonds =  df_ligand_bonds

        # ANGLES
        section_dict_itp =  sections_dict_itp['[ angles ]']
        colnames = ['ai', 'aj', 'ak', 'funct', 'c0', 'c1', 'c2', 'c3']
        df_ligand_angles = pd.DataFrame.from_dict(section_dict_itp, orient='index', columns=colnames)
        df_ligand_angles['ai'] = df_ligand_angles['ai'].astype(int)
        df_ligand_angles['aj'] = df_ligand_angles['aj'].astype(int)
        df_ligand_angles['ak'] = df_ligand_angles['ak'].astype(int)
        df_ligand_angles['ai'] = df_ligand_angles['ai'].map(ligand_sbtype_number_dict)
        df_ligand_angles['aj'] = df_ligand_angles['aj'].map(ligand_sbtype_number_dict)
        df_ligand_angles['ak'] = df_ligand_angles['ak'].map(ligand_sbtype_number_dict)
        self.ligand_angles =  df_ligand_angles

        # DIHEDRALS
        section_dict_itp = sections_dict_itp['[ dihedrals ]']
        colnames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
        df_ligand_dihedrals = pd.DataFrame.from_dict(section_dict_itp, orient='index', columns=colnames)
        df_ligand_dihedrals['ai'] = df_ligand_dihedrals['ai'].astype(int)
        df_ligand_dihedrals['aj'] = df_ligand_dihedrals['aj'].astype(int)
        df_ligand_dihedrals['ak'] = df_ligand_dihedrals['ak'].astype(int)
        df_ligand_dihedrals['al'] = df_ligand_dihedrals['al'].astype(int)
        df_ligand_dihedrals['ai'] = df_ligand_dihedrals['ai'].map(ligand_sbtype_number_dict)
        df_ligand_dihedrals['aj'] = df_ligand_dihedrals['aj'].map(ligand_sbtype_number_dict)
        df_ligand_dihedrals['ak'] = df_ligand_dihedrals['ak'].map(ligand_sbtype_number_dict)
        df_ligand_dihedrals['al'] = df_ligand_dihedrals['al'].map(ligand_sbtype_number_dict)
        self.ligand_dihedrals =  df_ligand_dihedrals

        # PAIRS
        # This is used when when want to read ligand pairs from the original topology
        # We might want to remove this part
        section_dict_itp = sections_dict_itp['[ pairs ]']
        colnames = ['ai', 'aj', 'funct']
        df_ligand_pairs = pd.DataFrame.from_dict(section_dict_itp, orient='index', columns=colnames)
        # TODO check if they want a value of pairs
        df_ligand_pairs['ai'] = df_ligand_pairs['ai'].astype(int)
        df_ligand_pairs['aj'] = df_ligand_pairs['aj'].astype(int)
        df_ligand_pairs['ai'] = df_ligand_pairs['ai'].map(ligand_sbtype_number_dict)
        df_ligand_pairs['aj'] = df_ligand_pairs['aj'].map(ligand_sbtype_number_dict)
        df_ligand_pairs['c6'] = 0.000000e+00
        df_ligand_pairs['c12_ai'] = df_ligand_pairs['ai'].map(ligand_sbtype_c12_dict)
        df_ligand_pairs['c12_aj'] = df_ligand_pairs['aj'].map(ligand_sbtype_c12_dict)
        # Now we keep it plain
        df_ligand_pairs['c12'] = (np.sqrt(df_ligand_pairs['c12_ai'] * df_ligand_pairs['c12_aj']))#*parameters['lj_reduction']
        df_ligand_pairs.drop(columns=['c12_ai', 'c12_aj'], inplace=True)
        self.ligand_pairs = df_ligand_pairs
        
        ## EXCLUSIONS
        # This one is obtained during the writing part
        #df_ligand_exclusions = df_ligand_pairs[['ai', 'aj']].copy()
        #self.ligand_exclusions = df_ligand_exclusions


class topology_ligand_bonds:
    def __init__(self, sections_dict):
        colnames = ['ai', 'aj', 'func', 'c0', 'c1']
        section_dict = sections_dict['[ bonds ]']
        bonds_df = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        bonds_df.reset_index(inplace=True, drop=True)

        self.df_topology_bonds = bonds_df
        bonds_pairs = list([(ai, aj) for ai, aj in zip(bonds_df['ai'].to_list(), bonds_df['aj'].to_list())])
        self.bond_pairs = bonds_pairs


class topology_ligand_angles:
    def __init__(self, sections_dict):
        colnames = ['ai', 'aj', 'ak', 'funct', 'c0', 'c1', 'c2', 'c3']
        section_dict = sections_dict['[ angles ]']
        angles_df = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        angles_df.reset_index(inplace=True, drop=True) 
        self.topology_angles = angles_df


class topology_ligand_dihedrals:
    def __init__(self, sections_dict):
        colnames = ['ai', 'aj', 'ak', 'al', 'c0', 'c1', 'c2', 'c3']
        section_dict = sections_dict['[ dihedrals ]']
        dihedrals_df = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        dihedrals_df.reset_index(inplace=True, drop=True)
        self.topology_dihedrals = dihedrals_df
