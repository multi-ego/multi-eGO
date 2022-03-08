import pandas as pd
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


class topology_moleculetype:
    def __init__(self, sections_dict):
        colnames = ['; Name', 'nrexcl']
        section_dict = sections_dict['[ moleculetype ]']
        moleculetype_df = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        self.topology_moleculetype = moleculetype_df


class topology_atoms:
    '''
    Cose
    '''
    def __init__(self, sections_dict):
        pd.options.mode.chained_assignment = None  # default='warn'
        colnames = ['atom_number', 'atom_type', 'residue_number', 'residue', 'atom', 'cgnr', 'charge', 'mass', 'typeB', 'chargeB', 'massB']

        section_dict = sections_dict['[ atoms ]']
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

class topology_bonds:
    def __init__(self, sections_dict):
        colnames = ['ai', 'aj', 'func', 'func_type']
        section_dict = sections_dict['[ bonds ]']
        bonds_df = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        bonds_df.reset_index(inplace=True, drop=True)

        self.df_topology_bonds = bonds_df
        bonds_pairs = list([(ai, aj) for ai, aj in zip(bonds_df['ai'].to_list(), bonds_df['aj'].to_list())])
        self.bond_pairs = bonds_pairs


class topology_angles:
    def __init__(self, sections_dict):
        colnames = ['ai', 'aj', 'ak', 'funct', 'c0']
        section_dict = sections_dict['[ angles ]']
        angles_df = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        angles_df.reset_index(inplace=True, drop=True) 
        self.topology_angles = angles_df


class topology_dihedrals:
    def __init__(self, sections_dict):
        colnames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0']
        section_dict = sections_dict['[ dihedrals ]']
        dihedrals_df = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        dihedrals_df.reset_index(inplace=True, drop=True)
        self.topology_dihedrals = dihedrals_df


class topology_impropers:
    def __init__(self, sections_dict):
        colnames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0']
        section_dict = sections_dict['[ impropers ]']
        impropers_df = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        impropers_df.reset_index(inplace=True, drop=True)
        self.topology_impropers = impropers_df


class topology_system:
    def __init__(self, sections_dict):
        section_dict = sections_dict['[ system ]']
        system_df = pd.DataFrame.from_dict(section_dict, orient='index')
        system_df.reset_index(inplace=True, drop=True)
        system_df['; Name'] = system_df[0] + system_df[1] + system_df[2] + system_df[3] + system_df[4] + system_df[5] + system_df[6]
        system_df = system_df['; Name']
        self.topology_system = system_df


class topology_molecules:
    def __init__(self, sections_dict):
        colnames = ['; Compound', '#mols']
        section_dict = sections_dict['[ molecules ]']
        molecules_df = pd.DataFrame.from_dict(section_dict, orient='index', columns=colnames)
        molecules_df.reset_index(inplace=True, drop=True)
        self.topology_molecules = molecules_df
