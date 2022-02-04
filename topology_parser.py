import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

def read_sections(topology_path, *topology_sections):
    sections_dict = {}
    for section in topology_sections:
        sections_dict[section] = line_parser(topology_path, section)
    return sections_dict


def line_parser(topology_path, section):
    section_dict = {}
    with open(topology_path) as f:
        for line_number, line in enumerate(f):
                if section in line:
                    section_dict['section_start'] = line_number
    section_dict = get_section_end(topology_path, section_dict)
    section_dict['lines_toskip'] = commented_lines(topology_path, section_dict)
    return section_dict


def get_section_end(topology_path, section_dict):
    with open(topology_path) as f:
        for line_number, line in enumerate(f):
            if not line.strip() and line_number > section_dict['section_start']:
                section_dict['section_end'] = line_number
                break
    
    return section_dict


def commented_lines(topology_path, section_dict):
    lines_toskip = list(range(0, section_dict['section_start']+1))
    with open(topology_path) as f:
        for line_number, line in enumerate(f):
            if line_number >= section_dict['section_start'] and line_number <= section_dict['section_end']:
                if line.startswith(';'):
                    lines_toskip.append(line_number)
    return lines_toskip

class topology_atoms:
    '''
    Cose
    '''
    def __init__(self, path):
        pd.options.mode.chained_assignment = None  # default='warn'
        # We read everything except charge, typeB, chargeB, massB
        colnames = ['atom_number', 'atom_type', 'residue_number', 'residue', 'atom', 'cgnr', 'mass']
        section_dict = read_sections(path, 'atoms')
        df_topology_atoms = pd.read_csv(path, names= colnames, usecols=[0,1,2,3,4,5,7], sep='\s+', header=None, nrows=section_dict['atoms']['section_end']-(len(section_dict['atoms']['lines_toskip'])), skiprows=section_dict['atoms']['lines_toskip'])

        # Changing the mass of the atoms section by adding the H
        df_topology_atoms['mass'].astype(float)
        # Adding H to backbone N

        # Adding an extra H to the N terminal
        mask = ((df_topology_atoms['residue'] == "PRO") & (df_topology_atoms['atom_type'] == 'N'))
        df_topology_atoms['mass'][mask] = df_topology_atoms['mass'][mask].astype(float).sub(1)
        # Removing an extra H to PRO 
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
        self.first_resid = 'N_'+str(df_topology_atoms['residue_number'][0])
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
    def __init__(self, path):
        colnames = ['ai', 'aj', 'func', 'func_type']
        section_dict = read_sections(path, 'bonds')
        topology_bonds = pd.read_csv(path, names=colnames, usecols=[0,1,2,3], sep='\s+', header=None, nrows=section_dict['bonds']['section_end']-(len(section_dict['bonds']['lines_toskip'])), skiprows=section_dict['bonds']['lines_toskip'])

        #for ai, aj in zip(topology_bonds['ai'].to_list(), topology_bonds['aj'].to_list()):
        #    print(ai, aj)

        bonds_pairs = list([(ai, aj) for ai, aj in zip(topology_bonds['ai'].to_list(), topology_bonds['aj'].to_list())])
        self.bond_pairs = bonds_pairs

