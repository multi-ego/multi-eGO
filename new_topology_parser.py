import pandas as pd

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
        # We read everything except charge, typeB, chargeB, massB
        colnames = ['atom_number', 'atom_type', 'residue_number', 'residue', 'atom', 'cgnr', 'mass']
        section_dict = read_sections(path, 'atoms')
        df_topology_atoms = pd.read_csv(path, names= colnames, usecols=[0,1,2,3,4,5,7], sep='\s+', header=None, nrows=section_dict['atoms']['section_end']-(len(section_dict['atoms']['lines_toskip'])), skiprows=section_dict['atoms']['lines_toskip'])
        self.df_topology_atoms = df_topology_atoms

class topology_bonds:
    def __init__(self, path):
        colnames = ['ai', 'aj', 'func', 'func_type']
        section_dict = read_sections(path, 'bonds')
        topology_bonds = pd.read_csv(path, names=colnames, usecols=[0,1,2,3], sep='\s+', header=None, nrows=section_dict['bonds']['section_end']-(len(section_dict['bonds']['lines_toskip'])), skiprows=section_dict['bonds']['lines_toskip'])

        #for ai, aj in zip(topology_bonds['ai'].to_list(), topology_bonds['aj'].to_list()):
        #    print(ai, aj)

        bonds_pairs = list([(ai, aj) for ai, aj in zip(topology_bonds['ai'].to_list(), topology_bonds['aj'].to_list())])
        self.bond_pairs = bonds_pairs


prova = topology_bonds(f'inputs/native_ABeta/topol.top')
print(len(prova.bond_pairs))