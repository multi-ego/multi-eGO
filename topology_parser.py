# Questo mi serve per leggere le topologie della fibrilla e della nativa
# Confrontarle per vedere la numerazione degli atomi e le cose sopra
# Sistemarle di conseguenza

from os import EX_SOFTWARE, linesep
import pandas as pd
#from read_input import topology_path

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
    
#topology_path = 'inputs/native_ABeta/topol.top'
#section_dict = read_sections(topology_path, 'atoms')#, 'bonds')



class topology_atoms:
    '''
    From .top file, read the atom section.
    '''
    def __init__(self, path):
    # We don't keep everything
        colnames = ['atom_number', 'atom_type', 'residue_number', 'residue', 'atom', 'cgnr', 'mass']
        section_dict = read_sections(path, 'atoms')
        df_topology_atoms = pd.read_csv(path, names= colnames, usecols=[0,1,2,3,4,5,7], sep='\s+', header=None, nrows=section_dict['atoms']['section_end']-(len(section_dict['atoms']['lines_toskip'])), skiprows=section_dict['atoms']['lines_toskip'])

        # Multi-eGO uses implicit Hydrogens which has to be added in the mass section.
        mask = ((df_topology_atoms['atom_type'] == 'N') | (df_topology_atoms['atom_type'] == 'OA')) | ((df_topology_atoms['residue'] == 'TYR') & ((df_topology_atoms['atom'] == 'CD1') | (df_topology_atoms['atom'] == 'CD2') | (df_topology_atoms['atom'] == 'CE1') | (df_topology_atoms['atom'] == 'CE2')))
        #df_topology_atoms['mass'][mask] = df_topology_atoms['mass'][mask].add(1)
        df_topology_atoms.loc[mask, 'mass'] = df_topology_atoms.loc[mask, 'mass'].add(1)

        # Same thing here with the N terminal which is charged
        mask = (df_topology_atoms['residue_number'] == 1) & (df_topology_atoms['atom'] == 'N')
        #df_topology_atoms['mass'][mask] = df_topology_atoms['mass'][mask].add(2)
        df_topology_atoms.loc[mask, 'mass'] = df_topology_atoms.loc[mask, 'mass'].add(2)

        # Multi-eGO atomtype definition
        df_topology_atoms['sb_type'] = df_topology_atoms['atom'] + '_' + df_topology_atoms['residue_number'].astype('string')
        self.df_topology_atoms = df_topology_atoms
        self.list_atom_numbers = df_topology_atoms['atom_number'].to_list()
        self.list_atom_type = df_topology_atoms['atom_type'].to_list()
        self.list_residue_number = df_topology_atoms['residue_number'].to_list()
        self.list_residue = df_topology_atoms['residue'].to_list()
        self.list_atoms = df_topology_atoms['atom'].to_list()
        self.list_atoms_mass = df_topology_atoms['mass'].to_list()

        # This is needed when we want to do some stuff only to the N terminal
        # TODO has to be resid == 1
        #first_resid = 'N_'+str(df_topology_atoms['residue'][0])

        # ACID pH
        # Selection of the aminoacids and the charged atoms (used for B2m)
        # TODO add some options for precise pH setting
        acid_ASP = df_topology_atoms[(df_topology_atoms['residue'] == "ASP") & ((df_topology_atoms['atom'] == "OD1") | (df_topology_atoms['atom'] == "OD2") | (df_topology_atoms['atom'] == "CG"))]
        acid_GLU = df_topology_atoms[(df_topology_atoms['residue'] == "GLU") & ((df_topology_atoms['atom'] == "OE1") | (df_topology_atoms['atom'] == "OE2") | (df_topology_atoms['atom'] == "CD"))]
        acid_HIS = df_topology_atoms[(df_topology_atoms['residue'] == "HIS") & ((df_topology_atoms['atom'] == "ND1") | (df_topology_atoms['atom'] == "CE1") | (df_topology_atoms['atom'] == "NE2") | (df_topology_atoms['atom'] == "CD2") | (df_topology_atoms['atom'] == "CG"))]
        list_acid_residues = pd.concat([acid_ASP, acid_GLU, acid_HIS], ignore_index = True)
        list_acid_residues = list_acid_residues['sb_type'].tolist()


class topology_bonds():
    def __init__(self, path):
        colnames = ['ai', 'aj', 'func', 'func_type']
        section_dict = read_sections(self.path, 'bonds')
        topology_bonds = pd.read_csv(self.path, names=colnames, usecols=[0,1,2,3], sep='\s+', header=None, nrows=section_dict['bonds']['section_end']-(len(section_dict['bonds']['lines_toskip'])), skiprows=section_dict['bonds']['lines_toskip'])

        #for ai, aj in zip(topology_bonds['ai'].to_list(), topology_bonds['aj'].to_list()):
        #    print(ai, aj)

        bonds_pairs = list([(ai, aj) for ai, aj in zip(topology_bonds['ai'].to_list(), topology_bonds['aj'].to_list())])

#atoms = topology_atoms(f'inputs/native_ABeta/topol.top')
#print(atoms.df_topology_atoms)
#from topology_definitions import raw_topology_atoms
# Are the same
#print(atoms.df_topology_atoms)
#print(raw_topology_atoms)

