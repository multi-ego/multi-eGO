import os
import datetime

import gromologist as gml
from collections import OrderedDict


class Top:
    def __init__(self, filename, gmx_dir=None, pdb=None, ignore_ifdef=False, define=None):
        """
        A class to represent and contain the Gromacs topology file and provide
        tools for editing topology elements
        :param filename: str, path to the .top file
        :param gmx_dir: str, Gromacs FF directory
        :param pdb: str, path to a matching PDB file
        :param ignore_ifdef: bool, whether to ignore #include statements within #ifdef blocks (e.g. posre.itp)
        :param define: dict, key:value pairs with variables that will be defined in .mdp
        """
        # TODO maybe allow for construction of a blank top with a possibility to read data later?
        if not gmx_dir:
            self._gromacs_dir = self._find_gmx_dir()
        else:
            self._gromacs_dir = gmx_dir
        self.pdb = None
        self.pdb = None if pdb is None else gml.Pdb(pdb, top=self)
        self.fname = filename
        self.top = self.fname.split('/')[-1]
        self.dir = os.getcwd() + '/' + '/'.join(self.fname.split('/')[:-1])
        self._contents = open(self.fname).readlines()
        self.defines = {}
        if define is not None:
            self.defines.update(define)
        self._include_all(ignore_ifdef)
        self.sections = []
        self.header = []
        self._parse_sections()
        if self.top.endswith('top'):
            self._read_system_properties()
        else:
            self.system, self.charge, self.natoms = None, 0, 0
        
    def __repr__(self):
        return "Topology with {} atoms and total charge {:.3f}".format(self.natoms, self.charge)
    
    @classmethod
    def _from_text(cls, text, gmx_dir=None, pdb=None, ignore_ifdef=False):
        with open('tmp_topfile.gromo', 'w') as tmp:
            tmp.write(text)
        instance = cls('tmp_topfile.gromo', gmx_dir, pdb, ignore_ifdef)
        os.remove('tmp_topfile.gromo')
        return instance
    
    @staticmethod
    def _find_gmx_dir():
        """
        Attempts to find Gromacs internal files to fall back to
        when default .itp files are included using the
        #include statement
        :return: str, path to share/gromacs/top directory
        """
        gmx = os.popen('which gmx 2> /dev/null').read().strip()
        if not gmx:
            gmx = os.popen('which gmx_mpi 2> /dev/null').read().strip()
        if not gmx:
            gmx = os.popen('which gmx_d 2> /dev/null').read().strip()
        if not gmx:
            gmx = os.popen('which grompp 2> /dev/null').read().strip()
        if gmx:
            gmx = '/'.join(gmx.split('/')[:-2]) + '/share/gromacs/top'
            print('Gromacs files found in directory {}'.format(gmx))
            return gmx
        else:
            print('No working Gromacs compilation found, assuming all file dependencies are referred to locally')
            return ""
    
    @property
    def molecules(self):
        return [s for s in self.sections if isinstance(s, gml.SectionMol)]

    @property
    def alchemical_molecules(self):
        return [s for s in self.sections if isinstance(s, gml.SectionMol) and s.is_alchemical]
    
    @property
    def parameters(self):
        return [s for s in self.sections if isinstance(s, gml.SectionParam)][0]
    
    def list_molecules(self):
        """
        Prints out a list of molecules contained in the System
        :return: None
        """
        for mol in self.system.keys():
            print("{:20s}{:>10d}".format(mol, self.system[mol]))

    def clear_ff_params(self, section='all'):
        used_params = []  # TODO add one more layer to clean atomtypes/pairtypes/nonbondeds
        for mol in self.molecules:
            used_params.extend(mol.find_used_ff_params(section=section))
        self.parameters.clean_unused(used_params, section=section)
    
    def add_pdb(self, pdbfile):
        """
        Allows to pair a PDB file with the topology after the instance was initialized
        :param pdbfile: str, path to PDB file
        :return: None
        """
        self.pdb = gml.Pdb(pdbfile, top=self)

    def add_ff_params(self, section='all'):
        for mol in self.molecules:
            mol.add_ff_params(add_section=section)

    def find_missing_ff_params(self, section='all'):
        for mol in self.molecules:
            mol.find_missing_ff_params(add_section=section)
    
    def add_params_file(self, paramfile):
        prmtop = Top._from_text('#include {}\n'.format(paramfile))
        paramsect_own = self.parameters
        paramsect_other = prmtop.parameters
        paramsect_own.subsections.extend(paramsect_other.subsections)
        paramsect_own._merge()
    
    def _include_all(self, ign_ifdef):
        """
        includes all .itp files in the .top file to facilitate processing
        :return: None
        """
        # TODO optionally insert .mdp defines in constructor & pass here?
        ignore_lines = self._find_ifdef_lines() if ign_ifdef else set()
        lines = [i for i in range(len(self._contents)) if self._contents[i].strip().startswith("#include")
                 and i not in ignore_lines]
        while len(lines) > 0:
            lnum = lines[0]
            to_include, extra_prefix = self._find_in_path(self._contents.pop(lnum).split()[1].strip('"\''))
            contents = self._add_prefix_to_include(open(to_include).readlines(), extra_prefix)
            self._contents[lnum:lnum] = contents
            ignore_lines = self._find_ifdef_lines() if ign_ifdef else set()
            lines = [i for i in range(len(self._contents)) if self._contents[i].startswith("#include")
                     and i not in ignore_lines]
    
    def _find_ifdef_lines(self):
        """
        Finds #ifdef/#endif blocks in the topology if user
        explicitly asks to ignore them
        :return: set, int numbers of all lines that fall within an #ifdef/#endif block
        """
        ignore_set = set()
        counter = 0
        for n, line in enumerate(self._contents):
            if line.strip().startswith("#ifdef") or line.strip().startswith("#ifndef"):
                counter += 1
            elif line.strip().startswith("#endif"):
                counter -= 1
            if counter > 0:
                ignore_set.add(n)
        return ignore_set
    
    def clear_sections(self):
        """
        Removes all SectionMol instances that are not part
        of the system definition in [ system ]
        # TODO optionally we could also delete all unused params
        :return: None
        """
        if self.system is None:
            raise AttributeError("System properties have not been read, this is likely not a complete .top file")
        sections_to_delete = []
        for section_num, section in enumerate(self.sections):
            if isinstance(section, gml.SectionMol) and section.mol_name not in self.system.keys():
                sections_to_delete.append(section_num)
        self.sections = [s for n, s in enumerate(self.sections) if n not in sections_to_delete]
    
    def _find_in_path(self, filename):
        """
        looks for a file to be included in either the current directory
        or in Gromacs directories (as given by user), in order to
        include all .itp files in a single .top file
        :param filename: str, name of the file to be searched for
        :return: None
        """
        if filename.strip().startswith('./'):
            filename = filename.strip()[2:]
        pref = '/'.join(filename.split('/')[:-1])
        suff = filename.split('/')[-1]
        if os.path.isfile(self.dir.rstrip(os.sep) + os.sep + pref + os.sep + suff):
            return self.dir.rstrip(os.sep) + os.sep + pref + os.sep + suff, pref
        elif os.path.isfile(self._gromacs_dir.rstrip(os.sep) + os.sep + pref + os.sep + suff):
            return self._gromacs_dir.rstrip(os.sep) + os.sep + pref + os.sep + suff, pref
        else:
            raise FileNotFoundError('file {} not found in neither local nor Gromacs directory.\n'
                                    'If the file is included in an #ifdef/#ifndef block, please try setting'
                                    ' ignore_ifdef=True'.format(filename))
    
    @staticmethod
    def _add_prefix_to_include(content, prefix):
        """
        Modifies #include statements if nested #includes
        point to different directories
        :param content: list of str, content of the included file
        :param prefix: str, directory name to add in the nested include
        :return: list of str, modified content
        """
        if prefix:
            for nline, line in enumerate(content):
                if line.strip().startswith("#include"):
                    try:
                        index = line.index('"')
                    except ValueError:
                        index = line.index("'")
                    newline = line[:index+1] + prefix + '/' + line[index+1:]
                    content[nline] = newline
        return content
    
    def _parse_sections(self):
        """
        Cuts the content in sections as defined by the position
        of special headers, and builds the self.sections list
        :return: None
        """
        special_sections = {'defaults', 'moleculetype', 'system'}
        special_lines = [n for n, l in enumerate(self._contents)
                         if l.strip() and l.strip().strip('[]').strip().split()[0] in special_sections]
        special_lines.append(len(self._contents))
        for beg, end in zip(special_lines[:-1], special_lines[1:]):
            self.sections.append(self._yield_sec(self._contents[beg:end]))
        # in case there are #defines at the very beginning (e.g. CHARMM36):
        for lnum in range(0, special_lines[0]):
            if not self._contents[lnum].lstrip().startswith(';') and self._contents[lnum].strip():
                entry = gml.Entry(self._contents[lnum].strip(), self.sections[0].subsections[0])
                self.header.append(entry)

    def _yield_sec(self, content):
        """
        Chooses which class (Section or derived classes)
        should be used for the particular set of entries
        :param content: list of str, slice of self.content
        :return: Section (or its subclass) instance
        """
        if 'defaults' in content[0]:
            return gml.SectionParam(content, self)
        elif 'moleculetype' in content[0]:
            return gml.SectionMol(content, self)
        else:
            return gml.Section(content, self)
        
    def _read_system_properties(self):
        """
        Reads in system composition based on the [ molecules ] section
        and calculates the number of atoms and charge of the system
        :return: OrderedDict, contains molecule_name:number_of_molecules pairs and preserves the order read from file
                 float, total charge of the system
                 int, number of atoms in the system
        """
        system_subsection = [s.get_subsection('molecules') for s in self.sections
                             if 'molecules' in [ss.header for ss in s.subsections]]
        molecules = OrderedDict()  # we want to preserve the order of molecules in the system for e.g. PDB checking
        natoms, charge = 0, 0
        if len(system_subsection) > 1:
            raise RuntimeError("Multiple 'molecules' subsection found in the topology, this is not allowed")
        elif len(system_subsection) == 0:
            print("Section 'molecules' not present in the topology, assuming this is an isolated .itp")
            self.system, self.charge, self.natoms = None, 0, 0
            return
        for e in system_subsection[0]:
            if e.content:
                molecules[e.content[0]] = int(e.content[1])
        for mol in molecules.keys():
            sect_mol = self.get_molecule(mol)
            natoms += molecules[mol] * sect_mol.natoms
            charge += molecules[mol] * sect_mol.charge
        self.system, self.charge, self.natoms = molecules, charge, natoms
    
    def recalc_sys_params(self):
        """
        Wrapper around read_system_properties that also
        recalculates relevant molecule properties after changes
        :return: None
        """
        sub_mol = [sub for sect in self.sections for sub in sect.subsections if isinstance(sub, gml.SubsectionAtom)]
        for sub in sub_mol:
            sub.calc_properties()
        self._read_system_properties()

    def explicit_defines(self):
        """
        Changes pre-defined keywords in parameter sets
        according to #define entries in FF params
        :return: None
        """
        self.parameters._get_defines()
        for m in self.molecules:
            for s in m.subsections:
                if isinstance(s, gml.SubsectionBonded):
                    s.explicit_defines()
        
    def get_molecule(self, mol_name):
        """
        Finds a molecule (SectionMol instance) whose mol_name
        matches the query name
        :param mol_name: name of the molecule
        :return: SectionMol instance
        """
        mol = [s for s in self.sections if isinstance(s, gml.SectionMol) and s.mol_name == mol_name]
        if len(mol) == 0:
            raise KeyError("Molecule {} is not defined in topology".format(mol_name))
        elif len(mol) > 1:
            raise RuntimeError("Molecule {} is duplicated in topology".format(mol_name))
        return mol[0]
    
    def check_pdb(self):
        """
        Compares the topology with a PDB object to check
        for consistency, just as gmx grompp does;
        if inconsistencies are found, prints a report
        :return: None
        """
        if self.pdb:
            self.pdb.check_top()
        else:
            raise AttributeError("No PDB file has been bound to this topology")
    
    def save_top(self, outname='merged.top', split=False):
        """
        Saves the combined topology to the specified file
        :param outname: str, file name for output
        :param split: bool, whether to split into individual .top files
        :return: None
        """
        outfile = open(outname, 'w')
        self._write_header(outfile)
        if not split:
            for section in self.sections:
                self._write_section(outfile, section)
        else:
            for section in self.sections:
                if isinstance(section, gml.SectionParam):
                    with open('ffparams.itp', 'w') as out_itp:
                        self._write_section(out_itp, section)
                    outfile.write('\n; Include ff parameters\n#include "ffparams.itp"\n')
                elif isinstance(section, gml.SectionMol):
                    with open('{}.itp'.format(section.mol_name), 'w') as out_itp:
                        self._write_section(out_itp, section)
                    outfile.write('\n; Include {mn} topology\n#include "{mn}.itp"\n'.format(mn=section.mol_name))
                else:
                    self._write_section(outfile, section)
        outfile.close()

    @staticmethod
    def _write_section(outfile, section):
        for subsection in section.subsections:
            outfile.write('\n[ {} ]\n'.format(subsection.write_header))
            for entry in subsection:
                str_entry = str(entry).rstrip() + '\n'
                outfile.write(str_entry)
                
    def _write_header(self, outfile):
        outname = outfile.name.split('/')[-1]
        outfile.write(";\n;  File {} was generated with the gromologist library\n"
                      ";  by user: {}\n;  on host: {}\n;  at date: {} \n;\n".format(outname,
                                                                                    os.getenv("USER"),
                                                                                    os.uname()[1],
                                                                                    datetime.datetime.now()))
        for entry in self.header:
            str_entry = str(entry).rstrip() + '\n'
            outfile.write(str_entry)
