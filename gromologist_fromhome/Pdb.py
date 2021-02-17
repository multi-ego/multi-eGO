import gromologist as gml
from collections import defaultdict


class Pdb:  # TODO optionally save as gro? & think of trajectories
    def __init__(self, filename=None, top=None, altloc='A', **kwargs):
        self.fname = filename
        if self.fname:
            if self.fname.endswith('gro'):
                self.atoms, self.box, self.remarks = self._parse_contents_gro([line.rstrip()
                                                                               for line in open(self.fname)])
            else:
                self.atoms, self.box, self.remarks = self._parse_contents([line.strip() for line in open(self.fname)])
        else:
            self.atoms, self.box, self.remarks = [], 3 * [100] + 3 * [90], []
        self.top = top if not isinstance(top, str) else gml.Top(top, **kwargs)
        if self.top and not self.top.pdb:
            self.top.pdb = self
        if not self.atoms and self.top:
            self.atoms = [Atom.from_top_entry(entry) for mol in self.top.molecules for entry in mol.atoms]
        self.altloc = altloc
        self._atom_format = "ATOM  {:>5d} {:4s}{:1s}{:4s}{:1s}{:>4d}{:1s}   " \
                            "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
        self._cryst_format = "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           1\n"

    def __repr__(self):
        return "PDB file {} with {} atoms".format(self.fname, len(self.atoms))

    @classmethod
    def from_selection(cls, pdb, selection):
        selected_indices = pdb.select_atoms(selection)
        new_pdb = Pdb()
        new_pdb.atoms = [atom for n, atom in enumerate(pdb.atoms) if n in selected_indices]
        new_pdb.box = pdb.box
        new_pdb.remarks = pdb.remarks
        new_pdb.altloc = pdb.altloc
        return new_pdb

    @classmethod
    def from_text(cls, text, orig_pdb):
        new_pdb = Pdb()
        try:
            ftype = orig_pdb.fname[-3:]
        except:
            ftype = orig_pdb.sfname[-3:]
        if ftype == 'pdb':
            new_pdb.atoms = cls._parse_contents([line.strip() for line in text.split('\n')])
        elif ftype == 'gro':
            new_pdb.atoms = cls._parse_contents_gro([line.strip() for line in text.split('\n')])
        new_pdb.box = orig_pdb.box
        new_pdb.remarks = orig_pdb.remarks
        new_pdb.altloc = orig_pdb.altloc
        return new_pdb

    def print_protein_sequence(self):
        mapping = {'ALA': 'A', 'CYS': 'C', 'CYX': 'C', 'CYM': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G',
                   'HIS': 'H', 'HIE': 'H', 'HID': 'H', 'HSD': 'H', 'HSE': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
                   'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
                   'TRP': 'W', 'TYR': 'Y', "GLUP": "E", "ASPP": "D"}
        chains = list({a.chain for a in self.atoms if a.atomname == 'CA'})
        if chains == [' ']:
            print("No chains in the molecule, run Pdb.add_chains() first")
            return
        for ch in sorted(chains):
            cas = set(self.select_atoms(f'name CA and chain {ch}'))
            atoms = [a for n, a in enumerate(self.atoms) if n in cas]
            print(''.join([mapping[i.resname] for i in atoms]))

    def print_nucleic_sequence(self):
        mapping = defaultdict(lambda: 'X')
        mapping.update({'DA': "A", 'DG': "G", 'DC': "C", 'DT': "T", 'DA5': "A", 'DG5': "G", 'DC5': "C", 'DT5': "T",
                        'DA3': "A", 'DG3': "G", 'DC3': "C", 'DT3': "T", 'RA': "A", 'RG': "G", 'RC': "C", 'RU': "U",
                        'RA5': "A", 'RG5': "G", 'RC5': "C", 'RU5': "U", 'RA3': "A", 'RG3': "G", 'RC3': "C", 'RU3': "U",
                        'A': "A", 'G': "G", 'C': "C", 'U': "U", 'A5': "A", 'G5': "G", 'C5': "C", 'U5': "U",
                        'A3': "A", 'G3': "G", 'C3': "C", 'U3': "U"})
        chains = list({a.chain for a in self.atoms if a.atomname == "O4'"})
        for ch in sorted(chains):
            cas = set(self.select_atoms(f"name O4' and chain {ch}"))
            atoms = [a for n, a in enumerate(self.atoms) if n in cas]
            print(''.join([mapping[i.resname] for i in atoms]))

    def find_missing(self):
        map_pro = {'ALA': 'A', 'CYS': 'C', 'CYX': 'C', 'CYM': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G',
                   'HIS': 'H', 'HIE': 'H', 'HID': 'H', 'HSD': 'H', 'HSE': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
                   'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
                   'TRP': 'W', 'TYR': 'Y', "GLUP": "E", "ASPP": "D"}
        map_nuc = {'DA': "A", 'DG': "G", 'DC': "C", 'DT': "T", 'DA5': "A", 'DG5': "G", 'DC5': "C", 'DT5': "T",
                   'DA3': "A", 'DG3': "G", 'DC3': "C", 'DT3': "T", 'RA': "A", 'RG': "G", 'RC': "C", 'RU': "U",
                   'RA5': "A", 'RG5': "G", 'RC5': "C", 'RU5': "U", 'RA3': "A", 'RG3': "G", 'RC3': "C", 'RU3': "U",
                   'A': "A", 'G': "G", 'C': "C", 'U': "U", 'A5': "A", 'G5': "G", 'C5': "C", 'U5': "U",
                   'A3': "A", 'G3': "G", 'C3': "C", 'U3': "U"}
        pro_bb = ['N', 'O', 'C', 'CA']
        pro_sc = {'A': ['CB'], 'C': ['CB', 'SG'], 'D': ['CB', 'CG', 'OD1', 'OD2'], 'E': ['CB', 'CG', 'CD', 'OE1',
                  'OE2'], 'F': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'], 'G': [], 'H': ['CB', 'CG', 'ND1', 'CE1',
                  'CD2', 'NE2'], 'I': ['CB', 'CG1', 'CG2', 'CD'], 'K': ['CB', 'CG', 'CD', 'CE', 'NZ'], 'L': ['CB', 'CG',
                  'CD1', 'CD2'], 'M': ['CB', 'CG', 'SD', 'CE'], 'N': ['CB', 'CG', 'OD1', 'ND2'], 'P': ['CB', 'CG',
                  'CD'], 'Q': ['CB', 'CG', 'CD', 'OE1', 'NE2'], 'R': ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
                  'S': ['CB', 'OG'], 'T': ['CB', 'CG2', 'OG1'], 'V': ['CB', 'CG1', 'CG2'], 'W': ['CB', 'CG', 'CD1',
                  'CD2', 'NE1', 'CE2', 'CE3', 'CH2', 'CZ2', 'CZ3'], 'Y': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ',
                  'OH']}
        alt = {'I': ['CD', 'CD1']}  # pretty temporary and non-extensible, need to rethink
        curr_res = ('X', 0)
        atomlist = []
        for at in self.atoms:
            if (at.resname, at.resnum) != curr_res:
                if 'CA' in atomlist:
                    full = set(pro_bb + pro_sc[map_pro[curr_res[0]]])
                    if not full.issubset(set(atomlist)):
                        if map_pro[curr_res[0]] in alt.keys():
                            modfull = set([alt[map_pro[curr_res[0]]][1] if x == alt[map_pro[curr_res[0]]][0] else x
                                           for x in full])
                            if not modfull.issubset(set(atomlist)):
                                print(f"atoms {modfull.difference(set(atomlist))} missing from residue {curr_res}")
                        else:
                            print(f"atoms {full.difference(set(atomlist))} missing from residue {curr_res}")

                # TODO implement for nucleic
                curr_res = (at.resname, at.resnum)
                atomlist = []
            atomlist.append(at.atomname)

    def add_chains(self, serials=None, chain=None, offset=0, maxwarn=100):
        """
        Given a matching Top instance, adds chain identifiers to atoms
        based on the (previously verified) matching between invididual
        molecules defined in Top and consecutive atoms in this Pdb instance.
        Solvent molecules are ommited, and only characters from A to Z
        are supported as valid chain identifiers.
        Optionally, given `serials` and `chain`, one can set all atoms
        with atom numbers in `serials` as belonging to chain `chain`.
        :param serials: iterable of ints, atom numbers to set chain for (default is all)
        :param chain: chain to set for serials (default is use consecutive letters)
        :param offset: int, start chain ordering from letter other than A
        :param maxwarn: int, max number of warnings before an error shows up
        :return: None
        """
        if not self.top:
            ch = chain if chain is not None else 'X'
            for atom in self.atoms:
                if atom.chain == ' ':
                    atom.chain = ch
        else:
            self.check_top(maxwarn=maxwarn)
            if (serials is None and chain is not None) or (serials is not None and chain is None):
                raise ValueError("Both serials and chain have to be specified simultaneously")
            excluded = {'SOL', 'HOH', 'TIP3', 'K', 'NA', 'CL', 'POT'}
            base_char = 65 + offset  # 65 is ASCII for "A"
            index = 0
            for mol_name in self.top.system.keys():
                if mol_name.upper() not in excluded:
                    n_mols = self.top.system[mol_name]
                    mol = self.top.get_molecule(mol_name)
                    n_atoms = mol.natoms
                    for m in range(n_mols):
                        for a in range(n_atoms):
                            if serials is None and chain is None:
                                self.atoms[index].chain = chr(base_char)
                            else:
                                if self.atoms[index].serial in serials:
                                    self.atoms[index].chain = chain
                            index += 1
                        base_char += 1
                        if base_char > 90:
                            return

    def check_top(self, maxwarn=20, fix_pdb=False, fix_top=False):
        if fix_pdb and fix_top:
            raise ValueError("You either want to fix topology or pdb naming")
        if self.top is None:
            raise ValueError("a Top object has not been assigned; molecule info missing")
        index, err = 0, 0
        self._remove_altloc()
        for mol_name in self.top.system.keys():
            mol = self.top.get_molecule(mol_name)
            n_mols = self.top.system[mol_name]
            atom_subsection = mol.get_subsection('atoms')
            atom_entries = [e for e in atom_subsection if isinstance(e, gml.EntryAtom)]
            for m in range(n_mols):
                for n, a in enumerate(atom_entries):
                    try:
                        rtrn = self._check_mismatch(atom_entries[n], self.atoms[index], mol_name)
                    except IndexError:
                        raise RuntimeError("Mismatch encountered: PDB has {} atoms while topology "
                                           "has {}".format(n + 1, index + 1, len(self.atoms), self.top.natoms))
                    if rtrn:
                        if fix_pdb:
                            self.atoms[index].atomname = atom_entries[n].atomname
                        elif fix_top:
                            atom_entries[n].atomname = self.atoms[index].atomname
                    err += rtrn
                    index += 1
                    if err > maxwarn > -1:
                        raise RuntimeError("Error: too many warnings; use maxwarn=N to allow for up to N exceptions,"
                                           "or -1 to allow for any number of them")
        print("Check passed, all names match")
    
    @staticmethod
    def _check_mismatch(atom_entry, atom_instance, mol_name):
        if atom_entry.atomname != atom_instance.atomname:
            print("Atoms {} ({}) in molecule {} topology and {} ({}) in .pdb have "
                  "non-matching names".format(atom_entry.num, atom_entry.atomname, mol_name,
                                              atom_instance.serial, atom_instance.atomname))
            return 1
        return 0
        
    def _remove_altloc(self):
        """
        We only keep one of the alternative locations in case
        there is more (by default, A is kept)
        :return: None
        """
        self.atoms = [a for a in self.atoms if a.altloc in [' ', self.altloc]]
    
    def _write_atom(self, atom):
        atom.serial %= 100000
        atom.resnum %= 10000
        if atom.occ > 1000:
            atom.occ %= 1000
        if atom.beta > 1000:
            atom.beta %= 1000
        return self._atom_format.format(atom.serial, atom.atomname, atom.altloc, atom.resname, atom.chain, atom.resnum,
                                        atom.insert, atom.x, atom.y, atom.z, atom.occ, atom.beta, atom.element)
    
    def renumber_atoms(self):
        for n, atom in enumerate(self.atoms, 1):
            atom.serial = n
    
    def renumber_residues(self):
        count = 1
        for n in range(len(self.atoms)):
            temp = count
            try:
                if self.atoms[n].resnum != self.atoms[n+1].resnum or self.atoms[n].chain != self.atoms[n+1].chain:
                    temp = count + 1
            except IndexError:
                pass
            self.atoms[n].resnum = count
            count = temp

    def reposition_atom_from_hook(self, atomsel, hooksel, bondlength, p1_sel, p2_sel):
        assert 1 == len(self.select_atoms(atomsel)) == len(self.select_atoms(hooksel))
        assert 1 == len(self.select_atoms(p1_sel)) == len(self.select_atoms(p2_sel))
        movable = self.atoms[self.select_atoms(atomsel)[0]]
        hook_xyz = self.get_coords()[self.select_atoms(hooksel)[0]]
        p1_xyz = self.get_coords()[self.select_atoms(p1_sel)[0]]
        p2_xyz = self.get_coords()[self.select_atoms(p2_sel)[0]]
        vec = [x2 - x1 for x1, x2 in zip(p1_xyz, p2_xyz)]
        vec_len = sum([x**2 for x in vec])**0.5
        scale = bondlength/vec_len
        movable.set_coords([h + scale*v for h, v in zip(hook_xyz, vec)])
    
    def match_order_by_top_names(self, arange=None):
        """
        Whenever PDB atoms have different ordering than .top ones
        but naming is consistent, we can use the ordering from .top
        to reorder PDB atoms.
        This can be done for the entire system if molecules are unique;
        otherwise, range has to be specified (using the Pythonic convention)
        to avoid ambiguities. In that case, matching molecules will only be
        looked for in the specified range.
        :param arange: tuple, start and end point of the modification (end is excluded)
        :return:
        """
        if self.top is None:
            raise ValueError("a Top object has not been assigned; molecule info missing")
        new_atoms = []
        index = 0
        for mol_name in self.top.system.keys():
            mol = self.top.get_molecule(mol_name)
            n_mols = self.top.system[mol_name]
            atom_subsection = mol.get_subsection('atoms')
            atom_entries = [e for e in atom_subsection if isinstance(e, gml.EntryAtom)]
            for _ in arange(n_mols):
                for a in atom_entries:
                    if not arange or arange[0] <= index < arange[1]:
                        pdb_loc = self.select_atoms("resname {} and resid {} and name {}".format(a.resname, a.resid,
                                                                                                 a.atomname))
                        pdb_loc = [loc for loc in pdb_loc if not arange or arange[0] <= index < arange[1]]
                        if len(pdb_loc) != 1:
                            raise ValueError("Could not proceed; for match-based renumbering, residue numberings "
                                             "have to be consistent between PDB and .top, atom names need to match, "
                                             "and molecules cannot be repeated.\nError encountered when processing"
                                             "residue {} with resid {}, atom name {}".format(a.resname, a.resid,
                                                                                             a.atomname))
                        new_atoms.append(self.atoms[list(pdb_loc)[0]])
                    index += 1
        self.atoms = new_atoms
    
    def insert_atom(self, index, base_atom, **kwargs):
        """
        Inserts an atom into the atomlist. The atom is defined by
        providing a base Atom instance and a number of keyworded modifiers,
        e.g. atomname="CA", serial=15, ...
        :param index: int, the new index of the inserted Atom in the Atoms list
        :param base_atom: Atom, instance based on which the new atom will be based
        :param kwargs: Atom attributes to be held by the new atom
        :return: None
        """
        new_atom = Atom(self._write_atom(base_atom))
        for kw in kwargs.keys():
            if kw not in {"serial", "atomname", "resname", "chain", "resnum", "x", "y", "z", "occ", "beta", "element"}:
                raise ValueError("{} is not a valid Atom attribute")
            new_atom.__setattr__(kw, kwargs[kw])
        self.atoms.insert(index, new_atom)
    
    def delete_atom(self, serial, renumber=False):
        num = [n for n, a in enumerate(self.atoms) if a.serial == serial]
        if len(num) == 0:
            raise ValueError('No atoms with serial number {}'.format(serial))
        elif len(num) > 1:
            raise ValueError('Multiple atoms with serial number {}; consider renumbering'.format(serial))
        atom = self.atoms.pop(num[0])
        print('Entry {} deleted from PDB'.format(str(atom)))
        if renumber:
            self.renumber_atoms()

    def select_atoms(self, selection_string):
        sel = gml.SelectionParser(self)
        return sel(selection_string)
    
    def same_residue_as(self, query_iter):
        new_list = []
        for atom in query_iter:
            residue, resid = self.atoms[atom].resname, self.atoms[atom].resnum
            matching = [n for n, a in enumerate(self.atoms) if a.resname == residue and a.resnum == resid]
            new_list.extend(matching)
        return set(new_list)
    
    def within(self, query_iter, threshold):
        new_list = []
        for n, atom in enumerate(self.atoms):
            if any([self._atoms_dist(atom, self.atoms[query]) <= threshold for query in query_iter]):
                new_list.append(n)
        return set(new_list)
    
    @staticmethod
    def _atoms_dist(at1, at2):
        return ((at2.x - at1.x)**2 + (at2.y - at1.y)**2 + (at2.z - at1.z)**2)**0.5
        
    @staticmethod
    def _parse_contents(contents):
        atoms, remarks = [], []
        box = [7.5, 7.5, 7.5, 90, 90, 90]  # generic default, will be overwritten if present
        for line in contents:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atoms.append(Atom(line))
            elif line.startswith("CRYST1"):
                box = [float(line[6+9*a:6+9*(a+1)]) for a in range(3)] + \
                      [float(line[33+7*a:33+7*(a+1)]) for a in range(3)]
            elif not line.startswith('TER') and not line.startswith('END'):
                remarks.append(line)
            if line.startswith('END') or line.startswith('ENDMDL'):
                break
        return atoms, tuple(box), remarks

    @staticmethod
    def _parse_contents_gro(contents):
        import math
        contents = [x for x in contents if x.strip()]
        atoms, remarks = [], []
        header = contents[0]
        remarks.append("TITLE     {}".format(header))
        natoms = int(contents[1])
        for line in contents[2:2+natoms]:
            atoms.append(Atom.from_gro(line))
        if len(contents[-1].split()) == 3:
            box = [float(10*x) for x in contents[-1]] + [90., 90., 90.]
        elif len(contents[-1].split()) == 9:
            boxline = [float(x) for x in contents[-1].split()]
            assert boxline[3] == boxline[4] == boxline[6] == 0
            box = [0.0] * 6
            box[0] = boxline[0]
            box[-1] = math.atan(boxline[1]/boxline[5]) if boxline[5] != 0 else 90.0  # TODO check twice
            box[1] = boxline[1]/math.sin(box[-1])
            box[2] = math.sqrt(boxline[7]**2 + boxline[8]**2 + boxline[2]**2)
            box[-2] = math.acos(boxline[7]/box[2])
            box[-3] = math.acos((boxline[8]*math.sin(box[-1]))/box[2] + math.cos(box[-2]) * math.cos(box[-1]))
            box[0], box[1], box[2] = box[0]*10, box[1]*10, box[2]*10
            box[3], box[4], box[5] = box[3]*180/math.pi, box[4]*180/math.pi, box[5]*180/math.pi
        else:
            raise RuntimeError('Can\'t read box properties')
        return atoms, tuple(box), remarks
        
    def fill_beta(self, values, serials=None, smooth=False, ignore_mem=False):
        """
        Enables user to write arbitrary values to the beta field
        of the PDB entry
        :param values: iterable; values to fill
        :param serials: iterable, optional; can be used to specify a subset of atoms
        :param smooth: if float, defines sigma (in Angstrom) for beta-value smoothing
        :param ignore_mem: bool, allows to ignore memory warnings
        :return: None
        """
        if any([v > 999 for v in values]):
            print("At least one value is too large to fit into the `beta` field, consider division "
                  "to make values smaller")
        if not serials:
            serials = [a.serial for a in self.atoms]
        if not all([i == j for i, j in zip(serials, sorted(serials))]):
            raise ValueError("Atoms' serial numbers are not sorted, consider running .renumber_all() first")
        if len(values) != len(serials):
            raise ValueError("Lists 'value' and 'serials' have inconsistent sizes: {} and {}".format(len(values),
                                                                                                     len(serials)))
        index = 0
        if not smooth:
            for atom in self.atoms:
                if atom.serial in serials:
                    atom.beta = values[index]
                    index += 1
        else:
            if len(serials) > 10000 and not ignore_mem:
                raise RuntimeError("Try to restrict the number of atoms (e.g. selecting CA only), or you're risking "
                                   "running out of memory. To proceed anyway, run again with ignore_mem=True")
            try:
                import numpy as np
            except ImportError:
                print("Needs numpy for smoothing, try installing it")
            else:
                atomnums = [n for n, atom in enumerate(self.atoms) if atom.serial in serials]
                coords = np.array(self.get_coords())[atomnums]
                for atom in self.atoms:
                    dists = np.linalg.norm(coords - np.array(atom.coords), axis=1)
                    weights = np.exp(-(dists**2/(2*smooth)))
                    weights /= np.sum(weights)
                    atom.beta = np.sum(values * weights)

    def interpolate_struct(self, other, num_inter, write=False):
        inter = []
        try:
            import numpy as np
            from copy import deepcopy
        except ImportError:
            raise RuntimeError("Needs numpy & deepcopy for interpolating, try installing it")
        else:
            self_atoms = np.array(self.get_coords())
            other_atoms = np.array(other.get_coords())
            for i in range(1, num_inter+1):
                incr = (other_atoms - self_atoms)/(num_inter + 1)
                pdb = deepcopy(self)
                pdb.set_coords(self_atoms + i * incr)
                inter.append(pdb)
        if write:
            for n, struct in enumerate([self] + inter + [other]):
                struct.save_pdb(f'interpolated_structure_{n}.pdb')
            return
        else:
            return [self] + inter + [other]

    def save_pdb(self, outname='out.pdb'):
        with open(outname, 'w') as outfile:
            outfile.write(self._cryst_format.format(*self.box))
            for line in self.remarks:
                outfile.write(line.strip() + '\n')
            for atom in self.atoms:
                outfile.write(self._write_atom(atom))
            outfile.write('ENDMDL\n')
    
    def get_coords(self, subset=None):
        if subset:
            return [[a.x, a.y, a.z] for a in [self.atoms[q] for q in subset]]
        else:
            return [[a.x, a.y, a.z] for a in self.atoms]
    
    def set_coords(self, new_coords):
        assert len(new_coords) == len(self.atoms)
        for atom, coords in zip(self.atoms, new_coords):
            assert len(coords) == 3
            atom.x, atom.y, atom.z = coords
            

class Atom:
    def __init__(self, line):
        self.label = line[:6].strip()
        self.serial = int(line[6:11].strip())
        self.atomname = line[12:16].strip()
        self.altloc = line[16:17]
        self.resname = line[17:21].strip()
        self.chain = line[21:22]
        self.resnum = int(line[22:26].strip())
        self.insert = line[26:27]
        self.x, self.y, self.z = [float(line[30+8*a:30+8*(a+1)]) for a in range(3)]
        self.occ = float(line[54:60].strip())
        self.beta = float(line[60:66].strip())
        try:
            self.element = line[76:78].strip()
        except IndexError:
            name = self.atomname.strip('1234567890')
            if name in 'CHONSP':
                self.element = name[:1]
            else:
                self.element = name[:2]

    @classmethod
    def from_gro(cls, line):
        data = "ATOM  {:>5d} {:4s} {:4s} {:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00"
        resnum = int(line[:5].strip()) % 10000
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atomnum = int(line[15:20].strip())
        x, y, z = [float(line[20+8*i:20+8*(i+1)].strip())*10 for i in range(3)]
        return cls(data.format(atomnum, atomname, resname, resnum, x, y, z))

    @classmethod
    def from_top_entry(cls, entry):
        data = "ATOM  {:>5d} {:4s} {:4s} {:>4d}    {:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00"
        resnum = entry.resid
        resname = entry.resname
        atomname = entry.atomname
        atomnum = entry.num
        x, y, z = [0, 0, 0]
        return cls(data.format(atomnum, atomname, resname, resnum, x, y, z))

    @property
    def coords(self):
        return [self.x, self.y, self.z]
    
    def set_coords(self, coords):
        self.x, self.y, self.z = coords
        
    def __repr__(self):
        chain = self.chain if self.chain != " " else "unspecified"
        return "Atom {} in residue {}{} of chain {}".format(self.atomname, self.resname, self.resnum, chain)


class Traj(Pdb):
    def __init__(self, structfile=None, compressed_traj=None, top=None, array=None, altloc='A', **kwargs):
        super().__init__(structfile, top, altloc, **kwargs)
        self.sfname = structfile
        self.ctfname = compressed_traj
        self.array = array
        self.box_traj = []
        if self.sfname is None and self.ctfname is None:
            self.coords = [[[]]]
        elif self.sfname is None and self.sfname.endswith('.pdb'):
            self.coords = self.get_coords_from_pdb()
        elif self.sfname is None and self.sfname.endswith('.gro'):
            self.coords = self.get_coords_from_gro()
        if self.ctfname is not None and self.ctfname.endswith('.xtc') and self.sfname is not None:
            import mdtraj as md
            traj = md.load(self.ctfname, top=self.sfname)
            self.coords = [frame.xyz[0] for frame in traj]
        if self.array and self.coords:
            self.coords = self.array
        self.top = top if not isinstance(top, str) else gml.Top(top, **kwargs)
        if self.top and not self.top.pdb:
            self.top.pdb = self
        self.altloc = altloc
        self.nframes = len(self.coords)
        self._atom_format = "ATOM  {:>5d} {:4s}{:1s}{:4s}{:1s}{:>4d}{:1s}   " \
                            "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n"
        self._cryst_format = "CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           1\n"

    def __repr__(self):
        return "trajectory file {} with {} frames and {} atoms".format(self.fname, self.nframes, len(self.atoms))

    def get_coords_from_pdb(self):
        coords = []
        content = [line for line in open(self.sfname)]
        term_lines = [0] + [n for n, line in enumerate(content)
                            if line.startswith('END') or line.startswith('ENDMDL')]
        for i, j in zip(term_lines[:-1], term_lines[1:]):
            frame = ''.join(content[i:j])
            struct = Pdb.from_text(frame, self)
            coords.append(struct.get_coords())
            # TODO extract box size & append to self.box_traj
        return coords

    def get_coords_from_gro(self):
        coords = []
        content = [line for line in open(self.sfname)]
        natoms = int(content[1].strip())
        per_frame = natoms + 3
        term_lines = range(0, len(content)+1, per_frame)
        for i, j in zip(term_lines[:-1], term_lines[1:]):
            struct = Pdb.from_text('\n'.join(content[i:j]), self)
            coords.append(struct.get_coords())
            # TODO extract box size & append to self.box_traj
        return coords

    def __str__(self):
        if not self.box_traj:
            text = self._cryst_format.format(*self.box)
        else:
            raise NotImplementedError("Variable box size not yet implemented")
        for frame in range(len(self.coords)):
            text = text + 'MODEL     {:>4d}\n'.format(frame+1)
            self.set_coords(self.coords[frame])
            for atom in self.atoms:
                text = text + self._write_atom(atom)
            text = text + 'ENDMDL\n'
        return text

    def save_traj_as_pdb(self, filename='traj.pdb'):
        with open(filename, 'w') as outfile:
            outfile.write(str(self))