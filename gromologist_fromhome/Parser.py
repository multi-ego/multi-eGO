from .Entries import EntryAtom


class SelectionParser:
    def __init__(self, master):
        self.master = master
        try:
            self.nat = len(self.master.atoms)
        except AttributeError:
            try:
                self.nat = self.master.natoms
            except AttributeError:
                raise TypeError("Can only parse selections with PDB or topology data")
        
    @staticmethod
    def ispdb(obj):
        try:
            _ = obj._cryst_format
        except AttributeError:
            return False
        else:
            return True
    
    def __call__(self, selection_string):
        protein_selection = "resname ALA ACE CYS ASP ASPP GLU GLUP PHE GLY HIS HID HIE HSD HSE ILE LYS LEU MET " \
                            "NME NMA ASN PRO GLN ARG SER THR VAL TRP"
        dna_selection = "resname DA DG DC DT DA5 DG5 DC5 DT5 DA3 DG3 DC3 DT3"
        rna_selection = "resname RA RG RC RT RA5 RG5 RC5 RT5 RA3 RG3 RC3 RT3"
        solvent_selection = "resname HOH TIP3 SOL K CL NA"
        selection_string = selection_string.replace('solvent', solvent_selection)
        selection_string = selection_string.replace('water', 'resname HOH TIP3 SOL')
        selection_string = selection_string.replace('protein', protein_selection)
        selection_string = selection_string.replace('nucleic', 'dna or rna')
        selection_string = selection_string.replace('dna', dna_selection)
        selection_string = selection_string.replace('rna', rna_selection)
        return sorted(list(self._select_set_atoms(selection_string)))
    
    def _select_set_atoms(self, selection_string, ):
        """
        :param selection_string: str
        :return:
        """
        assert isinstance(selection_string, str)
        parenth_ranges, operators = self._parse_sel_string(selection_string)
        if not parenth_ranges and not operators:
            if selection_string.strip().startswith("not "):
                return self._find_atoms(selection_string.lstrip()[4:], rev=True)
            else:
                return self._find_atoms(selection_string)
        elif parenth_ranges and not operators:
            if selection_string.strip().startswith("not "):
                set_all = {n for n in range(self.nat)}
                return set_all.difference(self._select_set_atoms(selection_string.strip()[4:].strip('()')))
            else:
                return self._select_set_atoms(selection_string.strip().strip('()'))
        else:
            first_op = selection_string[operators[0][0]:operators[0][1]].strip()
            if first_op == "and":
                return self._select_set_atoms(selection_string[:operators[0][0]]) \
                    .intersection(self._select_set_atoms(selection_string[operators[0][1]:]))
            elif first_op == "or":
                return self._select_set_atoms(selection_string[:operators[0][0]]) \
                    .union(self._select_set_atoms(selection_string[operators[0][1]:]))
    
    @staticmethod
    def _parse_sel_string(selection_string):
        parenth_ranges = []
        operators = []
        opened_parenth = 0
        beginning = 0
        for nc, char in enumerate(selection_string):
            if char == '(':
                opened_parenth += 1
                if beginning == 0:
                    beginning = nc
            elif char == ')':
                opened_parenth -= 1
                end = nc
                if opened_parenth == 0:
                    parenth_ranges.append((beginning, end))
                    beginning = 0
            if opened_parenth < 0:
                raise ValueError("Improper use of parentheses in selection string {}".format(selection_string))
            if selection_string[nc:nc + 5] == " and " and opened_parenth == 0:
                operators.append((nc, nc + 5))
            if selection_string[nc:nc + 4] == " or " and opened_parenth == 0:
                operators.append((nc, nc + 4))
        if opened_parenth != 0:
            raise ValueError("Improper use of parentheses in selection string {}".format(selection_string))
        return parenth_ranges, operators
    
    def _find_atoms(self, sel_string, rev=False):
        chosen = []
        same = False
        within = False
        if ' '.join(sel_string.split()[:3]) == "same residue as":
            sel_string = ' '.join(sel_string.split()[3:])
            same = True
        if sel_string.split()[0] == 'within' and sel_string.split()[2] == 'of':
            within = float(sel_string.split()[1])
            sel_string = ' '.join(sel_string.split()[3:])
        keyw = sel_string.split()[0]
        if self.ispdb(self.master):
            matchings = {"name": "atomname", "resid": "resnum", "resnum": "resnum", "element": "element",
                         "chain": "chain", "resname": "resname", "serial": "serial"}
        else:
            matchings = {"name": "atomname", "resid": "resid", "resnum": "resid", "mass": "mass",
                         "resname": "resname", "serial": "num", "type": "type"}
        try:
            vals = {int(x) for x in sel_string.split()[1:]}
        except ValueError:
            if " to " in sel_string:
                beg = int(sel_string.split()[1])
                end = int(sel_string.split()[3])
                vals = set(range(beg, end + 1))
            else:
                vals = set(sel_string.split()[1:])
        try:
            atomlist = self.master.atoms
        except AttributeError:
            atomlist = [e for e in self.master.get_subsection('atoms').entries if isinstance(e, EntryAtom)]
        for n, a in enumerate(atomlist):
            if not rev:
                if a.__getattribute__(matchings[keyw]) in vals:
                    chosen.append(n)
            else:
                if a.__getattribute__(matchings[keyw]) not in vals:
                    chosen.append(n)
        
        if same:
            chosen = self.master.same_residue_as(chosen)  # TODO implement for SectMol
        if within and self.ispdb(self.master):
            chosen = self.master.within(chosen, within)
        if within and not self.ispdb(self.master):
            raise ValueError("the within keyword only works for structural data, not topology")
        return set(chosen)
