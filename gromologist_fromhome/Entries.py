

class Entry:
    """
    A generic class representing a single line in the topology.
    In an entry, the actual content and the comments are kept
    in two separate variables.
    """
    def __init__(self, content, subsection):
        self.subsection = subsection
        semicol_index = content.find(';')
        if semicol_index >= 0:
            self.content = content[:semicol_index].strip().split()
            self.comment = ' ' + content[semicol_index:]
        else:
            self.content = content.strip().split()
            self.comment = ''
    
    @staticmethod
    def float_fmt(flt, fields=11, dpmax=8):
        """
        When a float of unknown precision is read, we do not want
        to clip off significant digits, but neither do we want
        to add too many decimal places. This function calculates
        how many dps we need to keep not to lose precision when
        handling ff params (by default, we clip to 8).
        :param flt: float, the number to be formatted
        :param fields: how many fields do we need overall in the fmt specifier
        :param dpmax: default limit on the number of decimal places
        :return: str, format specifier
        """
        try:
            nf = len(str(flt).split('.')[1])
        except IndexError:
            nf = 3
        if nf > dpmax:
            nf = dpmax
        return "{:>" + str(fields) + "." + str(nf) + "f}"

    @staticmethod
    def infer_type(val):
        try:
            _ = int(val)
        except ValueError:
            try:
                _ = float(val)
            except ValueError:
                return str
            else:
                return float
        else:
            return int
        
    def __bool__(self):
        if not self.content and not self.comment:
            return False
        return True
    
    def __getitem__(self, item):
        return self.content[item]
    
    def __str__(self):
        """
        Fallback if no explicit formatting is implemented
        :return: str
        """
        return ' '.join(self.content) + ' ' + self.comment + '\n'
    
    
class EntryBonded(Entry):
    """
    This Entry subclass is intended for entries that correspond
    to bonded interaction (bonds, pairs, angles, dihedrals)
    between specific atoms in the topology
    """
    fstr_suff = {('bonds', '1'): (float, float),
                 ('bonds', '2'): (float, float),
                 ('angles', '1'): (float, float),
                 ('angles', '2'): (float, float),
                 ('angles', '5'): (float, float, float, float),
                 ('dihedrals', '9'): (float, float, int),
                 ('dihedrals', '4'): (float, float, int),
                 ('impropers', '4'): (float, float, int),
                 ('dihedrals', '1'): (float, float, int),
                 ('dihedrals', '3'): (float, float, float, float, float, float),
                 ('dihedrals', '2'): (float, float),
                 ('impropers', '2'): (float, float),
                 ('cmap', '1'): (float,),
                 ('position_restraints', '1'): (float, float, float),
                 ('settles', '1'): (float, float)}

    def __init__(self, content, subsection):
        super().__init__(content, subsection)
        if subsection.header == 'exclusions':
            self.atoms_per_entry = len(self.content) - 1
        else:
            self.atoms_per_entry = type(self.subsection).n_atoms[self.subsection.header]
        self.atom_numbers = tuple([int(x) for x in self.content[:self.atoms_per_entry]])
        self.atom_resid = tuple([int(x) for x in self.content[:self.atom_resid]]) #REPORT
        self.interaction_type = self.content[self.atoms_per_entry]
        try:
            self.params_per_entry = len(EntryBonded.fstr_suff[(subsection.header, str(self.interaction_type))])
        except KeyError:
            self.params_per_entry = 0
            # type assignment should only be performed when asked to, i.e. outside of constructor, with read_types
        self.types_state_a = None
        self.types_state_b = None
        self.atom_names = None
        self.params_state_a = []
        self.params_state_a_entry = []
        self.params_state_b = []
        self.params_state_b_entry = []
        self.fstr_mod = []
        if len(self.content) > self.atoms_per_entry + 1:
            try:
                self.parse_bonded_params(self.content[self.atoms_per_entry + 1:])
            except Exception as e:
                print("While trying to process line {}, subsection {}:".format(content, self.subsection.header))
                raise e
        self.fstring = " ".join("{:>5d}" for _ in range(self.atoms_per_entry)) + " {:>5s}"

    def _fstr_suff(self, query):
        if self.fstr_mod:
            return self.fstr_mod
        else:
            return EntryBonded.fstr_suff[query]

    def explicit_defines(self):
        if self.params_state_a and isinstance(self.params_state_a[0], str):
            try:
                self.params_state_a = self.subsection.section.top.defines[self.params_state_a[0]]
                self.fstr_mod = [self.infer_type(x) for x in self.params_state_a]
            except:
                pass

    def read_types(self):
        atoms_sub = self.subsection.section.get_subsection('atoms') # REPORT queste sono tutte le linee di atoms
        atoms_sub.get_dicts()
        num_to_type_a = atoms_sub.num_to_type
        num_to_type_b = atoms_sub.num_to_type_b
        num_to_name = atoms_sub.num_to_name
        self.types_state_a = tuple(num_to_type_a[num] for num in self.atom_numbers)
        types_state_b = tuple(num_to_type_b[num] for num in self.atom_numbers)
        self.types_state_b = types_state_b if types_state_b != self.types_state_a else None
        self.atom_names = tuple(num_to_name[num] for num in self.atom_numbers)
        #self.atom_resid = tuple #REPORT

    def parse_bonded_params(self, excess_params):
        try:
            _ = EntryBonded.fstr_suff[(self.subsection.header, self.interaction_type)]
        except KeyError:
            print((self.subsection.header, self.interaction_type))
            raise RuntimeError("Line '{}' contains unrecognized parameters".format(self.content))
        else:
            # if len(excess_params) == 1 and len(types) > 1:
            #     try:
            #         params = self.subsection.section.top.defines[excess_params[0]]
            #     except KeyError:
            #         raise RuntimeError("Cannot process: ", excess_params)
            #     else:
            #         self.params_state_a = [types[n](prm) for n, prm in enumerate(params)]
            try:
                _ = [float(x) for x in excess_params]
            except ValueError:
                # self.fstr_mod = list(EntryBonded.fstr_suff[(self.subsection.header, self.interaction_type)])
                for n in range(len(excess_params)):
                    try:
                        _ = float(excess_params[n])
                    except ValueError:
                        if not excess_params[n] in self.subsection.section.top.defines.keys():
                            print(f'undefined parameter {excess_params} was found, try adding define={"KEY": value}'
                                  f'when initializing Top()')
                self.fstr_mod = [self.infer_type(x) for x in excess_params]
            types = self._fstr_suff((self.subsection.header, self.interaction_type))
            if len(excess_params) == len(types):
                self.params_state_a = [types[n](prm) for n, prm in enumerate(excess_params[:len(types)])]
            elif len(excess_params) == 2 * len(types):
                self.params_state_a = [types[n](prm) for n, prm in enumerate(excess_params[:len(types)])]
                self.params_state_b = [types[n](prm) for n, prm in enumerate(excess_params[len(types):])]
            else:
                raise RuntimeError("Cannot process: ", excess_params)

    def __str__(self):
        fmt_suff = ""
        for params in [self.params_state_a, self.params_state_b]:
            for parm in params:
                if isinstance(parm, int):
                    fmt_suff = fmt_suff + "{:>6d} "
                elif isinstance(parm, float):
                    fmt_suff = fmt_suff + self.float_fmt(parm)
                elif isinstance(parm, str):
                    if len(parm) > 14:
                        fmt_suff = fmt_suff + "{{:>{}s}} ".format(len(parm)+2)
                    else:
                        fmt_suff = fmt_suff + "{:>15s} "
        fstring = self.fstring + fmt_suff
        return fstring.format(*self.atom_numbers, self.interaction_type, *self.params_state_a, *self.params_state_b) \
            + ' ' + self.comment

        
class EntryParam(Entry):
    """
    This Entry subclass represents a line containing force field
    parameters, e.g. bondtypes, angletypes, cmaptypes, pairtypes etc.
    that map a set of atom types to a set of FF-specific values
    """
    def __init__(self, content, subsection, processed=False):
        super().__init__(content, subsection)
        self.atoms_per_entry = type(self.subsection).n_atoms[self.subsection.header]
        self.types = tuple(self.content[:self.atoms_per_entry])
        if self.subsection.header == 'cmaptypes' and processed:
            self.modifiers = self.content[self.atoms_per_entry + 1:self.atoms_per_entry + 3]
            self.params = [float(x) for x in self.content[self.atoms_per_entry + 3:]]
            self.interaction_type = self.content[self.atoms_per_entry]
        elif self.subsection.header == 'cmaptypes' and not processed:
            self.modifiers = []
            self.params = self.content[self.atoms_per_entry + 1:]
            self.interaction_type = self.content[self.atoms_per_entry]
        elif self.subsection.header == 'defaults':
            self.modifiers = self.content
            self.params = []
            self.interaction_type = ''
        elif self.subsection.header == 'atomtypes':
            self.modifiers = self.content[self.atoms_per_entry:self.atoms_per_entry + 4]
            self.params = [float(x) for x in self.content[self.atoms_per_entry + 4:]]
            self.interaction_type = ''
        else:
            self.params = [float(x) for x in self.content[self.atoms_per_entry + 1:]]
            self.modifiers = []
            self.interaction_type = self.content[self.atoms_per_entry]
        if self.subsection.header == 'dihedraltypes':
            if any([self.infer_type(x) == float for x in self.types]):
                self.types = self.content[:2]
                self.interaction_type = self.content[2]
                self.params = [float(x) for x in self.content[3:]]
        if self.subsection.header == 'dihedraltypes' and self.interaction_type in ('9', '4', '1'):
            self.params[-1] = int(self.params[-1])
        self.identifier = self.subsection.header + '-' + '-'.join(self.types) + '-' + self.interaction_type
            
    def format(self):
        fmt = {('bondtypes', '1'): "{:>8s} {:>8s}{:>6s}{:>13.8f}{:>13.2f} ",
               ('angletypes', '5'): "{:>8s} {:>8s} {:>8s}{:>6s}{:>13.6f}{:>13.6f}{:>13.8f}{:>13.2f} ",
               ('angletypes', '1'): "{:>8s} {:>8s} {:>8s}{:>6s}{:>13.8f}{:>13.2f} ",
               ('dihedraltypes', '9'): "{:>8s} {:>8s} {:>8s} {:>8s}{:>6s}{:>13.6f}{:>13.6f}{:>6d} ",
               ('dihedraltypes', '4'): "{:>8s} {:>8s} {:>8s} {:>8s}{:>6s}{:>13.6f}{:>13.6f}{:>6d} ",
               ('dihedraltypes', '3'): "{:>8s} {:>8s} {:>8s} {:>8s}{:>6s}{:>13.6f}{:>13.6f}{:>13.6f}{:>13.6f}"
                                       "{:>13.6f}{:>13.6f} ",
               ('dihedraltypes', '2'): "{:>8s} {:>8s} {:>8s} {:>8s}{:>6s}{:>13.6f}{:>13.6f} ",
               ('dihedraltypes', '1'): "{:>8s} {:>8s}{:>6s}{:>13.6f}{:>13.6f}{:>6d} ",
               ('atomtypes', ''): "{:>6s}{}{:>6s}{:>13s}{:>9s}{:>3s}{:>16.12f}{:>9.5f} ",
               ('pairtypes', '1'): "{:>8s} {:>8s}{:>3s}{:>16.12f}{:>16.12f} ",
               ('nonbond_params', '1'): "{:>8s} {:>8s}{:>3s}{:>20.16f}{:>20.16f} ",
               ('implicit_genborn_params', ''): " {:8s}{:8.4f}{:8.4f}{:8.4f}{:8.4f}{:8.4f} "}
        if (self.subsection.header, self.interaction_type) in fmt.keys():
            return fmt[(self.subsection.header, self.interaction_type)]
        else:
            return None
    
    def match(self, ext_typelist, int_type):
        if not ext_typelist or len(ext_typelist) != len(self.types):
            return False
        if self.interaction_type == int_type:
            if ext_typelist[0] == self.types[0] or ext_typelist[1] == self.types[1] \
                    or ext_typelist[-1] == self.types[-1] or ext_typelist[-2] == self.types[-2]:
                if all(ext_typelist[i] == self.types[i] for i in range(len(self.types)) if self.types[i] != 'X'):
                    return True
            if ext_typelist[0] == self.types[-1] or ext_typelist[1] == self.types[-2] \
                    or ext_typelist[-1] == self.types[0] or ext_typelist[-2] == self.types[1]:
                if all(ext_typelist[i] == self.types[len(self.types)-i-1] for i in range(len(self.types))
                       if self.types[i] != 'X'):
                    return True
        return False
    
    def __repr__(self):
        if len(self.params) <= 4:
            return "Parameters entry with atomtypes {}, interaction type {} " \
                   "and parameters {}".format(self.types,
                                              self.interaction_type,
                                              ', '.join([str(x) for x in self.params]))
        else:
            return "Parameters entry with atomtypes {}, interaction type {} " \
                   "and parameters {}...".format(self.types,
                                                 self.interaction_type,
                                                 ', '.join([str(x) for x in self.params[:4]]))
        
    def __str__(self):
        """
        For cmaptypes, we rearrange lines to retrieve the matrix
        format lost during read-in; for other entry types, we
        delegate formatting to Subsection.fmt
        :return:
        """
        if self.subsection.header == 'cmaptypes':
            first = ((8 * "{} ")[:-1] + "\\\n").format(*self.types, self.interaction_type, *self.modifiers)
            npar = len(self.params)
            last = '\\\n'.join([((10 * "{} ")[:-1]).format(*self.params[10*n:10*(n+1)]) for n in range(int(npar/10))])
            if 10 * int(npar/10) != npar:
                last = last + '\\\n' + \
                       (((npar-10*int(npar/10)) * "{} ")[:-1]).format(*self.params[10*int(npar/10):]) + '\n'
            return first + last
        elif self.format():
            try:
                return self.format().format(*self.types, self.interaction_type, *self.modifiers, *self.params) +\
                       ' ' + self.comment + '\n'
            except IndexError:
                print((*self.types, self.interaction_type, *self.modifiers, *self.params))
        else:
            return super().__str__()

        
class EntryAtom(Entry):
    """
    This Entry subclass corresponds to atoms defined in
    the [ atoms ] section of each molecule
    """
    def __init__(self, content, subsection):
        super().__init__(content, subsection)
        try:
            self.num, self.type, self.resid, self.resname, self.atomname, _, self.charge, self.mass = self.content[:8]
        except ValueError:
            self.num, self.type, self.resid, self.resname, self.atomname, _, self.charge = self.content[:7]
            atomtypes = self.subsection.section.top.parameters.get_subsection('atomtypes')
            matching = [atype for atype in atomtypes if isinstance(atype, EntryParam) and atype.types[0] == self.type]
            try:
                self.mass = float(matching[0].modifiers[1])
            except IndexError:
                self.mass = 0
        self.num, self.resid = int(self.num), int(self.resid)
        self.charge, self.mass = float(self.charge), float(self.mass)
        if len(self.content) == 11:
            self.type_b, self.charge_b, self.mass_b = self.content[8], float(self.content[9]), float(self.content[10])
        else:
            self.type_b, self.charge_b, self.mass_b = None, None, None
        self.fstring = "{:>6d} {:>11s} {:>7d}{:>7s}{:>7s}{:>7d}"
    
    def __str__(self):
        fstring = self.fstring + self.float_fmt(self.charge) + self.float_fmt(self.mass) + '   '
        if self.type_b:
            alch_fstring = "{:>11s}" + self.float_fmt(self.charge_b) + self.float_fmt(self.mass_b)
            fstring += alch_fstring
            return fstring.format(self.num, self.type, self.resid, self.resname, self.atomname, self.num,
                                  self.charge, self.mass, self.type_b, self.charge_b, self.mass_b) + self.comment
        else:
            return fstring.format(self.num, self.type, self.resid, self.resname, self.atomname, self.num,
                                  self.charge, self.mass) + ' ' + self.comment
