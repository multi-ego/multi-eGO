import re
import decimal
import sys

INPUT_FILE  = sys.argv[1]
OUTPUT_FILE = sys.argv[2]
EPSILON     = float(sys.argv[3])

atom_type_dict = {}
file_start = ""
file_end = ""
with open(INPUT_FILE, 'r') as f:
    atomtypes_found = False
    pairtypes_found = False

    for line in f.readlines():
        if ';' in line: line = line.split(';')[0]
        if '[ pairtypes ]' in line: pairtypes_found = True
        # if re.compile(r'\[\s*pairtypes\s*\]', line).search(): pairtypes_found = True
        # if re.compile(r'\[\s*atomtypes\s*\]').matched(): atomtypes_found = True

        if pairtypes_found:
            file_end += line
        else:
            line_data = re.split('\s+', line)
            if len(line_data) < 8: continue
            atom_type_dict[line_data[1]] = {
                'at.num' : line_data[2],
                'mass' : line_data[3],
                'charge' : line_data[4],
                'ptype' : line_data[5],
                'c6' : line_data[6],
                'c12' : line_data[7]
            }

for key in atom_type_dict.keys():
    c12 = float(atom_type_dict[key]['c12'])
    sigma = ( c12 / 4 / EPSILON ) ** ( 1 / 12 )
    atom_type_dict[key]['c6'] = '%.6e' % decimal.Decimal(str(4 * EPSILON * sigma ** 6))

with open(OUTPUT_FILE, 'w') as f:
    f.write('[ atomtypes ]\n')
    f.write(';name  at.num   mass      charge  ptype       c6           c12\n')
    for atdk in atom_type_dict.keys():
        atom_type = atom_type_dict[atdk]
        f.write(f'{atdk}')
        for k in atom_type.keys():
            f.write(f'\t{atom_type[k]}')
        f.write('\n')
    f.write('\n')
    f.write(file_end)
