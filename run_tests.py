import unittest
import subprocess
import shutil
import os

def read_outfile(path):
    out_string = ""
    with open(path, 'r') as f:
        for line in f.readlines():
            line = line.split(';')[0] if ';' in line else line
            out_string += line

    return out_string

def test_system(name, egos):
    out_egos = egos if egos == 'rc' else f'production_e{egos}_{egos}'
    topol_ref = read_outfile(f'test/test_outputs/{name}_{out_egos}/topol_GRETA.top')
    topol_test = read_outfile(f'outputs/{name}_{out_egos}/topol_GRETA.top')
    ffnonbonded_ref = read_outfile(f'test/test_outputs/{name}_{out_egos}/ffnonbonded.itp')
    ffnonbonded_test = read_outfile(f'outputs/{name}_{out_egos}/ffnonbonded.itp')
    return topol_ref, topol_test, ffnonbonded_ref, ffnonbonded_test


class TestOutputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        input_gprotein = 'inputs/gpref'
        input_abeta = 'inputs/abetaref'
        if not os.path.exists(input_gprotein): shutil.copytree('test/test_inputs/gpref', input_gprotein)
        if not os.path.exists(input_abeta): shutil.copytree('test/test_inputs/abetaref', input_abeta)
        subprocess.call(["python", "multiego.py", "--system=gpref", "--egos=rc"])
        subprocess.call(["python", "multiego.py", "--system=gpref", "--egos=production", "--epsilon=0.35", "--train_from=md_ensemble"])
        subprocess.call(["python", "multiego.py", "--system=abetaref", "--egos=production", "--epsilon=0.35", "--train_from=native_MD"])

    def test_gprotein_production(self):
        name = 'gpref'
        egos = 0.35

        topol_ref, topol_test, ffnonbonded_ref, ffnonbonded_test = test_system(name=name, egos=egos)
        self.assertEqual(topol_ref, topol_test, f"{name} :: {egos} topology not equal")
        self.assertEqual(ffnonbonded_ref, ffnonbonded_test, f"{name} :: {egos} nonbonded not equal")

    def test_gprotein_rc(self):
        name = 'gpref'
        egos = 'rc'

        topol_ref, topol_test, ffnonbonded_ref, ffnonbonded_test = test_system(name=name, egos=egos)
        self.assertEqual(topol_ref, topol_test, f"{name} :: {egos} topology not equal")
        self.assertEqual(ffnonbonded_ref, ffnonbonded_test, f"{name} :: {egos} nonbonded not equal")

    def test_abeta_production(self):
        name = 'abetaref'
        egos = 0.35
        
        topol_ref, topol_test, ffnonbonded_ref, ffnonbonded_test = test_system(name=name, egos=egos)
        self.assertEqual(topol_ref, topol_test, f"{name} :: {egos} topology not equal")
        self.assertEqual(ffnonbonded_ref, ffnonbonded_test, f"{name} :: {egos} nonbonded not equal")


if __name__ == '__main__':
    unittest.main()