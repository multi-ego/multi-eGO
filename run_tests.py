import unittest
import subprocess
import shutil
import os

def read_infile(path):
    '''
    Reads a test-case input file and parses the system name the multi-eGO
    command line parameters.

    Parameters
    ----------
    path : str
        The path to the test-case text file
    
    Returns
    -------
    input_list : list of list
        A list of the commands split at each whitespace
    test_systems : list
        A list containing the system names
    '''
    input_list = []
    test_systems = []
    with open(path, 'r') as f:
        for line in f.readlines():
            line = line.replace('\n', '').split(' ')
            system_index = line.index('--system') + 1

            input_list.append(line)
            test_systems.append(line[system_index])

    return input_list, test_systems


def read_outfile(path):
    '''
    Reads mutli-eGO putputfiles ignoring the comments

    Parameters
    ----------
    path : string
        A path to the multi-eGO output file be it ffnonbonded or topology

    Returns
    -------
    out_string : str
        The file contents 
    '''
    out_string = ""
    with open(path, 'r') as f:
        for line in f.readlines():
            line = line.split(';')[0] if ';' in line else line
            out_string += line

    return out_string


def prep_system_data(name, egos):
    '''
    Prepares system data to be compared and tested by reading all necessary files.

    Parameters
    ----------
    name : str
        Represents the system name (--system parameter in multi-eGO)
    egos : str or list
        Can take two types of values:
         - a string in the case of random coil
         - a list in case of production
        When egos is a list the list contains the two epsilon values for intra and inter.
    
    Returns
    -------
    topol_ref : str
        The contents of the reference topology which needs to be matched
    topol_test : str
        The contents of the newly created topology which needs to match topol_ref
    ffnonbonded_ref : str
        The contents of the reference ffnonbonded which needs to be matched
    ffnonbonded_test : str
        The contents of the newly created ffnonbonded which needs to match ffnonbonded_ref
    '''
    if egos == 'rc': out_egos = egos
    else: out_egos = f'production_e{egos[0]}_{egos[1]}'
    topol_ref = read_outfile(f'test/test_outputs/{name}_{out_egos}/topol_GRETA.top')
    topol_test = read_outfile(f'outputs/{name}_{out_egos}/topol_GRETA.top')
    ffnonbonded_ref = read_outfile(f'test/test_outputs/{name}_{out_egos}/ffnonbonded.itp')
    ffnonbonded_test = read_outfile(f'outputs/{name}_{out_egos}/ffnonbonded.itp')
    return topol_ref, topol_test, ffnonbonded_ref, ffnonbonded_test


def create_test_cases(test_case):
    '''
    Creates a test function based on the parameters. The metafunctions can be used with TestOutputs
    to automatically generate test cases.

    Parameters
    ----------
    test_case : list
        Contains the multi-eGO command line flags followed by arguments.

    Returns
    -------
    function_name : str
        The name of the function in format 'test_<system_name>_<egos>'
    function_template : function(self)
        A function taking only self as a parameter intended to be used as a unittest test case.
    '''
    # get system name
    system_index = test_case.index('--system') + 1
    system_name = test_case[system_index]

    # get egos type
    egos_index = test_case.index('--egos') + 1
    egos = test_case[egos_index]

    if egos == 'rc':
        name_suffix = 'rc'
    else:
        # figure out the name suffix and epsilons
        name_suffix = 'production'
        intra_epsilon_index = test_case.index('--epsilon') + 1
        intra_epsilon = test_case[intra_epsilon_index]
        if '--inter_epsilon' not in test_case:
            inter_epsilon = intra_epsilon
        else:
            inter_epsilon_index = egos_index if '--inter_epsilon' not in test_case else test_case.index('--inter_epsilon') + 1
            inter_epsilon = test_case[inter_epsilon_index]

    function_name = f'test_{system_name}_{name_suffix}'
    system_egos = 'rc' if name_suffix == 'rc' else [intra_epsilon, inter_epsilon]

    def function_template(self):
        name = system_name
        egos = system_egos

        topol_ref, topol_test, ffnonbonded_ref, ffnonbonded_test = prep_system_data(name=name, egos=egos)
        self.assertEqual(topol_ref, topol_test, f"{name} :: {egos} topology not equal")
        self.assertEqual(ffnonbonded_ref, ffnonbonded_test, f"{name} :: {egos} nonbonded not equal")


    return function_name, function_template


class TestOutputs(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        test_commands, test_systems = read_infile('test/test_cases.txt')
        input_directories = [ f'inputs/{system}' for system in test_systems ]
        input_gprotein = 'inputs/gpref'
        input_abeta = 'inputs/abetaref'
        input_ttr = 'inputs/ttrref'
        for system in test_systems:
            inputs_path = f'inputs/{system}'
            if os.path.exists(inputs_path): shutil.rmtree(inputs_path)
            shutil.copytree(f'test/test_inputs/{system}', inputs_path)

        error_codes = [ subprocess.call(["python", "multiego.py", *command]) for command in test_commands ]
        for e in error_codes: assert e == 0, "Test setup exited with non-zero error code"


if __name__ == '__main__':
    test_commands, test_systems = read_infile('test/test_cases.txt')
    for command in test_commands:
            function_name, new_method = create_test_cases(command)
            setattr(TestOutputs, function_name, new_method)

    unittest.main()