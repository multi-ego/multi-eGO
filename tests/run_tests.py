import unittest
import subprocess
import shutil
import os

TEST_ROOT = os.path.dirname(os.path.abspath(__file__))
MEGO_ROOT = os.path.abspath(os.path.join(TEST_ROOT, os.pardir))
# sys.path.append(MEGO_ROOT)


def read_infile(path):
    """
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
    """
    input_list = []
    test_systems = []
    with open(path, "r") as f:
        for line in f.readlines():
            if line[0] == "#":
                continue
            line = line.replace("\n", "").split(" ")
            system_index = line.index("--system") + 1

            input_list.append(line)
            test_systems.append(line[system_index])

    return input_list, test_systems


def read_outfile(path):
    """
    Reads multi-eGO output files ignoring the comments

    Parameters
    ----------
    path : string
        A path to the multi-eGO output file be it ffnonbonded or topology

    Returns
    -------
    out_string : str
        The file contents
    """
    out_string = ""
    with open(path, "r") as f:
        for line in f.readlines():
            line = line.split(";")[0] if ";" in line else line
            out_string += line

    return out_string


def prep_system_data(name, index):
    """
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
    """
    topol_ref = read_outfile(f"{TEST_ROOT}/test_outputs/{name}/case_{index}/topol_mego.top")
    topol_test = read_outfile(f"{MEGO_ROOT}/outputs/{name}/case_{index}/topol_mego.top")
    ffnonbonded_ref = read_outfile(f"{TEST_ROOT}/test_outputs/{name}/case_{index}/ffnonbonded.itp")
    ffnonbonded_test = read_outfile(f"{MEGO_ROOT}/outputs/{name}/case_{index}/ffnonbonded.itp")
    return topol_ref, topol_test, ffnonbonded_ref, ffnonbonded_test


def create_test_cases(test_case):
    """
    Creates a test function based on the parameters. The metafunctions can be used with TestOutputs
    to automatically generate test cases.

    Parameters
    ----------
    test_case : list
        Contains the multi-eGO command line flags followed by arguments.

    Returns
    -------
    function_name : str
        The name of the function in format 'test_<system_name>/case_<index>'
    function_template : function(self)
        A function taking only self as a parameter intended to be used as a unittest test case.
    """
    # get system name
    system_index = test_case.index("--system") + 1
    system_name = test_case[system_index]

    function_name_prefix = f"test_{system_name}/case_"

    # check if self alread has the function function_name
    idx = 1
    while hasattr(TestOutputs, f"{function_name_prefix}{idx}"):
        idx += 1
    function_name = f"{function_name_prefix}{idx}"

    def function_template(self):
        name = system_name

        topol_ref, topol_test, ffnonbonded_ref, ffnonbonded_test = prep_system_data(name=name, index=idx)
        self.assertEqual(topol_ref, topol_test, f"{name} :: {function_name} topology not equal")
        self.assertEqual(ffnonbonded_ref, ffnonbonded_test, f"{name} :: {function_name} nonbonded not equal")

    return function_name, function_template


class TestOutputs(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        test_commands, test_systems = read_infile(f"{TEST_ROOT}/test_cases.txt")
        # remake test_commands with only what comes before # if present
        test_commands = [[arg for arg in command if arg != "#"] for command in test_commands]
        # replace instances of TEST_ROOT in the commands with the actual path if TEST_ROOT is present
        test_commands = [[arg.replace("TEST_ROOT", TEST_ROOT) for arg in command] for command in test_commands]
        for system in test_systems:
            inputs_path = f"{MEGO_ROOT}/inputs/{system}"
            outputs_path = f"{MEGO_ROOT}/outputs/{system}"
            if os.path.exists(inputs_path):
                shutil.rmtree(inputs_path)
            if os.path.exists(outputs_path):
                shutil.rmtree(outputs_path)
            shutil.copytree(f"{TEST_ROOT}/test_inputs/{system}", inputs_path)

        error_codes = [subprocess.call(["python", f"{MEGO_ROOT}/multiego.py", *command]) for command in test_commands]
        for e in error_codes:
            assert e == 0, "Test setup exited with non-zero error code"


if __name__ == "__main__":
    test_commands, test_systems = read_infile(f"{TEST_ROOT}/test_cases.txt")
    for command in test_commands:
        function_name, new_method = create_test_cases(command)
        setattr(TestOutputs, function_name, new_method)

    unittest.main()
