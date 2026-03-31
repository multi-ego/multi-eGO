import sys
import unittest
import shutil
import subprocess
import os
import yaml

TEST_ROOT = os.path.dirname(os.path.abspath(__file__))
MEGO_ROOT = os.path.abspath(os.path.join(TEST_ROOT, os.pardir))


def _system_from_config(config_path):
    """
    Read the system name from a multi-eGO YAML config file.

    The config is a YAML sequence whose items are either plain strings
    (flags like ``no_header``) or single-key mappings (like ``system: gpref``).
    """
    with open(config_path) as f:
        data = yaml.safe_load(f)
    for item in data:
        if isinstance(item, dict) and "system" in item:
            return item["system"]
    raise ValueError(f"No 'system' key found in {config_path}")


def read_infile(path):
    """
    Reads a test-case input file and returns the command-line argument lists
    and corresponding system names.

    Lines starting with ``#`` are skipped.  Each active line is split into
    tokens and ``TEST_ROOT`` is expanded to the absolute tests directory.

    The system name is resolved in order:
    1. From a ``--system <name>`` token on the line itself.
    2. From the ``system:`` key inside the file given by ``--config <path>``.

    Lines with neither are skipped.

    Parameters
    ----------
    path : str
        Path to the test-case text file.

    Returns
    -------
    input_list : list of list of str
        Each entry is the argument list for one test run.
    test_systems : list of str
        The system name for each test run.
    """
    input_list = []
    test_systems = []
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = [arg.replace("TEST_ROOT", TEST_ROOT) for arg in line.split()]

            if "--system" in parts:
                system_name = parts[parts.index("--system") + 1]
            elif "--config" in parts:
                config_path = parts[parts.index("--config") + 1]
                system_name = _system_from_config(config_path)
            else:
                continue

            input_list.append(parts)
            test_systems.append(system_name)
    return input_list, test_systems


def read_outfile(path):
    """
    Reads a multi-eGO output file, stripping inline comments (everything after
    a ``;`` on each line).

    Parameters
    ----------
    path : str
        Path to the multi-eGO output file (ffnonbonded or topology).

    Returns
    -------
    out_string : str
        The file contents with comments removed.
    """
    out_string = ""
    with open(path) as f:
        for line in f:
            out_string += line.split(";")[0] if ";" in line else line
    return out_string


def prep_system_data(name, index):
    """
    Reads the reference and freshly generated output files for one test case.

    Parameters
    ----------
    name : str
        System name.
    index : int
        Test-case index (1-based).

    Returns
    -------
    topol_ref, topol_test, ffnonbonded_ref, ffnonbonded_test : str
        File contents ready for comparison.
    """
    topol_ref = read_outfile(f"{TEST_ROOT}/test_outputs/{name}/case_{index}/topol_mego.top")
    topol_test = read_outfile(f"{MEGO_ROOT}/outputs/{name}/case_{index}/topol_mego.top")
    ffnonbonded_ref = read_outfile(f"{TEST_ROOT}/test_outputs/{name}/case_{index}/ffnonbonded.itp")
    ffnonbonded_test = read_outfile(f"{MEGO_ROOT}/outputs/{name}/case_{index}/ffnonbonded.itp")
    return topol_ref, topol_test, ffnonbonded_ref, ffnonbonded_test


def create_test_cases(command, system_name):
    """
    Creates a test method for a single multi-eGO command-line invocation.

    The generated test method runs the command, then immediately compares
    the output against the reference files.  Each test is fully self-contained:
    it cleans its own output directory, executes multi-eGO, and performs the
    diff — so a failure in one case does not affect any other.

    Parameters
    ----------
    command : list of str
        The argument list for one test run.
    system_name : str
        The system name used to locate output files and name the test.

    Returns
    -------
    function_name : str
        Name in the format ``test_<system>/case_<n>``.
    function_template : callable
        A method for use with ``unittest.TestCase``.
    """
    prefix = f"test_{system_name}/case_"
    idx = 1
    while hasattr(TestOutputs, f"{prefix}{idx}"):
        idx += 1
    function_name = f"{prefix}{idx}"

    # Resolve the full command once at definition time.
    cmd = command
    if "--system" in cmd and "--inputs_dir" not in cmd:
        cmd = cmd + ["--inputs_dir", f"{TEST_ROOT}/test_inputs"]

    def function_template(self):
        # Remove only this case's output directory so previous cases are
        # not disturbed and the run starts from a clean state.
        output_path = f"{MEGO_ROOT}/outputs/{system_name}/case_{idx}"
        if os.path.exists(output_path):
            shutil.rmtree(output_path)

        ret = subprocess.call([sys.executable, f"{MEGO_ROOT}/multiego.py", *cmd])
        self.assertEqual(ret, 0, f"{system_name} :: {function_name} exited with non-zero error code")

        topol_ref, topol_test, ffnonbonded_ref, ffnonbonded_test = prep_system_data(system_name, idx)
        self.assertEqual(topol_ref, topol_test, f"{system_name} :: {function_name} topology not equal")
        self.assertEqual(ffnonbonded_ref, ffnonbonded_test, f"{system_name} :: {function_name} nonbonded not equal")

    return function_name, function_template


# ---------------------------------------------------------------------------
# Build test methods at *module* level so both pytest and unittest discover them
# ---------------------------------------------------------------------------

test_commands, test_systems = read_infile(f"{TEST_ROOT}/test_cases.txt")


class TestOutputs(unittest.TestCase):
    pass


for _cmd, _sys in zip(test_commands, test_systems):
    _name, _method = create_test_cases(_cmd, _sys)
    setattr(TestOutputs, _name, _method)


if __name__ == "__main__":
    unittest.main()
