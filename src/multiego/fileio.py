import importlib.metadata
import numpy as np
import pandas as pd
import glob
import os
import git
import time
import re
import json

from . import _term


def strip_gz_h5_suffix(filename):
    """
    Remove the '.gz' suffix from a filename if it ends with '.gz'.

    This function checks if the provided filename ends with the '.gz' suffix.
    If it does, the suffix is stripped (removed), and the modified filename is returned.
    If the filename does not end with '.gz', it is returned unchanged.

    Parameters:
    - filename (str): The filename to process.

    Returns:
    - str: The filename without the '.gz' suffix, if it was originally present.
           Otherwise, the original filename is returned.
    """
    if filename.endswith(".gz"):
        return filename[:-3]

    if filename.endswith(".h5"):
        return filename[:-3]

    return filename


def check_matrix_compatibility(input_path):
    """
    Check for matrix file compatibility by identifying any overlapping files
    that exist in both uncompressed ('.ndx') and compressed ('.ndx.gz') formats
    within a specified directory.

    This function searches for files with the patterns 'int??mat_?_?.ndx' and
    'int??mat_?_?.ndx.gz' in the provided input directory. It then checks for any
    common files that appear in both uncompressed and compressed forms.
    If such overlaps are found, a ValueError is raised indicating an issue
    with file compatibility, highlighting the names of the conflicting files.

    Parameters:
    - input_path (str): The path to the directory where the files will be checked.

    Raises:
    - ValueError: If files with both '.ndx' and '.ndx.gz' versions are found.

    Returns:
    - None: The function returns None but raises an error if incompatible files are found.
    """
    matrix_paths = glob.glob(f"{input_path}.ndx")
    matrix_paths_gz = glob.glob(f"{input_path}.ndx.gz")
    matrix_paths_h5 = glob.glob(f"{input_path}.ndx.h5")
    stripped_matrix_paths_gz_set = set(map(strip_gz_h5_suffix, matrix_paths_gz))
    stripped_matrix_paths_h5_set = set(map(strip_gz_h5_suffix, matrix_paths_h5))
    matrix_paths_set = set(matrix_paths)

    # Find intersection of the two sets
    common_files = matrix_paths_set.intersection(stripped_matrix_paths_gz_set)

    # Check if there are any common elements and raise an error if there are
    if common_files:
        raise ValueError(f"Error: Some files have both text and gz versions: {common_files}")

    # Find intersection of the two sets
    common_files = matrix_paths_set.intersection(stripped_matrix_paths_h5_set)

    # Check if there are any common elements and raise an error if there are
    if common_files:
        raise ValueError(f"Error: Some files have both text and hdf5 versions: {common_files}")

    # Find intersection of the two sets
    common_files = stripped_matrix_paths_gz_set.intersection(stripped_matrix_paths_h5_set)

    # Check if there are any common elements and raise an error if there are
    if common_files:
        raise ValueError(f"Error: Some files have both gz and hdf5 versions: {common_files}")


def check_mat_name(mat_name, ref):
    # Check name of matrix is either intramat_X_X or intermat_X_Y
    pattern = r"^(intra|inter)mat_\d+_\d+$"
    if not re.match(pattern, mat_name):
        raise ValueError(
            f"Wrong input matrix format {mat_name} in reference {ref}. \nContact matrix file(s) must be named as intramat_X_X.ndx(.gz/.h5) or intermat_X_Y.ndx(.gz/.h5)"
        )


def check_mat_extension(extension, ref):
    # pattern = r"^ndx(\.gz)?(\.h5)?$"
    # checks the extension of matrix name is either none or, .ndx(.gz/.h5)
    pattern = r"^(|\.ndx(\.gz|\.h5)?)$"
    if not re.match(pattern, extension):
        raise ValueError(
            f"Wrong input matrix format extension: {extension} in reference {ref}. \nContact matrix file(s) must be named as intramat_X_X.ndx(.gz/.h5) or intermat_X_Y.ndx(.gz/.h5)"
        )


def check_matrix_format(args):
    """
    Check the format of matrix files across multiple directories to ensure consistency
    and compatibility. This function specifically checks that there are no overlapping files
    in uncompressed ('.ndx') and compressed ('.ndx.gz') formats within the reference directory,
    training simulations, and check simulations directories.

    This function iterates through directories specified in the provided 'args' object. It starts
    by checking the reference directory for matrix file compatibility, then proceeds to check
    each training and checking simulation directory for similar issues.

    Parameters:
    - args (Namespace): An argparse.Namespace or similar object containing configuration settings.
      Expected keys include:
      - root_dir (str): The root directory under which all other directories are organized.
      - system (str): The specific system folder under 'root_dir' to use.
      - reference (str): The subdirectory within 'system' that contains the reference files.
      - train (list of str): A list of subdirectories within 'system' for training simulations.
      - check (list of str): A list of subdirectories within 'system' for checking simulations.

    Raises:
    - ValueError: If files with both '.ndx' and '.ndx.gz' versions are found in any checked directory.

    Returns:
    - None: The function returns None but raises an error if incompatible files are found in any directory.
    """
    for ref in args.input_refs:
        mat_appo = ref["matrix"].split(".")
        if len(mat_appo) > 1:
            extension = ref["matrix"].split(mat_appo[0])[1]
        else:
            extension = ""

        # Set name of matrix without extension
        ref["matrix"] = mat_appo[0]

        matrix_ref_path = f"{args.inputs_dir}/{args.system}/{ref['reference']}/{ref['matrix']}"
        check_mat_name(ref["matrix"], ref)
        check_mat_extension(extension, ref)
        check_matrix_compatibility(matrix_ref_path)
        for train in ref["train"]:
            matrix_train_path = f"{args.inputs_dir}/{args.system}/{train}/{ref['matrix']}"
            check_matrix_compatibility(matrix_train_path)


def read_symmetry_file(path):
    """
    Reads the symmetry file and returns a dictionary of the symmetry parameters.

        Parameters
        ----------
        path : str
            The path to the symmetry file

        Returns
        -------
        symmetry : dict
            The symmetry parameters as a dictionary
    """
    with open(path, "r") as file:
        lines = file.readlines()
    symmetry = parse_symmetry_list(lines)
    return symmetry


def parse_symmetry_list(symmetry_list):
    """
    Parse a symmetry string into a list of tuples.

    This function takes a string containing symmetry information and parses it into a list of tuples.
    Each tuple contains the symmetry information for a single interaction. The input string is expected
    to be formatted as a series of space-separated values, with each line representing a separate interaction.
    The values in each line are expected to be in the following order:
    - Name of the residue or molecule type
    - Name of the first atom
    - Name of the second atom

    Parameters
    ----------
    - symmetry_string : str
        A string containing symmetry information for interactions.

    Returns
    -------
    symmetry : list of tuple
        A list of tuples, with each tuple containing the symmetry information for a single interaction.
    """
    symmetry = []

    for line in symmetry_list:
        if "#" in line:
            line = line[: line.index("#")]
        line = line.replace("\n", "")
        line = line.strip()
        if not line:
            continue
        line = line.split(" ")
        line = [x for x in line if x]
        if len(line) < 3:
            continue

        symmetry.append(line)

    return symmetry


def read_molecular_contacts(path, ensemble_molecules_idx_sbtype_dictionary, simulation, h5=False):
    """
    Reads intra-/intermat files to determine molecular contact statistics.
    """
    _term.path(path)
    # Define column names and data types directly during read
    col_names = ["molecule_name_ai", "ai", "molecule_name_aj", "aj", "distance", "probability", "cutoff", "learned"]
    col_types = {
        "molecule_name_ai": "category",
        "ai": "category",
        "molecule_name_aj": "category",
        "aj": "category",
        "distance": np.float64,
        "probability": np.float64,
        "cutoff": np.float64,
        "learned": "Int64",  # Allows for integer with NaNs, which can be cast later
    }

    if not h5:
        contact_matrix = pd.read_csv(path, header=None, sep=r"\s+", names=col_names, dtype=col_types)
        contact_matrix["learned"] = contact_matrix["learned"].fillna(1).astype(bool)
    else:
        contact_matrix = pd.read_hdf(path, key="data", dtype=col_types)

    contact_matrix["learned"] = contact_matrix["learned"].astype(bool)

    # Validation checks using `query` for more efficient conditional filtering
    if contact_matrix.query("probability < 0 or probability > 1").shape[0] > 0:
        raise ValueError("ERROR: Probabilities should be between 0 and 1.")

    if contact_matrix.query("distance < 0 or distance > cutoff").shape[0] > 0:
        raise ValueError("ERROR: Distances should be between 0 and cutoff.")

    if contact_matrix.query("cutoff < 0").shape[0] > 0:
        raise ValueError("ERROR: Cutoff values cannot be negative.")

    # Check for NaN or infinite values in critical columns
    if contact_matrix[["probability", "distance", "cutoff"]].isnull().any().any():
        raise ValueError("ERROR: The matrix contains NaN values.")

    if np.isinf(contact_matrix[["probability", "distance", "cutoff"]].values).any():
        raise ValueError("ERROR: The matrix contains INF values.")

    molecule_names_dictionary = {
        name.split("_", 1)[0]: name.split("_", 1)[1] for name in ensemble_molecules_idx_sbtype_dictionary
    }

    # Access the first element and use it as a key in the dictionary
    name_mol_ai = "_" + molecule_names_dictionary[contact_matrix["molecule_name_ai"].iloc[0]]
    contact_matrix["molecule_name_ai"] = contact_matrix["molecule_name_ai"].cat.rename_categories(
        [category + name_mol_ai for category in contact_matrix["molecule_name_ai"].cat.categories]
    )

    name_mol_aj = "_" + molecule_names_dictionary[contact_matrix["molecule_name_aj"].iloc[0]]
    contact_matrix["molecule_name_aj"] = contact_matrix["molecule_name_aj"].cat.rename_categories(
        [category + name_mol_aj for category in contact_matrix["molecule_name_aj"].cat.categories]
    )

    contact_matrix["ai"] = contact_matrix["ai"].map(
        ensemble_molecules_idx_sbtype_dictionary[contact_matrix["molecule_name_ai"][0]]
    )
    contact_matrix["aj"] = contact_matrix["aj"].map(
        ensemble_molecules_idx_sbtype_dictionary[contact_matrix["molecule_name_aj"][0]]
    )

    name = path.split("/")[-1].split("_")
    len_ai = len(ensemble_molecules_idx_sbtype_dictionary[contact_matrix["molecule_name_ai"][0]])
    len_aj = len(ensemble_molecules_idx_sbtype_dictionary[contact_matrix["molecule_name_aj"][0]])
    if len_ai * len_aj != len(contact_matrix):
        raise Exception("The " + simulation + " topology and " + name[0] + " files are inconsistent")

    # Define a function to check the atom part
    # Vectorized split to extract the atom part
    ai_atoms = contact_matrix["ai"].str.split("_").str[0]
    aj_atoms = contact_matrix["aj"].str.split("_").str[0]

    # Create a mask for valid rows
    valid_rows = ~(
        (ai_atoms.str.startswith("H") & (ai_atoms != "H")) | (aj_atoms.str.startswith("H") & (aj_atoms != "H"))
    )

    contact_matrix = contact_matrix[valid_rows]

    contact_matrix = contact_matrix.assign(
        same_chain=name[0] == "intramat",
        source=pd.Categorical([simulation] * len(contact_matrix)),  # Convert to category
    )

    contact_matrix[["idx_ai", "idx_aj"]] = contact_matrix[["ai", "aj"]]
    contact_matrix.set_index(["idx_ai", "idx_aj"], inplace=True)

    return contact_matrix


def write_nonbonded(topology_dataframe, meGO_LJ, parameters, output_folder):
    """
    Writes the non-bonded parameter file ffnonbonded.itp.

    Parameters
    ----------
    topology_dataframe : pd.DataFrame
        The topology of the system as a dataframe
    meGO_LJ : pd.DataFrame
        The LJ c6 and c12 values which make up the nonbonded potential
    parameters : dict
        Contains the input parameters set from the terminal
    output_folder : str
        The path to the output directory
    """
    write_header = not parameters.no_header
    header = make_header(vars(parameters), write_header)
    with open(f"{output_folder}/ffnonbonded.itp", "w") as file:
        file.write(header)

        # write the defaults section
        file.write("\n[ defaults ]\n")
        file.write("; Include forcefield parameters\n")
        file.write("; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n")
        file.write("  1             1               no              1.0     1.0\n\n")

        file.write("[ atomtypes ]\n")
        if parameters.egos == "mg":
            atomtypes = topology_dataframe[
                ["sb_type", "atomic_number", "mass", "charge", "ptype", "mg_c6", "mg_c12"]
            ].copy()
            atomtypes.rename(columns={"mg_c6": "c6", "mg_c12": "c12"}, inplace=True)
        else:
            atomtypes = topology_dataframe[["sb_type", "atomic_number", "mass", "charge", "ptype", "c6", "c12"]].copy()

        # Only c6/c12 use scientific notation; mass and charge must keep their
        # default representation (e.g. "12.011" and "0"), so pre-format just
        # those two columns and call dataframe_to_write without float_format.
        # atomtypes is O(unique atom types) so the .map() overhead is negligible.
        atomtypes["c6"] = atomtypes["c6"].map("{:.6e}".format)
        atomtypes["c12"] = atomtypes["c12"].map("{:.6e}".format)
        file.write(dataframe_to_write(atomtypes))

        if not meGO_LJ.empty:
            file.write("\n\n[ nonbond_params ]\n")
            meGO_LJ.insert(5, ";", ";")
            # float_format is passed through to to_string(), so c6/c12 are formatted
            # at the C level without the Python-loop overhead of a prior .map() call.
            file.write(dataframe_to_write(meGO_LJ, float_format="%.6e"))


def write_model(meGO_ensemble, meGO_LJ, meGO_LJ_14, parameters, stat_str):
    """
    Takes care of the final print-out and the file writing of topology and ffnonbonded

    Parameters
    ----------
    meGO_ensemble : dict
        The meGO_ensemble object which contains all the system information
    meGO_LJ : pd.DataFrame
        Contains the c6 and c12 LJ parameters of the nonbonded potential
    meGO_LJ_14 : pd.DataFrame
        Contains the c6 and c12 LJ parameters of the pairs and exclusions
    parameters : dict
        A dictionaty of the command-line parsed parameters
    """
    output_dir = os.path.normpath(
        get_outdir_name(f"{parameters.outputs_dir}/{parameters.system}", parameters.explicit_name, parameters.egos)
    )
    create_output_directories(parameters, output_dir)
    write_topology(
        meGO_ensemble.topology_dataframe,
        meGO_ensemble.molecule_type_dict,
        meGO_ensemble.meGO_bonded_interactions,
        meGO_LJ_14,
        parameters,
        output_dir,
    )
    meGO_LJ = sort_LJ(meGO_ensemble, meGO_LJ)
    write_nonbonded(meGO_ensemble.topology_dataframe, meGO_LJ, parameters, output_dir)
    write_output_readme(meGO_LJ, parameters, output_dir, stat_str)
    return output_dir


def write_output_readme(meGO_LJ, parameters, output_dir, stat_str):
    """
    Writes a README file with the parameters used to generate the multi-eGO topology.

    Parameters
    ----------
    parameters : dict
        Contains the command-line parsed parameters
    """
    with open(f"{output_dir}/meGO.log", "w") as f:
        f.write(
            f"multi-eGO topology generated on {time.strftime('%d-%m-%Y %H:%M', time.localtime())} "
            f"using multi-eGO {_mego_version()}\n"
        )
        f.write("Parameters used to generate the topology:\n")
        for key, value in vars(parameters).items():
            f.write(f" - {key}: {value}\n")

        if parameters.egos == "production":
            # write contents of the symmetry file
            if parameters.symmetry:
                f.write("\nSymmetry file contents:\n")
                # symmetry = read_symmetry_file(parameters.symmetry)
                for line in parameters.symmetry:
                    f.write(f" - {' '.join(line)}\n")

            f.write("\nContact parameters:\n")
            f.write(stat_str)


def print_stats(meGO_LJ):
    # it would be nice to cycle over molecule types and print an half matrix with all the relevant information
    intrad_contacts = len(meGO_LJ.loc[(meGO_LJ["same_chain"])])
    interm_contacts = len(meGO_LJ.loc[~(meGO_LJ["same_chain"])])
    intrad_a_contacts = len(meGO_LJ.loc[(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)])
    interm_a_contacts = len(meGO_LJ.loc[~(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)])
    intrad_r_contacts = intrad_contacts - intrad_a_contacts
    interm_r_contacts = interm_contacts - interm_a_contacts
    intrad_a_ave_contacts = 0.000
    intrad_a_min_contacts = 0.000
    intrad_a_max_contacts = 0.000
    intrad_a_s_min_contacts = 0.000
    intrad_a_s_max_contacts = 0.000
    interm_a_ave_contacts = 0.000
    interm_a_min_contacts = 0.000
    interm_a_max_contacts = 0.000
    interm_a_s_min_contacts = 0.000
    interm_a_s_max_contacts = 0.000

    if intrad_a_contacts > 0:
        intrad_a_ave_contacts = meGO_LJ["epsilon"].loc[(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)].mean()
        intrad_a_min_contacts = meGO_LJ["epsilon"].loc[(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)].min()
        intrad_a_max_contacts = meGO_LJ["epsilon"].loc[(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)].max()
        intrad_a_s_min_contacts = meGO_LJ["sigma"].loc[(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)].min()
        intrad_a_s_max_contacts = meGO_LJ["sigma"].loc[(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)].max()

    if interm_a_contacts > 0:
        interm_a_ave_contacts = meGO_LJ["epsilon"].loc[~(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)].mean()
        interm_a_min_contacts = meGO_LJ["epsilon"].loc[~(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)].min()
        interm_a_max_contacts = meGO_LJ["epsilon"].loc[~(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)].max()
        interm_a_s_min_contacts = meGO_LJ["sigma"].loc[~(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)].min()
        interm_a_s_max_contacts = meGO_LJ["sigma"].loc[~(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)].max()

    stat_str = f"""
\t- LJ parameterization completed for a total of {len(meGO_LJ)} contacts.
\t- Attractive: intra-molecular: {intrad_a_contacts}, inter-molecular: {interm_a_contacts}
\t- Repulsive: intra-molecular: {intrad_r_contacts}, inter-molecular: {interm_r_contacts}
\t- The average epsilon is: {intrad_a_ave_contacts:5.3f} {interm_a_ave_contacts:5.3f} kJ/mol
\t- Epsilon range is: [{intrad_a_min_contacts:5.3f}:{intrad_a_max_contacts:5.3f}] [{interm_a_min_contacts:5.3f}:{interm_a_max_contacts:5.3f}] kJ/mol
\t- Sigma range is: [{intrad_a_s_min_contacts:5.3f}:{intrad_a_s_max_contacts:5.3f}] [{interm_a_s_min_contacts:5.3f}:{interm_a_s_max_contacts:5.3f}] nm

\t- RELEVANT MDP PARAMETERS:
\t- Suggested rlist value: {1.1*2.5*max(meGO_LJ['sigma'].max(), 1.00):4.2f} nm
\t- Suggested cut-off value: {2.5*max(meGO_LJ['sigma'].max(), 1.00):4.2f} nm
    """

    return stat_str


def get_outdir_name(output_dir, explicit_name, egos):
    """
    Returns the output directory name.

    Parameters
    ----------
    output_dir : str
        The path to the output directory
    explicit_name : str
        The name of the output directory

    Returns
    -------
    output_dir : str
        The path to the output directory
    """
    out = explicit_name
    if out == "":
        out = egos

    index = 1
    while os.path.exists(f"{output_dir}/{out}_{index}"):
        index += 1
        if index > 100:
            raise RuntimeError(f"too many directories in {output_dir}")
    output_dir = f"{output_dir}/{out}_{index}"

    return output_dir


def dataframe_to_write(df, float_format=None):
    """
    Returns a stringified and formatted dataframe and a message if the dataframe is empty.

    Parameters
    ----------
    df : pd.DataFrame
        The input dataframe
    float_format : str or callable, optional
        Format string (e.g. ``"%.6e"``) or callable passed to
        ``DataFrame.to_string`` for float columns.  A plain format string is
        automatically wrapped in a callable for compatibility with pandas >= 2.x.

    Returns
    -------
    str
        The stringified dataframe, or a warning comment string if *df* is empty.
    """
    if df.empty:
        # TODO insert and improve the following warning
        _term.warn("A topology parameter is empty. Check the reference topology.")
        return "; The following parameters where not parametrized on multi-eGO.\n; If this is not expected, check the reference topology."
    else:
        df.rename(columns={df.columns[0]: f"; {df.columns[0]}"}, inplace=True)
        # pandas >= 2.x requires float_format to be a callable, not a format string.
        if isinstance(float_format, str):
            fmt_str = float_format
            float_format = lambda x: fmt_str % x  # noqa: E731
        return df.to_string(index=False, float_format=float_format)


def _mego_version() -> str:
    """Return a human-readable multi-eGO version string.

    Combines the installed package version (from importlib.metadata) with the
    current git commit hash when available.  Degrades gracefully when running
    outside a git repository or from an uninstalled source tree.
    """
    try:
        version = f"v{importlib.metadata.version('mego')}"
    except importlib.metadata.PackageNotFoundError:
        version = "(uninstalled)"
    try:
        commit = git.Repo(search_parent_directories=True).head.object.hexsha[:8]
        return f"{version} (commit {commit})"
    except git.exc.InvalidGitRepositoryError:
        return version


def make_header(parameters, write_header):
    now = time.strftime("%d-%m-%Y %H:%M", time.localtime())

    header = f"""; Multi-eGO force field {_mego_version()}
; https://github.com/multi-ego/multi-eGO
; Please read and cite:
; Scalone, E. et al. PNAS 119, e2203181119 (2022) 10.1073/pnas.2203181119
; Bacic Toplek, F., Scalone, E. et al. JCTC 20, 459-468 (2024) 10.1021/acs.jctc.3c01182
"""
    if write_header:
        header += f"""
; Created on the {now}
; With the following parameters:
"""
        for parameter, value in parameters.items():
            if parameter == "no_header":
                continue
            elif parameter == "symmetry":
                header += ";\t- {:<26} :\n".format(parameter)
                for line in value:
                    header += f";\t  - {' '.join(line)}\n"
            elif isinstance(value, list):
                value = np.array(value, dtype=str)
                header += ";\t- {:<26} = {:<20}\n".format(parameter, ", ".join(value))
            elif type(value) is np.ndarray:
                value = np.array(value, dtype=str)
                header += ";\t- {:<26} = {:<20}\n".format(parameter, ", ".join(value))
            elif isinstance(value, dict):
                for key, val in value.items():
                    header += f";\t- {key} = {val}\n"
            elif not value:
                value = ""
                header += ";\t- {:<26} = {:<20}\n".format(parameter, ", ".join(value))
            else:
                header += ";\t- {:<26} = {:<20}\n".format(parameter, value)
    header += "\n"

    return header


def write_topology(
    topology_dataframe,
    molecule_type_dict,
    bonded_interactions_dict,
    meGO_LJ_14,
    parameters,
    output_folder,
):
    """
    Writes the topology output content into topol_mego.top

    Parameters
    ----------
    topology_dataframe : pd.DataFrame
        The topology of the multi-eGO system in dataframe format
    molecule_type_dict : dict
        not used yet
    bonded_interactions_dict : dict
        Contains the bonded interactions
    meGO_LJ_14 : pd.DataFrame
        Contains the c6 and c12 LJ parameters of the pairs and exclusions interactions
    parameters : dict
        Contains the command-line parsed parameters
    output_folder : str
        Path to the ouput directory
    """
    write_header = not parameters.no_header
    molecule_footer = []
    header = make_header(vars(parameters), write_header)

    with open(f"{output_folder}/topol_mego.top", "w") as file:
        header += """
; Include forcefield parameters
#include "ffnonbonded.itp"
"""

        file.write(header)
        for molecule, bonded_interactions in bonded_interactions_dict.items():
            exclusions = pd.DataFrame(columns=["ai", "aj"])
            pairs = meGO_LJ_14[molecule]
            if not pairs.empty:
                pairs.insert(5, ";", ";")
                bonded_interactions_dict[molecule]["pairs"] = pairs
                exclusions = pairs[["ai", "aj"]].copy()

            molecule_footer.append(molecule)
            molecule_header = f"""\n[ moleculetype ]
; Name\tnrexcl
{molecule}\t\t\t3

"""

            file.write(molecule_header)
            file.write("[ atoms ]\n")
            atom_selection_dataframe = topology_dataframe.loc[topology_dataframe["molecule_name"] == molecule][
                ["number", "sb_type", "resnum", "resname", "name", "cgnr"]
            ].copy()
            file.write(f"{dataframe_to_write(atom_selection_dataframe)}\n\n")
            # Here are written bonds, angles, dihedrals, impropers, and pairs
            for bonded_type, interactions in bonded_interactions.items():
                if interactions.empty:
                    continue
                label = "dihedrals" if bonded_type == "impropers" else bonded_type
                file.write(f"[ {label} ]\n")
                if bonded_type == "pairs":
                    file.write(dataframe_to_write(interactions, float_format="%.6e"))
                else:
                    file.write(dataframe_to_write(interactions))
                file.write("\n\n")
            file.write("[ exclusions ]\n")
            file.write(dataframe_to_write(exclusions))

        footer = f"""

; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif

[ system ]
{parameters.system}

[ molecules ]
; Compound #mols
"""

        file.write(footer)
        for molecule in molecule_footer:
            file.write(f"{molecule}\t\t\t1\n")


def create_output_directories(parameters, out_dir):
    """
    Creates the output directory

    Parameters
    ----------
    parameters : dict
        Contains the command-line parsed parameters

    Returns
    -------
    output_folder : str
        The path to the output directory
    """
    if not os.path.exists(parameters.outputs_dir) and not os.path.isdir(parameters.outputs_dir):
        os.mkdir(parameters.outputs_dir)
    if not os.path.exists(f"{parameters.outputs_dir}/{parameters.system}") and not os.path.isdir(
        f"{parameters.outputs_dir}/{parameters.system}"
    ):
        os.mkdir(f"{parameters.outputs_dir}/{parameters.system}")
    if not os.path.isdir(out_dir) and not os.path.exists(out_dir):
        os.mkdir(out_dir)


def check_files_existence(args):
    """
    Checks if relevant multi-eGO input files exist.

    Parameters
    ----------
    egos : str
        The egos mode of multi-eGO either 'rc' or 'production'
    system : str
        The system passed by terminal with the --system flag
    md_ensembles : list or list-like
        A list of ensembles to learn interactions from

    Raises
    ------
    FileNotFoundError
        If any of the files or directories does not exist
    """
    for ref in args.input_refs:
        md_ensembles = [ref["reference"]] + ref["train"] if args.egos == "production" else []

        if not os.path.exists(f"{args.inputs_dir}/{args.system}"):
            raise FileNotFoundError(f"Folder {args.inputs_dir}/{args.system}/ does not exist.")
        if not os.path.exists(f"{args.inputs_dir}/{args.system}/topol.top"):
            raise FileNotFoundError(f"File {args.inputs_dir}/{args.system}/topol.top does not exist.")

        for ensemble in md_ensembles:
            ensemble = f"{args.inputs_dir}/{args.system}/{ensemble}"
            if not os.path.exists(ensemble):
                raise FileNotFoundError(f"Folder {ensemble}/ does not exist.")
            else:
                top_files = glob.glob(f"{ensemble}/*.top")
                if not top_files:
                    raise FileNotFoundError(f"No .top files found in {ensemble}/")
                ndx_files = glob.glob(f"{ensemble}/*.ndx")
                ndx_files += glob.glob(f"{ensemble}/*.ndx.gz")
                ndx_files += glob.glob(f"{ensemble}/*.h5")
                if not ndx_files and not args.egos == "mg":
                    raise FileNotFoundError(
                        f"contact matrix input file(s) (e.g., intramat_1_1.ndx, etc.) were not found in {ensemble}/"
                    )


def read_intra_file(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} does not exist.")

    names = []
    epsilons = []
    with open(file_path, "r") as file:
        for line in file:
            name, param = line.strip().split(maxsplit=1)
            names.append(name)
            epsilons.append(float(param))
    epsilons = np.array(epsilons)
    return names, epsilons


def read_inter_file(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} does not exist.")

    with open(file_path, "r") as file:
        lines = file.readlines()

    # Extracting names and parameters
    names_col = np.array([line.split()[0] for line in lines[1:]])
    names_row = np.array(lines[0].split())

    # Check that the names are consistent on rows and columns (avoid mistakes)
    if np.any(names_row != names_col):
        raise ValueError(
            f"the names are inconsistent in the inter epsilon matrix:\n"
            f"  Rows:    {names_row}\n"
            f"  Columns: {names_col}\n"
            f"  Please fix to avoid silly mistakes."
        )

    epsilons = [line.split()[1:] for line in lines[1:]]
    epsilons = np.array(epsilons, dtype=float)
    if np.any(epsilons != epsilons.T):
        raise ValueError(f"the matrix of inter epsilon must be symmetric, check the input file {file_path}")
    return names_row, epsilons


def read_custom_c12_parameters(file):
    return pd.read_csv(file, names=["name", "at.num", "c12"], usecols=[0, 1, 6], header=0)


def parse_json(file_path):
    if not file_path:
        return {}
    try:
        with open(file_path, "r") as file:
            custom_dict = json.load(file)
        if not isinstance(custom_dict, dict):
            raise ValueError("Invalid dictionary format")
        return custom_dict
    except (json.JSONDecodeError, ValueError) as e:
        raise ValueError(f"Error reading custom dictionary {file_path}: {e}") from e


def sort_LJ(meGO_ensemble, meGO_LJ):
    """
    Sorts and deduplicates the final LJ DataFrame for output, ensuring ai <= aj
    ordering both within and across molecules.

    Parameters
    ----------
    meGO_ensemble : dict
        The initialized meGO ensemble.
    meGO_LJ : pd.DataFrame
        The LJ DataFrame from generate_LJ.

    Returns
    -------
    pd.DataFrame
        Sorted and deduplicated LJ DataFrame ready for writing.
    """
    meGO_LJ["type"] = 1
    meGO_LJ["number_ai"] = meGO_LJ["ai"].map(meGO_ensemble.sbtype_number_dict).astype(int)
    meGO_LJ["number_aj"] = meGO_LJ["aj"].map(meGO_ensemble.sbtype_number_dict).astype(int)

    meGO_LJ = meGO_LJ[(meGO_LJ["ai"].cat.codes <= meGO_LJ["aj"].cat.codes)].copy()

    (
        meGO_LJ["ai"],
        meGO_LJ["aj"],
        meGO_LJ["molecule_name_ai"],
        meGO_LJ["molecule_name_aj"],
        meGO_LJ["number_ai"],
        meGO_LJ["number_aj"],
    ) = np.where(
        (meGO_LJ["molecule_name_ai"].astype(str) <= meGO_LJ["molecule_name_aj"].astype(str)),
        [
            meGO_LJ["ai"],
            meGO_LJ["aj"],
            meGO_LJ["molecule_name_ai"],
            meGO_LJ["molecule_name_aj"],
            meGO_LJ["number_ai"],
            meGO_LJ["number_aj"],
        ],
        [
            meGO_LJ["aj"],
            meGO_LJ["ai"],
            meGO_LJ["molecule_name_aj"],
            meGO_LJ["molecule_name_ai"],
            meGO_LJ["number_aj"],
            meGO_LJ["number_ai"],
        ],
    )

    (
        meGO_LJ["ai"],
        meGO_LJ["aj"],
        meGO_LJ["molecule_name_ai"],
        meGO_LJ["molecule_name_aj"],
        meGO_LJ["number_ai"],
        meGO_LJ["number_aj"],
    ) = np.where(
        meGO_LJ["molecule_name_ai"] == meGO_LJ["molecule_name_aj"],
        np.where(
            meGO_LJ["number_ai"] <= meGO_LJ["number_aj"],
            [
                meGO_LJ["ai"],
                meGO_LJ["aj"],
                meGO_LJ["molecule_name_ai"],
                meGO_LJ["molecule_name_aj"],
                meGO_LJ["number_ai"],
                meGO_LJ["number_aj"],
            ],
            [
                meGO_LJ["aj"],
                meGO_LJ["ai"],
                meGO_LJ["molecule_name_aj"],
                meGO_LJ["molecule_name_ai"],
                meGO_LJ["number_aj"],
                meGO_LJ["number_ai"],
            ],
        ),
        [
            meGO_LJ["ai"],
            meGO_LJ["aj"],
            meGO_LJ["molecule_name_ai"],
            meGO_LJ["molecule_name_aj"],
            meGO_LJ["number_ai"],
            meGO_LJ["number_aj"],
        ],
    )

    meGO_LJ.sort_values(by=["molecule_name_ai", "molecule_name_aj", "number_ai", "number_aj"], inplace=True)

    final_fields = [
        "ai",
        "aj",
        "type",
        "c6",
        "c12",
        "same_chain",
        "source",
        "number_ai",
        "number_aj",
    ]

    return meGO_LJ[final_fields]
