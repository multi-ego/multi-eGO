import numpy as np
import pandas as pd
import time
import glob
import os
import yaml
import git


final_fields = [
    "ai",
    "aj",
    "type",
    "c6",
    "c12",
    "sigma",
    "epsilon",
    "probability",
    "rc_probability",
    "md_threshold",
    "rc_threshold",
    "rep",
    "cutoff",
    "molecule_name_ai",
    "molecule_name_aj",
    "same_chain",
    "source",
    "number_ai",
    "number_aj",
]


def read_config(file, args_dict):
    """
    Reads a YAML file and returns its content as a dictionary.

    Parameters
    ----------
    file : str
        The path to the YAML file

    Returns
    -------
    args_dict : dict
        The content of the YAML file as a dictionary
    """
    with open(file, "r") as f:
        yml = yaml.safe_load(f)
    # check if the keys in the yaml file are valid
    for element in yml:
        if type(element) is not dict:
            key = element
        else:
            key = list(element.keys())[0]
        if f"--{key}" not in args_dict:
            raise ValueError(f"ERROR: {key} in {file} is not a valid argument.")

    return yml


def combine_configurations(yml, args, args_dict):
    """
    Combines the configuration from a YAML file with the command-line arguments. By overwriting
    files from the YAML configuration with the command-line arguments, the function ensures that
    the command-line arguments take precedence over the configuration file. Overwriting is done
    directly on the args dictionary.

    Parameters
    ----------
    yml : dict
        The configuration from the YAML file
    args : dict
        The command-line arguments

    Returns
    -------
    dict
        The combined configuration
    """
    for element in yml:
        if type(element) is dict:
            key, value = list(element.items())[0]
            value = args_dict[f"--{key}"]["type"](value)
            parse_key = f"--{key}"
            default_value = args_dict[parse_key]["default"] if "default" in args_dict[parse_key] else None
            if hasattr(args, key) and getattr(args, key) is default_value:
                setattr(args, key, value)
        else:
            if hasattr(args, element):
                print(element, yml)
                setattr(args, element, True)

    return args


def strip_gz_suffix(filename):
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
    matrix_paths = glob.glob(f"{input_path}/int??mat_?_?.ndx")
    matrix_paths_gz = glob.glob(f"{input_path}/int??mat_?_?.ndx.gz")
    stripped_matrix_paths_gz_set = set(map(strip_gz_suffix, matrix_paths_gz))
    matrix_paths_set = set(matrix_paths)
    # Find intersection of the two sets
    common_files = matrix_paths_set.intersection(stripped_matrix_paths_gz_set)

    # Check if there are any common elements and raise an error if there are
    if common_files:
        raise ValueError(f"Error: Some files have both non-gz and gz versions: {common_files}")


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
    reference_path = f"{args.root_dir}/inputs/{args.system}/{args.reference}"
    check_matrix_compatibility(reference_path)
    for simulation in args.train:
        simulation_path = f"{args.root_dir}/inputs/{args.system}/{simulation}"
        check_matrix_compatibility(simulation_path)
    for simulation in args.check:
        simulation_path = f"{args.root_dir}/inputs/{args.system}/{simulation}"
        check_matrix_compatibility(simulation_path)


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


def read_molecular_contacts(path):
    """
    Reads intra-/intermat files to determine molecular contact statistics.
    """
    print("\t\t-", f"Reading {path}")

    # Define column names and data types directly during read
    col_names = ["molecule_number_ai", "ai", "molecule_number_aj", "aj", "distance", "probability", "cutoff", "intra_domain"]
    col_types = {
        "molecule_number_ai": str,
        "ai": str,
        "molecule_number_aj": str,
        "aj": str,
        "distance": np.float64,
        "probability": np.float64,
        "cutoff": np.float64,
    }
    contact_matrix = pd.read_csv(path, header=None, sep="\s+", names=col_names, dtype=col_types)
    contact_matrix["intra_domain"] = contact_matrix["intra_domain"].fillna(1).astype(bool)

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
    header = make_header(vars(parameters))
    with open(f"{output_folder}/ffnonbonded.itp", "w") as file:
        if write_header:
            file.write(header)

        # write the defaults section
        file.write("\n[ defaults ]\n")
        file.write("; Include forcefield parameters\n")
        file.write("#define _FF_MAGROS\n\n")
        file.write("; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n")
        file.write("  1             1               no              1.0     1.0\n\n")

        file.write("[ atomtypes ]\n")
        atomtypes = topology_dataframe[["sb_type", "atomic_number", "mass", "charge", "ptype", "c6", "c12"]].copy()
        atomtypes["c6"] = atomtypes["c6"].map(lambda x: "{:.6e}".format(x))
        atomtypes["c12"] = atomtypes["c12"].map(lambda x: "{:.6e}".format(x))
        file.write(dataframe_to_write(atomtypes))

        if not meGO_LJ.empty:
            file.write("\n\n[ nonbond_params ]\n")
            meGO_LJ["c6"] = meGO_LJ["c6"].map(lambda x: "{:.6e}".format(x))
            meGO_LJ["c12"] = meGO_LJ["c12"].map(lambda x: "{:.6e}".format(x))
            meGO_LJ.insert(5, ";", ";")
            meGO_LJ.drop(columns=["molecule_name_ai", "molecule_name_aj"], inplace=True)
            file.write(dataframe_to_write(meGO_LJ))


def write_model(meGO_ensemble, meGO_LJ, meGO_LJ_14, parameters):
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
    output_dir = get_outdir_name(
        f"{parameters.root_dir}/outputs/{parameters.system}", parameters.explicit_name, parameters.egos
    )
    create_output_directories(parameters, output_dir)
    meGO_LJ_out = meGO_LJ[final_fields].copy()
    write_topology(
        meGO_ensemble["topology_dataframe"],
        meGO_ensemble["molecule_type_dict"],
        meGO_ensemble["meGO_bonded_interactions"],
        meGO_LJ_14,
        parameters,
        output_dir,
    )
    write_nonbonded(meGO_ensemble["topology_dataframe"], meGO_LJ_out, parameters, output_dir)
    write_output_readme(meGO_LJ, parameters, output_dir)
    print(f"Output files written to {output_dir}")


def write_output_readme(meGO_LJ, parameters, output_dir):
    """
    Writes a README file with the parameters used to generate the multi-eGO topology.

    Parameters
    ----------
    parameters : dict
        Contains the command-line parsed parameters
    """
    repo = git.Repo(search_parent_directories=True)
    commit_hash = repo.head.object.hexsha
    with open(f"{output_dir}/meGO.log", "w") as f:
        f.write(
            f"multi-eGO topology generated on {time.strftime('%d-%m-%Y %H:%M', time.localtime())} using commit {commit_hash}\n"
        )
        f.write("Parameters used to generate the topology:\n")
        for key, value in vars(parameters).items():
            f.write(f" - {key}: {value}\n")

        # write contents of the symmetry file
        if parameters.symmetry:
            f.write("\nSymmetry file contents:\n")
            # symmetry = read_symmetry_file(parameters.symmetry)
            for line in parameters.symmetry:
                f.write(f" - {' '.join(line)}\n")

        f.write("\nContact parameters:\n")
        # write average data intra
        f.write("\n - Intermolecular contacts:\n")
        f.write(
            f"   - epsilon: {meGO_LJ.loc[(meGO_LJ['same_chain']) & (meGO_LJ['intra_domain']) & (meGO_LJ['epsilon'] > 0.0)]['epsilon'].mean():.3f} kJ/mol\n"
        )
        f.write(
            f"   - sigma: {meGO_LJ.loc[(meGO_LJ['same_chain']) & (meGO_LJ['intra_domain']) & (meGO_LJ['epsilon'] > 0.0)]['sigma'].mean():.3f} nm\n"
        )
        f.write(
            f"   - number of attractive contacts: {len(meGO_LJ.loc[(meGO_LJ['same_chain']) & (meGO_LJ['intra_domain']) & (meGO_LJ['epsilon'] > 0.0)])}\n"
        )
        f.write(
            f"   - number of repulsive contacts: {len(meGO_LJ.loc[(meGO_LJ['same_chain']) & (meGO_LJ['intra_domain']) & (meGO_LJ['epsilon'] < 0.0)])}\n"
        )

        # write average data interdomain
        f.write("\n - Interdomain contacts:\n")
        f.write(
            f"   - epsilon: {meGO_LJ.loc[(meGO_LJ['same_chain']) & (~meGO_LJ['intra_domain']) & (meGO_LJ['epsilon'] > 0.0)]['epsilon'].mean():.3f} kJ/mol\n"
        )
        f.write(
            f"   - sigma: {meGO_LJ.loc[(meGO_LJ['same_chain']) & (~meGO_LJ['intra_domain']) & (meGO_LJ['epsilon'] > 0.0)]['sigma'].mean():.3f} nm\n"
        )
        f.write(
            f"   - number of attractive contacts: {len(meGO_LJ.loc[(meGO_LJ['same_chain']) & (~meGO_LJ['intra_domain']) & (meGO_LJ['epsilon'] > 0.0)])}\n"
        )
        f.write(
            f"   - number of repulsive contacts: {len(meGO_LJ.loc[(meGO_LJ['same_chain']) & (~meGO_LJ['intra_domain']) & (meGO_LJ['epsilon'] < 0.0)])}\n"
        )

        # write average data inter
        f.write("\n - Intermolecular contacts:\n")
        f.write(
            f"   - epsilon: {meGO_LJ.loc[~(meGO_LJ['same_chain']) & (meGO_LJ['epsilon'] > 0.0)]['epsilon'].mean():.3f} kJ/mol\n"
        )
        f.write(f"   - sigma: {meGO_LJ.loc[~(meGO_LJ['same_chain']) & (meGO_LJ['epsilon'] > 0.0)]['sigma'].mean():.3f} nm\n")
        f.write(
            f"   - number of attractive contacts: {len(meGO_LJ.loc[~(meGO_LJ['same_chain']) & (meGO_LJ['epsilon'] > 0.0)])}\n"
        )
        f.write(
            f"   - number of repulsive contacts: {len(meGO_LJ.loc[~(meGO_LJ['same_chain']) & (meGO_LJ['epsilon'] < 0.0)])}\n"
        )

        # mdp parameters
        f.write("\nCutoff MDP parameters:\n")
        f.write(f" - Suggested rlist value: {1.1*2.5*meGO_LJ['sigma'].max():4.2f} nm\n")
        f.write(f" - Suggested cut-off value: {2.5*meGO_LJ['sigma'].max():4.2f} nm\n")


def print_stats(meGO_LJ):
    # it would be nice to cycle over molecule types and print an half matrix with all the relevant information
    intrad_contacts = len(meGO_LJ.loc[(meGO_LJ["same_chain"]) & (meGO_LJ["intra_domain"])])
    interd_contacts = len(meGO_LJ.loc[(meGO_LJ["same_chain"]) & (~meGO_LJ["intra_domain"])])
    interm_contacts = len(meGO_LJ.loc[~(meGO_LJ["same_chain"])])
    intrad_a_contacts = len(meGO_LJ.loc[(meGO_LJ["same_chain"]) & (meGO_LJ["intra_domain"]) & (meGO_LJ["epsilon"] > 0.0)])
    interd_a_contacts = len(meGO_LJ.loc[(meGO_LJ["same_chain"]) & (~meGO_LJ["intra_domain"]) & (meGO_LJ["epsilon"] > 0.0)])
    interm_a_contacts = len(meGO_LJ.loc[~(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)])
    intrad_r_contacts = intrad_contacts - intrad_a_contacts
    interd_r_contacts = interd_contacts - interd_a_contacts
    interm_r_contacts = interm_contacts - interm_a_contacts
    intrad_a_ave_contacts = 0.000
    intrad_a_min_contacts = 0.000
    intrad_a_max_contacts = 0.000
    intrad_a_s_min_contacts = 0.000
    intrad_a_s_max_contacts = 0.000
    interd_a_ave_contacts = 0.000
    interd_a_min_contacts = 0.000
    interd_a_max_contacts = 0.000
    interd_a_s_min_contacts = 0.000
    interd_a_s_max_contacts = 0.000
    interm_a_ave_contacts = 0.000
    interm_a_min_contacts = 0.000
    interm_a_max_contacts = 0.000
    interm_a_s_min_contacts = 0.000
    interm_a_s_max_contacts = 0.000

    if intrad_a_contacts > 0:
        intrad_a_ave_contacts = (
            meGO_LJ["epsilon"].loc[(meGO_LJ["same_chain"]) & (meGO_LJ["intra_domain"]) & (meGO_LJ["epsilon"] > 0.0)].mean()
        )
        intrad_a_min_contacts = (
            meGO_LJ["epsilon"].loc[(meGO_LJ["same_chain"]) & (meGO_LJ["intra_domain"]) & (meGO_LJ["epsilon"] > 0.0)].min()
        )
        intrad_a_max_contacts = (
            meGO_LJ["epsilon"].loc[(meGO_LJ["same_chain"]) & (meGO_LJ["intra_domain"]) & (meGO_LJ["epsilon"] > 0.0)].max()
        )
        intrad_a_s_min_contacts = (
            meGO_LJ["sigma"].loc[(meGO_LJ["same_chain"]) & (meGO_LJ["intra_domain"]) & (meGO_LJ["epsilon"] > 0.0)].min()
        )
        intrad_a_s_max_contacts = (
            meGO_LJ["sigma"].loc[(meGO_LJ["same_chain"]) & (meGO_LJ["intra_domain"]) & (meGO_LJ["epsilon"] > 0.0)].max()
        )

    if interd_a_contacts > 0:
        interd_a_ave_contacts = (
            meGO_LJ["epsilon"].loc[(meGO_LJ["same_chain"]) & (~meGO_LJ["intra_domain"]) & (meGO_LJ["epsilon"] > 0.0)].mean()
        )
        interd_a_min_contacts = (
            meGO_LJ["epsilon"].loc[(meGO_LJ["same_chain"]) & (~meGO_LJ["intra_domain"]) & (meGO_LJ["epsilon"] > 0.0)].min()
        )
        interd_a_max_contacts = (
            meGO_LJ["epsilon"].loc[(meGO_LJ["same_chain"]) & (~meGO_LJ["intra_domain"]) & (meGO_LJ["epsilon"] > 0.0)].max()
        )
        interd_a_s_min_contacts = (
            meGO_LJ["sigma"].loc[(meGO_LJ["same_chain"]) & (~meGO_LJ["intra_domain"]) & (meGO_LJ["epsilon"] > 0.0)].min()
        )
        interd_a_s_max_contacts = (
            meGO_LJ["sigma"].loc[(meGO_LJ["same_chain"]) & (~meGO_LJ["intra_domain"]) & (meGO_LJ["epsilon"] > 0.0)].max()
        )

    if interm_a_contacts > 0:
        interm_a_ave_contacts = meGO_LJ["epsilon"].loc[~(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)].mean()
        interm_a_min_contacts = meGO_LJ["epsilon"].loc[~(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)].min()
        interm_a_max_contacts = meGO_LJ["epsilon"].loc[~(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)].max()
        interm_a_s_min_contacts = meGO_LJ["sigma"].loc[~(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)].min()
        interm_a_s_max_contacts = meGO_LJ["sigma"].loc[~(meGO_LJ["same_chain"]) & (meGO_LJ["epsilon"] > 0.0)].max()

    print(
        f"""
    - LJ parameterization completed for a total of {len(meGO_LJ)} contacts.
    - Attractive: intra-domain: {intrad_a_contacts}, inter-domain: {interd_a_contacts}, inter-molecular: {interm_a_contacts}
    - Repulsive: intra-domain: {intrad_r_contacts}, inter-domain: {interd_r_contacts}, inter-molecular: {interm_r_contacts}
    - The average epsilon is: {intrad_a_ave_contacts:5.3f} {interd_a_ave_contacts:5.3f} {interm_a_ave_contacts:5.3f} kJ/mol
    - Epsilon range is: [{intrad_a_min_contacts:5.3f}:{intrad_a_max_contacts:5.3f}] [{interd_a_min_contacts:5.3f}:{interd_a_max_contacts:5.3f}] [{interm_a_min_contacts:5.3f}:{interm_a_max_contacts:5.3f}] kJ/mol
    - Sigma range is: [{intrad_a_s_min_contacts:5.3f}:{intrad_a_s_max_contacts:5.3f}] [{interd_a_s_min_contacts:5.3f}:{interd_a_s_max_contacts:5.3f}] [{interm_a_s_min_contacts:5.3f}:{interm_a_s_max_contacts:5.3f}] nm

    RELEVANT MDP PARAMETERS:
    - Suggested rlist value: {1.1*2.5*meGO_LJ['sigma'].max():4.2f} nm
    - Suggested cut-off value: {2.5*meGO_LJ['sigma'].max():4.2f} nm
    """
    )


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
            print(f"ERROR: too many directories in {output_dir}")
            exit()
    output_dir = f"{output_dir}/{out}_{index}"

    return output_dir


def dataframe_to_write(df):
    """
    Returns a stringified and formated dataframe and a message if the dataframe is empty.

    Parameters
    ----------
    df : pd.DataFrame
        The input dataframe

    Returns
    -------
    The stringified dataframe
    """
    if df.empty:
        # TODO insert and improve the following warning
        print("A topology parameter is empty. Check the reference topology.")
        return "; The following parameters where not parametrized on multi-eGO.\n; If this is not expected, check the reference topology."
    else:
        df.rename(columns={df.columns[0]: f"; {df.columns[0]}"}, inplace=True)
        return df.to_string(index=False)


def make_header(parameters):
    now = time.strftime("%d-%m-%Y %H:%M", time.localtime())

    header = f"""
; Multi-eGO force field version beta.1
; https://github.com/multi-ego/multi-eGO
; Please read and cite:
; Scalone, E. et al. PNAS 119, e2203181119 (2022) 10.1073/pnas.2203181119
; Bacic Toplek, F., Scalone, E. et al. JCTC 20, 459-468 (2024) 10.1021/acs.jctc.3c01182
; Created on the {now}
; With the following parameters:
"""
    for parameter, value in parameters.items():
        if parameter == "no_header":
            continue
        elif parameter == "symmetry":
            header += ";\t- {:<26}:\n".format(parameter)
            for line in value:
                header += f";\t  - {' '.join(line)}\n"
        elif parameter == "names_inter":
            n = value.size
            # indices_upper_tri = np.triu_indices(n)
            tuple_list = np.array([f"({value[i]}-{value[j]})" for i, j in zip(*np.triu_indices(n))], dtype=str)
            header += ";\t- {:<26} = {:<20}\n".format(parameter, ", ".join(tuple_list))
            continue
        elif type(value) is list:
            value = np.array(value, dtype=str)
            header += ";\t- {:<26} = {:<20}\n".format(parameter, ", ".join(value))
        elif type(value) is np.ndarray:
            value = np.array(value, dtype=str)
            header += ";\t- {:<26} = {:<20}\n".format(parameter, ", ".join(value))
        elif type(value) is dict:
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
    Writes the topology output content into GRETA_topol.top

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
    header = ""
    if write_header:
        header = make_header(vars(parameters))

    with open(f"{output_folder}/topol_GRETA.top", "w") as file:
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
                pairs["c6"] = pairs["c6"].map(lambda x: "{:.6e}".format(x))
                pairs["c12"] = pairs["c12"].map(lambda x: "{:.6e}".format(x))
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
            # Here are written bonds, angles, dihedrals and impropers
            for bonded_type, interactions in bonded_interactions.items():
                if interactions.empty:
                    continue
                else:
                    if bonded_type == "impropers":
                        file.write("[ dihedrals ]\n")
                    else:
                        file.write(f"[ {bonded_type} ]\n")
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


# TODO is it ever used?
def get_name(parameters):
    """
    Creates the output directory name.

    Parameters
    ----------
    parameters : dict
        Contains the parameters parsed from the terminal input

    Returns
    -------
    name : str
        The name of the output directory
    """
    if parameters.egos == "rc":
        name = f"{parameters.system}_{parameters.egos}"
    else:
        name = f"{parameters.system}_{parameters.egos}_epsis_intra{ '-'.join(np.array(parameters.multi_epsilon, dtype=str)) }_{parameters.inter_epsilon}"
    return name


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
    if not os.path.exists(f"{parameters.root_dir}/outputs") and not os.path.isdir(f"{parameters.root_dir}/outputs"):
        os.mkdir(f"{parameters.root_dir}/outputs")
    if not os.path.exists(f"{parameters.root_dir}/outputs/{parameters.system}") and not os.path.isdir(
        f"{parameters.root_dir}/outputs/{parameters.system}"
    ):
        os.mkdir(f"{parameters.root_dir}/outputs/{parameters.system}")
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
    md_ensembles = args.reference + args.train + args.check

    for ensemble in md_ensembles:
        ensemble = f"{args.root_dir}/inputs/{args.system}/{ensemble}"
        if not os.path.exists(ensemble):
            raise FileNotFoundError(f"Folder {ensemble}/ does not exist.")
        else:
            top_files = glob.glob(f"{ensemble}/*.top")
            if not top_files:
                raise FileNotFoundError(f"No .top files found in {ensemble}/")
            ndx_files = glob.glob(f"{ensemble}/*.ndx")
            ndx_files += glob.glob(f"{ensemble}/*.ndx.gz")
            if not ndx_files and not args.egos == "rc":
                raise FileNotFoundError(
                    f"contact matrix input file(s) (e.g., intramat_1_1.ndx, etc.) were not found in {ensemble}/"
                )


def read_intra_file(file_path):
    if not os.path.exists(file_path):
        print(f"File {file_path} does not exist.")
        exit()

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
        print(f"File {file_path} does not exist.")
        exit()

    with open(file_path, "r") as file:
        lines = file.readlines()

    # Extracting names and parameters
    names_col = np.array([line.split()[0] for line in lines[1:]])
    names_row = np.array(lines[0].split())

    # Check that the names are consistent on rows and columns (avoid mistakes)
    if np.any(names_row != names_col):
        print(
            f"""ERROR: the names are inconsistent in the inter epsilon matrix:
              Rows:{names_row}
              Columns:{names_col}
              Please fix to be sure to avoid silly mistakes
              """
        )
        exit()

    epsilons = [line.split()[1:] for line in lines[1:]]
    epsilons = np.array(epsilons, dtype=float)
    if np.any(epsilons != epsilons.T):
        print(f"ERROR: the matrix of inter epsilon must be symmetric, check the input file {file_path}")
        exit()
    return names_row, epsilons


def read_custom_c12_parameters(file):
    return pd.read_csv(file, names=["name", "at.num", "c12"], usecols=[0, 1, 6], header=0)
