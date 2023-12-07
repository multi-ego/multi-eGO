import pandas as pd
import time
import glob
import os


def read_molecular_contacts(path):
    """
    Reads intra-/intermat files to determine molecular contact statistics.

    Parameters
    ----------
    path : str
        The path to the file

    Returns
    -------
    contact_matrix : pd.DataFrame
        The content of the intra-/intermat file returned as a dataframe with columns
        ['molecule_number_ai', 'ai', 'molecule_number_aj', 'aj', 'distance', 'probability', 'cutoff']
    """

    print("\t-", f"Reading {path}")
    contact_matrix = pd.read_csv(path, header=None, sep="\s+")
    if contact_matrix.shape[1] == 7:
        contact_matrix.insert(7, 7, 1)
    contact_matrix.columns = [
        "molecule_number_ai",
        "ai",
        "molecule_number_aj",
        "aj",
        "distance",
        "probability",
        "cutoff",
        "intra_domain",
    ]
    contact_matrix["molecule_number_ai"] = contact_matrix["molecule_number_ai"].astype(str)
    contact_matrix["ai"] = contact_matrix["ai"].astype(str)
    contact_matrix["molecule_number_aj"] = contact_matrix["molecule_number_aj"].astype(str)
    contact_matrix["aj"] = contact_matrix["aj"].astype(str)

    return contact_matrix


def write_nonbonded(topology_dataframe, lj_potential, parameters, output_folder):
    """
    Writes the non-bonded parameter file ffnonbonded.itp.

    Parameters
    ----------
    topology_dataframe : pd.DataFrame
        The topology of the system as a dataframe
    lj_potential : pd.DataFrame
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
        file.write("[ atomtypes ]\n")
        atomtypes = topology_dataframe[["sb_type", "atomic_number", "mass", "charge", "ptype", "c6", "c12"]].copy()
        atomtypes["c6"] = atomtypes["c6"].map(lambda x: "{:.6e}".format(x))
        atomtypes["c12"] = atomtypes["c12"].map(lambda x: "{:.6e}".format(x))
        file.write(dataframe_to_write(atomtypes))

        if not lj_potential.empty:
            file.write("\n\n[ nonbond_params ]\n")
            lj_potential["c6"] = lj_potential["c6"].map(lambda x: "{:.6e}".format(x))
            lj_potential["c12"] = lj_potential["c12"].map(lambda x: "{:.6e}".format(x))
            lj_potential.insert(5, ";", ";")
            lj_potential.drop(columns=["molecule_name_ai", "molecule_name_aj"], inplace=True)
            file.write(dataframe_to_write(lj_potential))


def write_model(meGO_ensemble, meGO_LJ_potential, meGO_LJ_14, parameters, output_dir, suffix):
    """
    Takes care of the final print-out and the file writing of topology and ffnonbonded

    Parameters
    ----------
    meGO_ensemble : dict
        The meGO_ensemble object which contains all the system information
    meGO_LJ_potential : pd.DataFrame
        Contains the c6 and c12 LJ parameters of the nonbonded potential
    meGO_LJ_14 : pd.DataFrame
        Contains the c6 and c12 LJ parameters of the pairs and exclusions
    parameters : dict
        A dictionaty of the command-line parsed parameters
    output_dir : str
        Path to the output directory
    """
    print("- Writing Multi-eGO model")
    output_dir = f"{output_dir}"
    write_topology(
        meGO_ensemble["topology_dataframe"],
        meGO_ensemble["molecule_type_dict"],
        meGO_ensemble["meGO_bonded_interactions"],
        meGO_LJ_14,
        parameters,
        output_dir,
    )
    write_nonbonded(meGO_ensemble["topology_dataframe"], meGO_LJ_potential, parameters, output_dir)

    print("\n- The model is baked with the following parameters:\n")
    for argument, value in vars(parameters).items():
        if type(value) is list:
            print("\t- {:<20} = {:<20}".format(argument, ", ".join(value)))
        elif type(value) is not str:
            print("\t- {:<20} = {:<20}".format(argument, str(value)))
        else:
            print("\t- {:<20} = {:<20}".format(argument, value))
    if parameters.egos != "rc":
        print(
            f"""
        - LJ parameterization completed with a total of {len(meGO_LJ_potential)} contacts.
        - Attractive: {len(meGO_LJ_potential['epsilon'].loc[meGO_LJ_potential['epsilon']>0.])}
        - Repulsive: {len(meGO_LJ_potential['epsilon'].loc[meGO_LJ_potential['epsilon']<0.])}
        - The average epsilon is: {meGO_LJ_potential['epsilon'].loc[meGO_LJ_potential['epsilon']>0.].mean():{5}.{3}} kJ/mol
        - Epsilon range is: [{meGO_LJ_potential['epsilon'].loc[meGO_LJ_potential['epsilon']>0.].min():{5}.{3}}:{meGO_LJ_potential['epsilon'].max():{5}.{3}}] kJ/mol
        - Sigma range is: [{meGO_LJ_potential['sigma'].loc[meGO_LJ_potential['epsilon']>0.].min():{5}.{3}}:{meGO_LJ_potential['sigma'].loc[meGO_LJ_potential['epsilon']>0.].max():{5}.{3}}] nm

        RELEVANT MDP PARAMETERS:
        - Suggested rlist value: {1.1*2.5*meGO_LJ_potential['sigma'].loc[meGO_LJ_potential['epsilon']>0.].max():{4}.{3}} nm
        - Suggested cut-off value: {2.5*meGO_LJ_potential['sigma'].loc[meGO_LJ_potential['epsilon']>0.].max():{4}.{3}} nm
        """
        )
    print(f"\n- And it can be found in the following folder:\n{output_dir}")


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
; Multi-eGO force field version beta1
; https://github.com/multi-ego/multi-eGO
; Please read and cite:
; Scalone, E. et al. PNAS 119, e2203181119 (2022) 10.1073/pnas.2203181119
; Bacic Toplek, F., Scalone, E. et al. ChemRxiv (2023) 10.26434/chemrxiv-2023-67255-v2
; Created on the {now}
; With the following parameters:
"""
    for parameter, value in parameters.items():
        if parameter == "no_header":
            continue
        if type(value) is list:
            header += ";\t- {:<15} = {:<20}\n".format(parameter, ", ".join(value))
        elif not value:
            value = ""
            header += ";\t- {:<15} = {:<20}\n".format(parameter, ", ".join(value))
        else:
            header += ";\t- {:<15} = {:<20}\n".format(parameter, value)
    header += "\n"

    return header


def write_topology(
    topology_dataframe,
    molecule_type_dict,
    bonded_interactions_dict,
    lj_14,
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
    lj_14 : pd.DataFrame
        Contains the c6 and c12 LJ parameters of the pairs and exclusions interactions
    parameters : dict
        Contains the command-line parsed parameters
    output_folder : str
        Path to the ouput directory
    """
    write_header = not parameters.no_header
    molecule_footer = []
    header = make_header(vars(parameters))
    with open(f"{output_folder}/topol_GRETA.top", "w") as file:
        header += """
; Include forcefield parameters
#include "multi-ego-basic.ff/forcefield.itp"
"""

        if write_header:
            file.write(header)
        for molecule, bonded_interactions in bonded_interactions_dict.items():
            exclusions = pd.DataFrame(columns=["ai", "aj"])
            pairs = lj_14[molecule]
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
        name = f"{parameters.system}_{parameters.egos}_e{parameters.epsilon}_{parameters.inter_epsilon}"
    return name


def create_output_directories(parameters):
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
    if parameters.egos == "rc":
        name = f"{parameters.system}_{parameters.egos}"
        if parameters.out:
            name = f"{parameters.system}_{parameters.egos}_{parameters.out}"
    else:
        name = f"{parameters.system}_{parameters.egos}_e{parameters.epsilon}_{parameters.inter_epsilon}"
        if parameters.out:
            name = f"{parameters.system}_{parameters.egos}_e{parameters.epsilon}_{parameters.inter_epsilon}_{parameters.out}"
    output_folder = f"{parameters.root_dir}/outputs/{name}"

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    if os.path.isfile(f"{parameters.root_dir}/{output_folder}/ffnonbonded.itp"):
        os.remove(f"{parameters.root_dir}/{output_folder}/ffnonbonded.itp")
    if os.path.isfile(f"{parameters.root_dir}/{output_folder}/topol_GRETA.top"):
        os.remove(f"{parameters.root_dir}/{output_folder}/topol_GRETA.top")

    return output_folder


def check_files_existence(egos, system, root_dir, md_ensembles):
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
    for ensemble in md_ensembles:
        ensemble = f"{root_dir}/inputs/{system}/{ensemble}"
        if not os.path.exists(ensemble):
            raise FileNotFoundError(f"Folder {ensemble}/ does not exist.")
        else:
            top_files = glob.glob(f"{ensemble}/*.top")
            if not top_files:
                raise FileNotFoundError(f"No .top files found in {ensemble}/")
            ndx_files = glob.glob(f"{ensemble}/*.ndx")
            if not ndx_files and not egos == "rc":
                raise FileNotFoundError(f"No .ndx files found in {ensemble}/")
