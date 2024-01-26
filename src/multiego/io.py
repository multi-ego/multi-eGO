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
    contact_matrix["intra_domain"] = contact_matrix["intra_domain"].astype(bool)

    if len(contact_matrix.loc[(contact_matrix["probability"] < 0) | (contact_matrix["probability"] > 1)].values) > 0:
        print("ERROR: check your matrix, probabilities should be between 0 and 1.")
        exit()
    if (
        len(
            contact_matrix.loc[
                (contact_matrix["distance"] < 0) | (contact_matrix["distance"] > contact_matrix["cutoff"])
            ].values
        )
        > 0
    ):
        print("ERROR: check your matrix, distances should be between 0 and cutoff (last column)")
        exit()
    if len(contact_matrix.loc[(contact_matrix["cutoff"] < 0)].values) > 0:
        print("ERROR: check your matrix, cutoff values cannot be negative")
        exit()

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
    output_dir = create_output_directories(parameters)
    write_topology(
        meGO_ensemble["topology_dataframe"],
        meGO_ensemble["molecule_type_dict"],
        meGO_ensemble["meGO_bonded_interactions"],
        meGO_LJ_14,
        parameters,
        output_dir,
    )
    write_nonbonded(meGO_ensemble["topology_dataframe"], meGO_LJ, parameters, output_dir)

    print(f"{output_dir}")


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
    if not os.path.exists(f"{parameters.root_dir}/outputs"):
        os.mkdir(f"{parameters.root_dir}/outputs")

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
    md_ensembles = [args.reference] + args.train + args.check

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
