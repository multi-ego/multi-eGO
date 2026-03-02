import os
import re
import sys
import tempfile
import itertools

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from multiego.resources import type_definitions
from multiego.util import masking
from multiego import io

import argparse
import multiprocessing
import numpy as np
import pandas as pd
import parmed as pmd
import time
import warnings
import tarfile
import h5py
from scipy.special import logsumexp

d = {
    type_definitions.gromos_atp.name[i]: type_definitions.gromos_atp.rc_c12[i]
    for i in range(len(type_definitions.gromos_atp.name))
}

COLUMNS = ["mi", "ai", "mj", "aj", "c12dist", "p", "cutoff"]


def write_mat(df, output_file):
    if df.empty:  # Check if the DataFrame is empty
        print(f"Warning: The DataFrame is empty. No file will be written to {output_file}.")
        return

    df = df.rename(
        columns={
            "mi": "molecule_name_ai",
            "ai": "ai",
            "mj": "molecule_name_aj",
            "aj": "aj",
            "c12dist": "distance",
            "p": "probability",
        }
    )

    df["molecule_name_ai"] = df["molecule_name_ai"].astype("category")
    df["ai"] = df["ai"].astype("category")
    df["molecule_name_aj"] = df["molecule_name_aj"].astype("category")
    df["aj"] = df["aj"].astype("category")
    df["distance"] = df["distance"].astype("float64")
    df["probability"] = df["probability"].astype("float64")
    df["cutoff"] = df["cutoff"].astype("float64")

    # Force the column order
    ordered_columns = ["molecule_name_ai", "ai", "molecule_name_aj", "aj", "distance", "probability", "cutoff", "learned"]
    df = df[ordered_columns]

    # Save the data as HDF5 with compression
    df.to_hdf(output_file, key="data", mode="w", format="table", complib="blosc:lz4", complevel=9)


def read_topologies(mego_top, target_top):
    """
    Reads the input topologies using parmed. Ignores warnings to prevent printing
    of GromacsWarnings regarding 1-4 interactions commonly seen when using
    parmed in combination with multi-eGO topologies. In the case of the reference
    topology, the last atom number is changed to 1 to prevent parmed from allocating
    unnecessary memory.

    Parameters
    ----------
    mego_top : str
        Path to the multi-eGO topology obtained from gmx pdb2gmx with multi-ego-basic force fields
    target_top : str
        Path to the toplogy of the system on which the analysis is to be performed
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            topology_mego = pmd.load_file(mego_top)
        except Exception as e:
            print(f"ERROR {e} in read_topologies while reading {mego_top}")
            exit(1)
        try:
            dirname, basename = os.path.split(target_top)
            temp_ref = tempfile.NamedTemporaryFile(prefix=basename, dir=dirname)
            temp_ref.write(open(target_top, "rb").read())
            temp_ref.seek(0)
            molecules_tag = False
            with open(temp_ref.name, "r") as f:
                lines = f.readlines()
                lines = [x for x in lines if x.strip()]

            for i, line in enumerate(lines):
                if line.strip() == "" or line[0] == ";":
                    continue
                if line.strip() == "[ molecules ]":
                    molecules_tag = True
                    continue
                if line.strip().startswith("["):
                    molecules_tag = False
                if molecules_tag and re.match(r"\s*.+\s+\d+", lines[i]):
                    print(f"Changing molecule number in line {i} that is {lines[i].strip()} to 1")
                    lines[i] = re.sub(r"(\s*.+\s+)(\d+)", r"\g<1>1", lines[i])

            with open(temp_ref.name, "w") as f:
                f.writelines(lines)
            topology_ref = pmd.load_file(temp_ref.name)

        except Exception as e:
            print(f"ERROR {e} in read_topologies while reading {target_top}")
            exit(2)

    n_mol = len(list(topology_mego.molecules.keys()))
    mol_names = list(topology_mego.molecules.keys())
    mol_list = np.arange(1, n_mol + 1, 1)

    return topology_mego, topology_ref, n_mol, mol_names, mol_list


def map_if_exists(atom_name):
    """
    Maps an atom name to a multi-eGO atom name if possible

    Parameters
    ----------
    atom_name : str
        The atom name with which to attempt the mapping

    Return
    ------
    atom_name : str
        Mapped atom name. Equal to the input if mapping was not possible
    """
    if atom_name in type_definitions.from_ff_to_multiego.keys():
        return type_definitions.from_ff_to_multiego[atom_name]
    else:
        return atom_name


def generate_c12_values(df, types, combinations, molecule_type):
    """
    TODO
    ----
    Change symmetric to be a variable
    """
    all_c12 = np.sqrt(df["c12"].to_numpy() * df["c12"].to_numpy()[:, np.newaxis])
    c12_map = np.full(all_c12.shape, None)
    resnums = df["resnum"].to_numpy()

    if molecule_type == "protein":
        for combination in combinations:
            (name_1, name_2, factor, constant, shift) = combination
            # if factor is not None and constant is not None or factor == constant:
            #    raise RuntimeError("constant and error should be defined and mutualy exclusive")
            if factor is not None and constant is not None:
                operation = lambda x: np.minimum(factor * x, constant)
            elif factor is not None:
                operation = lambda x: factor * x
            elif constant is not None:
                operation = lambda _: constant
            else:
                raise ValueError("Either factor or constant must be specified.")

            combined_map = (types[name_1] & types[name_2][:, np.newaxis]) & (resnums + shift == resnums[:, np.newaxis])
            combined_map = combined_map | combined_map.T
            c12_map = np.where(combined_map, operation(all_c12), c12_map)

    c12_map = np.where(c12_map == None, all_c12, c12_map)

    return c12_map

def warning_cutoff_histo(cutoff, max_adaptive_cutoff):
    """
    Prints warning if the histogram cutoff is smaller as the maximum adaptive cutoff.

    Parameters
    ----------
    cutoff : float
        The cutoff of the histogram calculations. Parsed from the command-line in the standard programm.
    max_adaptive_cutoff : float
        The maximum adaptive cutoff calculated from the LJ c12 parameters.
    """
    print(
        f"""
    #############################

    -------------------
    WARNING
    -------------------

    Found an adaptive cutoff greater then the cutoff used to generate the histogram:
    histogram cutoff = {cutoff}
    maximum adaptive cutoff = {max_adaptive_cutoff}

    Be careful!. This could create errors.
    If this is not wanted, please recalculate the histograms setting the cutoff to at least cutoff={max_adaptive_cutoff}

    #############################
    """
    )



def calculate_matrices(args, mapping_dict):
    """
    Starts the main routine for calculating the intermat by:
     - reading the topologies
     - figuring out all the interacting molecules
     - calculating the cutoffs
     - and caclulating the probabilities
    The operation is finalized by writing out a csv with the name pattern intermat<_name>_{mol_i}_{mol_j}.ndx

    Parameters
    ----------
    args : dict
        The command-line parsed parameters
    """
    topology_mego, topology_ref, N_species, molecules_name, mol_list = read_topologies(args.mego_top, args.target_top)

    chain_list = []
    chains = [x for x in topology_mego.molecules]

    for i in chains:
        chain_list.append(
            (
                i,
                len(topology_mego.molecules[i][0].atoms),
                len(topology_mego.split()[list(topology_mego.molecules.keys()).index(i)][1]),
            )
        )

    # number of molecules per species
    N_mols = []
    for chain in chain_list:
        N_mols.append(chain[2])
    N_mols = np.array(N_mols)

    print(
        f"""
    Topology contains {N_species} molecules species. Namely {molecules_name}.
    Calculating intermat for all species\n\n
    """
    )
    for mol_i in mol_list:
        if args.intra:
            prefix = f"intra_mol_{mol_i}_{mol_i}"
            main_routine(mol_i, mol_i, topology_mego, topology_ref, molecules_name, prefix, mapping_dict)
        # fmt: off
        for mol_j in mol_list[mol_i - 1:]:
            # fmt: on
            if mol_i == mol_j and not args.same:
                continue
            if mol_i != mol_j and not args.cross:
                continue

            prefix = f"inter_mol_{mol_i}_{mol_j}"
            main_routine(mol_i, mol_j, topology_mego, topology_ref, molecules_name, prefix, mapping_dict)


def main_routine(mol_i, mol_j, topology_mego, topology_ref, molecules_name, prefix, mapping_dict):

    df = pd.DataFrame()

    topology_df_i = pd.DataFrame()
    topology_df_j = pd.DataFrame()

    # define matrix type (intra o inter)
    mat_type = prefix.split("_")[0]
    print(
        f"\nCalculating {mat_type} between molecule {mol_i} and {mol_j}: {molecules_name[mol_i-1]} and {molecules_name[mol_j-1]}"
    )

    protein_mego_i = topology_mego.molecules[list(topology_mego.molecules.keys())[mol_i - 1]][0]
    protein_mego_j = topology_mego.molecules[list(topology_mego.molecules.keys())[mol_j - 1]][0]

    protein_ref_i = topology_ref.molecules[list(topology_ref.molecules.keys())[mol_i - 1]][0]
    protein_ref_j = topology_ref.molecules[list(topology_ref.molecules.keys())[mol_j - 1]][0]

    original_size_i = len(protein_ref_i.atoms)
    original_size_j = len(protein_ref_j.atoms)

    d_protein_ref_indices_i = np.array([i + 1 for i in range(len(protein_ref_i.atoms))])
    protein_ref_indices_i = np.array(
        [
            i + 1
            for i in range(len(protein_ref_i.atoms))
            if (protein_ref_i[i].element_name != "H" or protein_ref_i[i].name == args.bkbn_H)
        ]
    )
    protein_ref_indices_j = np.array(
        [
            i + 1
            for i in range(len(protein_ref_j.atoms))
            if (protein_ref_j[i].element_name != "H" or protein_ref_j[i].name == args.bkbn_H)
        ]
    )

    # the following lines are needed to generate a dictionary atoms - residue includining hydrongen
    # that we use for neighbour contact
    d_protein_ref_i = [a for a in protein_ref_i.atoms]
    d_sorter_i = [str(x.residue.number) + map_if_exists(x.name) for x in d_protein_ref_i]
    # Use list comprehension with regex substitution
    d_sorter_i = [re.sub(r"[a-zA-Z]+[0-9]*", "", s) for s in d_sorter_i]
    # Convert to a NumPy array and then to integers
    d_sorter_i = np.array(d_sorter_i, dtype=int)
    # create full dictionary with ai to ri
    d_ref_ai_to_ri_i = dict(zip(d_protein_ref_indices_i, d_sorter_i))

    protein_ref_i = [a for a in protein_ref_i.atoms if (a.element_name != "H" or a.name == args.bkbn_H)]
    protein_ref_j = [a for a in protein_ref_j.atoms if (a.element_name != "H" or a.name == args.bkbn_H)]

    sorter_i = [str(x.residue.number) + map_if_exists(x.name) for x in protein_ref_i]
    sorter_mego_i = [str(x.residue.number) + x.name for x in protein_mego_i]

    sorter_j = [str(x.residue.number) + map_if_exists(x.name) for x in protein_ref_j]
    sorter_mego_j = [str(x.residue.number) + x.name for x in protein_mego_j]

    # this is needed to build dictionaries that include hydrogens

    sorter_i = [str(x.residue.number) + map_if_exists(x.name) for x in protein_ref_i]
    sorter_mego_i = [str(x.residue.number) + x.name for x in protein_mego_i]

    sorter_j = [str(x.residue.number) + map_if_exists(x.name) for x in protein_ref_j]
    sorter_mego_j = [str(x.residue.number) + x.name for x in protein_mego_j]

    # preparing topology of molecule i
    topology_df_i["ref_ai"] = protein_ref_indices_i
    topology_df_i["ref_type"] = [a.name for a in protein_ref_i]
    topology_df_i["resname"] = [a.residue.name for a in protein_ref_i]
    topology_df_i["resnum"] = [a.residue.idx for a in protein_ref_i]
    topology_df_i["sorter"] = sorter_i
    topology_df_i["ref_ri"] = topology_df_i["sorter"].str.replace("[a-zA-Z]+[0-9]*", "", regex=True).astype(int)
    topology_df_i.sort_values(by="sorter", inplace=True)
    topology_df_i["mego_ai"] = [a[0].idx for a in sorted(zip(protein_mego_i, sorter_mego_i), key=lambda x: x[1])]
    topology_df_i["mego_type"] = [a[0].type for a in sorted(zip(protein_mego_i, sorter_mego_i), key=lambda x: x[1])]
    topology_df_i["mego_name"] = [a[0].name for a in sorted(zip(protein_mego_i, sorter_mego_i), key=lambda x: x[1])]
    topology_df_i["name"] = topology_df_i["mego_name"]
    topology_df_i["type"] = topology_df_i["mego_type"]
    # need to sort back otherwise c12_cutoff are all wrong
    topology_df_i.sort_values(by="ref_ai", inplace=True)
    if args.custom_c12 is not None:
        custom_c12_dict = io.read_custom_c12_parameters(args.custom_c12)
        d_appo = {key: val for key, val in zip(custom_c12_dict.name, custom_c12_dict.rc_c12)}
        d.update(d_appo)

    topology_df_i["c12"] = topology_df_i["mego_type"].map(d)

    # preparing topology of molecule j
    topology_df_j["ref_ai"] = protein_ref_indices_j
    topology_df_j["ref_type"] = [a.name for a in protein_ref_j]
    topology_df_j["sorter"] = sorter_j
    topology_df_j["resname"] = [a.residue.name for a in protein_ref_j]
    topology_df_j["resnum"] = [a.residue.idx for a in protein_ref_j]
    topology_df_j["ref_ri"] = topology_df_j["sorter"].str.replace("[a-zA-Z]+[0-9]*", "", regex=True).astype(int)
    topology_df_j.sort_values(by="sorter", inplace=True)
    topology_df_j["mego_type"] = [a[0].type for a in sorted(zip(protein_mego_j, sorter_mego_j), key=lambda x: x[1])]
    topology_df_j["mego_name"] = [a[0].name for a in sorted(zip(protein_mego_j, sorter_mego_j), key=lambda x: x[1])]
    topology_df_j["name"] = topology_df_j["mego_name"]
    topology_df_j["type"] = topology_df_j["mego_type"]
    # need to sort back otherwise c12_cutoff are all wrong
    topology_df_j.sort_values(by="ref_ai", inplace=True)
    if args.custom_c12 is not None:
        custom_c12_dict = io.read_custom_c12_parameters(args.custom_c12)
        d_appo = {key: val for key, val in zip(custom_c12_dict.name, custom_c12_dict.rc_c12)}
        d.update(d_appo)

    topology_df_j["c12"] = topology_df_j["mego_type"].map(d)
    OO_mask = masking.create_matrix_mask(
        topology_df_i["mego_type"].to_numpy(),
        topology_df_j["mego_type"].to_numpy(),
        [("O", "O"), ("OM", "O")],
        symmetrize=True,
    )
    OMOM_mask = masking.create_matrix_mask(
        topology_df_i["mego_type"].to_numpy(),
        topology_df_j["mego_type"].to_numpy(),
        [("OM", "OM")],
        symmetrize=True,
    )
    HH_mask = masking.create_matrix_mask(
        topology_df_i["mego_type"].to_numpy(),
        topology_df_j["mego_type"].to_numpy(),
        [("H", "H")],
        symmetrize=True,
    )
    NN_mask = masking.create_matrix_mask(
        topology_df_i["mego_type"].to_numpy(),
        topology_df_j["mego_type"].to_numpy(),
        [("NL", "NL")],
        symmetrize=True,
    )

    if mat_type == "intra":
        first_aminoacid = topology_mego.residues[0].name
        if first_aminoacid in type_definitions.aminoacids_list:
            molecule_type = "protein"
        elif first_aminoacid in type_definitions.nucleic_acid_list:
            molecule_type = "nucleic_acid"
        else:
            molecule_type = "other"

        types = type_definitions.lj14_generator(topology_df_i)

        if molecule_type == "other":
            # read user pairs
            molecule_keys = list(topology_mego.molecules.keys())
            user_pairs = [
                (pair.atom1.idx, pair.atom2.idx, pair.type.epsilon * 4.184)
                for pair in topology_mego.molecules[molecule_keys[mol_i - 1]][0].adjusts
            ]
            user_pairs = [
                (
                    topology_df_i[topology_df_i["mego_ai"] == ai].index[0],
                    topology_df_i[topology_df_i["mego_ai"] == aj].index[0],
                    c12,
                )
                for ai, aj, c12 in user_pairs
            ]

        c12_values = generate_c12_values(topology_df_i, types, type_definitions.atom_type_combinations, molecule_type)

        # define all cutoff using combination rule values and OO_mask
        c12_cutoff = CUTOFF_FACTOR * np.power(np.where(OO_mask, type_definitions.mg_OO_c12_rep, c12_values), 1.0 / 12.0)
        # apply OMOM correction
        c12_cutoff = np.where(OMOM_mask, CUTOFF_FACTOR * np.power(type_definitions.mg_OMOM_c12_rep, 1.0 / 12.0), c12_cutoff)
        # apply HH correction
        c12_cutoff = np.where(HH_mask, CUTOFF_FACTOR * np.power(type_definitions.mg_HH_c12_rep, 1.0 / 12.0), c12_cutoff)
        # apply NN correction
        c12_cutoff = np.where(NN_mask, CUTOFF_FACTOR * np.power(type_definitions.mg_NN_c12_rep, 1.0 / 12.0), c12_cutoff)

        # apply the user pairs (overwrite all other rules)
        if molecule_type == "other":
            for ai, aj, c12 in user_pairs:
                ai = int(ai)
                aj = int(aj)
                if c12 > 0.0:
                    c12_cutoff[ai][aj] = CUTOFF_FACTOR * np.power(c12, 1.0 / 12.0)
                    c12_cutoff[aj][ai] = CUTOFF_FACTOR * np.power(c12, 1.0 / 12.0)

    if mat_type == "inter":
        # define all cutoff
        c12_cutoff = CUTOFF_FACTOR * np.where(
            OO_mask,
            np.power(type_definitions.mg_OO_c12_rep, 1.0 / 12.0),
            np.power(
                np.sqrt(topology_df_j["c12"].values * topology_df_i["c12"].values[:, np.newaxis]),
                1.0 / 12.0,
            ),
        )
        c12_cutoff = np.where(OMOM_mask, CUTOFF_FACTOR * np.power(type_definitions.mg_OMOM_c12_rep, 1.0 / 12.0), c12_cutoff)
        c12_cutoff = np.where(HH_mask, CUTOFF_FACTOR * np.power(type_definitions.mg_HH_c12_rep, 1.0 / 12.0), c12_cutoff)
        c12_cutoff = np.where(NN_mask, CUTOFF_FACTOR * np.power(type_definitions.mg_NN_c12_rep, 1.0 / 12.0), c12_cutoff)

    mismatched = topology_df_i.loc[topology_df_i["ref_type"].str[0] != topology_df_i["mego_name"].str[0]]
    if not mismatched.empty:
        raise ValueError(f"Mismatch found:\n{mismatched}, target and mego topology are not compatible")
    mismatched = topology_df_j.loc[topology_df_j["ref_type"].str[0] != topology_df_j["mego_name"].str[0]]
    if not mismatched.empty:
        raise ValueError(f"Mismatch found:\n{mismatched}, target and mego topology are not compatible")

    if np.any(c12_cutoff > args.cutoff):
        warning_cutoff_histo(args.cutoff, np.max(c12_cutoff))
    if np.isnan(c12_cutoff.astype(float)).any():
        warning_cutoff_histo(args.cutoff, np.max(c12_cutoff))

    print(c12_cutoff.flatten()[-15:])
    ########################
    # PARALLEL PROCESS START
    ########################
    print("WAAAAAAAAAAAA")
    # print(topology_df_i["mego_type"])
    combinations = [f"{a1}_{a2}" for a1, a2 in itertools.product(topology_df_i["mego_type"], topology_df_j["mego_type"])]
    # for each atom_pair in mapping_dict, create a copy with inverted atom-pair
    print(combinations)
    mapping_dict2 = mapping_dict.copy()
    mapping_dict2["atom_pair"] = mapping_dict2["atom_pair"].apply(lambda x: x.split("_")[1] + "_" + x.split("_")[0])
    mapping_dict = pd.concat([mapping_dict, mapping_dict2], ignore_index=True)
    mapping_dict.drop_duplicates(subset=["atom_pair"], inplace=True)

    dists = []
    cutoffs = []
    probs = []
    # for combination in combinations:
    #     print(combination)
    #     exp_avg = mapping_dict.loc[(mapping_dict["atom_pair"] == combination) | (mapping_dict["atom_pair"] == "_".join(combination.split("_")[::-1])), "exp_aver"].values
    #     prob = mapping_dict.loc[(mapping_dict["atom_pair"] == combination) | (mapping_dict["atom_pair"] == "_".join(combination.split("_")[::-1])), "probability"].values
    #     cutoff = mapping_dict.loc[(mapping_dict["atom_pair"] == combination) | (mapping_dict["atom_pair"] == "_".join(combination.split("_")[::-1])), "cutoff"].values
    #     if len(exp_avg) != 1 or len(prob) != 1 or len(cutoff) != 1:
    #          raise ValueError(f"Expected exactly one entry for combination {combination}, but got {len(exp_avg)} exp_avg, {len(prob)} prob and {len(cutoff)} cutoffs")
    #     dists.append(exp_avg[0])
    #     probs.append(prob[0])
    #     cutoffs.append(cutoff[0])
    #     print(f"combination {combination} has exp_avg {exp_avg}, prob {prob} and cutoff {cutoff}")
    df = pd.DataFrame()
    all_ai = [i for i in range(1, original_size_i + 1)]
    all_aj = [j for j in range(1, original_size_j + 1)]
    df["combination"] = combinations    
    df["mi"] = [mol_i for _ in range((original_size_i) * (original_size_j))]
    df["mj"] = [mol_j for _ in range((original_size_i) * (original_size_j))]
    df["ai"] = np.repeat(all_ai, (original_size_j))
    df["aj"] = np.tile(all_aj, original_size_i)
    df["c12dist"] = df["combination"].map(mapping_dict.set_index("atom_pair")["exp_aver"])
    df["p"] = df["combination"].map(mapping_dict.set_index("atom_pair")["probability"])
    cutoff_temp = df["combination"].map(mapping_dict.set_index("atom_pair")["cutoff"])
    cuts = []
    for ai, aj in zip(df["ai"], df["aj"]):
        # Check if both ai and aj are in reference indices
        if ai in protein_ref_indices_i and aj in protein_ref_indices_j:
            # Get matrix indices
            i = np.where(protein_ref_indices_i == ai)[0][0]
            j = np.where(protein_ref_indices_j == aj)[0][0]
            cuts.append(float(c12_cutoff[i, j]))
        else:
            cuts.append(0.0)
    # cuts = np.array([float(f"{cut:6f}") for cut in cuts])
    # if np.any(np.abs(np.array(cuts) - np.array(cutoff_temp)) > 1e-6):
    #     idx = np.where(np.abs(np.array(cuts) - np.array(cutoff_temp)) > 1e-6)[0]
        
    #     print(np.array(df["combination"].iloc[idx]))
    #     print(idx)
    #     print(f"Calculated cutoffs: {cuts[idx]}")
    #     print(f"Cutoffs from mapping_dict: {np.array(df['cutoff'].iloc[idx])}")
    #     raise ValueError(f"Cutoffs from mapping_dict do not match calculated cutoffs for molecule pair {mol_i} and {mol_j}. Mismatched indices: {idx}")
    df["cutoff"] = cuts
    df.loc[df["c12dist"] > df["cutoff"], "c12dist"] = df["cutoff"]
   

    df["learned"] = 1

    if mol_i == mol_j and args.nb:
        df["ai"] = df["ai"].astype(int)
        df["aj"] = df["aj"].astype(int)
        df["residue_diff"] = df.apply(lambda row: abs(d_ref_ai_to_ri_i[row["ai"]] - d_ref_ai_to_ri_i[row["aj"]]), axis=1)
        df.loc[df["residue_diff"] > 2, ["p", "c12dist", "learned"]] = 0
        df.drop(columns=["residue_diff"], inplace=True)

    if mol_i == mol_j and args.nonb:
        df["ai"] = df["ai"].astype(int)
        df["aj"] = df["aj"].astype(int)
        df["residue_diff"] = df.apply(lambda row: abs(d_ref_ai_to_ri_i[row["ai"]] - d_ref_ai_to_ri_i[row["aj"]]), axis=1)
        df.loc[df["residue_diff"] < 3, ["p", "c12dist", "learned"]] = 0
        df.drop(columns=["residue_diff"], inplace=True)

    df["mi"] = df["mi"].map("{:}".format)
    df["mj"] = df["mj"].map("{:}".format)
    df["ai"] = df["ai"].map("{:}".format)
    df["aj"] = df["aj"].map("{:}".format)
    df["c12dist"] = df["c12dist"].map("{:,.6f}".format)
    df["p"] = df["p"].map("{:,.6e}".format)
    df["cutoff"] = df["cutoff"].map("{:,.6f}".format)
    df.index = range(len(df.index))

    out_name = args.out_name + "_" if args.out_name else ""
    output_file = f"{args.out}/{mat_type}mat_{out_name}{mol_i}_{mol_j}.ndx.h5"
    print(f"Saving output for molecule {mol_i} and {mol_j} in {output_file}")
    write_mat(df, output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--map", required=True, help="Path to the mapping dictionary")
    parser.add_argument(
        "--mego_top",
        required=True,
        help="Path to the topology file of the system on which the histograms were calculated on",
    )
    parser.add_argument(
        "--target_top",
        required=True,
        help="Path to the topology file of the system on which the histograms were calculated on",
    )
    parser.add_argument("--bkbn_H", help="Name of backbone hydrogen (default H, charmm HN)", default="H")
    parser.add_argument(
        "--custom_c12",
        type=str,
        help="Custom dictionary of c12 for special molecules",
    )
    parser.add_argument(
        "--mode", help="Sets the caculation to be intra/same/cross for histograms processing", default="intra+same+cross"
    )
    parser.add_argument(
        "--out", help="Path to the output directory where the matrices will be saved", default="./"
    )
    parser.add_argument(
        "--cutoff",
        default=0.75,
        type=float,
        help="To be set to the max cutoff used for the accumulation of the histograms",
    )
    parser.add_argument(
        "--num_threads",
        default=1,
        type=int,
        help="Sets the number of processes to perform the calculation",
    )
    parser.add_argument(
        "--nb",
        action="store_true",
        help="consider contacts only between neighbour aminoacids",
    )
    parser.add_argument(
        "--nonb",
        action="store_true",
        help="consider contacts not between neighbour aminoacids",
    )
    parser.add_argument(
        "--out_name",
        help="""Sets the output name of files to be added to the default one: intermat_<out_name>_mi_mj.ndx or intramat_<out_name>_mi_mj.ndx""",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.map):
        print(f"The file '{args.map}' does not exist.")
        sys.exit()
    # check if output file exists
    else:
        mapping_dict = pd.read_csv(args.map, sep="\t")

    if not os.path.exists(args.out):
        print(f"The path '{args.out}' does not exist.")
        sys.exit()

    # Sets mode
    modes = np.array(args.mode.split("+"), dtype=str)
    modes_possible = np.array(["intra", "same", "cross"])
    args.intra = False
    args.same = False
    args.cross = False

    if not np.any(np.isin(modes, modes_possible)):
        raise ValueError(
            f"inserted mode {args.mode} is not correct and got evaluated to {modes}. Choose intra,same and or cross separated by '+', e.g.: intra+same or same+cross"
        )

    if "intra" in modes:
        args.intra = True
    if "same" in modes:
        args.same = True
    if "cross" in modes:
        args.cross = True

    N_BINS = args.cutoff / (0.01 / 4)
    DX = args.cutoff / N_BINS
    CUTOFF_FACTOR = 1.45
    print(
        f"""
    Starting with cutoff = {args.cutoff},
                  n_bins = {N_BINS},
                  dx     = {DX}
                  on {args.num_threads} threads
    """
    )

    calculate_matrices(args, mapping_dict)
