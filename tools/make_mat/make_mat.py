import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from multiego.resources import type_definitions
from multiego.util import masking
from multiego import io

import argparse
import itertools
import multiprocessing
import numpy as np
import pandas as pd
import parmed as pmd
import time
import warnings
import gzip
import tarfile

d = {
    type_definitions.gromos_atp.name[i]: type_definitions.gromos_atp.c12[i]
    for i in range(len(type_definitions.gromos_atp.name))
}

COLUMNS = ["mi", "ai", "mj", "aj", "c12dist", "p", "cutoff"]


def write_mat(df, output_file):
    out_content = df.to_string(index=False, header=False, columns=COLUMNS)
    out_content = out_content.replace("\n", "<")
    out_content = " ".join(out_content.split())
    out_content = out_content.replace("<", "\n")
    out_content += "\n"
    with gzip.open(output_file, "wt") as f:
        f.write(out_content)


def read_mat(name, protein_ref_indices, args, cumulative=False):
    path_prefix = f"{args.histo}"
    if args.tar:
        with tarfile.open(args.histo, "r:*") as tar:
            ref_df = pd.read_csv(tar.extractfile(name), header=None, sep="\s+", usecols=[0, *protein_ref_indices])
    else:
        ref_df = pd.read_csv(f"{path_prefix}/{name}", header=None, sep="\s+", usecols=[0, *protein_ref_indices])
    ref_df_columns = ["distance", *[str(x) for x in protein_ref_indices]]
    ref_df.columns = ref_df_columns
    ref_df.set_index("distance", inplace=True)

    return ref_df


def run_intra_(arguments):
    """
    Preforms the main routine of the histogram analysis to obtain the intra- and intermat files.
    Is used in combination with multiprocessing to speed up the calculations.

    Parameters
    ----------
    arguments : dict
        Contains all the command-line parsed arguments

    Returns
    -------
    out_path : str
        Path to the temporary file which contains a partial pd.DataFrame with the analyzed data
    """
    (
        args,
        protein_ref_indices_i,
        protein_ref_indices_j,
        original_size_j,
        c12_cutoff,
        mi,
        mj,
        frac_target_list,
    ) = arguments
    process = multiprocessing.current_process()
    df = pd.DataFrame()
    # We do not consider old histograms
    frac_target_list = [x for x in frac_target_list if x[0] != "#" and x[-1] != "#"]
    for i, ref_f in enumerate(frac_target_list):
        results_df = pd.DataFrame()
        ai = ref_f.split(".")[-2].split("_")[-1]

        if True:
            all_ai = [ai for _ in range(1, original_size_j + 1)]
            range_list = [str(x) for x in range(1, original_size_j + 1)]

            results_df["ai"] = np.array(all_ai).astype(int)
            results_df["aj"] = np.array(range_list).astype(int)

            results_df["mi"] = mi
            results_df["mj"] = mj
            results_df["c12dist"] = 0.0
            results_df["p"] = 0.0
            results_df["cutoff"] = 0.0

        if np.isin(int(ai), protein_ref_indices_i):
            cut_i = np.where(protein_ref_indices_i == int(ai))[0][0]

            # column mapping
            ref_df = read_mat(ref_f, protein_ref_indices_j, args)
            ref_df.loc[len(ref_df)] = c12_cutoff[cut_i]

            # calculate data
            c12dist = ref_df.apply(lambda x: c12_avg(ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).values
            p = ref_df.apply(
                lambda x: calculate_probability(ref_df.index.to_numpy(), weights=x.to_numpy()),
                axis=0,
            ).values
            results_df.loc[results_df["aj"].isin(protein_ref_indices_j), "c12dist"] = c12dist
            results_df.loc[results_df["aj"].isin(protein_ref_indices_j), "p"] = p
            results_df.loc[results_df["aj"].isin(protein_ref_indices_j), "cutoff"] = c12_cutoff[cut_i].astype(float)

        df = pd.concat([df, results_df])
        df = df.sort_values(by=["p", "c12dist"], ascending=True)

    df.fillna(0)
    out_path = f"mat_{process.pid}_t{time.time()}.part"
    df.to_csv(out_path, index=False)

    return out_path


def run_inter_(arguments):
    """
    Preforms the main routine of the histogram analysis to obtain the intra- and intermat files.
    Is used in combination with multiprocessing to speed up the calculations.

    Parameters
    ----------
    arguments : dict
        Contains all the command-line parsed arguments

    Returns
    -------
    out_path : str
        Path to the temporary file which contains a partial pd.DataFrame with the analyzed data
    """
    (
        args,
        protein_ref_indices_i,
        protein_ref_indices_j,
        original_size_j,
        c12_cutoff,
        mi,
        mj,
        frac_target_list,
    ) = arguments
    process = multiprocessing.current_process()
    df = pd.DataFrame()
    # We do not consider old histograms
    frac_target_list = [x for x in frac_target_list if x[0] != "#" and x[-1] != "#"]
    for i, ref_f in enumerate(frac_target_list):
        results_df = pd.DataFrame()
        ai = ref_f.split(".")[-2].split("_")[-1]

        if True:
            all_ai = [ai for _ in range(1, original_size_j + 1)]
            range_list = [str(x) for x in range(1, original_size_j + 1)]

            results_df["ai"] = np.array(all_ai).astype(int)
            results_df["aj"] = np.array(range_list).astype(int)

            results_df["mi"] = mi
            results_df["mj"] = mj
            results_df["c12dist"] = 0.0
            results_df["p"] = 0.0
            results_df["cutoff"] = 0.0

        if np.isin(int(ai), protein_ref_indices_i):
            cut_i = np.where(protein_ref_indices_i == int(ai))[0][0]

            # column mapping
            ref_df = read_mat(ref_f, protein_ref_indices_j, args)
            ref_df.loc[len(ref_df)] = c12_cutoff[cut_i]

            # repeat for cumulative
            c_ref_f = ref_f.replace("inter_mol_", "inter_mol_c_")
            c_ref_df = read_mat(c_ref_f, protein_ref_indices_j, args, True)
            c_ref_df.loc[len(c_ref_df)] = c12_cutoff[cut_i]

            # calculate data
            c12dist = ref_df.apply(lambda x: c12_avg(ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).values
            p = c_ref_df.apply(
                lambda x: get_cumulative_probability(c_ref_df.index.to_numpy(), weights=x.to_numpy()),
                axis=0,
            ).values

            results_df.loc[results_df["aj"].isin(protein_ref_indices_j), "c12dist"] = c12dist
            results_df.loc[results_df["aj"].isin(protein_ref_indices_j), "p"] = p
            results_df.loc[results_df["aj"].isin(protein_ref_indices_j), "cutoff"] = c12_cutoff[cut_i].astype(float)

        df = pd.concat([df, results_df])

        df = df.sort_values(by=["p", "c12dist"], ascending=True)

    df.fillna(0)
    out_path = f"mat_{process.pid}_t{time.time()}.part"
    df.to_csv(out_path, index=False)

    return out_path


def run_residue_inter_(arguments):
    """
    Preforms the main routine of the histogram analysis to obtain the intra- and intermat files.
    Is used in combination with multiprocessing to speed up the calculations.

    Parameters
    ----------
    arguments : dict
        Contains all the command-line parsed arguments

    Returns
    -------
    out_path : str
        Path to the temporary file which contains a partial pd.DataFrame with the analyzed data
    """
    (
        args,
        protein_ref_indices_i,
        protein_ref_indices_j,
        num_res_j,
        c12_cutoff,
        mi,
        mj,
        (ref_ai_to_ri_i, index_ai_to_ri_j),
        frac_target_list,
    ) = arguments
    process = multiprocessing.current_process()
    df = pd.DataFrame()
    # We do not consider old histograms
    for res in frac_target_list:
        p = 0.0
        c12dist = 0.0

        for ref_f in res:
            results_df = pd.DataFrame()
            ai = int(ref_f.split(".")[-2].split("_")[-1])
            ri = ref_ai_to_ri_i[ai]

            all_ai = [ri for _ in range(1, num_res_j + 1)]
            range_list = [str(x) for x in range(1, num_res_j + 1)]

            results_df["ai"] = np.array(all_ai).astype(int)
            results_df["aj"] = np.array(range_list).astype(int)

            results_df["mi"] = mi
            results_df["mj"] = mj
            results_df["c12dist"] = 0.0
            results_df["p"] = 0.0
            results_df["cutoff"] = 0.0

            if np.isin(int(ai), protein_ref_indices_i):
                cut_i = np.where(protein_ref_indices_i == int(ai))[0][0]

                # column mapping
                ref_df = read_mat(ref_f, protein_ref_indices_j, args)
                ref_df.loc[len(ref_df)] = c12_cutoff[cut_i]

                # repeat for cumulative
                c_ref_f = ref_f.replace("inter_mol_", "inter_mol_c_")
                c_ref_df = read_mat(c_ref_f, protein_ref_indices_j, args, True)
                c_ref_df.loc[len(c_ref_df)] = c12_cutoff[cut_i]

                # calculate data
                new_p = c_ref_df.apply(
                    lambda x: get_cumulative_probability(c_ref_df.index.to_numpy(), weights=x.to_numpy()),
                    axis=0,
                ).to_numpy()
                ridx = np.array([index_ai_to_ri_j[aj] for aj in range(len(new_p))])
                new_p = np.array([np.max(new_p[ridx == i]) for i in set(ridx)])
                greater_p = new_p > p
                p = np.where(greater_p, new_p, p)
                new_c12dist = ref_df.apply(lambda x: c12_avg(ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).to_numpy()
                new_c12dist = np.array([np.mean(new_c12dist[ridx == i]) for i in set(ridx)])
                c12dist = np.where(greater_p, new_c12dist, c12dist)

            results_df["c12dist"] = c12dist
            results_df["p"] = p

        df = pd.concat([df, results_df])

        df = df.sort_values(by=["p", "c12dist"], ascending=True)

    df.fillna(0)
    out_path = f"mat_{process.pid}_t{time.time()}.part"
    df.to_csv(out_path, index=False)

    return out_path


def read_topologies(mego_top, target_top):
    """
    Reads the input topologies using parmed. Ignores warnings to prevent printing
    of GromacsWarnings regarding 1-4 interactions commonly seen when using
    parmed in combination with multi-eGO topologies.

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
            topology_ref = pmd.load_file(target_top)
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


def hallfunction(values, weights):
    """
    A substitute to the allfunction without considering the cutoff.
    Does not change the data except for the removal of the cutoff from the array.

    Parameters
    ----------
    values : np.array
        The array of the histograms x values
    weights : np.array
        The array with the respective weights

    Return
    ------
    norm : float
        The normalization constant
    v : np.array
        The unchanged x values of the histogram
    w : np.array
        The unchanged weights of the histogram
    """
    v = values[:-1]
    w = weights[:-1]
    norm = np.sum(w)
    return norm, v, w


def allfunction(values, weights):
    """
    TODO rename pls

    Preprocesses arrays (histograms) to allow for proper analysis. Last values are removed from the arrays
    and should correspond to the respective cutoff for the histogram. The histograms are truncated
    according to the cutoff.

    Parameters
    ----------
    values : np.array
        The array of the histograms x values
    weights : np.array
        The array with the respective weights

    Returns
    -------
    cutoff : float
        The cutoff which is deduced by reading the last value of the weights array
    i : int
        The index at which the cutoff is greter or equal than the values array
    norm : float
        The new normalization constant after truncation
    v : np.array
        The truncated x values of the histogram according to the cutoff
    w : np.array
        The truncated weights of the histogram according to the cutoff
    """
    v = values[:-1]
    cutoff = weights[len(weights) - 1]
    w = weights[:-1]
    i = np.where(v <= cutoff)
    if not np.any(i):
        return 0, 0, 0, 0, 0  # check if empty
    i = i[0]
    w = w[i]
    v = v[i]
    norm = np.sum(w)
    i = i[-1]
    return cutoff, i, norm, v, w


def zero_callback(values, weights):
    """
    A zero callback doing nothing but returning the normalization constant, values and weights.
    The first two return values so the function can be switched with allfunction.

    Parameters
    ----------
    values : np.array
        The array of the histograms x values
    weights : np.array
        The array with the respective weights

    Returns
    -------
    None : None
        None
    None : None
        None
    np.sum(weights) : float
        The normalization constant
    values : np.array
        The array of the histograms x values
    weights : np.array
        The array with the respective weights
    """
    return None, None, np.sum(weights), values, weights


def get_cumulative_probability(values, weights, callback=allfunction):
    cutoff, i, norm, v, w = callback(values, weights)
    return weights[i]


def c12_avg(values, weights, callback=allfunction):
    """
    Calculates the c12 averaging of a histogram as 1 / ( (\sum_i^n w[i] * (1 / x[i])^12 ) / norm )^(1/12)

    Parameters
    ----------
    values : np.array
        The array of the histograms x values
    weights : np.array
        The array with the respective weights
    callback : function
        Preprocesses the data before going in to the analysis

    Returns
    -------
    The c12 average
    """
    cutoff, i, norm, v, w = callback(values, weights)
    if np.sum(w) == 0:
        return 0
    r = np.where(w > 0.0)
    v = v[r[0][0]:v.size]
    w = w[r[0][0]:w.size]

    res = np.maximum(cutoff / 4.5, 0.1)
    exp_aver = (1.0 / res) / np.log(np.sum(w * np.exp(1.0 / v / res)) / norm)

    return exp_aver


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
            if factor is not None and constant is not None or factor == constant:
                raise RuntimeError("constant and error should be defined and mutualy exclusive")
            if factor:
                operation = lambda x: factor * x
            if constant:
                operation = lambda _: constant
            combined_map = (types[name_1] & types[name_2][:, np.newaxis]) & (resnums + shift == resnums[:, np.newaxis])
            combined_map = combined_map | combined_map.T
            c12_map = np.where(combined_map, operation(all_c12), c12_map)

    c12_map = np.where(c12_map == None, all_c12, c12_map)

    return c12_map


def calculate_intra_probabilities(args):
    """
    Starts the main routine for calculating the intramat by:
     - reading the topologies
     - calculating the cutoffs
     - and caclulating the probabilities
    The operation is finalized by writing out a csv with the name pattern intramat<_name>_{mol_i}_{mol_j}.ndx

    Parameters
    ----------
    args : dict
        The command-line parsed parameters
    """
    (
        topology_mego,
        topology_ref,
        N_molecules,
        molecules_name,
        mol_list,
    ) = read_topologies(args.mego_top, args.target_top)

    print(
        f"""
    Topology contains {N_molecules} molecules species. Namely {molecules_name}.
    Calculating intramat for all species
    """
    )
    for i in range(N_molecules):
        print(f"\n Calculating intramat for molecule {mol_list[i]}: {molecules_name[i]}")
        df = pd.DataFrame()
        topology_df = pd.DataFrame()

        prefix = f"intra_mol_{mol_list[i]}_{mol_list[i]}"
        if args.tar:
            with tarfile.open(args.histo, "r:*") as tar:
                target_list = [x.name for x in tar.getmembers() if prefix in x.name]
        else:
            target_list = [x for x in os.listdir(args.histo) if prefix in x]

        protein_mego = topology_mego.molecules[list(topology_mego.molecules.keys())[i]][0]
        protein_ref = topology_ref.molecules[list(topology_ref.molecules.keys())[i]][0]
        original_size = len(protein_ref.atoms)
        protein_ref_indices = np.array([i + 1 for i in range(len(protein_ref.atoms)) if protein_ref[i].element_name != "H"])
        protein_ref = [a for a in protein_ref.atoms if a.element_name != "H"]

        sorter = [str(x.residue.number) + map_if_exists(x.name) for x in protein_ref]
        sorter_mego = [str(x.residue.number) + x.name for x in protein_mego]

        topology_df["ref_ai"] = protein_ref_indices
        topology_df["ref_type"] = [a.name for a in protein_ref]
        topology_df["resname"] = [a.residue.name for a in protein_ref]
        topology_df["resnum"] = [a.residue.idx for a in protein_ref]
        topology_df["sorter"] = sorter
        topology_df["ref_ri"] = topology_df["sorter"].str.replace("[a-zA-Z]+[0-9]*", "", regex=True).astype(int)
        topology_df.sort_values(by="sorter", inplace=True)
        topology_df["mego_ai"] = [a[0].idx for a in sorted(zip(protein_mego, sorter_mego), key=lambda x: x[1])]
        topology_df["mego_type"] = [a[0].type for a in sorted(zip(protein_mego, sorter_mego), key=lambda x: x[1])]
        topology_df["mego_name"] = [a[0].name for a in sorted(zip(protein_mego, sorter_mego), key=lambda x: x[1])]
        topology_df["name"] = topology_df["mego_name"]
        topology_df["type"] = topology_df["mego_type"]
        # need to sort back otherwise c12_cutoff are all wrong
        topology_df.sort_values(by="ref_ai", inplace=True)
        if args.custom_c12 is not None:
            custom_c12_dict = io.read_custom_c12_parameters(args.custom_c12)
            d_appo = {key: val for key, val in zip(custom_c12_dict.name, custom_c12_dict.c12)}
            d.update(d_appo)

        topology_df["c12"] = topology_df["mego_type"].map(d)
        first_aminoacid = topology_mego.residues[0].name
        if first_aminoacid in type_definitions.aminoacids_list:
            molecule_type = "protein"
        elif first_aminoacid in type_definitions.nucleic_acid_list:
            molecule_type = "nucleic_acid"
        else:
            molecule_type = "other"

        types = type_definitions.lj14_generator(topology_df)

        if molecule_type == "other":
            # read user pairs
            molecule_keys = list(topology_mego.molecules.keys())
            user_pairs = [
                (pair.atom1.idx, pair.atom2.idx, pair.type.epsilon * 4.184)
                for pair in topology_mego.molecules[molecule_keys[i]][0].adjusts
            ]
            user_pairs = [
                (topology_df[topology_df["mego_ai"] == ai].index[0], topology_df[topology_df["mego_ai"] == aj].index[0], c12)
                for ai, aj, c12 in user_pairs
            ]

        # create Datarame with the pairs and the c12 values
        c12_values = generate_c12_values(topology_df, types, type_definitions.atom_type_combinations, molecule_type)

        # consider special cases
        oxygen_mask = masking.create_matrix_mask(
            topology_df["mego_type"].to_numpy(),
            topology_df["mego_type"].to_numpy(),
            [("OM", "OM"), ("O", "O"), ("OM", "O")],
            symmetrize=True,
        )

        # define all cutoff
        c12_cutoff = CUTOFF_FACTOR * np.power(np.where(oxygen_mask, 11.4 * c12_values, c12_values), 1.0 / 12.0)

        # apply the user pairs (overwrite all other rules)
        if molecule_type == "other":
            for ai, aj, c12 in user_pairs:
                ai = int(ai)
                aj = int(aj)
                if c12 > 0.0:
                    c12_cutoff[ai][aj] = CUTOFF_FACTOR * np.power(c12, 1.0 / 12.0)
                    c12_cutoff[aj][ai] = CUTOFF_FACTOR * np.power(c12, 1.0 / 12.0)

        mismatched = topology_df.loc[topology_df["ref_type"].str[0] != topology_df["mego_name"].str[0]]
        if not mismatched.empty:
            raise ValueError(f"Mismatch found:\n{mismatched}, target and mego topology are not compatible")

        if np.any(c12_cutoff > args.cutoff):
            warning_cutoff_histo(args.cutoff, np.max(c12_cutoff))
        if np.isnan(c12_cutoff.astype(float)).any():
            warning_cutoff_histo(args.cutoff, np.max(c12_cutoff))

        ########################
        # PARALLEL PROCESS START
        ########################

        chunks = np.array_split(target_list, args.proc)
        pool = multiprocessing.Pool(args.proc)
        results = pool.map(
            run_intra_,
            [
                (
                    args,
                    protein_ref_indices,
                    protein_ref_indices,
                    original_size,
                    c12_cutoff,
                    mol_list[i],
                    mol_list[i],
                    x,
                )
                for x in chunks
            ],
        )
        pool.close()
        pool.join()

        ########################
        # PARALLEL PROCESS END
        ########################

        # concatenate and remove partial dataframes
        for name in results:
            try:
                part_df = pd.read_csv(name)
                df = pd.concat([df, part_df])
            except pd.errors.EmptyDataError:
                print(f"Ignoring partial dataframe in {name} as csv is empty")
        [os.remove(name) for name in results]

        df = df.astype(
            {
                "mi": "int32",
                "mj": "int32",
                "ai": "int32",
                "aj": "int32",
            }
        )

        df = df.sort_values(by=["mi", "mj", "ai", "aj"])
        df.drop_duplicates(subset=["mi", "ai", "mj", "aj"], inplace=True)

        df["mi"] = df["mi"].map("{:}".format)
        df["mj"] = df["mj"].map("{:}".format)
        df["ai"] = df["ai"].map("{:}".format)
        df["aj"] = df["aj"].map("{:}".format)
        df["c12dist"] = df["c12dist"].map("{:,.6f}".format)
        df["p"] = df["p"].map("{:,.6e}".format)
        df["cutoff"] = df["cutoff"].map("{:,.6f}".format)

        df.index = range(len(df.index))
        out_name = args.out_name + "_" if args.out_name else ""
        output_file = f"{args.out}/intramat_{out_name}{mol_list[i]}_{mol_list[i]}.ndx.gz"
        print(f"Saving output for molecule {mol_list[i]} in {output_file}")
        write_mat(df, output_file)


def calculate_inter_probabilities(args):
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
    pairs = list(itertools.combinations_with_replacement(mol_list, 2))

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
    columns = ["mi", "ai", "mj", "aj", "c12dist", "p", "cutoff"]
    for pair in pairs:
        df = pd.DataFrame()

        topology_df_i = pd.DataFrame()
        topology_df_j = pd.DataFrame()

        mol_i = pair[0]
        mol_j = pair[1]

        print(
            f"\nCalculating intermat between molecule {mol_i} and {mol_j}: {molecules_name[mol_i-1]} and {molecules_name[mol_j-1]}"
        )
        prefix = f"inter_mol_{mol_i}_{mol_j}"
        if args.tar:
            with tarfile.open(args.histo, "r:*") as tar:
                target_list = [x.name for x in tar.getmembers() if prefix in x.name]
        else:
            target_list = [x for x in os.listdir(args.histo) if prefix in x]

        protein_mego_i = topology_mego.molecules[list(topology_mego.molecules.keys())[mol_i - 1]][0]
        protein_mego_j = topology_mego.molecules[list(topology_mego.molecules.keys())[mol_j - 1]][0]

        protein_ref_i = topology_ref.molecules[list(topology_ref.molecules.keys())[mol_i - 1]][0]
        protein_ref_j = topology_ref.molecules[list(topology_ref.molecules.keys())[mol_j - 1]][0]

        original_size_i = len(protein_ref_i.atoms)
        original_size_j = len(protein_ref_j.atoms)

        if mol_i == mol_j:
            if N_mols[mol_i - 1] == 0:
                print(
                    f"Skipping intermolecular calculation between {mol_i} and {mol_j} cause the number of molecules of this species is only {N_mols[mol_i-1]}"
                )
                matrix_index = pd.MultiIndex.from_product(
                    [range(1, original_size_i + 1), range(1, original_size_j + 1)],
                    names=["ai", "aj"],
                )
                indeces_ai = np.array(list(matrix_index)).T[0]
                indeces_aj = np.array(list(matrix_index)).T[1]
                df = pd.DataFrame(columns=columns)
                df["mi"] = [mol_i for x in range(1, original_size_i * original_size_j + 1)]
                df["mj"] = [mol_j for x in range(1, original_size_i * original_size_j + 1)]
                df["ai"] = indeces_ai
                df["aj"] = indeces_aj
                df["c12dist"] = 0.0
                df["p"] = 0.0
                df["cutoff"] = 0.0
                df["mi"] = df["mi"].map("{:}".format)
                df["mj"] = df["mj"].map("{:}".format)
                df["ai"] = df["ai"].map("{:}".format)
                df["aj"] = df["aj"].map("{:}".format)
                df["c12dist"] = df["c12dist"].map("{:,.6f}".format)
                df["p"] = df["p"].map("{:,.6e}".format)
                df["cutoff"] = df["cutoff"].map("{:,.6f}".format)

                df.index = range(len(df.index))
                out_name = args.out_name + "_" if args.out_name else ""
                output_file = f"{args.out}/intermat_{out_name}{mol_i}_{mol_j}.ndx"

                df.to_csv(output_file, index=False, sep=" ", header=False, columns=columns)
                continue

        protein_ref_indices_i = np.array(
            [i + 1 for i in range(len(protein_ref_i.atoms)) if protein_ref_i[i].element_name != "H"]
        )
        protein_ref_indices_j = np.array(
            [i + 1 for i in range(len(protein_ref_j.atoms)) if protein_ref_j[i].element_name != "H"]
        )

        protein_ref_i = [a for a in protein_ref_i.atoms if a.element_name != "H"]
        protein_ref_j = [a for a in protein_ref_j.atoms if a.element_name != "H"]

        sorter_i = [str(x.residue.number) + map_if_exists(x.name) for x in protein_ref_i]
        sorter_mego_i = [str(x.residue.number) + x.name for x in protein_mego_i]

        sorter_j = [str(x.residue.number) + map_if_exists(x.name) for x in protein_ref_j]
        sorter_mego_j = [str(x.residue.number) + x.name for x in protein_mego_j]

        # preparing topology of molecule i
        topology_df_i["ref_ai"] = protein_ref_indices_i
        topology_df_i["ref_type"] = [a.name for a in protein_ref_i]
        topology_df_i["sorter"] = sorter_i
        topology_df_i["ref_ri"] = topology_df_i["sorter"].str.replace("[a-zA-Z]+[0-9]*", "", regex=True).astype(int)
        topology_df_i.sort_values(by="sorter", inplace=True)
        topology_df_i["mego_type"] = [a[0].type for a in sorted(zip(protein_mego_i, sorter_mego_i), key=lambda x: x[1])]
        topology_df_i["mego_name"] = [a[0].name for a in sorted(zip(protein_mego_i, sorter_mego_i), key=lambda x: x[1])]
        # need to sort back otherwise c12_cutoff are all wrong
        topology_df_i.sort_values(by="ref_ai", inplace=True)
        if args.custom_c12 is not None:
            custom_c12_dict = io.read_custom_c12_parameters(args.custom_c12)
            d_appo = {key: val for key, val in zip(custom_c12_dict.name, custom_c12_dict.c12)}
            d.update(d_appo)

        topology_df_i["c12"] = topology_df_i["mego_type"].map(d)

        # preparing topology of molecule j
        topology_df_j["ref_ai"] = protein_ref_indices_j
        topology_df_j["ref_type"] = [a.name for a in protein_ref_j]
        topology_df_j["sorter"] = sorter_j
        topology_df_j["ref_ri"] = topology_df_j["sorter"].str.replace("[a-zA-Z]+[0-9]*", "", regex=True).astype(int)
        topology_df_j.sort_values(by="sorter", inplace=True)
        topology_df_j["mego_type"] = [a[0].type for a in sorted(zip(protein_mego_j, sorter_mego_j), key=lambda x: x[1])]
        topology_df_j["mego_name"] = [a[0].name for a in sorted(zip(protein_mego_j, sorter_mego_j), key=lambda x: x[1])]
        # need to sort back otherwise c12_cutoff are all wrong
        topology_df_j.sort_values(by="ref_ai", inplace=True)
        if args.custom_c12 is not None:
            custom_c12_dict = io.read_custom_c12_parameters(args.custom_c12)
            d_appo = {key: val for key, val in zip(custom_c12_dict.name, custom_c12_dict.c12)}
            d.update(d_appo)

        topology_df_j["c12"] = topology_df_j["mego_type"].map(d)

        oxygen_mask = masking.create_matrix_mask(
            topology_df_i["mego_type"].to_numpy(),
            topology_df_j["mego_type"].to_numpy(),
            [("OM", "OM"), ("O", "O"), ("OM", "O")],
            symmetrize=True,
        )

        # define all cutoff
        c12_cutoff = CUTOFF_FACTOR * np.where(
            oxygen_mask,
            np.power(
                11.4 * np.sqrt(topology_df_j["c12"].values * topology_df_i["c12"].values[:, np.newaxis]),
                1.0 / 12.0,
            ),
            np.power(
                np.sqrt(topology_df_j["c12"].values * topology_df_i["c12"].values[:, np.newaxis]),
                1.0 / 12.0,
            ),
        )

        mismatched = topology_df_i.loc[topology_df_i["ref_type"].str[0] != topology_df_i["mego_name"].str[0]]
        if not mismatched.empty:
            raise ValueError(f"Mismatch found:\n{mismatched}, target and mego topology are not compatible")
        mismatched = topology_df_j.loc[topology_df_j["ref_type"].str[0] != topology_df_j["mego_name"].str[0]]
        if not mismatched.empty:
            raise ValueError(f"Mismatch found:\n{mismatched}, target and mego topology are not compatible")

        if args.residue:
            c12_cutoff = args.cutoff * np.ones(c12_cutoff.shape)
        if np.any(c12_cutoff > args.cutoff):
            warning_cutoff_histo(args.cutoff, np.max(c12_cutoff))
        if np.isnan(c12_cutoff.astype(float)).any():
            warning_cutoff_histo(args.cutoff, np.max(c12_cutoff))

        # create dictionary with ref_ai to ri
        ref_ai_to_ri_i = dict(zip(topology_df_i["ref_ai"], topology_df_i["ref_ri"]))
        ref_ai_to_ri_j = dict(zip(topology_df_j["ref_ai"], topology_df_j["ref_ri"]))
        # index_ai_to_ri_i = {k: v for k, v in enumerate(topology_df_i["ref_ri"])}
        index_ai_to_ri_j = {k: v for k, v in enumerate(topology_df_j["ref_ri"])}
        # create a dictionary with ref_ri to ai as a list of ai
        ref_ri_to_ai_i = {f"{mol_i}_{ri}": [] for ri in topology_df_i["ref_ri"]}
        ref_ri_to_ai_j = {f"{mol_j}_{ri}": [] for ri in topology_df_j["ref_ri"]}
        for ai, ri in ref_ai_to_ri_i.items():
            ref_ri_to_ai_i[f"{mol_i}_{ri}"].append(ai)
        for ai, ri in ref_ai_to_ri_j.items():
            ref_ri_to_ai_j[f"{mol_j}_{ri}"].append(ai)

        dict_m_m_r = {}
        for target in target_list:
            target_fields = target.replace(".dat", "").split("_")
            mi = int(target_fields[-4])
            mj = int(target_fields[-3])
            ai = int(target_fields[-1])
            if ai not in protein_ref_indices_i:
                continue
            ri = ref_ai_to_ri_i[ai]
            if (mi, mj, ri) in dict_m_m_r:
                dict_m_m_r[(mi, mj, ri)].append(target)
            else:
                dict_m_m_r[(mi, mj, ri)] = [target]

        ########################
        # PARALLEL PROCESS START
        ########################

        if not args.residue:
            chunks = np.array_split(target_list, args.proc)
        else:
            chunks = []
            n_threshold = sum([len(v) for v in dict_m_m_r.values()]) // args.proc
            chunk = []
            n = 0
            for k, v in dict_m_m_r.items():
                chunk.append(v)
                n += len(v)
                if n > n_threshold:
                    chunks.append(chunk)
                    chunk = []
                    n = 0
            chunks.append(chunk)
        pool = multiprocessing.Pool(args.proc)
        if args.residue:
            results = pool.map(
                run_residue_inter_,
                [
                    (
                        args,
                        protein_ref_indices_i,
                        protein_ref_indices_j,
                        len(ref_ri_to_ai_j),
                        c12_cutoff,
                        mol_i,
                        mol_j,
                        (ref_ai_to_ri_i, index_ai_to_ri_j),
                        x,
                    )
                    for x in chunks
                ],
            )
        else:
            results = pool.map(
                run_inter_,
                [
                    (
                        args,
                        protein_ref_indices_i,
                        protein_ref_indices_j,
                        original_size_j,
                        c12_cutoff,
                        mol_i,
                        mol_j,
                        x,
                    )
                    for x in chunks
                ],
            )
        pool.close()
        pool.join()

        ########################
        # PARALLEL PROCESS END
        ########################

        # concatenate and remove partial dataframes
        for name in results:
            try:
                part_df = pd.read_csv(name)
                df = pd.concat([df, part_df])
            except pd.errors.EmptyDataError:
                print(f"Ignoring partial dataframe in {name} as csv is empty")
        [os.remove(name) for name in results]
        df = df.astype({"mi": "int32", "mj": "int32", "ai": "int32", "aj": "int32"})

        df = df.sort_values(by=["mi", "mj", "ai", "aj"])
        df.drop_duplicates(subset=["mi", "ai", "mj", "aj"], inplace=True)

        df["mi"] = df["mi"].map("{:}".format)
        df["mj"] = df["mj"].map("{:}".format)
        df["ai"] = df["ai"].map("{:}".format)
        df["aj"] = df["aj"].map("{:}".format)
        df["c12dist"] = df["c12dist"].map("{:,.6f}".format)
        df["p"] = df["p"].map("{:,.6e}".format)
        df["cutoff"] = df["cutoff"].map("{:,.6f}".format)

        df.index = range(len(df.index))
        out_name = args.out_name + "_" if args.out_name else ""
        output_file = f"{args.out}/intermat_{out_name}{mol_i}_{mol_j}.ndx.gz"
        print(f"Saving output for molecule {mol_i} and {mol_j} in {output_file}")
        write_mat(df, output_file)


def calculate_probability(values, weights, callback=allfunction):
    """
    Calculates a plain probability accoring to \sum_x x * dx

    Parameters
    ----------
    values : np.array
        The array of the histograms x values
    weights : np.array
        The array with the respective weights
    callback : function
        Preprocesses the data before going in to the analysis

    Returns
    -------
    The probability of the histogram
    """
    dx = values[1] - values[0]
    cutoff, i, norm, v, w = callback(values, weights)
    return np.minimum(np.sum(w * dx), 1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--histo",
        type=str,
        required=True,
        help='Path to the directory containing the histograms. The histogram files should contain the prefix "intra_" for intra molecular contact descriptions and "inter_" for  inter molecular.',
    )
    parser.add_argument(
        "--target_top",
        required=True,
        help="Path to the topology file of the system on which the histograms were calculated on",
    )
    parser.add_argument(
        "--mego_top",
        required=True,
        help="""Path to the standard multi-eGO topology of the system generated by pdb2gmx""",
    )
    parser.add_argument(
        "--inter",
        action="store_true",
        help="Sets the caculation to be adapted to intermolecular calculations of histograms",
    )
    parser.add_argument("--out", default="./", help="""Sets the output path""")
    parser.add_argument(
        "--out_name",
        help="""Sets the output name of files to be added to the default one: intermat_<out_name>_mi_mj.ndx or intramat_<out_name>_mi_mj.ndx""",
    )
    parser.add_argument(
        "--proc",
        default=1,
        type=int,
        help="Sets the number of processes to perform the calculation",
    )
    parser.add_argument(
        "--cutoff",
        default=0.75,
        type=float,
        help="To be set to the max cutoff used for the accumulation of the histograms",
    )
    parser.add_argument(
        "--tar",
        action="store_true",
        help="Read from tar file instead of directory",
    )
    parser.add_argument(
        "--custom_c12",
        type=str,
        help="Custom dictionary of c12 for special molecules",
    )
    parser.add_argument(
        "--residue",
        action="store_true",
    )
    args = parser.parse_args()

    # check if output file exists
    if not os.path.exists(args.out):
        print(f"The path '{args.out}' does not exist.")
        sys.exit()

    if not args.tar and not os.path.isdir(args.histo):
        print(f"The path '{args.histo}' is not a directory.")
        sys.exit()

    if args.tar and not tarfile.is_tarfile(args.histo):
        print(f"The path '{args.histo}' is not a tar file.")
        sys.exit()

    if args.residue and not args.inter:
        print("Residue calculation is only possible for intermolecular calculations (not implemented yet for intramolecular).")
        sys.exit()

    N_BINS = args.cutoff / (0.01 / 4)
    DX = args.cutoff / N_BINS
    PREFIX = "inter_mol_" if args.inter else "intra_mol_"
    CUTOFF_FACTOR = 1.45
    print(
        f"""
    Starting with cutoff = {args.cutoff},
                  n_bins = {N_BINS},
                  dx     = {DX}
                  on {args.proc} threads
    """
    )

    if not args.inter:
        calculate_intra_probabilities(args)
    else:
        calculate_inter_probabilities(args)
