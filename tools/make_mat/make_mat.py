import os
import re
import sys
import tempfile
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

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from multiego import bonded
from multiego import type_definitions
from multiego import io

_DEFAULT_C12 = {
    type_definitions.gromos_atp.name[i]: type_definitions.gromos_atp.rc_c12[i]
    for i in range(len(type_definitions.gromos_atp.name))
}

COLUMNS = ["mi", "ai", "mj", "aj", "c12dist", "p", "cutoff"]
CUTOFF_FACTOR = 1.45


def create_matrix_mask(
    set1,
    set2,
    types,
    symmetrize=False,
    inner_op=lambda x, y: x == y,
    outer_op=lambda x, y: x * y,
):
    """
    Creates a 2-D boolean matrix mask based on comparison operations between two sets.

    For each ``(type1, type2)`` pair in ``types``, ``inner_op`` is applied
    element-wise to each set and the results are combined with ``outer_op``
    into a 2-D outer-product mask.  The final mask is the logical OR across
    all type pairs.  When ``symmetrize=True`` each ``(type1, type2)`` entry
    is automatically complemented by ``(type2, type1)``.

    Parameters
    ----------
    set1 : numpy.ndarray
        First 1-D array of values (e.g. atom-type labels), shape ``(N,)``.
    set2 : numpy.ndarray
        Second 1-D array of values (e.g. atom-type labels), shape ``(M,)``.
    types : list of tuple
        Pairs ``(type1, type2)`` used for element-wise comparisons.
    symmetrize : bool, optional
        If ``True``, automatically add the reversed pair ``(type2, type1)``
        for every entry in ``types`` (default ``False``).
    inner_op : callable, optional
        Element-wise comparison applied to each array against a type value
        (default ``lambda x, y: x == y``).
    outer_op : callable, optional
        Combines the two per-set boolean vectors into a 2-D boolean matrix
        (default ``lambda x, y: x * y``, i.e. outer logical AND via
        broadcasting).

    Returns
    -------
    numpy.ndarray
        2-D boolean matrix of shape ``(N, M)``.
    """
    # symmetrize type selection
    if symmetrize:
        types = list(set(types + [(t[1], t[0]) for t in types]))

    mask = np.full((set1.shape[0], set2.shape[0]), False)
    for type1, type2 in types:
        mask |= outer_op(inner_op(set1, type1), inner_op(set2, type2)[:, np.newaxis]).T

    return mask


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
    ordered_columns = [
        "molecule_name_ai",
        "ai",
        "molecule_name_aj",
        "aj",
        "distance",
        "probability",
        "cutoff",
        "learned",
    ]
    df = df[ordered_columns]

    # Save the data as HDF5 with compression
    df.to_hdf(output_file, key="data", mode="w", format="table", complib="blosc:lz4", complevel=9)


def read_mat(name, protein_ref_indices, args, cumulative=False):
    path_prefix = f"{args.histo}"
    if args.tar:
        with tarfile.open(args.histo, "r:*") as tar:
            ref_df = pd.read_csv(tar.extractfile(name), header=None, sep=r"\s+", usecols=[0, *protein_ref_indices])
            ref_df_columns = ["distance", *[str(x) for x in protein_ref_indices]]
            ref_df.columns = ref_df_columns
            ref_df.set_index("distance", inplace=True)
    else:
        if args.noh5:
            ref_df = pd.read_csv(f"{path_prefix}/{name}", header=None, sep=r"\s+", usecols=[0, *protein_ref_indices])
            ref_df_columns = ["distance", *[str(x) for x in protein_ref_indices]]
            ref_df.columns = ref_df_columns
            ref_df.set_index("distance", inplace=True)
        else:
            with h5py.File(f"{path_prefix}/{name}", "r") as f:
                if "density" not in f:
                    raise KeyError(f"Dataset 'density' not found in {name}")

                data = f["density"][:]  # Read full dataset
                # Extract the first column (distance) and the relevant protein_ref_indices columns
                distances = data[:, 0]  # First column is distance
                protein_data = data[:, protein_ref_indices]  # Select the relevant protein reference indices
                # Create a DataFrame
                ref_df = pd.DataFrame(protein_data, columns=[str(i) for i in protein_ref_indices])
                ref_df["distance"] = distances
                # Set 'distance' as the index
                ref_df.set_index("distance", inplace=True)

    return ref_df


def run_mat_(arguments):
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
        mat_type,
    ) = arguments
    process = multiprocessing.current_process()
    df = pd.DataFrame(columns=COLUMNS)
    # We do not consider old histograms
    frac_target_list = [x for x in frac_target_list if x[0] != "#" and x[-1] != "#"]
    for i, ref_f in enumerate(frac_target_list):
        print(f"\rProgress: [{i + 1}/{len(frac_target_list)}] {ref_f} ", end="", flush=True)
        results_df = pd.DataFrame()
        ai = ref_f.split(".")[-2].split("_")[-1]

        all_ai = [ai for _ in range(1, original_size_j + 1)]
        range_list = [str(x) for x in range(1, original_size_j + 1)]

        results_df["ai"] = np.array(all_ai).astype(int)
        results_df["mi"] = mi
        results_df["aj"] = np.array(range_list).astype(int)
        results_df["mj"] = mj
        results_df["c12dist"] = 0.0
        results_df["p"] = 0.0
        results_df["cutoff"] = 0.0

        if np.isin(int(ai), protein_ref_indices_i):
            cut_i = np.where(protein_ref_indices_i == int(ai))[0][0]

            # column mapping
            ref_df = read_mat(ref_f, protein_ref_indices_j, args)
            ref_df.loc[len(ref_df)] = c12_cutoff[cut_i]

            c12dist = ref_df.apply(lambda x: c12_avg(ref_df.index.to_numpy(), weights=x.to_numpy()), axis=0).values
            if mat_type == "intra":
                p = ref_df.apply(
                    lambda x: calculate_probability(ref_df.index.to_numpy(), weights=x.to_numpy()),
                    axis=0,
                ).values
            elif mat_type == "inter":
                # repeat for cumulative
                c_ref_f = ref_f.replace("inter_mol_", "inter_mol_c_")
                c_ref_df = read_mat(c_ref_f, protein_ref_indices_j, args, True)
                c_ref_df.loc[len(c_ref_df)] = c12_cutoff[cut_i]
                p = c_ref_df.apply(
                    lambda x: get_cumulative_probability(c_ref_df.index.to_numpy(), weights=x.to_numpy()),
                    axis=0,
                ).values

            results_df.loc[results_df["aj"].isin(protein_ref_indices_j), "c12dist"] = c12dist
            results_df.loc[results_df["aj"].isin(protein_ref_indices_j), "p"] = p
            results_df.loc[results_df["aj"].isin(protein_ref_indices_j), "cutoff"] = c12_cutoff[cut_i].astype(float)

        if df.empty:
            df = results_df.copy()
        else:
            if not results_df.empty:
                df = pd.concat([df, results_df])

    print("done.", flush=True)
    df = df.fillna(0)
    out_path = f"mat_{process.pid}_t{time.time()}.part"
    df.to_csv(out_path, index=False)

    return out_path


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
        Path to the multi-eGO topology obtained from gmx pdb2gmx with multi-ego-basic force fields.
    target_top : str
        Path to the topology of the system on which the analysis is to be performed.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            topology_mego = pmd.load_file(mego_top)
        except Exception as e:
            print(f"ERROR {e} in read_topologies while reading {mego_top}")
            sys.exit(1)
        try:
            dirname, basename = os.path.split(target_top)
            temp_ref = tempfile.NamedTemporaryFile(prefix=basename, dir=dirname)
            with open(target_top, "rb") as src:
                temp_ref.write(src.read())
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
            sys.exit(1)

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

    Returns
    -------
    atom_name : str
        Mapped atom name. Equal to the input if mapping was not possible.
    """
    if atom_name in type_definitions.from_ff_to_multiego.keys():
        return type_definitions.from_ff_to_multiego[atom_name]
    else:
        return atom_name


def get_col_params(values, weights):
    """
    Preprocesses histogram arrays for analysis. The last element of ``weights``
    encodes the cutoff; both arrays are truncated to only include bins at or
    below that cutoff.

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
    if i[0].size == 0:
        return 0, 0, 0, 0, 0  # check if empty
    i = i[0]
    w = w[i]
    v = v[i]
    norm = np.sum(w)
    i = i[-1]
    return cutoff, i, norm, v, w


def calculate_probability(values, weights):
    """
    Calculates a plain probability according to sum_x x * dx

    Parameters
    ----------
    values : np.array
        The array of the histograms x values
    weights : np.array
        The array with the respective weights

    Returns
    -------
    The probability of the histogram
    """
    dx = values[1] - values[0]
    cutoff, i, norm, v, w = get_col_params(values, weights)
    return np.minimum(np.sum(w * dx), 1)


def get_cumulative_probability(values, weights):
    """
    Returns the cumulative probability at the cutoff index.

    Parameters
    ----------
    values : numpy.ndarray
        The array of histogram x values.
    weights : numpy.ndarray
        The array of histogram weights; the last element encodes the cutoff.

    Returns
    -------
    float
        The cumulative weight at the bin index corresponding to the cutoff.
    """
    cutoff, i, norm, v, w = get_col_params(values, weights)
    return weights[i]


def c12_avg(values, weights):
    """
    Calculates the c12 exp averaging of a histogram

    Parameters
    ----------
    values : np.array
        The array of the histograms x values
    weights : np.array
        The array with the respective weights

    Returns
    -------
    The c12 average
    """
    cutoff, i, norm, v, w = get_col_params(values, weights)
    if np.sum(w) == 0:
        return 0
    r = np.where(w > 0.0)
    # fmt: off
    v = v[r[0][0]:v.size]
    w = w[r[0][0]:w.size]
    # fmt: on
    res = np.maximum(cutoff / 4.5, 0.1)
    log_exp_sum = logsumexp(1.0 / v / res, b=w) - np.log(norm)
    exp_aver = (1.0 / res) / log_exp_sum
    if exp_aver < 0.01:
        exp_aver = 0

    return exp_aver


def warning_cutoff_histo(cutoff, max_adaptive_cutoff):
    """
    Prints a warning when the adaptive cutoff exceeds the histogram cutoff.

    Parameters
    ----------
    cutoff : float
        The cutoff used during histogram accumulation.
    max_adaptive_cutoff : float
        The maximum adaptive c12 cutoff derived from LJ parameters.
    """
    print(
        f"\nWARNING: adaptive cutoff ({max_adaptive_cutoff:.4f} nm) exceeds histogram cutoff ({cutoff} nm).\n"
        f"  This may cause errors. Consider re-running histogram accumulation with --cutoff {max_adaptive_cutoff:.4f}\n"
    )


def generate_c12_values(df, types, combinations, molecule_type):
    """
    Compute the pairwise c12 cutoff matrix for a molecule.

    For protein molecules, per-combination rules from ``combinations`` are
    applied on top of the default geometric-mean c12 values.  Each rule
    specifies two atom-group names, an optional ``factor`` (scales the
    geometric mean), an optional ``constant`` (fixed ceiling), and a residue
    index ``shift``.  Non-protein molecules use the geometric mean directly.

    Parameters
    ----------
    df : pd.DataFrame
        Per-atom topology slice.  Must contain columns ``c12`` and ``resnum``.
    types : dict
        Mapping from group name to a 1-D boolean numpy array (as returned by
        ``bonded.proteins_atoms_mask``).
    combinations : list of tuple
        Each entry is ``(name_1, name_2, factor, constant, shift)`` defining
        a pairwise c12 override rule.
    molecule_type : str
        ``'protein'``, ``'nucleic_acid'``, or ``'other'``.  Rules in
        ``combinations`` are only applied when this is ``'protein'``.

    Returns
    -------
    numpy.ndarray
        2-D array of shape ``(N, N)`` containing the pairwise c12 values,
        where ``N`` is the number of atoms in ``df``.
    """
    all_c12 = np.sqrt(df["c12"].to_numpy() * df["c12"].to_numpy()[:, np.newaxis])
    c12_map = np.full(all_c12.shape, np.nan, dtype=float)
    resnums = df["resnum"].to_numpy()

    if molecule_type == "protein":
        for combination in combinations:
            name_1, name_2, factor, constant, shift = combination
            if factor is not None and constant is not None:

                def operation(x, _f=factor, _c=constant):
                    return np.minimum(_f * x, _c)

            elif factor is not None:

                def operation(x, _f=factor):
                    return _f * x

            elif constant is not None:

                def operation(_, _c=constant):
                    return _c

            else:
                raise ValueError("Either factor or constant must be specified.")

            combined_map = (types[name_1] & types[name_2][:, np.newaxis]) & (resnums + shift == resnums[:, np.newaxis])
            combined_map = combined_map | combined_map.T
            c12_map = np.where(combined_map, operation(all_c12), c12_map)

    c12_map = np.where(np.isnan(c12_map), all_c12, c12_map)

    return c12_map


def calculate_matrices(args):
    """
    Starts the main routine for calculating the intermat by:
     - reading the topologies
     - figuring out all the interacting molecules
     - calculating the cutoffs
     - and calculating the probabilities
    The operation is finalized by writing out a csv with the name pattern intermat<_name>_{mol_i}_{mol_j}.ndx

    Parameters
    ----------
    args : dict
        The command-line parsed parameters
    """
    topology_mego, topology_ref, N_species, molecules_name, mol_list = read_topologies(args.mego_top, args.target_top)

    print(f"Topology contains {N_species} molecule species: {molecules_name}")
    for mol_i in mol_list:
        if args.intra:
            prefix = f"intra_mol_{mol_i}_{mol_i}"
            main_routine(mol_i, mol_i, topology_mego, topology_ref, molecules_name, prefix, args)
        # fmt: off
        for mol_j in mol_list[mol_i - 1:]:
            # fmt: on
            if mol_i == mol_j and not args.same:
                continue
            if mol_i != mol_j and not args.cross:
                continue

            prefix = f"inter_mol_{mol_i}_{mol_j}"
            main_routine(mol_i, mol_j, topology_mego, topology_ref, molecules_name, prefix, args)


def main_routine(mol_i, mol_j, topology_mego, topology_ref, molecules_name, prefix, args):
    """
    Compute and write the intra- or inter-molecular contact matrix for a
    pair of molecule species.

    Builds per-atom topology DataFrames for both molecule species, computes
    adaptive c12 cutoffs, dispatches histogram analysis to worker processes
    via ``run_mat_``, reassembles the partial results, and writes the output
    as a compressed HDF5 file.

    Parameters
    ----------
    mol_i : int
        1-based index of the first molecule species.
    mol_j : int
        1-based index of the second molecule species.
    topology_mego : parmed.Structure
        multi-eGO topology (defines atom types and bonded terms).
    topology_ref : parmed.Structure
        Reference all-atom topology (defines the histogram atom indices).
    molecules_name : list of str
        Ordered list of molecule species names from the topology.
    prefix : str
        File-name prefix for histogram look-up and output naming,
        e.g. ``"intra_mol_1_1"`` or ``"inter_mol_1_2"``.
    args : argparse.Namespace
        Parsed command-line arguments.
    """

    topology_df_i = pd.DataFrame()
    topology_df_j = pd.DataFrame()

    # define matrix type (intra o inter)
    mat_type = prefix.split("_")[0]
    print(
        f"\nCalculating {mat_type} between molecule {mol_i} and {mol_j}: {molecules_name[mol_i-1]} and {molecules_name[mol_j-1]}"
    )
    if args.tar:
        with tarfile.open(args.histo, "r:*") as tar:
            target_list = [x.name for x in tar.getmembers() if prefix in x.name and x.name.endswith(".dat")]
    else:
        if args.noh5:
            target_list = [x for x in os.listdir(args.histo) if prefix in x and x.endswith(".dat")]
        else:
            target_list = [x for x in os.listdir(args.histo) if prefix in x and x.endswith(".h5")]

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
            if (protein_ref_i[i].element_name != "H" or protein_ref_i[i].name in {"H", "HN"} or protein_ref_i[i].name == args.bkbn_H)
        ]
    )
    protein_ref_indices_j = np.array(
        [
            i + 1
            for i in range(len(protein_ref_j.atoms))
            if (protein_ref_j[i].element_name != "H" or protein_ref_j[i].name in {"H", "HN"} or protein_ref_j[i].name == args.bkbn_H)
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

    protein_ref_i = [a for a in protein_ref_i.atoms if (a.element_name != "H" or a.name in {"H", "HN"} or (args.bkbn_H and a.name == args.bkbn_H))]
    protein_ref_j = [a for a in protein_ref_j.atoms if (a.element_name != "H" or a.name in {"H", "HN"} or (args.bkbn_H and a.name == args.bkbn_H))]

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
    c12_lookup = dict(_DEFAULT_C12)
    if args.custom_c12 is not None:
        custom_c12_dict = io.read_custom_c12_parameters(args.custom_c12)
        c12_lookup.update(zip(custom_c12_dict.name, custom_c12_dict.rc_c12))

    topology_df_i["c12"] = topology_df_i["mego_type"].map(c12_lookup)

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
    topology_df_j["c12"] = topology_df_j["mego_type"].map(c12_lookup)

    type_i = topology_df_i["mego_type"].to_numpy()
    type_j = topology_df_j["mego_type"].to_numpy()

    if mat_type == "intra":
        first_aminoacid = topology_mego.residues[0].name
        if first_aminoacid in type_definitions.aminoacids_list:
            molecule_type = "protein"
        elif first_aminoacid in type_definitions.nucleic_acid_list:
            molecule_type = "nucleic_acid"
        else:
            molecule_type = "other"

        types = bonded.proteins_atoms_mask(topology_df_i)

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
        c12_matrix = c12_values.copy()

        # special repulsive rules
        for rule in type_definitions.special_non_local:
            if rule["interaction"] != "rep":
                continue
            if rule["epsilon"] is None:
                continue

            types_i, types_j = rule["atomtypes"]
            pair_list = [(a, b) for a in types_i for b in types_j]
            mask = create_matrix_mask(
                type_i,
                type_j,
                pair_list,
                symmetrize=True,
            )
            c12_matrix = np.where(mask, rule["epsilon"], c12_matrix)

        c12_cutoff = CUTOFF_FACTOR * np.power(c12_matrix, 1.0 / 12.0)

        # apply the user pairs (overwrite all other rules)
        if molecule_type == "other":
            for ai, aj, c12 in user_pairs:
                ai = int(ai)
                aj = int(aj)
                if c12 > 0.0:
                    c12_cutoff[ai][aj] = CUTOFF_FACTOR * np.power(c12, 1.0 / 12.0)
                    c12_cutoff[aj][ai] = CUTOFF_FACTOR * np.power(c12, 1.0 / 12.0)

    if mat_type == "inter":
        c12_matrix = np.sqrt(topology_df_j["c12"].values * topology_df_i["c12"].values[:, np.newaxis])
        # special repulsive rules
        for rule in type_definitions.special_non_local:
            if rule["interaction"] != "rep":
                continue
            if rule["epsilon"] is None:
                continue

            types_i, types_j = rule["atomtypes"]
            pair_list = [(a, b) for a in types_i for b in types_j]
            mask = create_matrix_mask(
                type_i,
                type_j,
                pair_list,
                symmetrize=True,
            )
            c12_matrix = np.where(mask, rule["epsilon"], c12_matrix)

        c12_cutoff = CUTOFF_FACTOR * np.power(c12_matrix, 1.0 / 12.0)

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

    ########################
    # PARALLEL PROCESS START
    ########################

    if args.zero:
        df = pd.DataFrame()
        all_ai = [i for i in range(1, original_size_i + 1)]
        all_aj = [j for j in range(1, original_size_j + 1)]
        df["mi"] = [mol_i for _ in range((original_size_i) * (original_size_j))]
        df["mj"] = [mol_j for _ in range((original_size_i) * (original_size_j))]
        df["ai"] = np.repeat(all_ai, (original_size_j))
        df["aj"] = np.tile(all_aj, original_size_i)
        df["c12dist"] = 0.0
        df["p"] = 0.0
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
        df["cutoff"] = cuts
    else:
        chunks = np.array_split(target_list, args.num_threads)
        pool = multiprocessing.Pool(args.num_threads)
        results = pool.map(
            run_mat_,
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
                    mat_type,
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
        parts = []
        for name in results:
            try:
                parts.append(pd.read_csv(name))
            except pd.errors.EmptyDataError:
                print(f"Ignoring partial dataframe in {name} as csv is empty")
        df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
        for name in results:
            os.remove(name)
        df = df.astype({"mi": "int32", "mj": "int32", "ai": "int32", "aj": "int32"})

        df = df.sort_values(by=["mi", "mj", "ai", "aj"])
        df.drop_duplicates(subset=["mi", "ai", "mj", "aj"], inplace=True)

    df["learned"] = 1

    if mol_i == mol_j and args.nb:
        df["ai"] = df["ai"].astype(int)
        df["aj"] = df["aj"].astype(int)
        df["residue_diff"] = df.apply(
            lambda row: abs(d_ref_ai_to_ri_i[row["ai"]] - d_ref_ai_to_ri_i[row["aj"]]), axis=1
        )
        df.loc[df["residue_diff"] > 2, ["p", "c12dist", "learned"]] = 0
        df.drop(columns=["residue_diff"], inplace=True)

    if mol_i == mol_j and args.nonb:
        df["ai"] = df["ai"].astype(int)
        df["aj"] = df["aj"].astype(int)
        df["residue_diff"] = df.apply(
            lambda row: abs(d_ref_ai_to_ri_i[row["ai"]] - d_ref_ai_to_ri_i[row["aj"]]), axis=1
        )
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
    output_file = os.path.join(args.out, f"{mat_type}mat_{out_name}{mol_i}_{mol_j}.ndx.h5")
    print(f"Saving output for molecule {mol_i} and {mol_j} in {output_file}")
    write_mat(df, output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--histo",
        type=str,
        required=False,
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
        "--mode",
        help="Sets the calculation to be intra/same/cross for histograms processing",
        default="intra+same+cross",
    )
    parser.add_argument("--bkbn_H", help="Extra backbone H name to include beyond H and HN (e.g. a force-field-specific variant)", default="")
    parser.add_argument("--out", default="./", help="""Sets the output path""")
    parser.add_argument(
        "--out_name",
        help="""Sets the output name of files to be added to the default one: intermat_<out_name>_mi_mj.ndx or intramat_<out_name>_mi_mj.ndx""",
    )
    parser.add_argument(
        "--num_threads",
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
        "--tar",
        action="store_true",
        help="Read from tar file instead of directory",
    )
    parser.add_argument(
        "--noh5",
        action="store_true",
        help="Read from text file instead of hdf5",
    )
    parser.add_argument(
        "--custom_c12",
        type=str,
        help="Custom dictionary of c12 for special molecules",
    )
    parser.add_argument(
        "--zero",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()

    # check either histo or zero flag are set
    if not args.histo and not args.zero:
        raise ValueError("Either --histo or --zero flag must be set.")
    if args.histo and args.zero:
        raise ValueError("Both --histo and --zero flags cannot be set at the same time.")

    # check if output directory exists
    if not os.path.exists(args.out):
        print(f"ERROR: output path '{args.out}' does not exist.")
        sys.exit(1)

    if not args.zero and not args.tar:
        if not os.path.isdir(args.histo):
            print(f"ERROR: histogram path '{args.histo}' is not a directory.")
            sys.exit(1)

    if not args.zero and args.tar:
        if not tarfile.is_tarfile(args.histo):
            print(f"ERROR: histogram path '{args.histo}' is not a tar file.")
            sys.exit(1)

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

    print(f"\nStarting with cutoff = {args.cutoff} nm on {args.num_threads} thread(s)\n")

    calculate_matrices(args)
