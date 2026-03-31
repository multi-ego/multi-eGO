import numpy as np

# import sys
import argparse
import os
import parmed as pmd
import warnings
import pandas as pd


def find_atom_start(top, res_num):
    """
    Finds the starting atom associated to the residue
    """
    atom_idx = 0

    for i in range(res_num - 1):
        atom_idx += len(top.residues[i].atoms)

    return atom_idx


def find_atom_end(top, res_num):
    """
    Finds the ending atom associated to the residue
    """
    atom_idx = 0

    for i in range(res_num):
        atom_idx += len(top.residues[i].atoms)

    return atom_idx - 1


def dom_range(ranges_str):
    """
    Reads the ranges given in input as a string and puts them in output
    as a list of tuples checking that the ranges are non-decreasing and non-overlapping
    """

    print("\nReading domain ranges in which inserting intramats")
    doms = [(int(r.split("-")[0]), int(r.split("-")[1])) for r in ranges_str]

    if not all([x[0] <= x[1] for x in doms]):
        print("WARNING: Elements in each range should be non-decreasing e.g. dom_res 1-10 11-20 ...")

    if not all([x1[1] < x2[0] for x1, x2 in zip(doms[:-1], doms[1:])]):
        print("WARNING: Ranges should not overlap e.g. dom_res 1-10 11-20 ...")

    return doms


# TODO should re-use multiego reading topology function
def read_topologies(top):
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
        topology = pmd.load_file(top)

    # Return topology and a dataframe with:
    # molecule name, number of molecules?, residue list, atom_list_per_residue
    top_df = pd.DataFrame()
    # n_mol = len(list(topology.molecules.keys()))
    mol_names = list(topology.molecules.keys())
    top_df["name"] = mol_names
    # mol_list = np.arange(1, n_mol + 1, 1)
    res = []
    atoms = []
    atoms_name = []
    tot_num_atoms = []
    for name in mol_names:
        res.append([r.name for r in topology.molecules.values().mapping[name][0].residues])
        atoms.append([len(r.atoms) for r in topology.molecules.values().mapping[name][0].residues])
        atoms_name.append([[a.type for a in r.atoms] for r in topology.molecules.values().mapping[name][0].residues])
        tot_num_atoms.append(np.sum(np.array([len(r.atoms) for r in topology.molecules.values().mapping[name][0].residues])))
    top_df["residues"] = res
    top_df["atoms_per_res"] = atoms
    top_df["tot_atoms"] = tot_num_atoms
    top_df["atoms_name"] = atoms_name

    return topology, top_df


def readmat(intramat, h5=True):
    """
    Reads the intramat and checks that it has the correct format. Returns the intramat as a dataframe with correct column types.
    """
    if h5:
        intramat_md_df = pd.read_hdf(intramat, key="data")
    else:
        intramat_md = np.loadtxt(intramat, unpack=True)
        if intramat_md.shape[0] not in [7, 8]:
            raise ValueError(
                "Intramat should have 7 or 8 columns: molecule_name_ai, ai, molecule_name_aj, aj, distance, probability, cutoff, (optional) learned"
            )
        if intramat_md.shape[0] == 7:
            print("Intramat has 7 columns, domain mask will be added as an eigth column")
            intramat_md = np.concatenate((intramat_md, np.full(intramat_md.shape[1], False)[np.newaxis, :]), axis=0)

        col_types = {
            "molecule_name_ai": int,
            "ai": int,
            "molecule_name_aj": int,
            "aj": int,
            "distance": np.float64,
            "probability": np.float64,
            "cutoff": np.float64,
            "learned": int,  # Allows for integer with NaNs, which can be cast later
        }

        intramat_md_df = pd.DataFrame(columns=col_types.keys())
        for i, col in enumerate(col_types.keys()):
            intramat_md_df[col] = intramat_md[i].astype(col_types[col])

    return intramat_md_df


def write_intramat(contact_matrix, out_name, h5=False):
    """
    Writes the intramat in the correct format. If the intramat has 8 columns, it is assumed that the last column is the learned mask.
    """
    col_names = ["molecule_name_ai", "ai", "molecule_name_aj", "aj", "distance", "probability", "cutoff", "learned"]

    if not h5:
        col_types = {
            "molecule_name_ai": str,
            "ai": str,
            "molecule_name_aj": str,
            "aj": str,
            "distance": np.float64,
            "probability": np.float64,
            "cutoff": np.float64,
            "learned": "Int64",  # Allows for integer with NaNs, which can be cast later
        }
        contact_matrix = pd.DataFrame(contact_matrix.T, columns=col_names, dtype=col_types)

    # Read the input file with specified column names and data types
    contact_matrix["learned"] = contact_matrix["learned"].fillna(1).astype(bool)

    contact_matrix["molecule_name_ai"] = contact_matrix["molecule_name_ai"].astype("category")
    contact_matrix["ai"] = contact_matrix["ai"].astype("category")
    contact_matrix["molecule_name_aj"] = contact_matrix["molecule_name_aj"].astype("category")
    contact_matrix["aj"] = contact_matrix["aj"].astype("category")

    # Save the data as HDF5 with compression
    contact_matrix.to_hdf(out_name, key="data", mode="w", format="table", complib="blosc:lz4", complevel=9)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="TODO!")
    parser.add_argument("--intra", type=str, required=True, help="intramat to work on")
    parser.add_argument("--top", type=str, required=True)
    parser.add_argument(
        "--dom_res",
        nargs="+",
        type=str,
        default=[],
        help="list of residue indeces associated to the starting residues of new domains. Example: 15 44 ...",
        required=True,
    )
    parser.add_argument("--out", type=str, default="./", help="path for ouput")
    parser.add_argument("--invert", action="store_true", default=False, help="Inbert domain mask")

    args = parser.parse_args()

    if args.out:
        if not os.path.isdir(args.out):
            print(f"{args.out} does not exists. Insert an existing directory")
            exit()
        else:
            if args.out[-1] != "/":
                args.out = args.out + "/"

    # read topology
    topology_mego, top_df = read_topologies(args.top)

    # check if there is only one molecule. This code should modify only intramat of one molecule
    if len(top_df) > 1:
        raise ValueError("Only one molecule specie is allowed, topology contains more than one molecule")

    # define atom_num and res_num of the molecule
    n_atoms = top_df.tot_atoms[0]
    n_res = len(top_df.residues[0])

    ranges = dom_range(args.dom_res)

    print(f"\n Total number of residues {n_res} and total number of atoms {n_atoms} \n")

    # read intramat and check consistency
    intramat = args.intra
    if intramat.endswith(".h5"):
        print(f"Reading intramat {intramat} as HDF5")
        intra_md = readmat(intramat, h5=True)
    else:
        intra_md = readmat(intramat, h5=False)
    dim = int(np.sqrt(len(intra_md)))
    if dim != n_atoms:
        raise ValueError(f"ERROR: number of atoms in intramat ({dim}) does not correspond to that of topology ({n_atoms})")

    # define domain mask
    domain_mask_linear = np.full(dim**2, False)
    for r in ranges:
        start = find_atom_start(topology_mego, r[0])
        end = find_atom_end(topology_mego, r[1])
        if start >= end:
            appo_end = end
            end = start
            start = appo_end
            print(f"  Domain range: {r[0]}-{r[1]} INVERTED")
            print(f"     Atom index range start-end: {start+1} - {end+1}")
            print(f"     Number of atoms in domain range:  {end+1 - (start)}")
            print(f"     Atom and Residue of start-end {topology_mego.atoms[start]} - {topology_mego.atoms[end]}")
            print("\n")
            map_appo = np.invert(np.array([True if x >= start and x <= end else False for x in range(dim)]))
        else:
            print(f"  Domain range: {r[0]}-{r[1]}")
            print(f"     Atom index range start-end: {start+1} - {end+1}")
            print(f"     Number of atoms in domain range:  {end+1 - (start)}")
            print(f"     Atom and Residue of start-end {topology_mego.atoms[start]} - {topology_mego.atoms[end]}")
            print("\n")
            map_appo = np.array([True if x >= start and x <= end else False for x in range(dim)])
        domain_mask_linear = np.logical_or(domain_mask_linear, (map_appo * map_appo[:, np.newaxis]).reshape(dim**2))

    if args.invert:
        domain_mask_linear = np.logical_not(domain_mask_linear)
    print(domain_mask_linear)

    # set domain mask in intramat
    intra_md["learned"] = domain_mask_linear

    if "/" in intramat:
        intramat = intramat.split("/")[-1]
        if intramat.endswith(".h5"):
            intramat = intramat[:-3]
    if args.invert:
        out_name = f'{args.out}inverted_split_{"-".join(np.array(args.dom_res, dtype=str))}_{intramat}.h5'
    else:
        out_name = f'{args.out}split_{"-".join(np.array(args.dom_res, dtype=str))}_{intramat}.h5'

    # write output in h5 format
    write_intramat(intra_md, out_name, h5=True)

    print("Finished splitting")
