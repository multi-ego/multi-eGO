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
        raise ValueError("Elements in each range should be non-decreasing e.g. dom_res 1-10 11-20 ...")

    if not all([x1[1] < x2[0] for x1, x2 in zip(doms[:-1], doms[1:])]):
        raise ValueError("Ranges should not overlap e.g. dom_res 1-10 11-20 ...")

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
    intra_md = np.loadtxt(intramat, unpack=True)
    dim = int(np.sqrt(len(intra_md[0])))
    if dim != n_atoms:
        raise ValueError(f"ERROR: number of atoms in intramat ({dim}) does not correspond to that of topology ({n_atoms})")

    # define domain mask
    domain_mask_linear = np.full(dim**2, False)
    for r in ranges:
        start = find_atom_start(topology_mego, r[0])
        end = find_atom_end(topology_mego, r[1])
        print(f"  Domain range: {r[0]}-{r[1]}")
        print(f"     Atom index range start-end: {start+1} - {end+1}")
        print(f"     Number of atoms in domain range:  {end+1 - (start)}")
        print(f"     Atom and Residue of start-end {topology_mego.atoms[start]} - {topology_mego.atoms[end]}")
        print("\n")
        map_appo = np.array([True if x >= start and x <= end else False for x in range(dim)])
        domain_mask_linear = np.logical_or(domain_mask_linear, (map_appo * map_appo[:, np.newaxis]).reshape(dim**2))
    # domain_mask_linear = (domain_mask * domain_mask[:, np.newaxis]).reshape(dim**2)
    if args.invert:
        domain_mask_linear = np.logical_not(domain_mask_linear)
    print(domain_mask_linear)
    # add an eigth column with the domain_mask
    if intra_md.shape[0] == 7:
        intra_md = np.concatenate((intra_md, domain_mask_linear[np.newaxis, :]), axis=0)
    else:
        intra_md[7] = domain_mask_linear

    if "/" in intramat:
        intramat = intramat.split("/")[-1]

    if args.invert:
        out_name = f'{args.out}inverted_split_{"-".join(np.array(args.dom_res, dtype=str))}_{intramat}'
    else:
        out_name = f'{args.out}split_{"-".join(np.array(args.dom_res, dtype=str))}_{intramat}'
    np.savetxt(
        out_name,
        intra_md.T,
        delimiter=" ",
        fmt=["%i", "%i", "%i", "%i", "%2.6f", "%.6e", "%2.6f", "%1i"],
    )
    print("Finished splitting")
