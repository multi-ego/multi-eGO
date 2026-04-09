import os
import sys
import argparse

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Convert a plain-text (or .gz) contact matrix to a compressed HDF5 file."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=str,
        help="Path to the input contact matrix file (e.g., 'intramat_1_1.ndx' or 'intramat_1_1.ndx.gz')",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help=(
            "Path to the output HDF5 file. "
            "Defaults to the input path with '.gz' replaced by '.h5', or '.h5' appended."
        ),
    )
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"ERROR: input file '{args.input}' not found.")
        sys.exit(1)

    if args.output is not None:
        output_file = args.output
    elif args.input.endswith(".gz"):
        output_file = os.path.splitext(args.input)[0] + ".h5"
    else:
        # Covers .ndx and any other plain-text variant
        output_file = args.input + ".h5"

    col_names = [
        "molecule_name_ai",
        "ai",
        "molecule_name_aj",
        "aj",
        "distance",
        "probability",
        "cutoff",
        "learned",
    ]
    col_types = {
        "molecule_name_ai": str,
        "ai": str,
        "molecule_name_aj": str,
        "aj": str,
        "distance": np.float64,
        "probability": np.float64,
        "cutoff": np.float64,
        "learned": "Int64",
    }

    contact_matrix = pd.read_csv(args.input, header=None, sep=r"\s+", names=col_names, dtype=col_types)
    contact_matrix["learned"] = contact_matrix["learned"].fillna(1).astype(bool)
    contact_matrix["molecule_name_ai"] = contact_matrix["molecule_name_ai"].astype("category")
    contact_matrix["ai"] = contact_matrix["ai"].astype("category")
    contact_matrix["molecule_name_aj"] = contact_matrix["molecule_name_aj"].astype("category")
    contact_matrix["aj"] = contact_matrix["aj"].astype("category")

    contact_matrix.to_hdf(output_file, key="data", mode="w", format="table", complib="blosc:lz4", complevel=9)
    print(f"Data successfully saved to {output_file}")


if __name__ == "__main__":
    main()
