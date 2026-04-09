import os
import sys
import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Convert a compressed HDF5 contact matrix to a plain-text file.")
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=str,
        help="Path to the input HDF5 contact matrix file (e.g., 'intramat_1_1.ndx.h5')",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help=("Path to the output text file. " "Defaults to the input path with the '.h5' extension removed."),
    )
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"ERROR: input file '{args.input}' not found.")
        sys.exit(1)

    output_file = args.output if args.output is not None else os.path.splitext(args.input)[0]

    contact_matrix = pd.read_hdf(args.input, key="data")

    expected_columns = [
        "molecule_name_ai",
        "ai",
        "molecule_name_aj",
        "aj",
        "distance",
        "probability",
        "cutoff",
        "learned",
    ]
    for col in expected_columns:
        if col not in contact_matrix.columns:
            raise ValueError(f"Column '{col}' not found in the input HDF5 file.")

    contact_matrix["learned"] = contact_matrix["learned"].astype(int)
    contact_matrix.to_csv(output_file, sep=" ", index=False, header=False)
    print(f"Data successfully saved to {output_file}")


if __name__ == "__main__":
    main()
