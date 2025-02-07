import pandas as pd
import numpy as np
import argparse
import os

# Set up argument parser
parser = argparse.ArgumentParser(description="Read a molecular contact matrix and save it as HDF5.")
parser.add_argument(
    "-i",
    "--input_file",
    required=True,
    type=str,
    help="Path to the input contact matrix file (e.g., 'intramat_1_1.ndx.gz')",
)

# Add an optional argument for output file name
parser.add_argument(
    "-o",
    "--output_file",
    type=str,
    default=None,
    help="Path to the output text file (e.g., 'intramat_1_1.ndx.h5'). Defaults to input file name plus '.h5' extension.",
)

args = parser.parse_args()

if args.output_file is None:
    if args.input_file.endswith(".gz"):
        # Automatically generate the output filename by replacing .gz with .h5
        output_file = os.path.splitext(args.input_file)[0] + ".h5"
    elif args.input_file.endswith(".ndx"):
        output_file = os.path.splitext(args.input_file)[0] + ".ndx.h5"

else:
    output_file = args.output_file

# Define column names and data types
col_names = ["molecule_name_ai", "ai", "molecule_name_aj", "aj", "distance", "probability", "cutoff", "learned"]
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

# Read the input file with specified column names and data types
contact_matrix = pd.read_csv(args.input_file, header=None, sep=r"\s+", names=col_names, dtype=col_types)
contact_matrix["learned"] = contact_matrix["learned"].fillna(1).astype(bool)

contact_matrix["molecule_name_ai"] = contact_matrix["molecule_name_ai"].astype("category")
contact_matrix["ai"] = contact_matrix["ai"].astype("category")
contact_matrix["molecule_name_aj"] = contact_matrix["molecule_name_aj"].astype("category")
contact_matrix["aj"] = contact_matrix["aj"].astype("category")

# Save the data as HDF5 with compression
contact_matrix.to_hdf(output_file, key="data", mode="w", format="table", complib="blosc:lz4", complevel=9)

print(f"Data successfully saved to {output_file}")
