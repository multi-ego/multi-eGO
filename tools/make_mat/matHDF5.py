import pandas as pd
import numpy as np
import argparse
import os

# Set up argument parser
parser = argparse.ArgumentParser(description="Read a molecular contact matrix and save it as HDF5.")
parser.add_argument(
    "input_file",
    type=str,
    help="Path to the input contact matrix file (e.g., 'intramat_1_1.ndx.gz')",
)

args = parser.parse_args()

# Automatically generate the output filename by replacing .gz with .h5
output_file = os.path.splitext(args.input_file)[0] + ".h5"

# Define column names and data types
col_names = ["molecule_name_ai", "ai", "molecule_name_aj", "aj", "distance", "probability", "cutoff", "intra_domain"]
col_types = {
    "molecule_name_ai": str,
    "ai": str,
    "molecule_name_aj": str,
    "aj": str,
    "distance": np.float32,
    "probability": np.float32,
    "cutoff": np.float32,
    "intra_domain": "Int64"  # Allows for integer with NaNs, which can be cast later
}

# Read the input file with specified column names and data types
contact_matrix = pd.read_csv(args.input_file, header=None, sep="\s+", names=col_names, dtype=col_types)
contact_matrix["intra_domain"] = contact_matrix["intra_domain"].fillna(1).astype(bool)

# Save the data as HDF5 with compression
contact_matrix.to_hdf(
    output_file,
    key='data',
    mode='w',
    complib='blosc:lz4',
    complevel=9
)

print(f"Data successfully saved to {output_file}")
