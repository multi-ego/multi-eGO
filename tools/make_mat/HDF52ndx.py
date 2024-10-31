import pandas as pd
import numpy as np
import argparse
import os

# Set up argument parser
parser = argparse.ArgumentParser(description="Read a HDF5 molecular contact matrix and save it as text file.")
parser.add_argument(
    "-i",
    "--input_file",
    type=str,
    help="Path to the input contact matrix file (e.g., 'intramat_1_1.ndx.h5')",
)

# Add an optional argument for output file name
parser.add_argument(
    "-o",
    "--output_file",
    type=str,
    default=None,
    help="Path to the output text file (e.g., 'intramat_1_1.ndx'). Defaults to input file name without '.h5' extension.",
)

args = parser.parse_args()

if args.output_file is None:
    output_file = os.path.splitext(args.input_file)[0]
else:
    output_file = args.output_file

# Read the HDF5 file
contact_matrix = pd.read_hdf(args.input_file, key="data")

# Define the column names if needed (ensure they are consistent with HDF5)
col_names = ["molecule_name_ai", "ai", "molecule_name_aj", "aj", "distance", "probability", "cutoff", "intra_domain"]

# Check if the columns exist in the DataFrame
for col in col_names:
    if col not in contact_matrix.columns:
        raise ValueError(f"Column '{col}' not found in the input HDF5 file.")

contact_matrix["intra_domain"] = contact_matrix["intra_domain"].astype(int)

# Save the DataFrame as a text file with space-separated values
contact_matrix.to_csv(output_file, sep=" ", index=False, header=False)

print(f"Data successfully saved to {output_file}")
