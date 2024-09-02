import sys
import io
import os
import logging
import argparse
import pyarrow
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

logging.basicConfig(
    stream=sys.stdout,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
)

logging.info(f"[{os.path.basename(__file__)}]...")

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Analysis Script")
parser.add_argument('--local_tmp_dir', type=str,
                    required=True, help="Local temporary directory")
parser.add_argument('--local_output_dir', type=str,
                    required=True, help="Local output directory")
args = parser.parse_args()

# Use the provided local directories
local_tmp_dir = args.local_tmp_dir
input_dir = os.path.join(local_tmp_dir, "pyinput")
output_dir = args.local_output_dir
os.chdir(output_dir)

# ---------------------------------------------------------

file_path_feather = os.path.join(input_dir, "adjusted_expres_df.feather")
adjusted_expres_df = pd.read_feather(file_path_feather)
adjusted_expres_df.set_index('rownames', inplace=True)

logging.info("Data loaded.")

# ---------------------------------------------------------

logging.info("Convert the DataFrame to a sparse matrix...")

adjusted_expres_df[adjusted_expres_df < 0.4] = 0

sparse_matrix = sparse.csr_matrix(adjusted_expres_df.values)

# ---------------------------------------------------------

logging.info("Saving data...")

mmwrite(os.path.join(output_dir, 'sparse_matrix.mtx'), sparse_matrix)

adjusted_expres_df.index.to_series().to_csv(os.path.join(output_dir, 'rows.csv'),
                                            index=False, header=False)

adjusted_expres_df.columns.to_series().to_csv(os.path.join(output_dir, 'cols.csv'),
                                              index=False, header=False)

logging.info("Done!")
