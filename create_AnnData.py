import sys
import io
import os
import logging
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

logging.basicConfig(
    stream=sys.stdout,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
)

logging.info(f"[{os.path.basename(__file__)}]...")

# Parse command-line arguments
parser = argparse.ArgumentParser(description="PAGA Analysis Script")
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

file_path_meta_data_csv = os.path.join(input_dir, "OneK1K_meta_data.csv")
cell_metadata = (pd.read_csv(file_path_meta_data_csv, index_col=0)[['predicted.celltype.l2', 'cell_type']]
                 .rename(columns={'cell_type': 'celltype',
                                  'predicted.celltype.l2': 'celltype_pred'}))

file_path_feather = os.path.join(input_dir, "adjusted_expres_df.feather")
adjusted_expres_df = pd.read_feather(file_path_feather)
adjusted_expres_df.set_index('rownames', inplace=True)

file_path_gene_snp_df = os.path.join(input_dir, "gene_snp_df.csv")
gene_snp_df = pd.read_csv(file_path_gene_snp_df)

gene_metadata = (gene_snp_df[['GENE_ID', 'GENE']]
                 .drop_duplicates()
                 .rename(columns={'GENE': 'gene_short_name'}))
gene_metadata.set_index('GENE_ID', inplace=True)

file_path_cell_type_df = os.path.join(input_dir, "cell_type.csv")
cell_type_df = pd.read_csv(file_path_cell_type_df)

logging.info("Data loaded.")

# ---------------------------------------------------------

cell_metadata['original_index'] = cell_metadata.index

cell_metadata = pd.merge(
    cell_metadata,
    cell_type_df[['celltype_pred', 'category', 'category_abbr']],
    on='celltype_pred',
    how='left'
)

cell_metadata.set_index('original_index', inplace=True)
cell_metadata.index.name = None

cell_metadata = cell_metadata.loc[adjusted_expres_df.columns, :]

gene_metadata = gene_metadata.loc[adjusted_expres_df.index, :]
gene_metadata.index.name = "GENE_ID"

logging.info("Cell metadata and gene metadata reordered.")

# ---------------------------------------------------------

adata = ad.AnnData(X=adjusted_expres_df.T,
                   obs=cell_metadata,
                   var=gene_metadata)
logging.info(adata)
adata.write("onek1k_adata.h5ad")

logging.info("AnnData created.")
