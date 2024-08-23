import sys
import io
import os
import logging
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

working_dir = ""

input_dir = os.path.join(working_dir, "input")
output_dir = os.path.join(working_dir, "output")

# ---------------------------------------------------------

file_path_meta_data_csv = os.path.join(input_dir, "OneK1K_meta_data.csv")
cell_metadata = pd.read_csv(file_path_meta_data_csv, index_col=0)

file_path_feather = os.path.join(input_dir, "adjusted_expres_df.feather")
adjusted_expres_df = pd.read_feather(file_path_feather)
adjusted_expres_df.set_index('rownames', inplace=True)

file_path_gene_snp_df = os.path.join(input_dir, "gene_snp_df.csv")
gene_snp_df = pd.read_csv(file_path_gene_snp_df)

gene_metadata = (gene_snp_df[['GENE_ID', 'GENE']]
                 .drop_duplicates()
                 .rename(columns={'GENE': 'gene_short_name'}))
gene_metadata.set_index('GENE_ID', inplace=True)

logging.info("Data loaded.")

# ---------------------------------------------------------

cell_metadata = cell_metadata.loc[adjusted_expres_df.columns, :]

gene_metadata = gene_metadata.loc[adjusted_expres_df.index, :]
gene_metadata.index.name = "GENE_ID"

logging.info("Cell metadata and gene metadata reordered.")

# ---------------------------------------------------------

adata = ad.AnnData(X=adjusted_expres_df.T,
                   obs=cell_metadata,
                   var=gene_metadata)
logging.info(adata)
file_path_adata = os.path.join(input_dir, "onek1k_adata.h5ad")
adata.write(os.path.join(input_dir, "onek1k_adata.h5ad"))

logging.info("AnnData created.")