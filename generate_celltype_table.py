import os
import pandas as pd

working_dir = ""

input_dir = os.path.join(working_dir, "input")
output_dir = os.path.join(working_dir, "output")
os.chdir(output_dir)

file_path_meta_data_csv = os.path.join(input_dir, "OneK1K_meta_data.csv")
cell_metadata = pd.read_csv(file_path_meta_data_csv, index_col=0)

# ---------------------------------------------------------

cell_type_df = (cell_metadata[["cell_type", "predicted.celltype.l2"]]
                .drop_duplicates()
                .sort_values("cell_type")
                .reset_index(drop=True)
                .rename(columns={'cell_type': 'celltype',
                                 'predicted.celltype.l2': 'celltype_pred'}))

# ---------------------------------------------------------

category_mapping = {
    'CD16 Mono': 'Monocyte',
    'CD14 Mono': 'Monocyte',
    'NK_CD56bright': 'Natural Killer Cell',
    'NK': 'Natural Killer Cell',
    'NK Proliferating': 'Natural Killer Cell',
    'CD4 Proliferating': 'T Cell',
    'CD4 CTL': 'T Cell',
    'CD8 Proliferating': 'T Cell',
    'CD4 TCM': 'T Cell',
    'CD8 TCM': 'T Cell',
    'dnT': 'T Cell',
    'CD4 TEM': 'T Cell',
    'CD8 TEM': 'T Cell',
    'gdT': 'T Cell',
    'MAIT': 'T Cell',
    'CD4 Naive': 'T Cell',
    'CD8 Naive': 'T Cell',
    'Treg': 'T Cell',
    'cDC1': 'Dendritic Cell',
    'cDC2': 'Dendritic Cell',
    'ASDC': 'Dendritic Cell',
    'pDC': 'Dendritic Cell',
    'B memory': 'B Cell',
    'B naive': 'B Cell',
    'Plasmablast': 'B Cell',
    'B intermediate': 'B Cell',
    'HSPC': 'Hematopoietic Stem/Progenitor Cell',
    'Eryth': 'Erythrocyte',
    'Platelet': 'Platelet',
    'ILC': 'Innate Lymphoid Cell',
    'Doublet': 'Doublet'
}

abbreviation_mapping = {
    'Monocyte': 'Mono',
    'Natural Killer Cell': 'NK',
    'T Cell': 'T Cell',
    'Dendritic Cell': 'DC',
    'B Cell': 'B Cell',
    'Hematopoietic Stem/Progenitor Cell': 'HSPC',
    'Erythrocyte': 'Eryth',
    'Platelet': 'Platelet',
    'Innate Lymphoid Cell': 'ILC',
    'Doublet': 'Doublet'
}

cell_type_df['category'] = cell_type_df['celltype_pred'].map(category_mapping)
cell_type_df['category_abbr'] = cell_type_df['category'].map(abbreviation_mapping)
cell_type_df.to_csv("cell_type.csv", index=False)

print('Done!')

