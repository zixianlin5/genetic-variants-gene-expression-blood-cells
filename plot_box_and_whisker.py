import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

working_dir = ""

input_dir = os.path.join(working_dir, "input")
output_dir = os.path.join(working_dir, "output")
os.chdir(output_dir)

file_path_pca_data = os.path.join(input_dir, 'visual_data.csv')
visual_data = pd.read_csv(file_path_pca_data)

print("Data loaded.")

# ---------------------------------------------------------

cell_type_abbreviations = {
    "central memory CD4-positive, alpha-beta T cell": "CD4+ TCM",
    "natural killer cell": "NK",
    "naive thymus-derived CD4-positive, alpha-beta T cell": "CD4+ TN",
    "effector memory CD8-positive, alpha-beta T cell": "CD8+ TEM",
    "naive B cell": "B naive",
    "central memory CD8-positive, alpha-beta T cell": "CD8+ TCM",
    "effector memory CD4-positive, alpha-beta T cell": "CD4+ TEM",
    "memory B cell": "B mem",
    "erythrocyte": "Ery",
    "naive thymus-derived CD8-positive, alpha-beta T cell": "CD8+ TN",
    "CD16-negative, CD56-bright natural killer cell, human": "CD56bright NK",
    "transitional stage B cell": "B trans",
    "CD14-low, CD16-positive monocyte": "CD14low CD16+ Mono",
    "mucosal invariant T cell": "MAIT",
    "regulatory T cell": "Treg",
    "CD14-positive monocyte": "CD14+ Mono",
    "gamma-delta T cell": "Tγδ",
    "CD4-positive, alpha-beta cytotoxic T cell": "CD4+ CTL",
    "double negative thymocyte": "DN Thymo",
    "hematopoietic precursor cell": "HPC",
    "platelet": "Plt",
    "plasmablast": "PB",
    "innate lymphoid cell": "ILC",
    "conventional dendritic cell": "cDC",
    "plasmacytoid dendritic cell": "pDC",
    "CD8-positive, alpha-beta T cell": "CD8+ T",
    "CD4-positive, alpha-beta T cell": "CD4+ T",
    "dendritic cell": "DC"
}

visual_data['cell_type'] = visual_data['cell_type'].map(cell_type_abbreviations)

# ---------------------------------------------------------

for pc in [f'PC_{i}' for i in range(1, 15)]:
    plt.figure(figsize=(10, 6))

    order = visual_data.groupby('cell_type')[pc].mean().sort_values(ascending=False).index

    sns.boxplot(x='cell_type', y=pc, data=visual_data, order=order)

    plt.title(f'Distribution of {pc.replace("_", "")} across cell types', fontsize=16)
    plt.xticks(rotation=45, ha='right', fontsize=10)
    plt.xlabel('Cell Type', fontsize=14)
    plt.ylabel(pc.replace("_", ""), fontsize=14)
    plt.tight_layout()

    output_file = f'distribution_of_{pc}_across_cell_types.png'
    plt.savefig(output_file, dpi=300)

    plt.close()
    print(f"Figure saved: {output_file}")
    
