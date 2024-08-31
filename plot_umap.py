import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
import numpy as np

working_dir = ""

input_dir = os.path.join(working_dir, "input")
output_dir = os.path.join(working_dir, "output")
os.chdir(output_dir)

file_path_pca_data = os.path.join(input_dir, 'visual_data.csv')
visual_data = pd.read_csv(file_path_pca_data)

print("Data loaded.")

# ---------------------------------------------------------

custom_palette = [
    "#ff9999",  # CD4+ TCM
    "#66b3ff",  # NK
    "#99ff99",  # CD4+ TN
    "#ffcc99",  # CD8+ TEM
    "#c2c2f0",  # B naive
    "#ffb3e6",  # CD8+ TCM
    "#c4e17f",  # CD4+ TEM
    "#ff6666",  # B mem
    "#ffcc00",  # Ery
    "#c6a300",  # CD8+ TN
    "#00cc99",  # CD56bright NK
    "#009999",  # B trans
    "#66ff66",  # CD14low CD16+ Mono
    "#ff9966",  # MAIT
    "#ff99cc",  # Treg
    "#6699ff",  # CD14+ Mono
    "#cc99ff",  # Tγδ
    "#66ccff",  # CD4+ CTL
    "#ff6699",  # DN Thymo
    "#ff66cc",  # HPC
    "#99ffcc",  # Plt
    "#ccff66",  # PB
    "#ff6600",  # ILC
    "#cc99cc",  # cDC
    "#ccff99",  # pDC
    "#669999",  # CD8+ T
    "#6666ff",  # CD4+ T
    "#66ffff"  # DC
]

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

plt.figure(figsize=(12, 10))
palette = sns.husl_palette(28)

sns.scatterplot(
    x='UMAP_1', y='UMAP_2',
    data=visual_data,
    hue='cell_type',
    palette=custom_palette,
    s=2, alpha=0.7,
    marker=".",
    edgecolor=None,
    legend=False
)

plt.xlabel('UMAP 1', fontsize=16)
plt.ylabel('UMAP 2', fontsize=16)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

offsets = {
    'ILC': (0, 0.4),
    'CD56bright NK': (-0.4, -0.2),
    'CD4+ T': (0, 0.1)
}

for cell_type, group in visual_data.groupby('cell_type'):
    counts, xedges, yedges = np.histogram2d(group['UMAP_1'], group['UMAP_2'], bins=50)

    max_bin_idx = np.unravel_index(np.argmax(counts), counts.shape)

    max_density_x = (xedges[max_bin_idx[0]] + xedges[max_bin_idx[0] + 1]) / 2
    max_density_y = (yedges[max_bin_idx[1]] + yedges[max_bin_idx[1] + 1]) / 2

    offset_x, offset_y = offsets.get(cell_type, (0, 0))

    plt.text(
        max_density_x + offset_x,
        max_density_y + offset_y,
        cell_type,
        fontsize=12,
        ha='center',
        va='center'
    )

plt.tight_layout()

output_file_path = os.path.join(output_dir, "UMAP_plot.png")
plt.savefig(output_file_path, dpi=300, bbox_inches='tight')
plt.close()

print("Done!")

