import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

np.float = np.float64

working_dir = ""

input_dir = os.path.join(working_dir, "input")
output_dir = os.path.join(working_dir, "output")
os.chdir(output_dir)

file_path_loadings = os.path.join(output_dir, "Top_Loadings.xlsx")
df = pd.read_excel(file_path_loadings, header=0, index_col=None)


# ---------------------------------------------------------


def plot_loadings(df, pcs, output_file, ncol=2, figsize=(15, 10)):
    num_pcs = len(pcs)

    fig, axes = plt.subplots(nrows=(num_pcs + ncol - 1) // ncol, ncols=ncol, figsize=figsize)
    axes = axes.flatten()

    for i, pc in enumerate(pcs):
        ax = axes[i]
        data = df[df['PC'] == pc]
        sns.barplot(x='Loading', y='Gene', data=data, ax=ax, palette="GnBu_r")
        ax.set_title(f"Top Genes for {pc}", fontsize=16)
        ax.set_xlabel("Loading", fontsize=16)
        ax.set_ylabel("Gene", fontsize=16)
        ax.tick_params(axis='both', which='major', labelsize=10)

    # Remove any unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


pcs = df['PC'].unique()

pcs_part1 = pcs[:8]
pcs_part2 = pcs[8:14]

plot_loadings(df, pcs_part1, os.path.join(output_dir, "Top_genes_PC_part1.png"), ncol=2, figsize=(13, 20))
plot_loadings(df, pcs_part2, os.path.join(output_dir, "Top_genes_PC_part2.png"), ncol=2, figsize=(13, 15))

print("Done!")

