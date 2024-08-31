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

sns.set(style="white")
sns.set_context("notebook")

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 14
plt.rcParams['font.weight'] = 'normal'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'normal'

palette = sns.husl_palette(28)

for j in range(2):
    fig, axes = plt.subplots(3, 3, figsize=(18, 25))
    axes = axes.flatten()
    scatter = None
    for i in range(9):
        idx = i + j * 9
        ax = axes[i]

        scatter = sns.scatterplot(
            x=visual_data[f'PC_{idx + 1}'],
            y=visual_data[f'PC_{(idx + 2) if (idx + 2) <= 19 else 1}'],
            hue=visual_data['cell_type'],
            ax=ax, palette=palette,
            s=5, marker=".",
            edgecolor=None
        )

        ax.set_xlabel(f'PC{idx + 1}', fontsize=14)
        ax.set_ylabel(f'PC{(idx + 2) if (idx + 2) <= 19 else 1}', fontsize=14)
        ax.set_title(
            f'PCA Plot of PC{idx + 1} and PC{(idx + 2) if (idx + 2) <= 19 else 1}',
            weight='bold', fontsize=16
        )

        ax.spines['bottom'].set_color('black')
        ax.spines['left'].set_color('black')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='x')
        ax.tick_params(axis='y')
        ax.legend_.remove()

    handles, labels = scatter.get_legend_handles_labels()

    new_labels = [cell_type_abbreviations.get(label, label) for label in labels]

    fig.legend(
        handles, new_labels,
        loc='lower center',
        bbox_to_anchor=(0.5, -0.09),
        ncol=6,
        frameon=False,
        handletextpad=0.1,
        fontsize=18,
        handlelength=2,
        markerscale=1.5,
        columnspacing=4
    )

    plt.tight_layout(rect=[0, 0, 1, 0.73])

    output_file = f"pca_plots_part_{j + 1}.png"
    fig.savefig(output_file, dpi=300, bbox_inches='tight')

    plt.close(fig)
    print(f"Figure saved: {output_file}")


