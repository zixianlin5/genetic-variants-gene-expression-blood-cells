import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

working_dir = ""

input_dir = os.path.join(working_dir, "input")
output_dir = os.path.join(working_dir, "output")
os.chdir(output_dir)

file_path_meta_data_csv = os.path.join(input_dir, "OneK1K_meta_data.csv")
cell_metadata = pd.read_csv(file_path_meta_data_csv, index_col=0)

# ---------------------------------------------------------

## Donor-related distributions

unique_donors = cell_metadata.drop_duplicates(subset='donor_id')
fig1, axes1 = plt.subplots(1, 3, figsize=(18, 6))

# Subplot 1: Donor Gender Distribution as a Pie Chart
sex_counts = unique_donors['sex'].value_counts()
# Adjusting the labels to be placed inside the pie chart
sex_labels = [f'{label}\n{count} ({count / len(unique_donors) * 100:.1f}%)'
              for label, count in zip(sex_counts.index, sex_counts.values)]

# Plotting the pie chart with labels inside
axes1[0].pie(sex_counts, labels=sex_labels, colors=['#7fbfbf', '#ffd17f'],
             startangle=90, textprops={'fontsize': 14}, labeldistance=0.25)
axes1[0].set_title('Gender Distribution of Donors', fontsize=16)

# Subplot 2: Donor Age Distribution
sns.histplot(unique_donors['age'], bins=20, kde=True, ax=axes1[1], color='teal')
axes1[1].set_title('Distribution of Donor Age', fontsize=16)
axes1[1].set_xlabel('Age', fontsize=14)
axes1[1].set_ylabel('Number of Donors', fontsize=14)

# Subplot 3: Distribution of Donor Cell Contributions
sns.histplot(cell_metadata['donor_id'].value_counts(), bins=20, kde=True, ax=axes1[2], color='orange')
axes1[2].set_title('Distribution of Cells Contributed by Donors', fontsize=16)
axes1[2].set_xlabel('Number of Cells', fontsize=14)
axes1[2].set_ylabel('Number of Donors', fontsize=14)

# Adjust tick labels font size
for ax in axes1:
    ax.tick_params(axis='both', which='major', labelsize=12)

fig1.tight_layout()
plt.show(block=True)

output_file1 = os.path.join(output_dir, "donor_related_distributions.png")
fig1.savefig(output_file1)

plt.close(fig1)

# ---------------------------------------------------------

## Cell quality metrics

fig2, axes2 = plt.subplots(1, 3, figsize=(18, 6))

# Subplot 1: nCount_RNA Distribution with adjusted annotations
sns.histplot(cell_metadata['nCount_RNA'], bins=20, kde=True, ax=axes2[0], color='teal')
axes2[0].set_title('Distribution of RNA Counts per Cell', fontsize=16)
axes2[0].set_xlabel('RNA Counts', fontsize=14)
axes2[0].set_ylabel('Number of Cells', fontsize=14)

min_nCount = cell_metadata['nCount_RNA'].min()
max_nCount = cell_metadata['nCount_RNA'].max()
axes2[0].annotate(f'Min: {min_nCount}', xy=(min_nCount, 0), xytext=(min_nCount + 15000, 100000),
                  arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.5), fontsize=12)
axes2[0].annotate(f'Max: {max_nCount}', xy=(max_nCount, 0), xytext=(max_nCount - 30000, 100000),
                  arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.5), fontsize=12)
# axes2[0].annotate(f'Min: {min_nCount}', xy=(min_nCount, 0), xytext=(min_nCount + 7000, 200000),
#                   arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.5), fontsize=12)
# axes2[0].annotate(f'Max: {max_nCount}', xy=(max_nCount, 0), xytext=(max_nCount - 15000, 120000),
#                   arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.5), fontsize=12)

# Subplot 2: nFeature_RNA Distribution with adjusted annotations
sns.histplot(cell_metadata['nFeature_RNA'], bins=20, kde=True, ax=axes2[1], color='purple')
axes2[1].set_title('Distribution of Features Detected per Cell', fontsize=16)
axes2[1].set_xlabel('Number of Features', fontsize=14)
axes2[1].set_ylabel('Number of Cells', fontsize=14)

min_nFeature = cell_metadata['nFeature_RNA'].min()
max_nFeature = cell_metadata['nFeature_RNA'].max()
axes2[1].annotate(f'Min: {min_nFeature}', xy=(min_nFeature, 0), xytext=(min_nFeature + 1500, 100000),
                  arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.5), fontsize=12)
axes2[1].annotate(f'Max: {max_nFeature}', xy=(max_nFeature, 0), xytext=(max_nFeature - 3000, 100000),
                  arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.5), fontsize=12)
# axes2[1].annotate(f'Min: {min_nFeature}', xy=(min_nFeature, 0), xytext=(min_nFeature + 1500, 150000),
#                   arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.5), fontsize=12)
# axes2[1].annotate(f'Max: {max_nFeature}', xy=(max_nFeature, 0), xytext=(max_nFeature - 1500, 120000),
#                   arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.5), fontsize=12)

# Subplot 3: percent.mt Distribution
sns.histplot(cell_metadata['percent.mt'], bins=20, kde=True, ax=axes2[2], color='orange')
axes2[2].set_title('Distribution of Mitochondrial Gene\nPercentage per Cell', fontsize=16)
axes2[2].set_xlabel('Percent Mitochondrial Genes', fontsize=14)
axes2[2].set_ylabel('Number of Cells', fontsize=14)

for ax in axes2:
    ax.tick_params(axis='both', which='major', labelsize=12)

fig2.tight_layout()
plt.show(block=True)

output_file2 = os.path.join(output_dir, "cell_quality_metrics_adjusted.png")
fig2.savefig(output_file2)

plt.close(fig2)

# ---------------------------------------------------------

## Distribution of Cell Type

celltype_count = (cell_metadata[["cell_type", "orig.ident"]]
                  .groupby("cell_type")
                  .count()
                  .sort_values(by="orig.ident"))

output_file = os.path.join(output_dir, "celltype_count.xlsx")
celltype_count.to_excel(output_file)

