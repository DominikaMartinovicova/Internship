import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

print('It works')
adata = sc.read_h5ad("/net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_complete/data/combined_all/adata_final_l012.h5ad")
print(adata)

# Exclude 'malignant' cell type
adata = adata[adata.obs['l2_celltype'] != 'Malignant']
print(adata)


# Create a contingency table of cell types versus diagnosis
contingency_table = pd.crosstab(adata.obs['l2_celltype'], adata.obs['Diagnosis_label'])
print(contingency_table)

# Calculate the percentage of each cell type in each disease
percentage_table = contingency_table.div(contingency_table.sum(axis=0), axis=1) * 100
print(percentage_table)

# Reset the index to use cell_label as a column and make the dataframe long
percentage_table_reset = percentage_table.reset_index().melt(id_vars='l2_celltype', var_name='Diagnosis_label', value_name='Percentage')
print(percentage_table_reset)

# Compute the standard deviation for each cell type and diagnosis
percentage_values = adata.obs.groupby(['l2_celltype', 'Diagnosis_label']).apply(lambda x: x['l2_celltype'].value_counts(normalize=True) * 100).reset_index['Percentage']
print('pct values worked')


std_table = percentage_values.groupby(['l2_celltype', 'Diagnosis_label'])['Percentage'].std()
print(std_table)

# Merge the standard deviation data with the percentage_table_reset
percentage_table_reset = percentage_table_reset.merge(std_table, on=['l2_celltype', 'Diagnosis_label'], how='left')

# Set the color palette to distinguish diseases
#palette = sns.color_palette("Set1", n_colors=len(percentage_table_reset['Diagnosis_label'].unique()))
palette=['darkorange', 'royalblue', 'red']

# Create the bar plot
plt.figure(figsize=(10, 6))
sns.barplot(x='l2_celltype', y='Percentage', hue='Diagnosis_label', data=percentage_table_reset, palette=palette, yerr=percentage_table_reset['std_percentage'])

# Customize plot
print('Plotting...')
plt.title('Cell Type Distribution by Disease')
plt.xlabel('Cell Type')
plt.ylabel('Percentage')
plt.xticks(rotation=70)
plt.legend(title='Diagnosis', bbox_to_anchor=(1.05, 1), loc='upper left')

# Display the plot
plt.tight_layout()

# Save the plot
plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_complete/composition_fig/cell_type_per_disease_w_errorb_wo_malignant.png')