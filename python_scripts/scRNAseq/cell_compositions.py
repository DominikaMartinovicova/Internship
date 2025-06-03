import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Load your AnnData object (assuming 'adata' is already loaded)
print('Loading adata...')
adata = sc.read_h5ad("/trinity/home/dmartinovicova/snakemake/data/combined_all/adata_final.h5ad")

# Check the first few rows of the metadata to understand its structure
print(adata.obs.head())

# Group by Diagnosis_label and cell type and count occurrences
composition = adata.obs.groupby(['Diagnosis_label', 'l1_celltype']).size().unstack(fill_value=0)

# Plot the bar graph
print('Plotting...')

# Normalize by row (i.e., by disease), so that each row adds up to 100%
composition_percent = composition.div(composition.sum(axis=1), axis=0) * 100

# Plot the bar graph
custom_colors = ['royalblue', 'darkorange', 'mediumseagreen', 'crimson', 'darkviolet', 'steelblue', 'sienna', 'lightcoral', 'darkkhaki','palevioletred', 'mediumturquoise', 'lightblue', 'dimgrey', 'mediumpurple', 'chocolate', 'firebrick' ]
ax = composition_percent.plot(kind='bar', stacked=True, figsize=(8, 6), color=custom_colors)

# Add labels and title
ax.set_ylabel('Percentage of Cells')
ax.set_title('Cell Type Composition as Percentage in Each Disease')

# Rotate x-axis labels for better readability
plt.xticks(rotation=0)
plt.tight_layout()

# Adjust the layout to add space for the legend
plt.subplots_adjust(right=0.7)  # Increase right margin to make space for the legend

# Move the legend outside the plot to the right
ax.legend(title="Cell Types", bbox_to_anchor=(1.4, 0.5), loc='center right')

# Save the figure
plt.savefig('/trinity/home/dmartinovicova/snakemake/composition_fig/cell_composition_percentage.png')

# Optionally, show the plot
plt.show()






# Filter out malignant cells (assuming 'cell_type' column has 'malignant' label)
# Adjust this condition based on your specific labeling for malignant cells
non_malignant_adata = adata[~adata.obs['l1_celltype'].isin(['Malignant']), :]

# Group by diagnosislabel and cell type, then count occurrences
composition = non_malignant_adata.obs.groupby(['Diagnosis_label', 'l1_celltype']).size().unstack(fill_value=0)

# Normalize by row (i.e., by disease), so that each row adds up to 100%
composition_percent = composition.div(composition.sum(axis=1), axis=0) * 100

# Plot the bar graph
custom_colors = ['royalblue', 'darkorange', 'mediumseagreen', 'crimson', 'darkviolet', 'steelblue', 'lightcoral', 'darkkhaki','palevioletred', 'mediumturquoise', 'lightblue', 'dimgrey', 'mediumpurple', 'chocolate', 'firebrick' ]
ax = composition_percent.plot(kind='bar', stacked=True, figsize=(8, 6), color=custom_colors)

# Add labels and title
ax.set_ylabel('Percentage of Cells')
ax.set_title('Cell Type Composition (Excluding Malignant Cells) in Each Disease')

plt.xticks(rotation=0)
plt.tight_layout()

# Adjust the layout to add space for the legend
plt.subplots_adjust(right=0.7)  # Increase right margin to make space for the legend

# Move the legend outside the plot to the right
ax.legend(title="Cell Types", bbox_to_anchor=(1.4, 0.5), loc='center right')

# Save the figure
plt.savefig('/trinity/home/dmartinovicova/snakemake/composition_fig/cell_composition_non_malignant_percentage.png')

# Optionally, show the plot
plt.show()