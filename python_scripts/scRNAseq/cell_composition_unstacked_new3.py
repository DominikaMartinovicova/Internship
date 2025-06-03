import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Assuming your AnnData object is named `adata` and has columns 'celltype', 'disease', and 'patient' in `adata.obs`
adata = sc.read("/net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_complete/data/combined_all/adata_final_l012.h5ad")
print(adata)
adata = adata[adata.obs['l2_celltype'] != 'Malignant']

# 1. Calculate the proportions of each cell type per disease-patient group
# First, create a contingency table for disease, patient, and celltype
contingency = pd.crosstab([adata.obs['Diagnosis_label'], adata.obs['case_id']], adata.obs['l2_celltype'])
print(contingency)

# Calculate the proportions for each disease-patient group (relative counts)
proportions = contingency.div(contingency.sum(axis=1), axis=0)
print(proportions)

# 2. Calculate the mean proportion per disease for each cell type
# Group by disease and calculate the mean across patients
mean_proportions_per_disease = proportions.groupby(level=0).mean()  # Mean across patients for each disease
std_proportions_per_disease = proportions.groupby(level=0).std()
print(mean_proportions_per_disease)

# Reset index to facilitate plotting
mean_proportions_per_disease = mean_proportions_per_disease.reset_index()
std_proportions_per_disease = std_proportions_per_disease.reset_index()
print(mean_proportions_per_disease)

# 3. Prepare data for plotting
# Melt the data to long format for easier plotting with seaborn
plot_data = pd.melt(mean_proportions_per_disease, id_vars=['Diagnosis_label'], var_name='l2_celltype', value_name='mean_proportion')
std_data = pd.melt(std_proportions_per_disease, id_vars=['Diagnosis_label'], var_name='l2_celltype', value_name='std_proportion')

# reshape the dataframe into a wide format for Values
vals = plot_data.pivot(index='l2_celltype', columns='Diagnosis_label', values='mean_proportion')
print(vals)

# reshape the dataframe into a wide format for Errors
yerr = std_data.pivot(index='l2_celltype', columns='Diagnosis_label', values='std_proportion')
print(yerr)

# plot vals with yerr
ax = vals.plot(kind='bar', yerr=yerr, figsize=(10, 6), 
     error_kw={'elinewidth': 1},
     color=['orange', 'royalblue', 'red'],  # Specify bar colors
     width=0.9)

# Adding title and labels
plt.title("Average Proportion of Cell Types per Disease with Error Bars (Standard Deviation)")
plt.ylabel("Mean Proportion")
plt.xlabel("Cell Type")
plt.xticks(rotation=70)
plt.tight_layout()
plt.legend(title='Disease', loc='upper right')
plt.show()

plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_complete/composition_fig/cell_type_per_disease_w_errorb.png')



