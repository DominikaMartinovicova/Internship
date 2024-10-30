#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# QC_scRNA_before_Norm.py 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# QC for raw single cell data before normalization
#
# Author: Shiva Najjary (s.najjary@amsterdamumc.nl)
#
# TODO:
# 1) Fix reading and saving path for data
#
# History:
#  30-7-2024: File creation, write code
#-------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
# Read libraries
#------------------------------------------------------------------------------------------------------------

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad 

#------------------------------------------------------------------------------------------------------
#Reading adata file
#-------------------------------------------------------------------------------------------------------
# Path to the directory containing merged adata file
input_dir = "/trinity/home/dmartinovicova/all_data/OMIX920/output/GBM4_Astro3.h5ad"
adata = ad.read_h5ad(input_dir)

# Remove duplicate gene names
adata.var_names_make_unique()

print(adata) 

#---------------------------------------------------------------------------------------------------------
#1-Perform QC to identify the threshold for min and max genes to be filtered
#--------------------------------------------------------------------------------------------------------
# Calculate the number of genes and counts per cell
adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1  # Number of genes expressed per cell
adata.obs['n_counts'] = adata.X.sum(axis=1).A1  # Total counts per cell

genes_per_cell = adata.obs['n_genes']
counts_per_cell = adata.obs['n_counts']

#------------------------------------------------------------------------------------------------------------
# Save the QC metrics to a CSV file
#------------------------------------------------------------------------------------------------------------

qc_metrics_df = pd.DataFrame({
    'genes_per_cell': genes_per_cell,
    'counts_per_cell': counts_per_cell
})

qc_metrics_df.to_csv("/trinity/home/dmartinovicova/all_data/OMIX920/output/qc_metrics.csv", index=False)

#------------------------------------------------------------------------------------------------------------
# Calculate QC metrics
#------------------------------------------------------------------------------------------------------------

adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1  # Number of genes expressed per cell
adata.obs['n_counts'] = adata.X.sum(axis=1).A1  # Total counts per cell

# Save the QC metrics to a CSV file
qc_metrics_df = pd.DataFrame({
    'genes_per_cell': adata.obs['n_genes'],
    'counts_per_cell': adata.obs['n_counts'],
    'pct_mito': adata.obs['pct_counts_mt']
})
qc_metrics_df.to_csv('/trinity/home/dmartinovicova/all_data/OMIX920/output/qc_metrics_2.csv', index=False)

# Plot number of genes per cell
plt.figure(figsize=(8, 6))
plt.plot(range(len(adata.obs['n_genes'])), np.sort(adata.obs['n_genes']), 'o', alpha=0.5)
plt.xlabel('Cell')
plt.ylabel('Number of Genes')
plt.yscale('linear')  # Change to 'log' if needed
plt.title('Number of Genes per Cell (Ordered)')
plt.savefig('/trinity/home/dmartinovicova/all_data/OMIX920/output/genes_p_cell.png')



#------------------------------------------------------------------------------------------------------------
# Plot mitochondrial gene percentage
#------------------------------------------------------------------------------------------------------------

# Identify mitochondrial genes (if they have 'MT-' in their gene names)
adata.var['mito'] = adata.var_names.str.startswith('MT-')

# Calculate the percentage of mitochondrial genes
mito_gene_counts = np.sum(adata[:, adata.var['mito']].X, axis=1).A1
counts_per_cell = adata.X.sum(axis=1).A1
pct_mito = mito_gene_counts / counts_per_cell * 100

# Add pct_mito to adata observations
adata.obs['percent_mt'] = pct_mito

# Plotting the percentage of mitochondrial counts
plt.figure(figsize=(10, 6))
plt.hist(pct_mito, bins=50, edgecolor='black')
plt.xlabel('Percentage of Mitochondrial Counts')
plt.ylabel('Number of Cells')
plt.title('Distribution of Mitochondrial Gene Percentages')

# Save the plot
plt.savefig('/trinity/home/dmartinovicova/all_data/OMIX920/output/percent_mitochondrial_counts_distribution.png')
adata.obs[['percent_mt']].to_csv('/trinity/home/dmartinovicova/all_data/OMIX920/output/mito_percentage.csv')

# Cumulative plot
plt.figure(figsize=(10, 6))
plt.plot(np.sort(pct_mito), 'o', alpha=0.5)
plt.xlabel('Cells sorted by percentage mitochondrial counts')
plt.ylabel('Percentage mitochondrial counts')
plt.title('Percentage of Mitochondrial Counts per Cell')

# Define the maximum percentage of mitochondrial counts
MAX_PCT_MITO = 10

# Add a horizontal line at the maximum percentage
plt.axhline(y=MAX_PCT_MITO, color='red', linestyle='--')

# Save the plot
plt.savefig('/trinity/home/dmartinovicova/all_data/OMIX920/output/percent_mitochondrial_counts_cumulative.png')

# Determine the threshold
print(f"Suggested threshold: {MAX_PCT_MITO}%")

#------------------------------------------------------------------------------------------------------------
# To Find a Suitable Minimum Cell Threshold
#-------------------------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Calculate the number of cells in which each gene is expressed
gene_expression_counts = (adata.X > 0).sum(axis=0).A1  # Convert sparse matrix to array

# Define a range of thresholds to test (for example from 1 to 100 cells)
thresholds = np.arange(1, 201, 1)

# Calculate the number of genes expressed in at least each threshold
genes_expressed_at_least = [np.sum(gene_expression_counts >= t) for t in thresholds]

# Plot the number of genes expressed vs. the minimum number of cells
plt.figure(figsize=(30, 15))
plt.plot(thresholds, genes_expressed_at_least, marker='o')
plt.xlabel('Minimum number of cells')
plt.ylabel('Number of genes expressed')
plt.title('Number of genes expressed vs. Minimum number of cells')
plt.grid(True)

# Set x-axis limits and ticks
plt.xlim(0, 201)  # Set the range for x-axis
plt.xticks(np.arange(0, 201, 3))  # Set ticks every 10 units

#save plot
plt.savefig('/trinity/home/dmartinovicova/all_data/OMIX920/output/min_cells_threshold_plot.png')

# Save the results to a CSV file
threshold_df = pd.DataFrame({
    'min_cells_threshold': thresholds,
    'num_genes_expressed': genes_expressed_at_least
})

threshold_df.to_csv('/trinity/home/dmartinovicova/all_data/OMIX920/output/min_cells_threshold.csv', index=False)
























