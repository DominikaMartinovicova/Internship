import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the AnnData object (replace with your actual file path)
print('Reading adata...')
adata = sc.read("/net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_complete/data/combined_all/adata_raw_w_labels_l012.h5ad")

# Check if the data is sparse (it likely is)
# If it is sparse, we can directly sum along axis 1 (cells)
if hasattr(adata.X, 'sum'):
    # Sum the expression across genes (columns) for each cell (row)
    total_counts_per_cell = adata.X.sum(axis=1).A1  # .A1 converts it to a 1D array

# Create a DataFrame with the total read counts and metadata
print('Creating a DataFrame...')
read_counts = pd.DataFrame({
    'Total_reads': total_counts_per_cell,
    'Cell_type': adata.obs['l2_celltype'],
    'Cohort': adata.obs['Cohort']
})
print(read_counts)
# Group by 'cell_type' and 'cohort' and calculate min, max, and range of read counts
read_counts_grouped = read_counts.groupby(['Cell_type', 'Cohort'])['Total_reads'].agg(['min', 'max', 'mean', 'std'])

# Display the results
print(read_counts_grouped)

# Plot a boxplot for read count distribution by cell type and cohort
plt.figure(figsize=(12, 8))
sns.boxplot(x='Cell_type', y='Total_reads', hue='Cohort', data=read_counts, palette='Paired')
plt.xticks(rotation=45, ha='right')
plt.title("Read count distribution per cell type in each cohort")
plt.show()
plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova/figures/figures/readcount_range_percelltype_percohort.png')