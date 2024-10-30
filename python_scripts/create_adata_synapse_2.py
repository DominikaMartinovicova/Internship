import anndata as ad
import scanpy as sc
import numpy as np

# Load the .h5ad file
adata = ad.read_h5ad("/trinity/home/dmartinovicova/all_data/Verhaak_synapse/analysis_scRNAseq_tumor_counts.h5ad")
print(adata.obs)
# Batch information is stored in adata.obs['batch']
batch_column = 'batch'

# Split the data by batch
batches = adata.obs[batch_column].unique()  # Unique batch identifiers
batch_data_list = []

for batch in batches:
    # Subset data for each batch
    batch_data = adata[adata.obs[batch_column] == batch].copy()
    
    # Remove unnecessary information, keeping only matrix, var (genes), and obs (cells)
    batch_data.obsm = None
    batch_data.varm = None
    batch_data.uns = {}
    batch_data.obsp = None
    batch_data.varp = None
    batch_data.raw = None
    
    # Add batch dataset to list
    batch_data_list.append(batch_data)

# Find common genes across all batches
# Get list of genes for each batch
gene_sets = [set(batch_data.var_names) for batch_data in batch_data_list]

# Find intersection of all gene sets
common_genes = set.intersection(*gene_sets)

# Convert to a sorted list for consistency
common_genes = sorted(list(common_genes))

#Filter each batch dataset to keep only common genes
filtered_batch_data_list = []
for batch_data in batch_data_list:
    # Filter the batch dataset to keep only the common genes
    batch_data = batch_data[:, common_genes].copy()
    filtered_batch_data_list.append(batch_data)


#Concatenate all batch datasets into one final AnnData object
final_adata = ad.concat(filtered_batch_data_list, join='inner', label=batch_column, keys=batches)

# Add the study identifier
final_adata.obs_names = final_adata.obs_names + '_synapse'

# Create gene_id column and move it to position [1]
final_adata.var['gene_id'] = final_adata.var.index
cols= list(final_adata.var.columns)
cols.insert(1, cols.pop(cols.index('gene_id')))
final_adata.var = final_adata.var[cols]  

final_adata.var = final_adata.var.sort_values(by='gene_id')

# Add mitochondrial
final_adata.var['mt'] = final_adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(final_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# Save the final combined dataset
final_adata.write("/trinity/home/dmartinovicova/all_data/Verhaak_synapse/output/Verhaak_synapse_2.h5ad")

# Optional: Check the shape of the final dataset (should be number of cells x number of common genes)
print(final_adata.shape)
