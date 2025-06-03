import scanpy as sc

study_id = "_OMIX920"

# Load the two h5ad files
adata_cancer1 = sc.read_h5ad("/trinity/home/dmartinovicova/all_data/OMIX920/output/Astro3.h5ad")
adata_cancer2 = sc.read_h5ad("/trinity/home/dmartinovicova/all_data/OMIX920/output/GBM4.h5ad")

# Optionally: Print basic information about the datasets
print(adata_cancer1)
print(adata_cancer2)

# Concatenate the AnnData objects based on common genes
common_genes = adata_cancer1.var_names.intersection(adata_cancer2.var_names)

# Subset the AnnData objects to keep only common genes
adata_cancer1_common = adata_cancer1[:, common_genes].copy()
adata_cancer2_common = adata_cancer2[:, common_genes].copy()

# Concatenate the AnnData objects
combined_adata = adata_cancer1_common.concatenate(adata_cancer2_common)

# Add the study identifier
combined_adata.obs_names = combined_adata.obs_names + study_id

# Create gene_id column and move it to position [1]
combined_adata.var['gene_id'] = combined_adata.var.index
cols= list(combined_adata.var.columns)
cols.insert(1, cols.pop(cols.index('gene_id')))
combined_adata.var = combined_adata.var[cols]  

combined_adata.var = combined_adata.var.sort_values(by='gene_id')

# Add mitochondrial
combined_adata.var['mt'] = combined_adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(combined_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Remove batch column
if 'batch' in combined_adata.obs.columns:
  combined_adata.obs = combined_adata.obs.drop(columns = ['batch'])


# Save the combined AnnData object
combined_adata.write("/trinity/home/dmartinovicova/all_data/OMIX920/output/GBM4_Astro3.h5ad")

print("Combined AnnData file saved")

