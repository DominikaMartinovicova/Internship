import scanpy as sc

adata = sc.read_h5ad("/trinity/home/dmartinovicova/snakemake/data/combined_all/adata_clustered_bc_ing_denovo_complete2.h5ad")

# List the columns you want to remove from adata.obs 
columns_to_remove = ['kept_labels', 
                      'kept_denovo_myeloid',
                      'de_novo',
                      'denovo_myeloid',
                      'denovo_immune',
                      'complete_denovo_myeloid',
                      'complete_denovo'] 

# Remove the specified columns from adata_copy.obs
adata.obs = adata.obs.drop(columns=columns_to_remove)


# Rename the columns
adata.obs = adata.obs.rename(columns={'kept_denovo': 'l1_celltype',
                                      'kept_denovo_TAM': 'l0_celltype'})

# Display the first few rows of the updated 'obs' to confirm the changes
print(adata.obs)
print(adata)

adata.write("/trinity/home/dmartinovicova/snakemake/data/combined_all/adata_final.h5ad")