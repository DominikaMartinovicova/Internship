import scanpy as sc

# Load your AnnData object (adata)
adata = sc.read_h5ad("/trinity/home/dmartinovicova/snakemake/data/QC_results/all_normalized.h5ad")

# Run PCA before UMAP
print('Running PCA...')
sc.tl.pca(adata, svd_solver='arpack')

# Compute neighbors (this is needed before running UMAP)
print('Computing neighbors...')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=25)

# Run UMAP
print('Running UMAP...')
sc.tl.umap(adata)

# Plot the UMAP
print('Plotting UMAP...')
sc.pl.umap(adata, color=['batch'], save='UMAP.png')  
