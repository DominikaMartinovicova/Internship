import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("/trinity/home/dmartinovicova/snakemake/data/combined_all/adata_final.h5ad")

# Extract the cell identifier and the cell type columns
# Assuming 'cell_id' is the index of the AnnData object, and 'l1_celltype' and 'l0_celltype' are columns in adata.obs
cell_data = adata.obs[['l1_celltype', 'l0_celltype']].copy()

# Optionally, reset the index to add 'cell_id' as a column (if needed)
cell_data['cells'] = adata.obs.index.str.replace(r'-(ref|target)$', '', regex=True)

# Reorder columns if necessary (to have cell_id first)
cell_data = cell_data[['cells', 'l1_celltype', 'l0_celltype']]

print(cell_data)

# Save the DataFrame to a CSV file
cell_data.to_csv('/trinity/home/dmartinovicova/snakemake/data/annotations_markers/l1l0_celltypes.csv', index=False)

# Print a message to confirm
print('Cell type data saved to "l1l0_celltypes.csv"')