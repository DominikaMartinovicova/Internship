import pandas as pd
import scanpy as sc
import scipy.sparse


annotations_file = "/trinity/home/dmartinovicova/all_data/GSE163120_RAW/GSM4972210_annot.Human.GBM.R1_2_3_4_4nc.csv.gz"

# Step 1: Load your CSV file
data = pd.read_csv("/trinity/home/dmartinovicova/all_data/GSE163120_RAW/GSM4972210_Human.GBM.R1_2_3_4_4nc.filtered.gene.bc.matrix.csv.gz", index_col=0, compression='gzip')

data = data.T

# Step 2: Extract the features (genes) and barcodes (cells)
# Rows (index) are the barcodes (cell names), columns are the features (gene names)
barcodes = data.index
genes = data.columns

print(barcodes)

# Step 3: Create a dense matrix from the data
matrix = data.values

# Step 4: Convert the matrix into a sparse format (CSR)
sparse_matrix = scipy.sparse.csr_matrix(matrix)

# Step 5: Create an AnnData object with the sparse matrix
adata = sc.AnnData(X=sparse_matrix)

# Step 6: Add the gene names and cell barcodes to the AnnData object
adata.var_names = genes  # Features/genes
adata.obs_names = barcodes  # Barcodes/cells

print(adata)

# Add the cell type column from annotation file
annotations_df = pd.read_csv(annotations_file, compression="gzip")
print(annotations_df.head())

annotations_df.set_index('cell', inplace=True)

# Merge annotations with adata
adata.obs = adata.obs.merge(
    annotations_df[['cluster']],  # Select only the 'cell_state' column
    left_index=True,  # Match this with 'cell_barcode' in adata.obs
    right_index=True,  # Use index from annotations_df for matching
    how='left'  # Keep all rows from adata.obs
)

adata.obs.rename(columns={'cluster': 'cell_type'}, inplace=True)
print(adata.obs)

print(adata.obs)
# Step 7: Save the AnnData object if desired
adata.write('/trinity/home/dmartinovicova/all_data/GSE163120_RAW/output/R.h5ad')





