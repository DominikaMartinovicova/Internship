import os
import scanpy as sc
import numpy as np
from scipy.sparse import csr_matrix

# Specify the folder containing your .h5ad files
folder_path = '/net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_complete/data/processed'

# Loop through each file in the folder
for filename in os.listdir(folder_path):
    if filename.endswith('.h5ad'):
        # Load the AnnData object from the file
        adata = sc.read(os.path.join(folder_path, filename))
        print(adata)
        # Ensure the data is sparse (AnnData.X can be a sparse matrix)
        if isinstance(adata.X, csr_matrix):  # Check if it's sparse
            # Get the minimum and maximum values directly from the sparse matrix
            min_val = adata.X.min()
            max_val = adata.X.max()
        else:
            # If it's dense, you can still get the min/max
            min_val = np.min(adata.X)
            max_val = np.max(adata.X)

        # Optionally print the result for each dataset
        print(f"{filename}: Min = {min_val}, Max = {max_val}")
