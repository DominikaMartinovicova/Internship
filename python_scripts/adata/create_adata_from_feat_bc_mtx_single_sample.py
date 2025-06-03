import os
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from anndata import AnnData


#-----------------------------------------------------------------------------
# 0. Prepare variables
#-----------------------------------------------------------------------------

study_id = "_GSE182109_broadSC"
data_dir = "/trinity/home/dmartinovicova/all_data/GSE182109_bSC/"
output_file = "/trinity/home/dmartinovicova/all_data/GSE182109_bSC/GSE182109_bSC.h5ad"

matrix_path = os.path.join(data_dir, 'Raw_matrix.mtx.gz') 
barcodes_path = os.path.join(data_dir, 'Raw_barcodes.tsv.gz') 
features_path = os.path.join(data_dir, 'Raw_genes.tsv.gz') 

#-----------------------------------------------------------------------------
# 1. Create adata object
#-----------------------------------------------------------------------------

# Load the sparse matrix from the .mtx.gz file
print("Loading matrix...")
matrix = mmread(matrix_path).tocsc().T

# Load the metadata (genes and barcodes) from .tsv files
print("Loading genes and barcodes...")
genes = pd.read_csv(features_path, sep="\t", header=None)[0].values
barcodes = pd.read_csv(barcodes_path, sep='\t', header=None)[0].values

# Create the AnnData object
adata = sc.AnnData(X=matrix, var=pd.DataFrame(index=genes), obs=pd.DataFrame(index=barcodes))
print("Adata from matrix created")

#-----------------------------------------------------------------------------
# 2. Add information
#-----------------------------------------------------------------------------

# Identify mit genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Add study identifier
adata.obs_names = adata.obs_names + study_id

print(adata.var)
print(adata.obs)

# Save the AnnData object as an .h5ad file
adata.write(output_file)

