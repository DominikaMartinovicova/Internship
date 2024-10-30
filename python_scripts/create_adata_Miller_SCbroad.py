import pandas as pd
import scanpy as sc
from scipy.io import mmread
from anndata import AnnData

study_id = "_Miller_broadSC"

# Load the sparse matrix from the .mtx.gz file
print("Loading matrix...")
matrix = mmread("/trinity/home/dmartinovicova/all_data/Miller_broadSC/matrix.mtx.gz").tocsc().transpose()

# Load the metadata (genes and barcodes) from .tsv files
print("Loading genes and barcodes...")
genes = pd.read_csv("/trinity/home/dmartinovicova/all_data/Miller_broadSC/genes.tsv", header=None, sep="\t")[0].values
barcodes = pd.read_csv("/trinity/home/dmartinovicova/all_data/Miller_broadSC/barcodes.tsv", sep='\t', header=None)[0].values
print(len(genes))
print(len(barcodes))

# Create the AnnData object
adata = sc.AnnData(X=matrix, var=pd.DataFrame(index=genes), obs=pd.DataFrame(index=barcodes))
print("adata from matrix created")

# Create gene_id column and move it to position [1]
adata.var['gene_id'] = adata.var.index
cols= list(adata.var.columns)
cols.insert(1, cols.pop(cols.index('gene_id')))
adata.var = adata.var[cols]  

adata.var = adata.var.sort_values(by='gene_id')

# Add mitochondrial
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata.obs_names = adata.obs_names + study_id

print(adata.var)
print(adata.obs)

# Save the AnnData object as an .h5ad file
adata.write("/trinity/home/dmartinovicova/all_data/Miller_broadSC/adata_Miller_broadSC.h5ad")


