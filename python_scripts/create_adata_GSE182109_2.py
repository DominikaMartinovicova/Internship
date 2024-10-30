import os  
import gzip
import scanpy as sc 
import scipy.io 
import pandas as pd 
import numpy as np 
from scipy.sparse import csr_matrix 

# Directory containing your 44 samples 
data_dir = '/trinity/home/dmartinovicova/all_data/GSE182109_RAW/' 
output_file = '/trinity/home/dmartinovicova/all_data/GSE182109_RAW/output/combined_adata_GSE182109.h5ad'  
study_id = '_GSE182109'

# Get a list of all files in the directory 
files = os.listdir(data_dir) 

# Initialize a set to hold the unique sample IDs 
sample_ids = set() 

for file in files: 
    if file.endswith('matrix.mtx.gz'): 
        # Extract the part of the filename before 'matrix.mtx' to get the sample ID 
        sample_id = file.replace('matrix.mtx.gz', '') 
        sample_ids.add(sample_id) 

# Convert the set to a sorted list for consistency 
sample_ids = sorted(list(sample_ids)) 

# Print the sample IDs or use them in your code 
print(sample_ids) 

# To store the overlapping genes across all samples 
overlapping_genes = None 

# To store the data for each sample 
adata_list = [] 

# Iterate over each sample 
for sample_id in sample_ids: 
    # Define paths for the matrix, barcodes, and features files 
    matrix_path = os.path.join(data_dir, f'{sample_id}matrix.mtx.gz') 
    barcodes_path = os.path.join(data_dir, f'{sample_id}barcodes.tsv.gz') 
    features_path = os.path.join(data_dir, f'{sample_id}features.tsv.gz') 

    # Load the matrix (gene expression data) 
    with gzip.open(matrix_path, 'rt') as f:
      matrix = scipy.io.mmread(f).tocsr() 

    # Load the barcodes (cell names) 
    with gzip.open(barcodes_path, 'rt') as f:
      barcodes = pd.read_csv(f, sep='\t', header=None)[0].values 

    # Load the features (gene names) 
    with gzip.open(features_path, 'rt') as f:
      features = pd.read_csv(f, sep='\t', header=None)[1].values 

    # Find overlapping genes 
    if overlapping_genes is None: 
        overlapping_genes = set(features) 
    else: 
        overlapping_genes = overlapping_genes.intersection(set(features))

# Convert overlapping_genes to list for consistent ordering 
overlapping_genes = sorted(list(overlapping_genes)) 

a=0
# Iterate again to extract the data only for overlapping genes and merge it 
for sample_id in sample_ids: 
    a+=1
    # Define paths for the matrix, barcodes, and features files 
    matrix_path = os.path.join(data_dir, f'{sample_id}matrix.mtx.gz') 
    barcodes_path = os.path.join(data_dir, f'{sample_id}barcodes.tsv.gz') 
    features_path = os.path.join(data_dir, f'{sample_id}features.tsv.gz')  

    # Load the matrix (gene expression data) 
    with gzip.open(matrix_path, 'rt') as f:
      matrix = scipy.io.mmread(f).tocsr() 
      matrix = matrix.T
      print(matrix.shape)
      print(a)
    # Load the barcodes (cell names) 
    with gzip.open(barcodes_path, 'rt') as f:
      barcodes = pd.read_csv(f, sep='\t', header=None)[0].values 

    # Load the features (gene names) 
    with gzip.open(features_path, 'rt') as f:
      features = pd.read_csv(f, sep='\t', header=None)[1].values 

    # Create a boolean mask to filter only the overlapping genes 
    gene_mask = np.isin(features, overlapping_genes) 
    matrix_filtered = matrix[:, gene_mask] 

    # Create an AnnData object for this sample 
    adata = sc.AnnData(X=matrix_filtered, var=pd.DataFrame(index=features[gene_mask]), obs=pd.DataFrame(index=barcodes)) 
    
    adata.var_names_make_unique()
    
    # Store the AnnData object 
    adata_list.append(adata) 

# Now concatenate all the AnnData objects into one 
combined_adata = adata_list[0].concatenate(adata_list[1:], join='inner') 

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

# Save the final combined AnnData object 
combined_adata.write(output_file) 

print(combined_adata)
print(combined_adata.var)
print(combined_adata.obs)
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 