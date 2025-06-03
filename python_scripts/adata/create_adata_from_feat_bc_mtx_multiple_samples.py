import os  
import gzip
import scanpy as sc 
import scipy.io 
import pandas as pd 
import numpy as np 
from scipy.sparse import csr_matrix 
import anndata as ad

data_dir =      # Path to the folder containing the files
study_id =      # Study identificator
#file_ending =             # File ending to recognize relevant files in folder
output_file =   # Directory and name of the output file

#-----------------------------------------------------------------------------
# 1. Get list of all the samples
#-----------------------------------------------------------------------------

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
print(sample_ids) 
print(len(sample_ids))

#-----------------------------------------------------------------------------
# 2. Create adata file for each sample
#-----------------------------------------------------------------------------

# To store the data for each sample 
adatas = [] 

# Iterate again to extract the data only for overlapping genes and concatenate it 
for sample_id in sample_ids: 
    print('Processing sample' + str(sample_id))
    # Define paths for the matrix, barcodes, and features files 
    matrix_path = os.path.join(data_dir, f'{sample_id}matrix.mtx.gz') 
    barcodes_path = os.path.join(data_dir, f'{sample_id}barcodes.tsv.gz') 
    features_path = os.path.join(data_dir, f'{sample_id}features.tsv.gz')  

    # Load the matrix (gene expression data) 
    with gzip.open(matrix_path, 'rt') as f:
      matrix = scipy.io.mmread(f).tocsr() 
      matrix = matrix.T
      print(matrix.shape)
    
    # Load the barcodes (cell names) 
    barcodes = pd.read_csv(barcodes_path, sep='\t', compression='gzip', header=None)[0].values
    
    # Load the features (gene names) 
    features = pd.read_csv(features_path, sep='\t', compression='gzip', header=None)[1].values

    # Create an AnnData object for this sample 
    adata = sc.AnnData(X=matrix, var=pd.DataFrame(index=features), obs=pd.DataFrame(index=barcodes)) 
    
    adata.var_names_make_unique()
    
    # Optional depending on the desired case_id 
    #adata.obs['case_id'] = sample_id.split('_')[1]
    
    # Store the AnnData object 
    adatas.append(adata) 


#-----------------------------------------------------------------------------
# 3. Identify common genes and subset for them
#-----------------------------------------------------------------------------

# Ensure common genes across all AnnData objects
print('Finding common genes...')
common_genes = set(adatas[0].var_names)
for adata in adatas[1:]:
    common_genes &= set(adata.var_names)
common_genes = sorted(list(common_genes))  # Sort for consistency
print('Number of common genes:' + str(len(common_genes)))

# Subset each AnnData to only include the common genes
adatas = [adata[:, common_genes] for adata in adatas]

# Combine all AnnData objects
print('Combining adata...')
combined_adata = ad.concat(adatas, axis=0, join='inner')  # Join on common genes


#-----------------------------------------------------------------------------
# 4. Add information
#-----------------------------------------------------------------------------

# Add the study identifier
print('Adding study_id')
combined_adata.obs_names = combined_adata.obs_names + study_id

# Identify mitochondrial genes and calculate qc
combined_adata.var['mt'] = combined_adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(combined_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Remove batch column
if 'batch' in combined_adata.obs.columns:
  combined_adata.obs = combined_adata.obs.drop(columns = ['batch'])

# Inspect the combined adata
print(combined_adata)
print(combined_adata.var)
print(combined_adata.obs)

# Save the final combined AnnData object 
combined_adata.write(output_file) 


## Adjust cell_id
#combined_adata.obs['case_id'] = combined_adata.obs['case_id'].apply(lambda x: '-'.join(x.split('-')[:2]))






 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 