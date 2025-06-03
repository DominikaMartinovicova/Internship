import os
import scanpy as sc
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix


folder_path = "/trinity/home/dmartinovicova/all_data/GSE163120_RAW/"   # Path to the folder containing the files
study_id =  '_GSE163120'    # Study identificator
file_ending = 'filtered.gene.bc.matrix.csv.gz'  # File ending to recognize relevant files in folder
output_file =  '/trinity/home/dmartinovicova/all_data/GSE163120_RAW/test_gse163120.h5ad' # Directory and name of the output file

#-----------------------------------------------------------------------------
# 1. Get list of all the samples
#-----------------------------------------------------------------------------

# List all files ending with file_ending
file_list = [f for f in os.listdir(folder_path) if f.endswith(file_ending)]
print(file_list)
print(len(file_list))

#-----------------------------------------------------------------------------
# 2. Create adata for each sample
#-----------------------------------------------------------------------------

# Initialize a list to store AnnData objects
adatas = []

# Process each file
print('Processing files...')
for file in file_list:
    print('Processing file:' + str(file))
    file_path = os.path.join(folder_path, file)
    
    # Load the data (adjust sep according to .csv or .tsv file)
    data = pd.read_csv(file_path, compression='gzip', index_col=0) 
    
    # Convert data to a sparse matrix if it's large and transpose
    sparse_data = csr_matrix(data.values).T
    
    # Create AnnData object
    adata = ad.AnnData(X=sparse_data, var=pd.DataFrame(index=data.index), obs=pd.DataFrame(index=data.columns))
    
    # Optional based on desired case_id
    # adata.obs['case_id'] = file.split('_')[1]  
    
    adatas.append(adata)

#-----------------------------------------------------------------------------
# 3. Find common genes and subset for them 
#-----------------------------------------------------------------------------

# Ensure common genes across all AnnData objects
print('Finding common genes...')
common_genes = set(adatas[0].var_names)
for adata in adatas[1:]:
    common_genes &= set(adata.var_names)
common_genes = sorted(list(common_genes))  # Sort for consistency

# Subset each AnnData to only include the common genes
adatas = [adata[:, common_genes] for adata in adatas]

# Combine all AnnData objects
print('Combining adata...')
combined_adata = ad.concat(adatas, axis=0, join='inner')  # Join on common genes

#-----------------------------------------------------------------------------
# 4. Add information 
#-----------------------------------------------------------------------------

# Identify mitochondrial genes and calculate qc
combined_adata.var['mt'] = combined_adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(combined_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Remove batch column
if 'batch' in combined_adata.obs.columns:
  combined_adata.obs = combined_adata.obs.drop(columns = ['batch'])

# Add study identifier to keep track of the original dataset of the cells
combined_adata.obs_names = combined_adata.obs_names + study_id

# Inspect the combined AnnData
print(combined_adata)
print(combined_adata.var)
print(combined_adata.obs)

# Save the combined AnnData
combined_adata.write(output_file)
