# All adata files should contain: 
# n_genes_by_counts 
# total_counts
# total_counts_mt
# pct_counts_mt
# cell_id - original obs from the downloaded dataset
# cell_type (only if annotated)
# IDH - Mut/WT
# Diagnosis 
# Cohort
# Primary_Recurrent
# case_id - sample id
#
# Adjust creating and transfering columns as necessary


import scanpy as sc
import pandas as pd

#-------------------------------------------------------------------
# 0. Prepare variables
#-------------------------------------------------------------------

metadata_file = "/trinity/home/dmartinovicova/all_data/GSE182109_bSC/Meta_GBM_2.txt"
metadata_file2 = "/trinity/home/dmartinovicova/all_data/GSE182109_RAW/clinical_info.csv"
adata_file = "/trinity/home/dmartinovicova/all_data/GSE182109_bSC/GSE182109_bSC.h5ad"
output_file = "/trinity/home/dmartinovicova/all_data/GSE182109_bSC/GSE182109_bSC_complete.h5ad"
study_id = '_GSE182109_broadSC'
common_col =  'cell_id'   # matching column in metadata file and adata to help the transfer
common_col2 =  'case_id'
cohort = 'GSE182109'

meta = pd.read_csv(metadata_file, sep=',')
meta2 = pd.read_csv(metadata_file2)
print(meta)
print(meta.columns)
adata = sc.read_h5ad(adata_file)

#-------------------------------------------------------------------
# 1. Create missing columns
#-------------------------------------------------------------------

# Create column that will later become index
adata.obs['later_index'] = adata.obs.index.tolist()

print(adata)
print(adata.obs)

#Create column cell_id (original obs)
adata.obs['cell_id'] = [name.replace(study_id,'') for name in adata.obs_names]

# Cohort
print('Adding Cohort...')
adata.obs['Cohort'] = cohort


#-------------------------------------------------------------------
# 2. Transfer columns from metadata/clinical data
#-------------------------------------------------------------------

# Add columns from metadata to adata
print('Adding metadata...')

adata.obs = adata.obs.merge(
    meta[[common_col, 'cell_type','cell_type2', 'case_id']],
    left_on=common_col,      # Column in adata.obs to match
    right_on=common_col,  # Column in df to match
    how='left'              # Retain all rows from adata.obs
)

print(adata.obs)

# Add columns from metadata2 to adata
print('Adding metadata2...')

adata.obs = adata.obs.merge(
    meta2[[common_col2, 'IDH', 'Primary_Recurrent','Diagnosis']],
    left_on=common_col2,      # Column in adata.obs to match
    right_on=common_col2,  # Column in df to match
    how='left'              # Retain all rows from adata.obs
)
print(adata)
print(adata.obs)

# Now set 'later_index' as the index
adata.obs.set_index('later_index', inplace=True)
adata.obs.index.name = None


print(adata)
print(adata.obs)

# Save complete adata file
adata.write(output_file)