import scanpy as sc
import pandas as pd

study_id = "_Miller_broadSC"


# Step 1: Load the TSV file
print("Loading annot...")
annotations_file = "/trinity/home/dmartinovicova/all_data/Miller_broadSC/GBM2_MetaData_Single_Cell_Portal.txt"
annotations_df = pd.read_csv(annotations_file, sep='\t', low_memory=False)
annotations_df = annotations_df.drop(0)

# Check the structure of your annotations
print(annotations_df.head())

# Step 2: Load the AnnData object
print("Loading adata...")
adata = sc.read_h5ad("/trinity/home/dmartinovicova/all_data/Miller_broadSC/adata_Miller_broadSC.h5ad")

# Check the structure of the AnnData object
print(adata)


# Remove study identifier and create cell_barcodes column
print("Creating cell_barcode...")
adata.obs['cell_barcode'] = adata.obs.index.str.replace(study_id, '', regex=False)
new_columns = ['cell_barcode'] + [col for col in adata.obs.columns if col != 'cell_barcode']
adata.obs = adata.obs[new_columns]
print(adata.obs)

new_columns = ['cell_barcode'] + [col for col in adata.obs.columns if col != 'cell_barcode']
adata.obs = adata.obs[new_columns]

print(adata.obs)

# Step 3: Ensure the cell barcodes are in the correct format
# Set the index of the annotations DataFrame to the 'cell_barcode' column
#annotations_df['cell_barcode'] = annotations_df['NAME'].apply(lambda x: x.split('_')[1])
#annotations_df.set_index('cell_barcode', inplace=True)
#print(annotations_df.head())

# Step 4: Merge annotations into AnnData
# The 'obs' attribute of adata should have a column that corresponds to the barcodes
# Assuming adata.obs contains a column named 'cell_barcode' for matching
print("Adding annot...")
adata.obs = adata.obs.merge(
    annotations_df[['NAME','Total_Tumor_Annotation']],  # Select only the 'cell_state' column
    left_on="cell_barcode",  # Match this with 'cell_barcode' in adata.obs
    right_on="NAME",  # Use "NAME" from annotations_df for matching
    how='left'  # Keep all rows from adata.obs
)
print(adata.obs)

# Remove batch column
if 'batch' in adata.obs.columns:
  adata.obs = adata.obs.drop(columns = ['batch'])

adata.obs.set_index('cell_barcode', inplace=True)
adata.obs_names = adata.obs_names + study_id

# Optionally, rename the new column to something descriptive
adata.obs.rename(columns={'Total_Tumor_Annotation': 'cell_type'}, inplace=True)
print(adata.obs)
# Step 5: Save the updated AnnData object
adata.write("/trinity/home/dmartinovicova/all_data/Miller_broadSC/adata_Miller_broadSC.h5ad")
