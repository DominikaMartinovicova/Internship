import scanpy as sc
import pandas as pd

study_id = "_synapse"


# Step 1: Load the TSV file
annotations_file = "/trinity/home/dmartinovicova/snakemake/data/annotations/Verhaak_synapse_metadata.tsv"
annotations_df = pd.read_csv(annotations_file, sep='\t')

# Check the structure of your annotations
print(annotations_df.head())

# Step 2: Load the AnnData object
adata = sc.read("/trinity/home/dmartinovicova/all_data/Verhaak_synapse/output/Verhaak_synapse_2.h5ad")

# Check the structure of the AnnData object
print(adata)

# Remove study identifier and create cell_barcodes column
adata.obs['cell_barcode'] = adata.obs.index.str.replace(study_id, '', regex=False)
new_columns = ['cell_barcode'] + [col for col in adata.obs.columns if col != 'cell_barcode']
adata.obs = adata.obs[new_columns]

print(adata.obs)

# Step 3: Ensure the cell barcodes are in the correct format
# Set the index of the annotations DataFrame to the 'cell_barcode' column
annotations_df.set_index('cell_barcode', inplace=True)

# Step 4: Merge annotations into AnnData
# The 'obs' attribute of adata should have a column that corresponds to the barcodes
# Assuming adata.obs contains a column named 'cell_barcode' for matching
adata.obs = adata.obs.merge(
    annotations_df[['cell_state']],  # Select only the 'cell_state' column
    left_on="cell_barcode",  # Match this with 'cell_barcode' in adata.obs
    right_index=True,  # Use index from annotations_df for matching
    how='left'  # Keep all rows from adata.obs
)
print(adata.obs)

# Remove batch column
if 'batch' in adata.obs.columns:
  adata.obs = adata.obs.drop(columns = ['batch'])

# Optionally, rename the new column to something descriptive
adata.obs.rename(columns={'cell_state': 'cell_type'}, inplace=True)
print(adata.obs)
# Step 5: Save the updated AnnData object
adata.write("/trinity/home/dmartinovicova/snakemake/data/adata/Verhaak_synapse.h5ad")
