import pandas as pd
import scanpy as sc

study_id = "_OMIX_A"

# Load the AnnData object (e.g., adata)
adata = sc.read_h5ad("/trinity/home/dmartinovicova/all_data/OMIX920/output/Astro3.h5ad")

# Load the cell annotations CSV
annotations = pd.read_csv("/trinity/home/dmartinovicova/all_data/OMIX920/output/Astro_annotations.csv")

# Ensure 'cell_barcodes' is set as the index to match the AnnData object
annotations.set_index("cell_barcodes", inplace=True)

# Add cell annotations to adata by matching cell barcodes
# Here we assume that `adata.obs_names` are the cell barcodes in `adata`
adata.obs = adata.obs.join(annotations, how="left")

# Add the study identifier
adata.obs_names = adata.obs_names + study_id

# Now `adata.obs` has an additional column 'cell_type' with annotations
print(adata.obs.head())

adata.write("/trinity/home/dmartinovicova/all_data/OMIX920/output/Astro3_ann.h5ad")