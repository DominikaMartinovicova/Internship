import scanpy as sc
import pandas as pd

adata_file = "/trinity/home/dmartinovicova/snakemake/data/combined_3ann_gse182109_2un/adata_combined.h5ad"
csv_file = "/trinity/home/dmartinovicova/snakemake/data/annotations/mapping.csv"

adata = sc.read_h5ad(adata_file)
mapping_df = pd.read_csv(csv_file)

# Create dictionaries to store the matched detailed phenotype with category
print('Creating dictionaries...')
category = dict(zip(mapping_df['cell_type2'], mapping_df['level1_celltype']))
category2 = dict(zip(mapping_df['cell_type2'], mapping_df['level0_celltype']))

# Create cell_type2 column
adata.obs['cell_type2'] = adata.obs['cell_type2'].astype('object')
adata.obs['cell_type'] = adata.obs['cell_type'].astype('object')
adata.obs['cell_type2'] = adata.obs['cell_type2'].fillna(adata.obs['cell_type'])


# Map the 'cell_type2' to the new, more general category and create a new column 'level1' and 'level0'
print('Mapping...')
adata.obs['level1_celltype'] = adata.obs['cell_type2'].map(category)
adata.obs['level0_celltype'] = adata.obs['cell_type2'].map(category2)


# Optionally, check the result
print(adata)
print(adata.obs)


# Save into file
adata.write("/trinity/home/dmartinovicova/snakemake/data/combined_3ann_gse182109_new/combined_adata_mapped.h5ad")