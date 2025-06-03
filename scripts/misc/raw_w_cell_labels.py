import scanpy as sc

# Load your AnnData objects
adata1 = sc.read("/net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_complete/data/combined_all/adata_final_l012.h5ad")  # This file contains the labels you want to transfer
adata2 = sc.read("/net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_complete/data/combined_all/adata_combined.h5ad")  # This is the target AnnData

print(adata1.obs)
print(adata2.obs)


# Remove '-ref' and '-target' from the obs_names in adata1
adata1.obs_names = adata1.obs_names.str.replace('-ref', '').str.replace('-target', '')
print(adata1.obs)
# Check if the cell barcodes in adata2 exist in adata1
matching_barcodes = adata2.obs_names.isin(adata1.obs_names)

if matching_barcodes.all():
    print("All cell barcodes in adata1 are present in adata2.")
else:
    print("Some cell barcodes in adata1 are missing from adata2.")


adata2_subset = adata2[matching_barcodes, :]

# Transfer the label column based on matching barcodes in adata2
adata2_subset.obs['l1_celltype'] = adata1.obs['l1_celltype'].loc[adata2_subset.obs_names]
adata2_subset.obs['l2_celltype'] = adata1.obs['l2_celltype'].loc[adata2_subset.obs_names]
adata2_subset.obs['Diagnosis_label'] = adata1.obs['Diagnosis_label'].loc[adata2_subset.obs_names]

print(adata1)
print(adata2)
print(adata2_subset)
print(adata2_subset.obs)

# Verify that labels were added correctly
print(adata2_subset.obs[['l2_celltype']].head())

adata2_subset.write("/net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_complete/data/combined_all/adata_raw_w_labels_l012.h5ad")
