import os
import scanpy as sc

#folder to adata files
folder_path = '/trinity/home/dmartinovicova/all_data/GSE182109_RAW/output'

#list all the adata files from the folder
adata_files = [f for f in os.listdir(folder_path) if f.endswith('.h5ad')]

adatas = []

#loop through the adata files and add adata file to the list

for adata_file in adata_files:
  adata_path = os.path.join(folder_path, adata_file)
  adata = sc.read_h5ad(adata_path)
  adatas.append(adata)
  print(adata)


#check for common genes throughout the adata files
common_genes = set(adatas[0].var_names)
for adata in adatas[1:]:
  common_genes.intersection_update(adata.var_names)

#sort the common_genes
common_genes = sorted(common_genes)

#take subset of the adata dataset containing only common genes info
for i in range(len(adatas)):
  adatas[i] = adatas[i][:, common_genes]
  print(adatas[i])

#merge the filtere adata files
merged_adata = adatas[0].concatenate(*adatas[1:], join='inner')

#add study identifier
merged_adata.obs_names = merged_adata.obs_names + '_GSE182109'

print(merged_adata.obs_names, merged_adata.var_names)

#create output file
output_file = os.path.join(folder_path, 'merged/merged_samples_GSE182109.h5ad') #output directory
merged_adata.write(output_file)