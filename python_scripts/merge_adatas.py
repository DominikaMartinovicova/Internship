import os
import scanpy as sc

#change according to the particular dataset
folder_path = '/trinity/home/dmartinovicova/all_data/GSE163120_RAW/output'
study_identifier = '_GSE163120'
output_folder_file = 'merged/merged_samples_GSE163120.h5ad'

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
combined_adata = adatas[0].concatenate(*adatas[1:], join='inner')

#add study identifier
combined_adata.obs_names = combined_adata.obs_names + study_identifier

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


print(combined_adata.var_names, combined_adata.obs_names)

#create output file
output_file = os.path.join(folder_path, output_folder_file) #output directory
combined_adata.write(output_file)