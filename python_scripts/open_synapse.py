import os
import scanpy as sc

#path to the adata file
folder_path = '/trinity/home/dmartinovicova/all_data/Verhaak_synapse'

#study identifier
study_id = '_synapse'

#read the adata file and add study identifier
adata_path = os.path.join(folder_path, 'analysis_scRNAseq_tumor_counts.h5ad')
adata = sc.read_h5ad(adata_path)
adata.obs_names = adata.obs_names + study_id

#path to output file - adata with study identifier
output_file = os.path.join(folder_path, 'output/Verhaak_synapse.h5ad')
adata.write(output_file)

print(adata.obs_names, adata.var_names)