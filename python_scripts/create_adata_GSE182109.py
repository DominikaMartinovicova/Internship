import os
import scanpy as sc
import gzip

#path to the folder with files
folder_path = '/trinity/home/dmartinovicova/all_data/GSE182109_RAW'
output_folder_path = '/trinity/home/dmartinovicova/all_data/GSE182109_RAW/output'

#list all the files from the folder with .gz
files = [f for f in os.listdir(folder_path) if f.endswith('.gz')]

#put files corresponding to the same sample under one key in dictionary
sample_files = {}

#loop over the files and get the filepaths
for f in files:
  
  #get sample name
  if f.endswith('_barcodes.tsv.gz'):
    sample_name = f.replace('barcodes.tsv.gz', '')
  elif f.endswith('_features.tsv.gz'):
    sample_name = f.replace('features.tsv.gz', '')
  elif f.endswith('_matrix.mtx.gz'):
    sample_name = f.replace('matrix.mtx.gz', '')
    
  #create a dictionary key with sample name and add keys for the corresponding files
  if sample_name not in sample_files:
    sample_files[sample_name] = {'barcodes':None, 'features':None, 'matrix':None}
  
  #assign filepath to the correct files
  if f.endswith('_barcodes.tsv.gz'):
    sample_files[sample_name]['barcodes'] = os.path.join(folder_path, f)
  elif f.endswith('_features.tsv.gz'):
    sample_files[sample_name]['features'] = os.path.join(folder_path, f)
  elif f.endswith('_matrix.mtx.gz'):
    sample_files[sample_name]['matrix'] = os.path.join(folder_path, f)
    
    
#get files and create adata file
for sample_name, file_path in sample_files.items():
  barcodes_path = file_path['barcodes']
  features_path = file_path['features']
  matrix_path = file_path['matrix']
  
  adata = sc.read_10x_mtx(folder_path, var_names='gene_symbols', prefix = sample_name,cache=False)
  
  #save adata file
  output_file = os.path.join(output_folder_path, f"{sample_name}.h5ad")
  adata.write(output_file)
  
  print(f"{sample_name} saved to {output_file}")

#print(adata.obs_names) #to check for added study identifier
