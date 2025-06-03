import os
import csv
import itertools
import scanpy as sc
import matplotlib.pyplot as plt

print('This is running')

# Load adata
print("Loading adata...")
adata = sc.read_h5ad("//net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_complete/data/combined_all/adata_final_l012.h5ad")
dataset_char = 'poster'
markers = "/net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_complete/data/annotations_markers/markers_poster.csv"

file = open(markers, "r")
markers_list = list(csv.reader(file, delimiter=","))
file.close()

# Flatten the list using itertools.chain
markers_list = list(itertools.chain.from_iterable(markers_list))

# Find markers that are in adata.var
print('Finding marker genes in adata.var...')
markers_in_adata = [gene for gene in markers_list if gene in adata.raw.var_names]
print(len(markers_in_adata))

# Set working directory (to save dotplots in desired folder)
print("Setting working directory...")
os.chdir("/net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_complete/dotplots")

# List different ways of groupby
groups = ['l2_celltype']#['Cohort', 'cell_type', 'cell_type2', 'leiden', 'level0_celltype','level1_celltype']

#Create dotplots
print("Creating dotplots...")
for group in groups:
  if group in adata.var_names or group in adata.obs.columns:
      sc.set_figure_params(scanpy=True, fontsize=14)
      sc.pl.dotplot(adata, markers_in_adata, show=False, groupby=group, save=f'_{dataset_char}_{group}.png', use_raw=True)
  else:
    print(f"Warning: {group} not found in adata.var_names or adata.obs.columns.")




## Subset adata for 'ref' and 'target' batches
#adata_ref = adata[adata.obs['batch'] == 'ref']
#adata_target = adata[adata.obs['batch'] == 'target']
#
## Create dot plot for 'ref' batch
#print("Creating dot plot for 'ref' batch...")
#sc.pl.dotplot(adata_ref, 
#              var_names=markers_in_adata, 
#              groupby='level1_celltype',  # You can replace 'cell_type' with any other categorical column
#              show=False, 
#              title='Marker Genes - ref', 
#              save=f'_{dataset_char}_ref_level1_celltype.png',
#              use_raw=True)
#
## Create dot plot for 'target' batch
#print("Creating dot plot for 'target' batch...")
#sc.pl.dotplot(adata_target, 
#              var_names=markers_in_adata, 
#              groupby='level1_celltype',  # You can replace 'cell_type' with any other categorical column
#              show=False, 
#              title='Marker Genes - target', 
#              save=f'_{dataset_char}_target_level1_celltype.png',
#              use_raw=True)








