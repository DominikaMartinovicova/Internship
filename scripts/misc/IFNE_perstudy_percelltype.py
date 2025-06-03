import os
import scanpy as sc
import pandas as pd
import glob

# Change working directory
os.chdir("/trinity/home/dmartinovicova/snakemake/dotplots/IFNE_perstudy_percelltype")

# Read cell labels file + get all the desired files from folder path
folder_path = "/trinity/home/dmartinovicova/snakemake/data/individual_for_IFNe(with_l1celltype)/"
files = glob.glob(os.path.join(folder_path, '*_clustered.h5ad'))
print(files)

cell_labels = pd.read_csv('/trinity/home/dmartinovicova/snakemake/data/annotations_markers/l1l0_celltypes.csv')
print(cell_labels)

# Loop over all datasets and create a new file containing cell labels + plot dotplot showing the expression of IFNE in cell types
for file in files:
  print(f'Processing {file}...')
  # Read the adata file
  adata = sc.read_h5ad(file)
  study = adata.obs['Cohort'][1]
  adata.obs['cells']=adata.obs.index
  # Add the cell labels
  adata.obs = adata.obs.merge(
      cell_labels[['cells','l1_celltype', 'l0_celltype']], 
      left_on='cells', 
      right_on='cells', 
      how='left')
  adata.obs.set_index('cells', inplace=True)
  adata.obs.index.name = None
  print(adata.obs.head())
  adata.write(f'/trinity/home/dmartinovicova/snakemake/data/individual_for_IFNe(with_l1celltype)/{study}_with_l1l0.h5ad')
  # Plot dotplot
  if 'IFNE' in adata.raw.var_names:
      sc.pl.dotplot(adata, 'IFNE', show=False, groupby='l1_celltype', save=f'_{study}_IFNE.png', use_raw =True)
  else:
    print(f"Warning: IFNE not found in {study} adata.var_names or adata.obs.columns.")

