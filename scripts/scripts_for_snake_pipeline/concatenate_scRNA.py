#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# concatenate_scRNA.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Concatenate adata files to one united adata object + cell type label harmonization
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
"""
        python3 scripts/concatenate_scRNA.py \
        -i {input.processed_adata_dir} \
        -phen_map {input.phenotype_map} \
        -o {output.adata_concatenated} \
        -cache_dir {params.cache_dir} 
"""

#================================================================================================================
# 0 Import libraries
#================================================================================================================
import argparse
import anndata as ad
import os
import scanpy as sc
import numpy as np
import argparse
import pandas as pd

#================================================================================================================
# 1 Parse arguments
#================================================================================================================

def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog='python3 concatenate_adatas.py',
                                     formatter_class=argparse.RawTextHelpFormatter, description='Concatenate scRNA adata files')
    parser.add_argument('-i', help='Input adata filtered directory', dest='input', type=str, required=True)
    parser.add_argument('-o', help='Path to output file', dest='output', type=str, required=True)
    parser.add_argument('-cache_dir', help='Directory to save cache', dest='cache_dir', type=str, required=True)
    parser.add_argument('-phen_map', help='File with cell types and mapping', dest='phen_map', type=str, required=True)
    args = parser.parse_args()
    return args

args = parse_args()

# set directory to save cache
sc.settings.verbosity = 3
sc.settings.cachedir = args.cache_dir


#================================================================================================================
# 2 Concatenate adatas
#================================================================================================================
# Get a list of all files in the directory 
data_dir = args.input
files = os.listdir(data_dir) 
print(files)

# Load each adata file
adata_list = [ad.read_h5ad(data_dir + '/' + f) for f in files]
print(adata_list)

# Find common genes across all AnnData files
common_genes = set(adata_list[0].var_names)

for adata in adata_list[1:]:
    common_genes.intersection_update(adata.var_names)
    
common_genes = sorted(list(common_genes))
print(len(common_genes))

# Filter each AnnData object to only include common genes
adata_list = [adata[:, list(common_genes)].copy() for adata in adata_list]

# Remove any 'batch' column in `obs` for each AnnData object
for adata in adata_list:
    if "batch" in adata.obs:
        adata.obs = adata.obs.drop(columns=["batch"])

# Concatenate along the `obs` axis
print('Concatenating...')
adata_combined = ad.concat(adata_list, join="outer")

# Change values to integers (float not accepted)
adata_combined.X = np.round(adata_combined.X).astype(int)

print(adata_combined)
print(adata_combined.obs)

#================================================================================================================
# 3 Map cell types (level0, level1, cell_type, cell_type2)
#================================================================================================================
if 'cell_type' in adata_combined.obs.columns:
  mapping_df = pd.read_csv(args.phen_map)
  
  # Create dictionaries to store the matched detailed phenotype with category
  print('Creating dictionaries...')
  category = dict(zip(mapping_df['cell_type2'], mapping_df['level1_celltype']))
  category2 = dict(zip(mapping_df['cell_type2'], mapping_df['level0_celltype']))
  
  # Create cell_type2 column
  adata_combined.obs['cell_type2'] = adata_combined.obs['cell_type2'].astype('object')
  adata_combined.obs['cell_type'] = adata_combined.obs['cell_type'].astype('object')
  adata_combined.obs['cell_type2'] = adata_combined.obs['cell_type2'].fillna(adata_combined.obs['cell_type'])
  
  # Map the 'cell_type2' to the new, more general category and create a new column 'level1' and 'level0'
  print('Mapping...')
  adata_combined.obs['level1_celltype'] = adata_combined.obs['cell_type2'].map(category)
  adata_combined.obs['level0_celltype'] = adata_combined.obs['cell_type2'].map(category2)
  
  # Remove indifferent cell types
  print(adata_combined)
  adata_combined = adata_combined[~adata_combined.obs['cell_type2'].isin(['Other', 'Proliferating']), :]
  
  # Optionally, check the result
  print(adata_combined)
  print(adata_combined.obs)
  print(adata_combined.obs['cell_type'].value_counts())
 
#================================================================================================================
# 3 Write adata file
#================================================================================================================
print('Writing file...')
adata_combined.write(args.output)

