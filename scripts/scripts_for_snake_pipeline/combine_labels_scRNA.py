#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# combine_labels_scRNA.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Combine the refined labels from myeloid, immune and all adatas and add l0 and l1 cell type labels to the united adata and create final adata file.
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl) 
#
# Usage:
     """
     python3 scripts/combine_labels_scRNA.py \
     -i_all {input.adata_clustered} \
     -i_i {input.immune_clustered} \
     -i_m {input.myeloid_clustered} \
     -o {output.adata_final} \
     -cache_dir {params.cache_dir} 
     """

#================================================================================
# 0 Import libraries
#================================================================================
import argparse
import scanpy as sc
import pandas as pd

#================================================================================
# 1 Parse arguments
#================================================================================
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 combine_labels_scRNAseq.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        'Combine refined labels in two levels and create final adata file')
    parser.add_argument('-i_all', help='path to input united adata file',
                        dest='input_all',
                        type=str)
    parser.add_argument('-i_i', help='path to input immune subset file',
                        dest='input_immune',
                        type=str)    
    parser.add_argument('-i_m', help='path to input myeloid subset file',
                        dest='input_myeloid',
                        type=str)
    parser.add_argument('-o', help='path to final adata file',
                        dest='output',
                        type=str)
    parser.add_argument('-cache_dir', help='Directory to save cache',
                        dest='cache_dir',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

# Set directory to save cache
sc.settings.verbosity = 3
sc.settings.cachedir = args.cache_dir

#================================================================================
# 2 Read files
#================================================================================
print("Reading the files...")
myeloid = sc.read_h5ad(args.input_myeloid)
immune = sc.read_h5ad(args.input_immune)
all = sc.read_h5ad(args.input_all)

#================================================================================
# 3 Map the phenotypes
#================================================================================
#-------------------------------------------------------------------------------
# 3.1 Map the new cell labels to their corresponding cells in the united adata 
#-------------------------------------------------------------------------------
# Create new columns with the labels from immune and myeloid
all.obs['denovo_myeloid'] = all.obs.index.map(myeloid.obs['de_novo'])
all.obs['denovo_immune'] = all.obs.index.map(immune.obs['de_novo'])
print(all.obs)

# Create column with malignant, fibro, peri, and endo labeled from the original studies
kept_labels = ['Malignant', 'Fibroblast', 'Pericyte', 'Endothelial']
all.obs['kept_labels'] = all.obs['level1_celltype'].apply(lambda x: x if x in kept_labels else pd.NA)
print(all.obs)

#-------------------------------------------------------------------------------
# 3.2 Combine labels together to create a complete column with all the cell labels 
#-------------------------------------------------------------------------------
# l1_celltype
all.obs['l1_celltype'] = all.obs['kept_labels'].combine_first(all.obs['denovo_immune']).combine_first(all.obs['denovo_myeloid']).combine_first(all.obs['de_novo'])

# Rename Microglia and Macrophages to TAM for l0_celltype
TAM_subtypes=['Microglia', 'Macrophage']
all.obs['l0_celltype'] = all.obs['l1_celltype'].apply(lambda x: 'TAM' if x in TAM_subtypes else x)
print(all.obs)

#================================================================================
# 4 Remove unwanted columns and NA values
#================================================================================
# List the columns you want to remove from adata.obs 
columns_to_remove = ['kept_labels', 
                      'de_novo',
                      'denovo_myeloid',
                      'denovo_immune'] 

# Remove 
adata.obs = adata.obs.drop(columns=columns_to_remove)

#================================================================================
# 5 Write file
#================================================================================
all.write(args.output)

