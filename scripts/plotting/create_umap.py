import os
import scanpy as sc
import matplotlib.pyplot as plt

print('This is running')
# Load adata
print("Loading adata...")
adata = sc.read_h5ad("/net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_test/data/phenotyping/adata_final.h5ad")
dataset_char = 'final'

# Set working directory (to save umaps in desired folder)
print("Setting working directory...")
os.chdir("/net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_test/")

## Define a color dictionary
#color_map = {'B_cell':'royalblue',
#              'CD4_Tcell':'darkorange',
#              'CD8_Tcell':'mediumseagreen',
#              'Endothelial':'crimson',
#              'Fibroblast':'darkviolet',
#              'Malignant':'sienna',
#              'Monocyte':'darkkhaki',
#              'NK_cell':'palevioletred',
#              'Oligodendrocyte':'mediumturquoise',
#              'Pericyte':'lightblue',
#              'Plasma_B':'dimgrey',
#              'TAM':'lightpink',
#              'T_reg':'mediumpurple',
#              'cDC':'chocolate',
#              'pDC':'firebrick',
#              'Microglia':'lightcoral',
#              'Macrophage':'steelblue'
#              }


# List desired attributes to plot umap
#attributes = ['Diagnosis','Cohort', 'cell_type', 'cell_type2', 'leiden', 'level0_celltype','level1_celltype', 'de_novo', 'denovo_myeloid','denovo_immune','complete_denovo','complete_denovo_myeloid', 'kept_labels','kept_denovo','kept_denovo_myeloid']
  
#attributes = ['l0_celltype',
#  'Cohort', 'cell_type', 'cell_type2', 'leiden', 'harmonized_celltype',
#  'CD3D', 'CD3E', 'CD3G','CD8A','CD4', 'CTLA4','FOXP3','IL2RA', 'IL7R', 'GZMK',
#  'CD79A','MS4A1', 'SDC1', 
#  'GNLY', 'NKG7', 'KLRC1', 'CD56', 'NCAM1',
#  'P2RY12', 'TMEM119', 'TGFBR1', 'IBA1','SAL11', 'CD45', 'CX3CR1', 'TREM2', 'CCR2', 'GPNMB',
#  'PLP1', 'MOG', 'MAG', 'MBP',
#  'LYZ', 'S100A9', 'S100A8', 'FCN1',
#  'CD68', 'CD14','CD11B',
#  'DCN', 'PDGFRB', 'VCAN', 'FBLN2', 'COL1A1',
#  'RGS5', 'TPM2',
#  'VWF','CLDN5',
#  'FCER1A', 'TCF4','LILRA4', 'CLEC4C', 'BDCA2', 'CD1C',
#  'CD86', 'C1QA', 'FCGR1A', 'TNF', 'IL6', 'IL12', 'IL-6', 'CD163', 'CD204', 'CD206', 'IRF5', 'IL10', 'NFKB','H2-EB1', 'H2EB1', 
#  'MS4A2','HDC','CTSG', 'FCER1B', 'CG','IGER','APY',
#  'SOX2', 'STMN1', 'MKI67',
#  'HLA-DPB1','ANPEP',
#  'IFNE', 'IFNB', 'IFNG','IFN-E', 'IFN-B', 'IFN-G']

attributes = ['l1_celltype', 'l2_celltype','Diagnosis_label','harmonized_celltype','l0_celltype']

#attributes = ['NDRG2','GFAP','S100B','A2B5','ALDH1L1','EAAT1','EAAT2','HES-1','HES1','SOX9','PEA-15','PEA15','VIM','CD44','CLU','ITGA6','ID3','GJA1','AQP4','SLC1A3','HEPACAM']

##Create umaps
#print("Creating umaps...")
#for attribute in attributes:
#  if attribute in adata.raw.var_names or attribute in adata.obs.columns:
#    if attribute == 'leiden':
#      sc.pl.umap(adata, color=attribute, palette=color_map, show=False, save=f'_{dataset_char}_{attribute}.png', legend_loc='on data', use_raw=True)
#    else:
#      sc.pl.umap(adata, color=attribute, palette=color_map, show=False, save=f'_{dataset_char}_{attribute}.png', use_raw=True)
#  else:
#    print(f"Warning: {attribute} not found in adata.var_names or adata.obs.columns.")
#  #plt.savefig(f'{dataset_char}_{attribute}.png', dpi=300)
  
  
#Create umaps
print("Creating umaps...")
for attribute in attributes:
  if attribute in adata.raw.var_names or attribute in adata.obs.columns:
    if attribute == 'leiden':
      sc.pl.umap(adata, color=attribute, show=False, save=f'_{dataset_char}_{attribute}.png', legend_loc='on data', use_raw=True)
    else:
      sc.pl.umap(adata, color=attribute, show=False, save=f'_{dataset_char}_{attribute}.png', use_raw=True)
  else:
    print(f"Warning: {attribute} not found in adata.var_names or adata.obs.columns.")