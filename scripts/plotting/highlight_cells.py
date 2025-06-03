import scanpy as sc
import os

adata_dir = "/trinity/home/dmartinovicova/snakemake/data/combined_all/adata_clustered_bc_ing_denovo_complete.h5ad"
study = "OMIX920"
#cell = ['AP-microglia', 'a-microglia', 'h-microglia','i-microglia','s-mac 1','s-mac 2'] #GSE182109 celltype2
#cell = ['Myeloid'] #GSE182109
#cell = ['Mast cells', 'Monocytes', 'TAM 1','TAM 2', 'prol. TAM'] #GSE163120
#cell = ['MONO_DC', 'Microglia', 'T_NK'] #OMIX920
#cell = ['Dendritic cell', 'Myeloid', 'Granulocyte'] #Synapse
col_name = "omix_monodc_micro_tnk"

os.chdir("/trinity/home/dmartinovicova/snakemake/umaps_perstudy_inall/")

#Read adata
print('Reading adata...')
adata = sc.read_h5ad(adata_dir)
print(adata.obs)

# Filter for the Synapse cohort and Myeloid label

highlight_cells = adata[(adata.obs['Cohort'] == study) & 
                        (adata.obs['cell_type'].isin(cell))]


## Create a new column in 'adata.obs' to mark the cells for highlighting
#print('Creating new column...')
#adata.obs[col_name] = 'Other'
#adata.obs.loc[highlight_cells.obs.index, col_name] = highlight_cells.obs['cell_type2']

# Generate UMAP plot with highlighting
print('Generating umap...')
sc.pl.umap(highlight_cells, color='cell_type', 
            size=2,
            save=f'_inall_{col_name}.png')
