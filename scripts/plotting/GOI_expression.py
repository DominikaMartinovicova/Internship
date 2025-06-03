import matplotlib.pyplot as plt
import scanpy as sc
import os

# Change working directory
os.chdir("/trinity/home/dmartinovicova/snakemake/umaps/CDKN2A_IFNE/")

# Identify files with clustered data
print("Idetifying clustered adatas...")
folder_path = "/trinity/home/dmartinovicova/snakemake/data/combined_all/"
all_files = os.listdir(folder_path)
clustered_adatas = [f for f in all_files if f.endswith('clustered_bc_ing.h5ad')]
print(clustered_adatas)

# Genes of interest
goi = ['CDKN2A', 'CDKN2B', 'IFNE']

#===========================================================
# Check expression of genes of interest and plot into umap
#===========================================================
#-----------------------------------------------------------
# Loop over all clustered files to check individually for goi expression
#-----------------------------------------------------------
print("Looping over clustered adatas...")
for file_name in clustered_adatas: 
    print("Checking in ", file_name)
    # Get the clustered file 
    file_path = os.path.join(folder_path, file_name)
    adata = sc.read_h5ad(file_path)
    
    # Save study name by omitting '_clustered.h5ad'
    study = file_name.replace('_clustered_bc_ing.h5ad','')
    
    # Loop over all genes of interest
    for g in goi:
        print("Checking for ", g)
        #Check if gene is in adata.raw.var_names
        if g in adata.raw.var_names:

            # Store values of gene expression
            g_exp = adata.raw[:,g].X.toarray()
            
            # Create new column with gene expression as dichotomous variable
            adata.obs[f'{g}_expressed'] = (g_exp>0).astype(int)
            
            # Plot UMAP
            sc.pl.umap(adata, color=f'{g}_expressed',save=f'_{study}_{g}.png', use_raw=True)
            print(f"Dichotomous variable {g}_expressed added to adata.obs.")
        
        else:
            print(f"Gene {g} not found in the {study}.")
           
     
    
