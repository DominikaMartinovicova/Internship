import pandas as pd
import scipy.io
import scanpy as sc
import numpy as np


#load matrices
matrix_1 = scipy.io.mmread('/trinity/home/dmartinovicova/all_data/Neftel_broadSC/IDHwtGBM.processed.10X.counts.mtx').tocsr().transpose()
matrix_2 = scipy.io.mmread('/trinity/home/dmartinovicova/all_data/Neftel_broadSC/IDHwtGBM.processed.10X.counts.2.mtx').tocsr().transpose()

#load barcodes
barcodes_1 = pd.read_csv('/trinity/home/dmartinovicova/all_data/Neftel_broadSC/cells1.proc.tsv', header = None, sep='\t')
barcodes_2 = pd.read_csv('/trinity/home/dmartinovicova/all_data/Neftel_broadSC/cells2.proc.tsv', header = None, sep='\t')

#load genes
features_1 = pd.read_csv('/trinity/home/dmartinovicova/all_data/Neftel_broadSC/genes1.tsv', header = None, sep='\t')
features_2 = pd.read_csv('/trinity/home/dmartinovicova/all_data/Neftel_broadSC/genes2.tsv', header = None, sep='\t')
print(features_1)
print(barcodes_1)
print(features_2)
print(barcodes_2)

#create 2 adata files
adata_1 = sc.AnnData(X=matrix_1)
adata_1.obs_names = barcodes_1[0].values
adata_1.var_names = features_1[1].values
print('first adata created')

adata_1.write('/trinity/home/dmartinovicova/all_data/Neftel_broadSC/output/adata_1.h5ad')

adata_2 = sc.AnnData(X=matrix_2)
adata_2.obs_names = barcodes_2[0].values
adata_2.var_names = features_2[1].values

adata_2.write('/trinity/home/dmartinovicova/all_data/Neftel_broadSC/output/adata_2.h5ad')

print('both adatas done')

