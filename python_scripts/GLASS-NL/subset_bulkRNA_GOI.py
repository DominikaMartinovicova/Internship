#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# subset_bulkRNA_IFN.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Subset bulkRNAseq data to only readcounts from IFN genes. 
#
#-------------------------------------------------------------------------------
# 0 Import packages and prepare variables
#-------------------------------------------------------------------------------
import pandas as pd

# Load the data (assuming your CSV has gene symbols as row names and sample names as columns)
bulk = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_RNA-seq/Readcounts/GLASS_only_expression_gSymbol_primary.csv', index_col=1)
bulk = bulk.drop(bulk.columns[0],axis=1)

# List of genes for which we want to extract read counts
genes_of_interest = ['IFNE', 'IFNA2', 'IFNB1', 'IFNG']

#-------------------------------------------------------------------------------
# 1 Subset and save 
#-------------------------------------------------------------------------------
# Extract the rows for these genes
subset_data = bulk.loc[genes_of_interest]
print(subset_data)

# Save the table to a new CSV file
subset_data.to_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_RNA-seq/Readcounts/IFN_readcounts.tsv', sep='\t')


