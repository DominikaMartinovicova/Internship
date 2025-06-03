#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Bland-Altman.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Bland-Altman (a.k.a Tukey mean diffrerence plot) to compare the differences in 
# expression between bulk and single cell data and identify genes whose expression 
# differs too much.
#
#-------------------------------------------------------------------------------
# 0 Import packages and read the data
#-------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

# Read data
print('Reading data...')
scRNA_pseudobulk_data = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova/Statescope/test/data/aryamaan_original/Psedudobulk_dominika_all.csv', index_col = 0)
bulkRNA_data = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_RNA-seq/Readcounts/GLASS_expression_filtered.csv', index_col = 0)

#-------------------------------------------------------------------------------
# 1 Preprocess the data
#-------------------------------------------------------------------------------
# Normalize per cell (columns) to counts per 10,000
norm_scRNA = scRNA_pseudobulk_data.divide(scRNA_pseudobulk_data.sum(axis=0), axis=1) * 1e4
scRNA_counts_df = np.log1p(norm_scRNA)
print(scRNA_counts_df)

# Sum duplicates
bulkRNA_df = bulkRNA_data.groupby(bulkRNA_data.index).sum()

# Log normalize
norm_bulkRNA = bulkRNA_df.divide(bulkRNA_df.sum(axis=0), axis=1) * 1e4
bulk_counts_df = np.log1p(norm_bulkRNA)
print(bulk_counts_df)

# Find genes in common
genes = np.intersect1d(scRNA_counts_df.index, bulk_counts_df.index)
print(len(genes))
scRNA = scRNA_counts_df.loc[genes].values
bulk = bulk_counts_df.loc[genes].values
print(scRNA.shape, bulk.shape)     # The length of genes should be the same

#-------------------------------------------------------------------------------
# 2 Calculate the differences and identify genes that are too different 
#-------------------------------------------------------------------------------
# Calculate means and differences
mean_scRNA = np.mean(scRNA, axis=1)  
mean_bulk = np.mean(bulk, axis=1)
means = np.mean([mean_scRNA, mean_bulk], axis=0)
print(len(means))

differences = mean_scRNA - mean_bulk
print(len(differences)) # Has to be the same as length means

# DataFrame for easier handling
data = pd.DataFrame({'gene': genes, 'mean': means, 'difference': differences})
print(data)

# Calculate the mean of the difference and the standard deviation
mean_diff = np.mean(data['difference'])
sd_diff = np.std(data['difference'])

# Determine cutoffs for significant differences
upper_cutoff = mean_diff + 1.96 * sd_diff
lower_cutoff = mean_diff - 1.96 * sd_diff

# Identify significantly different genes
data['significant'] = (data['difference'] > upper_cutoff) | (data['difference'] < lower_cutoff)
sig_dif = data[data['significant']]['gene']
print(sig_dif)

# Save to .csv
sig_dif.to_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_DNA-seq/equivalent_sc_bulk_genes.csv')

#-------------------------------------------------------------------------------
# 3 Plot Bland-Altman plot 
#-------------------------------------------------------------------------------
# Plotting the Bland-Altman plot
plt.figure(figsize=(10, 6))
plt.scatter(data['mean'], data['difference'], color=np.where(data['significant'], 'red', 'grey'), alpha=0.5)
plt.axhline(y=mean_diff, color='blue', linestyle='dashed', label='Mean difference')
plt.axhline(y=upper_cutoff, color='red', linestyle='dashed', label='Upper cutoff')
plt.axhline(y=lower_cutoff, color='red', linestyle='dashed', label='Lower cutoff')
plt.title('Bland-Altman Plot')
plt.xlabel('Mean Log Expression (scRNA and Bulk)')
plt.ylabel('Difference in Log Expression (scRNA - Bulk)')
plt.legend()
plt.show()

plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/Bland_Altman.png')



# # Check how many different genes are considered markers in signature
# signature = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova/OncoBLADE/glioma_signature_all.txt', sep='\t')
# print(signature)

# markers = signature[signature['IsMarker']==True]['Gene']
# print(markers)

# overlapping = np.intersect1d(markers, sig_dif)
# print(overlapping)
