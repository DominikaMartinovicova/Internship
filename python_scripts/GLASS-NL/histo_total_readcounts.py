#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# histo_total_readcount.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Calculate the total readcount of bulkRNAseq per sample
#
#-------------------------------------------------------------------------------
# 0 Import packages and prepare variables
#-------------------------------------------------------------------------------
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_RNA-seq/Readcounts/GLASS_only_expression.tsv', sep='\t')

#-------------------------------------------------------------------------------
# 1 Sum and plot
#-------------------------------------------------------------------------------
sum_data = data.sum().iloc[1:]      # Sum the readcounts

# Plot
plt.figure(figsize=(10, 6))
plt.hist(sum_data, color='skyblue')

# Add labels and title
plt.xlabel('Total Counts')
plt.ylabel('Frequency')
plt.title('Frequency of total read counts')

# Adjust layout for better spacing
plt.tight_layout()

# Save
plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/total_counts.png')

