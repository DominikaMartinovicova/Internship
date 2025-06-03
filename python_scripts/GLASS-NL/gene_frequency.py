#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# gene_frequency.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Inspect the frequency of a certain gene being expressed. 
#
#-------------------------------------------------------------------------------
# 0 Import packages and prepare variables
#-------------------------------------------------------------------------------
# Import packages
import pandas as pd
import matplotlib.pyplot as plt

# Load bulk and specifiy genes of interest
bulk = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_RNA-seq/Readcounts/GLASS_only_expression_gSymbol_primary.csv')
goi = ['DMRTA1', 'CDKN2B', 'CDKN2A', 'MTAP', 'IFNA2', 'IFNB1', 'IFNG', 'IFNE', 'HACD4']

#-------------------------------------------------------------------------------
# 1 Count the frequency
#-------------------------------------------------------------------------------
# Initialize lists to store the counts for 0 and >0 expression for each gene
count_values_0 = []  # Counts for expression = 0
count_values_nonzero = []  # Counts for expression > 0
genes = []  # To store the names of the genes

for gene in goi:
    print(gene)
    # Check which columns have expression ==0 or >0
    gene_rows = bulk[bulk['GeneSymbol'] == gene]
    num_columns_gene = (gene_rows.iloc[:, 2:] > 0).sum(axis=1)   # Save indices of samples with expression
    num_columns_nogene = (gene_rows.iloc[:, 2:] == 0).sum(axis=1)    # Save indices of samples with no expression

    # Store the sums of > 0 and = 0 expressions for the current gene
    count_values_nonzero.append(num_columns_gene.sum())
    count_values_0.append(num_columns_nogene.sum())
    genes.append(gene)

#-------------------------------------------------------------------------------
# 2 Plot a barplot with the calculated frequencies
#-------------------------------------------------------------------------------
# Create a bar plot for all genes
x = range(len(goi))  # Position of the genes on the x-axis
width = 0.35  # Width of the bars
plt.figure(figsize=(10, 6))

# Plot bars for expression = 0 and > 0 for each gene
plt.bar(x, count_values_0, width, label='Expression = 0', color='salmon')
plt.bar([p + width for p in x], count_values_nonzero, width, label='Expression > 0', color='skyblue')

# Add labels, title, and legend
plt.xlabel('Genes')
plt.ylabel('Frequency')
plt.title('Expression Frequency for Genes of Interest')
plt.xticks([p + width / 2 for p in x], genes)  # Place gene names at the center of the grouped bars
plt.legend()

# Save the plot as a PNG file
plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/genes_expression_frequency.png')
