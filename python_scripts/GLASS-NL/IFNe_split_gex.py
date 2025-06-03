#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# IFNe_split_gex.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Subset for samples with zero and nonzero IFNE counts and check the frequency 
# of other genes being/not being expressed
#
#-------------------------------------------------------------------------------
# 0 Import packages and prepare variables
#-------------------------------------------------------------------------------
# Import packages
import pandas as pd
import matplotlib.pyplot as plt

# Load the data (assuming your CSV has gene symbols as row names and sample names as columns)
bulk = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_RNA-seq/Readcounts/GLASS_only_expression_gSymbol_primary.csv', index_col=1)
bulk = bulk.drop(bulk.columns[0], axis=1)
# List of genes of interest (excluding IFNE)
goi = ['DMRTA1', 'CDKN2B', 'CDKN2A', 'MTAP', 'IFNA2', 'IFNB1', 'IFNG', 'HACD4', 'TNF']
# Extract the IFNE expression values (assuming IFNE is one of the rows)

#-------------------------------------------------------------------------------
# 1 Calculate the frequency of >0 and ==0 readcounts for other GOI when daa subset based on IFNE
#-------------------------------------------------------------------------------
ifne_expression = bulk.loc['IFNE'].values  # Get IFNE expression values across all samples

# Split the samples into two groups based on IFNE expression (zero or greater than zero)
group_0 = bulk.loc[:, ifne_expression == 0]  # Samples where IFNE = 0
group_nonzero = bulk.loc[:, ifne_expression > 0]  # Samples where IFNE > 0

# Initialize lists to store the counts for 0 and >0 expression for each gene in each group
count_values_0_group_0 = []  # Counts for expression = 0 in group 0
count_values_nonzero_group_0 = []  # Counts for expression > 0 in group 0
count_values_0_group_nonzero = []  # Counts for expression = 0 in group nonzero
count_values_nonzero_group_nonzero = []  # Counts for expression > 0 in group nonzero
genes = []  # To store the names of the genes

# Loop through each gene of interest (except IFNE)
for gene in goi:
    # Group 0: Samples with IFNE = 0
    gene_rows_group_0 = group_0.loc[gene]
    num_columns_gene_group_0 = (gene_rows_group_0 > 0).sum()  # Expression > 0
    num_columns_nogene_group_0 = (gene_rows_group_0 == 0).sum()  # Expression = 0
    
    # Group nonzero: Samples with IFNE > 0
    gene_rows_group_nonzero = group_nonzero.loc[gene]
    num_columns_gene_nonzero = (gene_rows_group_nonzero > 0).sum()  # Expression > 0
    num_columns_nogene_group_nonzero = (gene_rows_group_nonzero == 0).sum()  # Expression = 0
    
    # Store the results for the current gene in both groups
    count_values_nonzero_group_0.append(num_columns_gene_group_0)
    count_values_0_group_0.append(num_columns_nogene_group_0)
    count_values_nonzero_group_nonzero.append(num_columns_gene_nonzero)
    count_values_0_group_nonzero.append(num_columns_nogene_group_nonzero)
    genes.append(gene)

# Create a bar plot for all genes in Group 0 (IFNE = 0)
x = range(len(goi))  # Position of the genes on the x-axis
width = 0.35  # Width of the bars


#-------------------------------------------------------------------------------
# 2 Plot
#-------------------------------------------------------------------------------
plt.figure(figsize=(12, 6))

# Plot bars for expression = 0 and > 0 for each gene in Group 0
plt.bar(x, count_values_0_group_0, width, label='Expression = 0 (IFNE = 0)', color='salmon')
plt.bar([p + width for p in x], count_values_nonzero_group_0, width, label='Expression > 0 (IFNE = 0)', color='skyblue')

# Add labels, title, and legend
plt.xlabel('Genes', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.title('Expression Frequency for Genes of Interest (IFNE = 0)', fontsize=16)
plt.xticks([p + width / 2 for p in x], genes, ha='right', fontsize=12)  # Rotate x-axis labels
plt.legend(fontsize=12)

# Save the plot as a PNG file for Group 0
plt.tight_layout()  # Adjust layout to prevent overlap
plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/gene_frequency_in_bulk/genes_expression_frequency_IFNE_0.png')

# Show the plot
plt.show()

# Create a bar plot for all genes in Group Nonzero (IFNE > 0)
plt.figure(figsize=(12, 6))

# Plot bars for expression = 0 and > 0 for each gene in Group Nonzero
plt.bar(x, count_values_0_group_nonzero, width, label='Expression = 0 (IFNE > 0)', color='salmon')
plt.bar([p + width for p in x], count_values_nonzero_group_nonzero, width, label='Expression > 0 (IFNE > 0)', color='skyblue')

# Add labels, title, and legend
plt.xlabel('Genes', fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.title('Expression Frequency for Genes of Interest (IFNE > 0)', fontsize=16)
plt.xticks([p + width / 2 for p in x], genes, ha='right', fontsize=12)  # Rotate x-axis labels
plt.legend(fontsize=12)

# Save the plot as a PNG file for Group Nonzero
plt.tight_layout()  # Adjust layout to prevent overlap
plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/gene_frequency_in_bulk/genes_expression_frequency_IFNE_nonzero.png')

# Show the plot
plt.show()
