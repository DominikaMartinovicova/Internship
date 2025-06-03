#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# gene_counts.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Count the number of transcripts of a specific gene. 
#
#-------------------------------------------------------------------------------
# 0 Import packages and prepare variables
#-------------------------------------------------------------------------------
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the RNA readcounts from CSV
print('Reading data...')
df = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova/GLASS-NL/GLASS-NL_RNA-seq/Readcounts/GLASS_only_expression_gSymbol_primary.csv', index_col=1)
df=df.drop(df.columns[0], axis=1)
df = df.T  # Transpose to make the gene names as column names
print(df)

# List of other genes you want to correlate with IFNE
goi = ['DMRTA1', 'CDKN2B', 'CDKN2A', 'MTAP', 'IFNA2', 'IFNB1', 'IFNG', 'IFNE', 'HACD4', 'TNF']
print(goi)

#-------------------------------------------------------------------------------
# 1 Calculate correlation between each two genes
#-------------------------------------------------------------------------------
# Extract the genes of interest from the dataframe
df_subset = df[goi]

# Calculate the correlation matrix
correlation_matrix = df_subset.corr()

# Visualize the correlation matrix using a heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt='.2f', vmin=-1, vmax=1)
plt.title("Correlation between the GOI")
plt.show()

# Save
plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/correlation_gex.png')

#-------------------------------------------------------------------------------
# 2 Plot a scatterplot of IFNE and each GOI
#-------------------------------------------------------------------------------
for gene in goi:
    if gene == 'IFNE':
        continue
    else:
        plt.figure(figsize=(6, 4))
        sns.scatterplot(x=df['IFNE'], y=df[gene])
        plt.xlabel('IFNE Expression')
        plt.ylabel(f'{gene} Expression')
        plt.title(f'Scatter plot of IFNE vs {gene}')
        plt.savefig(f'/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/IFNe_correlation_gex/{gene}_scatter_corr.png')