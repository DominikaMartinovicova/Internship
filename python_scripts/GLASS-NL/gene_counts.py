#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# gene_counts.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Count the number of transcripts of a specific gene. 
#
#-------------------------------------------------------------------------------
# 0 Import packages and prepare variables
#-------------------------------------------------------------------------------
# Import packages
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read the data
data = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_RNA-seq/Readcounts/GLASS_only_expression_gSymbol_primary.csv', index_col=1) 
data=data.drop(data.columns[0], axis=1) # ignore the ENSEMBL IDs
print(data.head())

#-------------------------------------------------------------------------------
# 1 Count the gene transcripts
#-------------------------------------------------------------------------------
# Specify genes of interest 
goi = ['IFNA2', 'IFNB1', 'IFNG', 'IFNE','CDKN2A', 'CDKN2B', 'MTAP', 'DMRTA1', 'HACD4']

# Loop over all the genes in goi list and plot histogram, bargraph and check for correlation between total readcounts and readcounts of each gene
for gene in goi:
    if gene in data.index:
        print(f'Processing {gene} gene...')
        gene_data = data.loc[gene]
        # Get the frequency of each unique value in the gene row/column
        frequency = gene_data.value_counts()
        print(frequency)

        # Increase the font size
        plt.rcParams.update({'font.size': 12})
        
        # Make a histogram
        plt.figure(figsize=(10, 6))
        plt.hist(gene_data, bins=20, color='skyblue', edgecolor='white')
        plt.xlabel(f'RNA Read Counts for {gene} Gene')
        plt.ylabel('Frequency')
        plt.title(f'Histogram of {gene} Gene Read Counts Across Samples')
        #tick_values = range(0, max(frequency.index)+100, 100)
        #plt.xticks(tick_values)
        plt.savefig(f'/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/{gene}_counts_histo.png')

        # Make a barplot
        plt.figure(figsize=(12, 6))
        plt.bar(gene_data.index, gene_data.values, color='skyblue', edgecolor='white')
        plt.xlabel('Samples')
        plt.ylabel(f'RNA Read Counts for {gene} Gene')
        plt.title(f'RNA Read Counts for {gene} Gene Across Samples')
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(f'/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/{gene}_counts_bar.png')

        # Check for correlation between total readcount and the goi readcount 
        sum_data = data.sum()
        sum_data=sum_data.reset_index()
        gene_data=gene_data.reset_index()

        data_merged = pd.merge(gene_data, sum_data, on='index')
        print(data_merged)
        data_merged.columns = ['ID',gene, 'Total']
        print(data_merged)

        correlation = data_merged[gene].corr(data_merged['Total'])
        print(f"Correlation: {correlation}")

        # Plot the data
        plt.figure(figsize=(8, 6))
        sns.scatterplot(x=data_merged['Total'], y=data_merged[gene])
        plt.title(f"Correlation = {correlation:.2f}")
        plt.ylabel(f'{gene} Readcounts')
        plt.xlabel('Total Readcounts')

        plt.savefig(f'/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/correlation_rc_{gene}_total.png')
    else:
        print(f'Gene {gene} not found in data')
