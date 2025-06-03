import pandas as pd

# Open the file, skip the first line, and read it as a TSV
df = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_RNA-seq/Readcounts/GLASS.LGG.EMC.RNA.readcounts.deduplicated_s_2.txt', sep='\t', skiprows=1)

# Split columns starting from the 7th (index 6) on '/' and update column names with the 9th element (index 8)
new_columns = df.columns.tolist()

for i in range(6, len(df.columns)):
#    value = str(new_columns.iloc[i])  # Convert the value to string 
    new_column_name = new_columns[i].split('/')[9]  # Get the 9th element (index 8) after splitting by '/'
    new_columns[i] = new_column_name

# Set the new column names
df.columns = new_columns

# Print the first 5 rows and first 5 columns
print(df.iloc[:5, :10])

# Create a new df with only read counts
df_new = df.iloc[:, [0] + list(range(6, len(df.columns)))]
print(df_new.head())

# Save the new df
#df_new.to_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_RNA-seq/Readcounts/GLASS_only_expression.tsv',sep='\t', index = False)

import pybiomart as pbm

dataset = pbm.Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
results = dataset.query(attributes=['ensembl_gene_id', 'hgnc_symbol'])

# Assuming you now have a mapping of Ensembl Gene ID to Gene Symbol
gene_id_to_symbol = dict(zip(results['Gene stable ID'], results['HGNC symbol']))

# Remove version from the Ensemble Gene ID to find matching Gene Symbol
df_new['Geneid'] = df_new['Geneid'].str.split('.').str[0]

# Replace Ensembl Gene ID with the corresponding Gene Symbol in your DataFrame
df_new['GeneSymbol'] = df_new['Geneid'].map(gene_id_to_symbol).fillna(df_new['Geneid'])
print(df_new.head())

# Reorder columns so that 'GeneSymbol' is the second column
column_order = [df_new.columns[0], 'GeneSymbol'] + [col for col in df_new.columns if col not in ['GeneSymbol'] and col != df_new.columns[0]]
df_new = df_new[column_order]

# Display the updated dataframe
print(df_new.head())

# Save the updated DataFrame to a new file
df_new.to_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_RNA-seq/Readcounts/GLASS_only_expression_gSymbol.tsv', sep='\t', index=False)


