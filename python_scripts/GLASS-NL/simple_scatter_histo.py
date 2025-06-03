import pandas as pd
import matplotlib.pyplot as plt

# Load CNV file, metadata
file='/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_DNA-seq/CDKN2B_IFNe_CNV.txt'
metafile='/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_RNA-seq/metadata/primary_samples_metadata.csv'

data=pd.read_csv(file, sep='\t', index_col=0, header=0, skiprows=[1,2,3])
meta=pd.read_csv(metafile, index_col=0, header=0)
print(data)
print(meta)

# Subset for primary samples
primary_samples=meta['WES_ID'].tolist()
print(primary_samples)

primary_data = data.loc[data.index.isin(primary_samples)]
print(primary_data)

# Plot
plt.scatter(primary_data[primary_data.columns[0]], primary_data[primary_data.columns[1]])  # Using the first and second column as x and y
plt.xlabel(primary_data.columns[0])  # Label for the x-axis
plt.ylabel(primary_data.columns[1])  # Label for the y-axis
plt.title('Scatter Plot')  # Title for the plot
plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/DNA_seq/scatter_CDKN2B_IFNe.png')

first_column_freq = primary_data[primary_data.columns[0]].value_counts()  # Frequency of values in the first column
second_column_freq = primary_data[primary_data.columns[1]].value_counts()  # Frequency of values in the second column

# Step 2: Plot the histograms of the value frequencies for both columns
plt.figure(figsize=(12, 6))

# Plot histogram for the first column
plt.subplot(1, 2, 1)
first_column_freq.plot(kind='bar', color='skyblue')
plt.title(f'Frequency of Values in {primary_data.columns[0]}')
plt.xlabel('Value')
plt.ylabel('Frequency')

# Plot histogram for the second column
plt.subplot(1, 2, 2)
second_column_freq.plot(kind='bar', color='coral')
plt.title(f'Frequency of Values in {primary_data.columns[1]}')
plt.xlabel('Value')
plt.ylabel('Frequency')

plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/DNA_seq/histo_CNV_CDKN2B_IFNe.png')