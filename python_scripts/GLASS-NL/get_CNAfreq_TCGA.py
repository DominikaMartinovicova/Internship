#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# get_CNAfreq_TCGA.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Plot the frequency of IFNE deletion in TCGA data.
#
#-------------------------------------------------------------------------------
# 0 Import packages and prepare variables
#-------------------------------------------------------------------------------
# Import packages
import pandas as pd
import matplotlib.pyplot as plt

data = "/net/beegfs/cfg/tgac/dmartinovicova_new/TCGA/LGG_IFN_CNA.txt"

# Read the data
df = pd.read_csv(data, sep="\t", index_col=0)
print(df)

#-------------------------------------------------------------------------------
# 1 Counting samples with deletion
#-------------------------------------------------------------------------------
counts = df['V1'].value_counts()
print("Value counts in the 'V1' column:")
print(counts)

# Create a histogram of the 'value' column
plt.figure(figsize=(8, 6))
plt.hist(df['V1'], bins=range(-2, 3), edgecolor='black', align='left')  # Use bins for -2, -1, 0
plt.xticks(range(-2, 1))  # Set x-ticks to match the values
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram of Values')

# Save 
plt.savefig("/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/TCGA/IFNe_del_frequency.png")
