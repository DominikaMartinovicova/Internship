import anndata as ad 
import numpy as np 
import matplotlib.pyplot as plt 

  

# Step 1: Load the AnnData object from the .h5ad file 
file_path = "/trinity/home/dmartinovicova/all_data/Verhaak_synapse/output/Verhaak_synapse_2.h5ad"
adata = ad.read_h5ad(file_path) 

# Step 2: Calculate the column sums (summing across features for each cell) 
# If adata.X is sparse, convert it to dense before summing 
if not isinstance(adata.X, np.ndarray): 
    column_sums = np.array(adata.X.sum(axis=0)).flatten() 
else: 
    column_sums = adata.X.sum(axis=0) 


# Print with two decimal points
column_sums_rounded = np.round(column_sums, 2)
print(column_sums_rounded)

# Convert the column sums to integers
column_sums_int = column_sums.astype(int)
print(column_sums_int)

print(column_sums)

## Step 3: Create the histogram 
#plt.figure(figsize=(8, 6)) 
#plt.hist(column_sums, bins=1000, color='skyblue', edgecolor='black') 
#plt.title("Distribution of Column Sums (Total Expression per Cell)") 
#plt.xlim(0,1000000)
#plt.xlabel("Total Expression per Cell") 
#plt.ylabel("Frequency") 
#plt.grid(True) 
#
#  
#
## Step 4: Show the plot 
#plt.show() 

 