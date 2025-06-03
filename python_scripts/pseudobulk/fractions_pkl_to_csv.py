#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# fractions_pkl_to_csv.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Transform the .pkl dictionary with fraction (from Juriaan's old script) to DataFrame
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#-------------------------------------------------------------------------------
# 0 Import packages and prepare variables
#-------------------------------------------------------------------------------
import pickle
import pandas as pd
import numpy as np

#-------------------------------------------------------------------------------
# 1 Read the .pkl dictionary with fractions
#-------------------------------------------------------------------------------
with open('/net/beegfs/cfg/tgac/dmartinovicova/pseudobulk/pseudobulk_fractions_l1.pkl', 'rb') as file:
    frac = pickle.load(file)

rows = []
for sample, cell_types in frac.items():
    for cell_type, fraction in cell_types.items():
        rows.append({
            'Sample': sample,
            'Cell Type': cell_type,
            'Fraction': fraction
        })

#-------------------------------------------------------------------------------
# 2 Create a DataFrame
#-------------------------------------------------------------------------------
df = pd.DataFrame(rows)

# Pivot the DataFrame to make samples as rows and cell types as columns
df_pivoted = df.pivot(index='Sample', columns='Cell Type', values='Fraction')

# Remove the names of the index column and row
df_pivoted.columns.name = None
df_pivoted.index.name = None
print(df_pivoted)
df_pivoted.to_csv('/net/beegfs/cfg/tgac/dmartinovicova/pseudobulk/pseudobulk_fractions_l1.csv')     # Save as .csv

# Subset only the Malignant column 
df_malignant = df_pivoted[['Malignant']]
print(df_malignant)
df_malignant.to_csv('/net/beegfs/cfg/tgac/dmartinovicova/pseudobulk/pseudobulk_malignant_fractions_l1.csv')     # Save as .csv


#-------------------------------------------------------------------------------
# 3 Create Expectation as numpy array
#-------------------------------------------------------------------------------
# Intialize Expectation (Nsample x Ncell with None for non-tumor celltypes)
Expectation = np.zeros((len(df_pivoted.index), len(df_pivoted.columns))) + np.nan
print(Expectation)
# iterate over samples
for i in range(len(df_pivoted.index)):
    # iterate over celltypes
    for j,celltype in enumerate(df_pivoted.columns):
        if celltype in ['Cancer_cell', 'Tumor cell', 'Malignant']:
            # fetch true tumor purity and add to array
            Expectation[i,j] = df_pivoted.iloc[i,j]
        else:
            pass

print(Expectation)

if (Expectation == 0).any():
    print('The Expectation contains 0 which is not allowed. Consider giving a very small value (0.01)')
    Expectation[Expectation == 0] = 0.01  # Replace 0 with a small value like 0.01
    print("Zero values replaced with 0.01 to avoid issues.")
if (Expectation == 1).any():
    print('The Expectation contains 1 which is not allowed. Consider giving a very large value (0.99)')
    Expectation[Expectation == 1] = 0.99  # Replace 1 with a small value like 0.99
    print("1 values replaced with 0.99 to avoid issues.")

if (Expectation == 0).any():
    raise ValueError('The Expectation contains 0 which is not allowed. Consider giving a very small value (0.01)')
if (Expectation == 1).any():
    raise ValueError('The Expectation contains 1 which is not allowed. Consider giving a very large value (0.99)')


with open("/net/beegfs/cfg/tgac/dmartinovicova/pseudobulk/pseudobulk_fractions_malignant.pkl", 'wb') as file:
  pickle.dump(Expectation, file)