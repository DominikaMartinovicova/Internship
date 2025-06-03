#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# tumor_fraction_for_deconvolution.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Create suitable tumor fraction dataframe to provide prior expectation to Statescope deconvolution.
#
#-------------------------------------------------------------------------------
# 0 Import packages 
#-------------------------------------------------------------------------------
import pandas as pd

#-------------------------------------------------------------------------------
# 1 Reaad the data and shift the index
#-------------------------------------------------------------------------------
# Get sample names
glass_ACE = pd.read_csv("/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_DNA-seq/GLASSFULL_ACE_curatedFilled_001.txt", sep=" ", index_col=1)
print(glass_ACE)
glass_ACE = glass_ACE.drop(columns=['samplenames'])
samples = glass_ACE.index
print(samples)

# Get cell types in the signature to be deconvolved
pseudobulk = pd.read_csv("/net/beegfs/cfg/tgac/dmartinovicova/Statescope/test/data/aryamaan_original/GT_dominika_all.csv")
pseudobulk = pseudobulk.drop(columns=['case_id'])
cell_types = pseudobulk.columns
print(cell_types)

# Create a dataframe
tumor_fractions=pd.DataFrame(index=samples, columns=cell_types)
tumor_fractions['Malignant'] = glass_ACE['ManualPurity']    # Fill in the malignant column with the values from glass_ACE, leave all the other columns empty
tumor_fractions.index.name = 'case_id'
print(tumor_fractions)

# Replace 0 with 0.01 and 1 with 0.99
if (tumor_fractions == 0).any().any():
    print("Replacing 0 values with 0.01")
if (tumor_fractions == 1).any().any():
    print("Replacing 1 values with 0.99")

tumor_fractions = tumor_fractions.mask(tumor_fractions == 0, 0.01).mask(tumor_fractions == 1, 0.99)
print(tumor_fractions)

#-------------------------------------------------------------------------------
# 2 Save the tumor fraction file
#-------------------------------------------------------------------------------
tumor_fractions.to_csv("/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_DNA-seq/GLASS_tumorF_for_Statescope.csv")