#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# histo_tumorF.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Plot histogram of tumor fraction values from GLASSNL samples. 
#
#-------------------------------------------------------------------------------
# 0 Import packages and prepare variables
#-------------------------------------------------------------------------------
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

tumorF = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_DNA-seq/GLASS_tumorF_for_Statescope_filtered.csv', usecols=['Malignant'])
print(tumorF)

#-------------------------------------------------------------------------------
# 1 Plot histogram
#-------------------------------------------------------------------------------
# Plot histogram
sns.set(style="whitegrid")
plt.figure(figsize=(8, 6))

plt.hist(tumorF, bins=10, color="dodgerblue", alpha=0.6)

# Add labels and title
plt.title("Tumor fraction frequency", fontsize=15)
plt.xlabel("Tumor fraction", fontsize=12)
plt.ylabel("Frequency", fontsize=12)

# Save
plt.savefig('/net/beegfs/cfg/tgac/dmartinovicova_new/graphs/DNA_seq/tumorF_histo.png')