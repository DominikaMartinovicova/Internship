#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# filter_GLASS_samples.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Filter unsuitable GLASS samples which have normal CNA profile and do not have 
# matching RNA and DNA seq data. 
#
#-------------------------------------------------------------------------------
# 0 Import packages and prepare variables
#-------------------------------------------------------------------------------
# Import packages
import pandas as pd

# Read in the data
print("Reading data...")
metadata = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_RNA-seq/metadata/Master_Datasheet_ALL_METHODS_27012023.csv', usecols=['GS_ID', 'WES_ID'], index_col='GS_ID')
RNA_readcounts = pd.read_csv("/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_RNA-seq/Readcounts/GLASS_only_expression_gSymbol.tsv", sep="\t", index_col='GeneSymbol')
tumorF_data = pd.read_csv("/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_DNA-seq/GLASS_tumorF_for_Statescope.csv", index_col='case_id')
excluded_samples = pd.read_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_RNA-seq/metadata/excluded.csv', header=None)
print(excluded_samples)

print(len(tumorF_data.index))

#-------------------------------------------------------------------------------
# 1 Filter samples
#-------------------------------------------------------------------------------
metadata = metadata.reset_index().dropna().set_index('GS_ID')   # Remove rows with NaN for any of the columns. Excludes samples without matching RNA and DNA sequencing data.
RNA_readcounts = RNA_readcounts.drop(columns=['Geneid']).sort_index(axis=0).T   # Remove Geneid column with ENSEMBL ID

# Check for which samples there are both DNA and RNA sequencing
overlapping_samples = metadata.index.intersection(RNA_readcounts.index)
print(f"Overlapping samples: {len(overlapping_samples)}")

# Subset for only those samples in both dataframes -> Retain only those samples that have both DNA and RNA seq data available
RNA_readcounts = RNA_readcounts.loc[overlapping_samples]
RNA_readcounts.index = RNA_readcounts.index.map(metadata['WES_ID'])     # reindex to have the WES-ID as row index

# Check for samples that have both RNAseq and tumor fraction data available + remove samples that are in the excluded samples list (samples with normal CNA profile)
overlapping_samples2 = tumorF_data.index.intersection(RNA_readcounts.index).difference(excluded_samples[0])
print(overlapping_samples2)

RNA_readcounts = RNA_readcounts.loc[overlapping_samples2]       # Subset only for the intersected samples

RNA_readcounts = RNA_readcounts.T
RNA_readcounts.index.name = None
print(RNA_readcounts)

tumorF_data = tumorF_data.loc[overlapping_samples2]     # Subset only for the intersected samples
print(tumorF_data)

#-------------------------------------------------------------------------------
# 2 Save to .csv new bulk dataframe and new tumor fraction fro deconvolution
#-------------------------------------------------------------------------------
RNA_readcounts.to_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_RNA-seq/Readcounts/GLASS_expression_filtered.csv')
tumorF_data.to_csv('/net/beegfs/cfg/tgac/dmartinovicova_new/GLASS-NL_DNA-seq/GLASS_tumorF_for_Statescope_filtered.csv')