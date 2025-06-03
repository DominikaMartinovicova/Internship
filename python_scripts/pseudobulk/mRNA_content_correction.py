#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# mRNA_content_correction.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Correcting for diverse mRNA content in cell types. 
# Calculate average mRNA content per cell in a cell type and divide the deconvolved 
# cell type fraction by the average mRNA content.
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#-------------------------------------------------------------------------------
# 0 Import packages and prepare variables
#-------------------------------------------------------------------------------
import pandas as pd
import scanpy as sc
import numpy as np
import argparse

#-------------------------------------------------------------------------------
# 1 Parse command line arguments
#-------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 Transcriptome size correction',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Correct for varying transcritpome size  ')
    parser.add_argument('-adata', help='path to adata file',
                        dest='adata',
                        type=str)
    parser.add_argument('-f', help='path to deconvolved fractions file',
                        dest='fractions',
                        type=str)
    parser.add_argument('-o', help='path to corrected output file',
                        dest='output_corrected',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

#-------------------------------------------------------------------------------
# 2 Reading data
#-------------------------------------------------------------------------------
print("Reading data...")
adata_file = args.adata

fractions_df = pd.read_csv(args.fractions, index_col=0)
fractions_df = fractions_df[sorted(fractions_df.columns)]  # Sort columns alphabetically
print(fractions_df)

adata = sc.read_h5ad(adata_file)
print(adata)

#-------------------------------------------------------------------------------
# 3 Correction
#-------------------------------------------------------------------------------
# Convert expression matrix to a dense DataFrame (cells x genes)
if isinstance(adata.X, np.ndarray):
    expr_df = pd.DataFrame(adata.X, index=adata.obs.index)
else:
    expr_df = pd.DataFrame(adata.X.toarray(), index=adata.obs.index)  # Ensure dense array
print(expr_df.shape)

# Sum all gene counts for each cell
adata.obs['TotalReadCount'] = expr_df.sum(axis=1)

# Group by Cell Type (l2_celltype) and compute sum and count
grouped = adata.obs.groupby("l2_celltype").agg(
    TotalReadCount=("TotalReadCount", "sum"), 
    NumCells=("TotalReadCount", "count")
).reset_index()

grouped.loc[grouped["NumCells"] == 0, "NumCells"] = np.nan  # Avoid division by zero by replacing 0 with NaN before division
grouped["AvgReadCountPerCell"] = grouped["TotalReadCount"] / grouped["NumCells"]    # Compute average read count per cell per cell type

# Convert to dictionary for easy lookup
avg_counts_dict = grouped.set_index("l2_celltype")["AvgReadCountPerCell"].to_dict()

# Adjust fractions by dividing by the global avg read count per cell type
print('Adjusting fractions...')
adjusted_fractions = fractions_df.copy()

for cell_type in adjusted_fractions.columns:
    if cell_type in avg_counts_dict:
    # Divide fraction by mRNA content for each cell type
        adjusted_fractions[cell_type] = fractions_df[cell_type] / avg_counts_dict[cell_type]


# Scale the adjusted fractions so they sum to 1 for each sample
normalization_factor = adjusted_fractions.sum(axis=1)
normalized_adjusted_fractions = adjusted_fractions.div(normalization_factor, axis=0)

print(normalized_adjusted_fractions.head())

#-------------------------------------------------------------------------------
# 4 Save into a .csv
#-------------------------------------------------------------------------------
normalized_adjusted_fractions.to_csv(args.output_corrected)
