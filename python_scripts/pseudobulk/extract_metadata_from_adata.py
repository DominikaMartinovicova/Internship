#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# extract_metadata_from_adata.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Extract metadata from adata file
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#-------------------------------------------------------------------------------
# 0 Import packages and prepare variables
#-------------------------------------------------------------------------------
import scanpy as sc
import pandas as pd
import argparse

#-------------------------------------------------------------------------------
# 1 Parse command line arguments
#-------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 Metadata from adata',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Extracty metadata from adata file  ')
    parser.add_argument('-adata', help='path to adata file',
                        dest='adata',
                        type=str)
    parser.add_argument('-o', help='path to output file',
                        dest='output',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

#-------------------------------------------------------------------------------
# 2 Extract metadata info
#-------------------------------------------------------------------------------
adata = sc.read_h5ad(args.adata)  # Read adata
metadata = adata.obs[['case_id', 'Cohort', 'Diagnosis', 'Diagnosis_label', 'Primary_Recurrent', 'IDH']]   # Extract relevant metadata columns

# Remove duplicate values and set 'case_id' as the index
metadata = metadata.drop_duplicates('case_id')
metadata.set_index('case_id', inplace=True)

print(metadata)  # Display the extracted metadata

#-------------------------------------------------------------------------------
# 3 Save into a .csv
#-------------------------------------------------------------------------------
metadata.to_csv(args.output)