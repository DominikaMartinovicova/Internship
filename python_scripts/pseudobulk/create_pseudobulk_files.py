#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# create_pseudobulk_files.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# 1. Create pseudobulk data from adata object by accumulating the read counts for each gene per sample as DataFrame. Needed for Statescope.
# 2. Save fractions of all cell types per sample and subset for Malignant fraction and replace other with NaN. Needed as Expectation for deconvolution with Statescope.
# 1. Create pseudobulk data from adata object by accumulating the read counts for each gene per sample as .pkl. Needed for oncoBLADE.
# 3. Create an array with NaN values for all cell fractions but Malignant. Needed as Expectation for deconvolution with oncoBLADE.
# 4. Optionally, create minipseudobulk and corresponding minifiles for testing.
#
#-------------------------------------------------------------------------------
# 0 Import packages
#-------------------------------------------------------------------------------
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import pickle

#-------------------------------------------------------------------------------
# 1 Create pseudobulk as DataFrame
#-------------------------------------------------------------------------------
def create_pseudobulk(adata, patient_id, pseudobulk_file):
    print('Creating pseudobulk...')
    # Group by 'case_id' and sum the gene expression counts
    print('Saving adata.X into expression_matrix variable...')
    expression_matrix = adata.X

    # Check if matrix is sparse (if not and the file is too large the script might crash)
    print('Checking if sparse matrix...')
    if isinstance(expression_matrix, csr_matrix):
        # Create a DataFrame where rows are 'case_id' and columns are gene names, sum each 
        print('Creating data frame...')
        summed_expression = pd.DataFrame.sparse.from_spmatrix(expression_matrix, index=adata.obs[patient_id], columns=adata.var_names).groupby(level=0).sum()
        summed_expression.index.name = None
    else:
        # If it's not sparse (should not happen if adata.X is large), use regular dense handling
        print("Expression matrix is not sparse.")

    # Print resulting DataFrame
    print(summed_expression.head())
    print(summed_expression.index)

    # Save
    print('Saving to csv...')
    summed_expression.to_csv(pseudobulk_file)
    return summed_expression

#-------------------------------------------------------------------------------
# 2 Calculate True Cell Fractions and DataFrame Expectation
#-------------------------------------------------------------------------------
def calculate_fractions_expectation_csv(adata, patient_id, cell_type_column, cell_fractions_file, malignant_fraction_file):
    # Count the number of cells per cell type in each sample
    print("Calculating cell type fractions...")
    cell_type_counts = adata.obs.groupby([patient_id, cell_type_column]).size().unstack(fill_value=0)

    total_cells_per_sample = adata.obs.groupby(patient_id).size()   # Calculate the total number of cells per sample
    cell_type_fractions = cell_type_counts.div(total_cells_per_sample, axis=0)  # Compute the fraction of each cell type in each sample

    # Reorder the columns based on order of cell types in signature
    signature_order = ['Malignant','Macrophage','Oligodendrocyte','Monocyte','CD4_Tcell','Endothelial','CD8_Tcell',
                    'Microglia','Fibroblast','cDC','NK_cell','pDC','Pericyte','B_cell','T_reg','Plasma_B']
    cell_type_fractions = cell_type_fractions[signature_order]
    print(cell_type_fractions.head())   # Print the first few rows of cell type fractions

    # Save the dataframe with all the cell type fractions
    print('Saving cell type fractions to CSV...')
    cell_type_fractions.to_csv(cell_fractions_file)

    # Keep values only for Malignant celltype, replace all other values with nan
    print("Calculating Expectation...")
    Expectation=cell_type_fractions
    Expectation.loc[:, Expectation.columns != 'Malignant'] = np.nan
    print(Expectation)

    Expectation = Expectation.map(lambda x: 0.01 if x == 0 else (0.99 if x == 1 else x), na_action='ignore')    # Replace 0 with 0.1 and 1 with 0.99 (necessary for deconvolution)

    # Save the dataframe with only malignant fraction and others nan to .csv
    print('Saving Expectations to CSV...')
    Expectation.to_csv(malignant_fraction_file)

    return cell_type_fractions, Expectation

# -------------------------------------------------------------------------------
# 3 Create pseudobulk as dictionary (adapted from Juriaan's script simulate_bulk.py)
# -------------------------------------------------------------------------------
def create_pseudobulk_dict(adata, pseudobulk_output):
    # obtain raw gene expression
    raw_data = adata.X.toarray()
    # obtain patient IDs
    patient_ids = adata.obs['case_id'].tolist()
    # initialize dict
    simulated_bulk = dict()
    # iterate over patients
    for patient_id in sorted(set(patient_ids), key=patient_ids.index):
        print(patient_id)
        # obtain indices of 
        patient_indices = [index for index, element in enumerate(patient_ids) if element == patient_id]
        # subset raw data
        filtered_data = raw_data[patient_indices]
        # add sum of counts to dict
        summed_counts = np.sum(filtered_data, axis = 0)
        # store in dict
        simulated_bulk[patient_id] = summed_counts

    sorted_bulk = {key:simulated_bulk[key] for key in sorted(simulated_bulk)}
    output_dict = {'pseudobulk' : sorted_bulk, 'genelist': adata.var_names}

    with open(pseudobulk_output, "wb") as out:
        pickle.dump(output_dict, out)

# -------------------------------------------------------------------------------
# 4 Create an array with only malignant cell fractions
# -------------------------------------------------------------------------------
def create_Expectation_array(cell_type_fractions, malignant_fraction_file_pkl):
    # Intialize Expectation (Nsample x Ncell with NaN for non-tumor celltypes)
    print('Creating array with true malignant fractions...')
    Expectation = np.zeros((len(cell_type_fractions.index), len(cell_type_fractions.columns))) + np.nan
    print(Expectation.shape)
    # iterate over samples
    for i in range(len(cell_type_fractions.index)):
        # iterate over celltypes
        for j,celltype in enumerate(cell_type_fractions.columns):
            if celltype in ['Cancer_cell', 'Tumor cell', 'Malignant']:
                # fetch true tumor purity and add to array
                Expectation[i,j] = cell_type_fractions.iloc[i,j]
            else:
                pass

    print(Expectation)

    # Check if Expectation contains 0 or 1 values which cause issues with computation in deconvolution
    if (Expectation == 0).any():
        print('The Expectation contains 0 which is not allowed. Consider giving a very small value (0.01)')
        Expectation[Expectation == 0] = 0.01  # Replace 0 with a small value like 0.01
        print("Zero values replaced with 0.01 to avoid issues.")
    if (Expectation == 1).any():
        print('The Expectation contains 1 which is not allowed. Consider giving a very large value (0.99)')
        Expectation[Expectation == 1] = 0.99  # Replace 1 with a large value like 0.99
        print("1 values replaced with 0.99 to avoid issues.")

    # Save the numpy array as .pkl
    print('Saving Expectations to .pkl...')
    with open(malignant_fraction_file_pkl, 'wb') as file:
        pickle.dump(Expectation, file)

#-------------------------------------------------------------------------------
# 5 Create minifiles - pseudobulk, cell fractions, Expectation
#-------------------------------------------------------------------------------
def create_mini(pseudobulk, ct_fractions, Expectation):
    print('Creating mini files for testing deconvolution...')
    # Create minipseudobulk
    mini_pseudo = pseudobulk.iloc[0:20,]
    print(mini_pseudo)
    print('Saving to csv...')
    mini_pseudo.to_csv("/net/beegfs/cfg/tgac/dmartinovicova/Statescope/test/data/midipseudobulk_dominika.csv")

    # Create mini celltype fractions
    mini_cell_type_fraction = ct_fractions.iloc[0:20,]
    print(mini_cell_type_fraction)
    mini_cell_type_fraction.to_csv("/net/beegfs/cfg/tgac/dmartinovicova/Statescope/test/data/midi_true_fractions_dominika.csv")

    # Create mini Expectation
    mini_expectation = Expectation[:20]
    print(mini_expectation)
    mini_expectation.to_csv("/net/beegfs/cfg/tgac/dmartinovicova/Statescope/test/data/midi_malignant_fractions_dominika.csv")

#-------------------------------------------------------------------------------
# 6 Call functions
#-------------------------------------------------------------------------------
def main():
    # Specify variables
    patient_id = 'case_id' # Name of the column with patient ids in adata.obs
    cell_type_column = 'l2_celltype' # Name of the column with cell types in adata.obs
    adata_file = "/net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_complete/data/combined_all/adata_raw_w_labels_l012_synapse.h5ad"
    
    # output files
    pseudobulk_file = "/net/beegfs/cfg/tgac/dmartinovicova/Statescope/test/data/synapse_pseudobulk_dominika0.csv"
    dict_pseudobulk_file = ""
    cell_fractions_file = "/net/beegfs/cfg/tgac/dmartinovicova/Statescope/test/data/synapse_cell_fractions_dominika0.csv"
    malignant_fraction_file = "/net/beegfs/cfg/tgac/dmartinovicova/Statescope/test/data/synapse_malignant_fraction_dominika0.csv"
    malignant_fraction_file_pkl = "/net/beegfs/cfg/tgac/dmartinovicova/Statescope/test/data/synapse_malignant_fraction_dominika0.pkl"


    # Read adata
    print('Reading adata...')
    adata = sc.read_h5ad(adata_file)
    print(adata)

    # Call functions of choice
    pseudobulk = create_pseudobulk(adata, patient_id, pseudobulk_file)
    ct_fractions, Expectation = calculate_fractions_expectation_csv(adata, patient_id, cell_type_column, cell_fractions_file, malignant_fraction_file)
    create_pseudobulk_dict(adata, dict_pseudobulk_file)
    #create_Expectation_array(ct_fractions, malignant_fraction_file_pkl)
    #create_mini(pseudobulk, ct_fractions, Expectation)

main()
