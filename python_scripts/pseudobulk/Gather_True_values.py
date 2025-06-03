#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Gather_True_values.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Gather true fractions and profiles for evaluation
#
# Author: Jurriaan Janssen (j.janssen4@amsterdamumc.nl)
#
# Usage:
"""
python3 scripts/Gather_True_values.py \
        -i {input.adata_final} \
        -fold {wildcards.fold} \
        -celltypes {wildcards.celltype_level} \
        -o_profiles {output.True_profiles} \
        -o_fractions {output.True_fractions} \
	-o_matrices {output.True_matrices}
"""
#
# TODO:
# 1) Code is messy and repetitive
#
# History:
#  30-11-2021: File creation, write code
#  08-12-2021: Edits for adata_final
#  12-05-2022: Edits for final_celltype
#  01-08-2022: Edits for CV
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
import argparse
import scanpy as sc
import pickle
import numpy as np

#-------------------------------------------------------------------------------
# 1.1 Parse command line arguments
#-------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 Gather_True_values.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Gather true fractions and profiles for evaluation  ')
    parser.add_argument('-i', help='path to input files',
                        dest='input',
                        type=str)
    parser.add_argument('-i_matrix', help='path to input files',
                        dest='input_matrix',
                        type=str)
    parser.add_argument('-fold', help='Cross validation fold to create signature for',
                        dest='fold',
                        type=int,
                        default=1000)
    parser.add_argument('-seed', help='Seed for random shuffling of patients',
                        dest='seed',
                        type=int,
                        default = 1)
    parser.add_argument('-celltypes', help='Celltypes to create signature for',
                        dest='celltype',
                        choices = ["major_celltype","final_celltype"],
                        type=str)
    parser.add_argument('-o_profiles', help='path to profiles output file',
                        dest='output_profiles',
                        type=str)
    parser.add_argument('-o_fractions', help='path to fractions output file',
                        dest='output_fractions',
                        type=str)
    parser.add_argument('-o_matrices', help='path to matrices output file',           
                        dest='output_matrices',                                    
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

#-------------------------------------------------------------------------------
# 2.1 read Anndata object
#-------------------------------------------------------------------------------
print("Reading adata...")
adata = sc.read_h5ad("/net/beegfs/cfg/tgac/dmartinovicova/scRNA_snake_complete/data/combined_all/adata_raw_w_labels_l012_70cases_new.h5ad")

#-------------------------------------------------------------------------------
# 3.1 Obtain per patient fractions and profiles
#-------------------------------------------------------------------------------
# fetch major and minor celltypes
print("Creating dictionaries...")
celltypes = adata.obs['l1_celltype'].unique().tolist()

# intialize dict to store fractions
fractions_dict = dict()

# intialize dict to store profiles
profiles_dict = dict()

# initialize dict to store reconstructed matrices
matrix_dict = dict()
i=0
# iterate over patients
for patient_id in adata.obs['case_id'].unique():
    i += 1
    print(f'Processing patient: {i} {patient_id}')
    Profile_celltype = dict()
    Count_celltype = dict()
    Matrix_celltype = dict()
    # iterate over celltypes
    for celltype in celltypes:
        celltype_data = adata[adata.obs['l1_celltype'] == celltype]
        # obtain indices of celltype in patient
        patient_indices = [index for index, element in enumerate(celltype_data.obs['case_id']) if element == patient_id]
        # if no patient indices (i.e. a cell type not found for this patient)
        if not patient_indices:
            Ncells = 0
            Profile = np.empty(adata.raw.shape[1])
            Profile[:] = np.nan
        else:
            filtered_data = celltype_data[patient_indices]
            # get Ncells
            Ncells = filtered_data.shape[0]
            # obtain profile
            Profile = filtered_data.raw.X.toarray().mean(axis = 0)
        # store in count dict
        Count_celltype[celltype] = Ncells
        # store profile
        Profile_celltype[celltype] = Profile
    # Store profile in eventual dict
    profiles_dict[patient_id] =  Profile_celltype


    matrix_dict[patient_id] =  Matrix_celltype
    profiles_dict['celltype_list'] = celltypes
    #print(profiles_dict)

    # calculate total number of cells in patient
    Ntot_cells = sum(Count_celltype.values())
    # change counts to fractions
    for celltype,count in Count_celltype.items():
        Count_celltype[celltype] = count / Ntot_cells
        Matrix_celltype[celltype] = Profile_celltype[celltype] * Count_celltype[celltype]
    # store results
    fractions_dict[patient_id] = Count_celltype
    fractions_dict = {key:fractions_dict[key] for key in sorted(fractions_dict)}
    #print(fractions_dict)

    matrix_dict[patient_id] = Matrix_celltype
    matrix_dict = {key:matrix_dict[key] for key in sorted(matrix_dict)}
    #print(matrix_dict)

#-------------------------------------------------------------------------------
# 4.1 Write to file
#-------------------------------------------------------------------------------
print('Writing files...')
with open('/net/beegfs/cfg/tgac/dmartinovicova/pseudobulk/pseudobulk_profiles_l1.pkl', "wb") as out:
    pickle.dump(profiles_dict, out)
with open('/net/beegfs/cfg/tgac/dmartinovicova/pseudobulk/pseudobulk_fractions_l1.pkl', "wb") as out:
    pickle.dump(fractions_dict, out)
with open('/net/beegfs/cfg/tgac/dmartinovicova/pseudobulk/pseudobulk_matrices_l1.pkl', "wb") as out:
    pickle.dump(matrix_dict, out)