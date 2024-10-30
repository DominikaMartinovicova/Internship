#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# QC_scRNA.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# QC for scRNAseq counts of each dataset separately before Normalization
# *Adapted for cervical cancer data*
#
# Author: Shiva NAjjary (s.najjary@amsterdamumc.nl)
#
# Usage:
"""
        python3 scripts/QC_scRNA.py \
        -i {input.matrix} \
        -o {output.adata_processed} \
        -min_genes {params.min_genes} \
        -max_genes {params.max_genes} \
        -min_cells {params.min_cells} \
        -max_mt_pct {params.max_mt_pct} \
        -cache_dir {params.cache_dir} \
        -plot_dir {params.plot_dir}
"""
#expand(output_dir + "/merged_data_august/{dataset}_QC/{dataset}_QC.h5ad", dataset=datasets)
# TODO:
# 1) 
#
# History:
#  6-8-2024: File creation, write code
#-------------------------------------------------------------------------------
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import os
import numpy as np
import scrublet as scr

#-------------------------------------------------------------------------------
# 1.1 Parse command line arguments
#-------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog='python3 QC_scRNA.py',
                                     formatter_class=argparse.RawTextHelpFormatter, description='QC scRNAseq counts using scanpy')
    parser.add_argument('-i', help='Input matrix', dest='input', type=str, required=True)
    parser.add_argument('-o', help='Path to output file', dest='output', type=str, required=True)
    parser.add_argument('-min_genes', help='Minimal genes per cell', dest='min_genes', type=int, default=0)
    parser.add_argument('-max_genes', help='Maximal genes per cell', dest='max_genes', type=int, default=7000)
    parser.add_argument('-max_mt_pct', help='Maximal percentage of mitochondrial genes', dest='max_mt_pct', type=int, default=10)
    parser.add_argument('-min_cells', help='Minimal cells', dest='min_cells', type=int, default=30)
    parser.add_argument('-cache_dir', help='Directory to save cache', dest='cache_dir', type=str, required=True)
    parser.add_argument('-plot_dir', help='Directory to save plots', dest='plot_dir', type=str, required=True)
    args = parser.parse_args()
    return args

args = parse_args()

#-------------------------------------------------------------------------------
# 1.2 Configure scanpy and prepare working directory
#-------------------------------------------------------------------------------
# set directory to save cache
sc.settings.verbosity = 3
sc.settings.cachedir = args.cache_dir

# create plot directories
if not os.path.exists(args.plot_dir):
    os.makedirs(args.plot_dir)
if not os.path.exists(os.path.join(args.plot_dir, "before_filtering")):
    os.makedirs(os.path.join(args.plot_dir, "before_filtering"))
if not os.path.exists(os.path.join(args.plot_dir, "after_filtering")):
    os.makedirs(os.path.join(args.plot_dir, "after_filtering"))

#-------------------------------------------------------------------------------------------------------------------------------------------
# 2.1 read matrix data to adata object
#-------------------------------------------------------------------------------------------------------------------------------------------
# Read the 10x matrix data
print("Reading scRNAseq counts from directory")
adata = sc.read_h5ad(args.input)

# Remove duplicate gene names
adata.var_names_make_unique()

print(adata)

# Make a copy for initial filtering
adata_copy = adata.copy()

#---------------------------------------------------------------------------------------------------------------------------------------------
# 2.2 Plot top 50 expressed genes and QC metrics
#---------------------------------------------------------------------------------------------------------------------------------------------
print("Plotting and saving top 50 expressed genes before normalization in " + args.plot_dir + "Top50_expressed_genes_prenorm.png")
sc.pl.highest_expr_genes(adata, n_top=50,show = False)
plt.savefig(args.plot_dir + "Top50_expressed_genes_prenorm.png")
#--------------------------------------------------------------------------------------------------------------------------------------------
# Calculate/Plot QC metrics
#--------------------------------------------------------------------------------------------------------------------------------------------
# annotate mitochrondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
print("Plotting and saving QC metrics before filtering in "+ args.plot_dir + "before_filtering/")
#------------------------------------------------------------------------------------------------------------------------------------------
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, show = False)
plt.savefig(args.plot_dir + "before_filtering/" + "QC_metrics.png")

#--------------------------------------------------------------------------------------------------------------------------------------
# Plot percentage of mitochondrial counts vs total counts
print("Plotting percentage mitochondrial counts vs total counts")
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', show = False)
plt.savefig(args.plot_dir + "before_filtering/" + "total_counts_vs_pct_mt.png")
#-------------------------------------------------------------------------------------------------------------------------------------
print("Plotting and saving n_genes_by_counts vs total_counts")
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show = False)
plt.savefig(args.plot_dir + "before_filtering/" + "total_counts_vs_ngenes.png")
#-------------------------------------------------------------------------------------------------------------------------------------------
# 2.2 Filtering Steps
#-------------------------------------------------------------------------------------------------------------------------------------
# 2.3 Remove cells by min/max genes per cell
#-------------------------------------------------------------------------------------------------------------------------------------------
# Filter cells by min/max genes per cell
print(f"Filtering cells: min genes = {args.min_genes}, max genes = {args.max_genes}")
sc.pp.filter_cells(adata, min_genes=args.min_genes)
sc.pp.filter_cells(adata, max_genes=args.max_genes)
print(adata.n_obs, "cells remaining")

# Filter genes by min number of cells
print(f"Filtering genes: min cells = {args.min_cells}")
sc.pp.filter_genes(adata, min_cells=args.min_cells)
print(adata.n_obs, "cells remaining")


#-------------------------------------------------------------------------------------------------------------------------------------------
# 2.4 Recalculate QC metrics after filtering
#------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate mt_coutns after filtering
print("Recalculating QC metrics after filtering")
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


#------------------------------------------------------------------------------------------
print("Checking raw data before setting:")
print(adata.raw)

# Store raw counts
print("Storing raw counts...")
adata.raw = adata 
print(adata.raw)


#------------------------------------------------------------------------------------------------------------------------------------------------
# 2.5 Filter cells by mitochondrial gene expression and plot QC metrics
#------------------------------------------------------------------------------------------------------------------------------------------------
 # annotate the group of mitochondrial genes as 'mt'
#adata = adata[adata.obs.n_genes_by_counts < 7000, :]
print("Retain cells with a maximum precentage of mitochondrial genes of",str(args.max_mt_pct) + "%")
adata = adata[adata.obs.pct_counts_mt < args.max_mt_pct, :]
print(adata.n_obs, "cells remaining")

#---------------------------------------------------------------------------------------------------------------------------------------------------
# 2.3 Plot QC Metrics After Filtering
#--------------------------------------------------------------------------------------------------------------------------------------------------
print("Plotting and saving QC metrics after filtering in "+  args.plot_dir + "after_filtering/")
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, show = False)
plt.savefig(args.plot_dir + "after_filtering/" + "QC_metrics.png")

# Plot percentage of mitochondrial counts vs total counts after filtering
print("Plotting percentage mitochondrial counts vs total counts after filtering")
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', show = False)
plt.savefig(args.plot_dir + "after_filtering/" + "total_counts_vs_pct_mt.png")

# Plot number of genes vs total counts after filtering
print("Plotting and saving n_genes_by_counts vs total_counts")
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show = False)
plt.savefig(args.plot_dir + "after_filtering/" + "total_counts_vs_ngenes.png")

#--------------------------------------------------------------------------------------------------------------------------------------------------
# 2.4 Extra checking, not necessary for all
#--------------------------------------------------------------------------------------------------------------------------------------------------
print("Filtered data summary:")
print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")


def check_raw_data(adata):
    print("Checking raw data:")
    
    # Compare shapes
    print(f"Raw data matrix shape: {adata.raw.X.shape}")
    print(f"Data matrix shape: {adata.X.shape}")

    # Check if raw data and current data are the same
    if (adata.raw.X != adata.X).nnz == 0:
        print("The raw data and current data are identical.")
    else:
        print("The raw data and current data differ.")
    
    # Check first few values
    print("First few rows of raw data:")
    print(adata.raw.X[:5, :5].toarray() if hasattr(adata.raw.X, 'toarray') else adata.raw.X[:5, :5])
    
    print("First few rows of current data:")
    print(adata.X[:5, :5].toarray() if hasattr(adata.X, 'toarray') else adata.X[:5, :5])
    
    # Mean and variance
    print("Mean of raw data:", np.mean(adata.raw.X.toarray() if hasattr(adata.raw.X, 'toarray') else adata.raw.X))
    print("Mean of current data:", np.mean(adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X))
    
    print("Variance of raw data:", np.var(adata.raw.X.toarray() if hasattr(adata.raw.X, 'toarray') else adata.raw.X))
    print("Variance of current data:", np.var(adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X))

# Run this function after loading and processing your data
check_raw_data(adata)

#-------------------------------------------------------------------------------
# 2.5 Write filtered_adata.h5ad
#-------------------------------------------------------------------------------
print("Writing to file")
adata.write(args.output)





