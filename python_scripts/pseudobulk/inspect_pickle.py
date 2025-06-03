#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# inspect_pickle.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Inspect pickle files through debugging function 
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#-------------------------------------------------------------------------------
# 0 Import packages
#-------------------------------------------------------------------------------
import pickle
import os
import sys
import pandas as pd
# import anndata 
# import scanpy as sc 
# import numpy as np
# import glob
# import torch


os.chdir('/net/beegfs/cfg/tgac/dmartinovicova/Statescope/test/')
os.chdir('../')
sys.path.insert(1, 'Statescope/')
sys.path.insert(1, 'Statescope/StateDiscovery/lib/pymf/')
import Statescope
from Statescope import Initialize_Statescope

os.chdir('/net/beegfs/cfg/tgac/dmartinovicova_new/')
sys.path.insert(1, 'Statescope/Statescope/')
sys.path.insert(1, 'Statescope/Statescope/StateDiscovery/lib/pymf/')
sys.path.insert(1, 'Statescope/')
sys.path.insert(1, 'OncoBLADE/scripts/')
from BLADE import Framework_Iterative
print('Running...')

#-------------------------------------------------------------------------------
# 1 Read the .pkl files of interest
#-------------------------------------------------------------------------------
with open("/net/beegfs/cfg/tgac/dmartinovicova/GLASS-NL/Statescope/Statescope_model_initialized_GLASSNL2_fitsDominika.pkl", 'rb') as file:
    Statescope = pickle.load(file)

with open("/net/beegfs/cfg/tgac/dmartinovicova/GLASS-NL/OncoBLADE/Deconvolution/GLASSNL_oncoBLADE_output.pkl", 'rb') as file:
    onco = pickle.load(file)

print('Hello')
