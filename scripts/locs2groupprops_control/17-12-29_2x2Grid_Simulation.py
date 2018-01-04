# Template for creating iput variables
# for locs_groupprops.locs2groupprops(path,ac_lastframe,ignore_dark)
# 
#   path: List of full paths to grouped locs file ('*locs_picked.hdf5')
#   ac_lastframe: List of last frames up to which autocorreleation function is usdeed for fit
#   ignore_dark: List of how many pseudo dark frames are ignored equal to Picasso 
#%%
# Load modules
import os #platform independent paths
import numpy as np
# Load&Reload own modules
import locs_groupprops as l2grp
import importlib
importlib.reload(l2grp)

#%%
# Define folder of locs_picked.hdf5 file
dir_name='E:/Projects/qPAINT/data/18-01-02_OwnSimulation_1ds/0ignore_lin'
# Define names of locs_picked.hdf5 file
file_names=['01_locs_picked.hdf5']
file_names.extend(['02_locs_picked.hdf5'])
file_names.extend(['03_locs_picked.hdf5'])
file_names.extend(['04_locs_picked.hdf5'])
file_names.extend(['05_locs_picked.hdf5'])
file_names.extend(['06_locs_picked.hdf5'])
file_names.extend(['07_locs_picked.hdf5'])
file_names.extend(['08_locs_picked.hdf5'])
# Create full path list
path=[os.path.join(dir_name,file_name) for file_name in file_names]

## Define up to which frame ac is fitted
ac_lastframe=[ 5e3]*4+[6e2]*4

## Define ignore_dark equal to Picasso for files in path
ignore_dark=[0]*8

#%%
# Get group properties of all groups in locs by using module l2grp
for i in range(0,np.size(path)):
    groupprops=l2grp.locs2groupprops(path[i],ac_lastframe[i],ignore_dark[i])