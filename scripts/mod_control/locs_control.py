# Template for creating iput variables
# for locs_groupprops.locs2groupprops(path,ignore_dark)
# 
#   path: List of full paths to grouped locs file ('*locs_picked.hdf5')
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
dir_names=['/fs/pool/pool-schwille-paint/Data/D134/18-06-22/N12-3_10nM-P1modA_p30_T23_10MHz-g300_field1_1/18-06-22_FS/']
#dir_names.extend([''])

# Define names of locs_picked.hdf5 file
file_names=['N12-3_10nM-P1modA_p30_T23_10MHz-g300_field1_1_MMStack_Pos0.ome_locs_picked.hdf5']
#file_names.extend([''])

# Create full path list
path=[]
for i in range(0, len(file_names)):
    path.append(os.path.join(dir_names[i],file_names[i]))

# Ignore dark frames <=ignore_dark 
ignore_dark=[1]*4

# Get group properties of all groups in locs by using module l2grp
for i in range(0,np.size(path)):
    groupprops=l2grp.locs2groupprops(path[i],ignore_dark[i])    
     
