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
dir_names=[r'E:\Projects\qPAINT\data\18-01-12\sample01_p40_convampf_3x_T32_1\JS_18-01-16\crop4sum_ng=1000']

# Define names of locs_picked.hdf5 file
file_names=['sample01_p40_convampf_3x_T32_1_MMStack_Pos0.ome_locs_sum.hdf5']

# Create full path list
path=[]
for i in range(0, len(file_names)):
    path.append(os.path.join(dir_names[i],file_names[i]))

# Expected number of docking sites
NoDocks=[12]#*16
#%%
# Get group properties of all groups in locs by using module l2grp
for i in range(0,np.size(path)):
    groupprops=l2grp.sum2groupprops(path[i],NoDocks[i])