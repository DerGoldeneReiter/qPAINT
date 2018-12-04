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

# Define folder of locs_picked.hdf5 file
dir_names=[] 
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/z.PAINT-checks/18-11-02_SDS-P1-9nt/SDS9nt-80_P1-20nM_p100mW-8deg_flat_1/18-11-02_FS/'])
# Define names of locs_picked.hdf5 file
file_names=[] 
file_names.extend(['SDS9nt-80_P1-20nM_p100mW-8deg_flat_1_MMStack.ome_locs_picked_d1-0.hdf5'])
# Create full path list
path=[]
for i in range(0, len(file_names)):
    path.append(os.path.join(dir_names[i],file_names[i]))

# Ignore dark frames <=ignore_dark 
ignore_dark=[1]*5
#%%
# Get group properties of all groups in locs by using module l2grp
for i in range(0,np.size(path)):
    groupprops=l2grp.locs2groupprops(path[i],ignore_dark[i])
     
