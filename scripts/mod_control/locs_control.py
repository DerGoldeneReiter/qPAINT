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
dir_names=[r'E:\Projects\qPAINT\data\18-01-12\sample01_p40_convampf_3x_T25_1\JS_18-01-16\ng=400_single']

# Define names of locs_picked.hdf5 file
file_names=['sample01_p40_convampf_3x_T25_1_MMStack_Pos0.ome_locs_render_picked.hdf5']
#file_names.extend(['02_locs_picked.hdf5'])
#file_names.extend(['03_locs_picked.hdf5'])
#file_names.extend(['04_locs_picked.hdf5'])
#file_names.extend(['05_locs_picked.hdf5'])
#file_names.extend(['06_locs_picked.hdf5'])
#file_names.extend(['07_locs_picked.hdf5'])
#file_names.extend(['08_locs_picked.hdf5'])


# Create full path list
path=[]
for i in range(0, len(file_names)):
    path.append(os.path.join(dir_names[i],file_names[i]))

# Possible locs2groupprops input arguments
ignore_dark=[1]*8
NoDocks=[1]*8
# Load p_inf_1 calibration array
p_inf_1=np.load(r'E:\Flo\repos\qPAINT\scripts\analysis\18-01-04_OwnSimulation_Screening\p_inf_1.npy')


#%%
# Get group properties of all groups in locs by using module l2grp

for i in range(0,np.size(path)):
    groupprops=l2grp.locs2groupprops(path[i],ignore_dark[i],NoDocks=NoDocks[i])    

#for i in range(0,np.size(path)):
#    groupprops=l2grp.locs2groupprops(path[i],ignore_dark[i],p_inf_1=p_inf_1[i,:])       
