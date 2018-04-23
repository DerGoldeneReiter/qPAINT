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
dir_names=['/fs/pool/pool-schwille-paint/Data/Simulation/18-04-19_copasi_biexp/']

# Define names of locs_picked.hdf5 file
file_names=['taub1s&taub5s_kon2e6_n12_c10nM_locs_picked.hdf5']
#file_names.extend(['taub2s_kon2e6_n24_c10nM_locs_picked.hdf5'])
#file_names.extend(['taub2s_kon2e6_n24_c20nM_locs_picked.hdf5'])
#file_names.extend(['taub2s_kon2e6_n24_c50nM_locs_picked.hdf5'])
#file_names.extend(['taub2s_kon2e6_n24_c100nM_locs_picked.hdf5'])

# Create full path list
path=[]
for i in range(0, len(file_names)):
    path.append(os.path.join(dir_names[i],file_names[i]))

# Possible locs2groupprops input arguments
ignore_dark=[0]*5
NoDocks=[12]*5
# Load p_inf_1 calibration array
#p_inf_1=np.load(r'E:\Flo\repos\qPAINT\scripts\analysis\18-01-04_OwnSimulation_Screening\p_inf_1.npy')


#%%
# Get group properties of all groups in locs by using module l2grp

for i in range(0,np.size(path)):
    groupprops=l2grp.locs2groupprops(path[i],ignore_dark[i],NoDocks=NoDocks[i])    

#for i in range(0,np.size(path)):
#    groupprops=l2grp.locs2groupprops(path[i],ignore_dark[i],p_inf_1=p_inf_1[i,:])       
