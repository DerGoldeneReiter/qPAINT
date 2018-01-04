# Template for creating simulation locs
# for simulate_locs.simulate_locs(path,tau_b,tau_d,NoFrames,NoDocks,NoGroups)
# 
#   path: List of full paths to grouped locs files to be created ('*locs_picked.hdf5')
#   ac_lastframe: List of last frames up to which autocorreleation function is usdeed for fit
#   ignore_dark: List of how many pseudo dark frames are ignored equal to Picasso 
#%%
# Load modules
import os #platform independent paths
import numpy as np
import importlib

# Load&Reload own modules
import locs_groupprops as l2grp
importlib.reload(l2grp)
import simulate_locs
importlib.reload(simulate_locs)

#%%
# Define folder of locs_picked.hdf5 file
dir_name=['E:/Projects/qPAINT/data/18-01-03_OwnSimulation_1ds_ConcentrationScreen/20frames']*4+['E:/Projects/qPAINT/data/18-01-03_OwnSimulation_1ds_ConcentrationScreen/100frames']*4
# Define names of locs_picked.hdf5 file
file_names=['01']
file_names.extend(['02'])
file_names.extend(['03'])
file_names.extend(['04'])
file_names.extend(['01'])
file_names.extend(['02'])
file_names.extend(['03'])
file_names.extend(['04'])

# Create full path list
path=[]
for i in range(0,np.size(file_names)): 
    path.append(os.path.join(dir_name[i],file_names[i]))
    
## Define simulation parameters
tau_b=[20]*4+[100]*4 # Bright time [Frames]
tau_d=[10000,2000,1000,500]
tau_d=tau_d+tau_d # Dark time [Frames]
NoFrames=[100000]*8
NoDocks=[1]*8
NoGroups=[400]*8

#%%
# Simulation for all paths
for i in range(0,np.size(path)):
    locs=simulate_locs.simulate_locs(path[i],tau_b[i],tau_d[i],NoFrames[i],NoDocks[i],NoGroups[i])