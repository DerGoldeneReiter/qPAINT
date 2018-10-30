
# Load modules
import os #platform independent paths
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt #plotting

# Load&Reload own modules
import simulate_locs 
import importlib
importlib.reload(simulate_locs)

#%%
# Define folder of locs_picked.hdf5 file
dir_names=[]
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/Simulation/18-09-18_FS_copasi_bleaching/bright3-3s_kon1-5e6_c10nM/']*5)
# Define names of locs_picked.hdf5 file
file_names=[]
file_names.extend(['12.txt'])
file_names.extend(['10.txt'])
file_names.extend(['08.txt'])
file_names.extend(['06.txt'])
file_names.extend(['04.txt'])


# Create full path list
path=[]
for i in range(0, len(file_names)):
    path.append(os.path.join(dir_names[i],file_names[i]))


# Set number of intervals as defined in COPASI
intervals=5000
# Set interval_size as defined in COPASI
interval_size=0.2

for i in range(0,np.size(path)):
    locs=simulate_locs.copasi2locs(path[i],interval_size,intervals)
#    locs=simulate_locs.copasi2locs_double(path[i],interval_size,intervals)
