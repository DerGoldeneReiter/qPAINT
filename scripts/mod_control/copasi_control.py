
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
dir_names=['/fs/pool/pool-schwille-paint/Data/Simulation/18-04-16_copasi/']*5
# Define names of locs_picked.hdf5 file
file_names=['taub2s_kon2e6_n1_c5nM.txt']
file_names.extend(['taub2s_kon2e6_n1_c10nM.txt'])
file_names.extend(['taub2s_kon2e6_n1_c20nM.txt'])
file_names.extend(['taub2s_kon2e6_n1_c50nM.txt'])
file_names.extend(['taub2s_kon2e6_n1_c100nM.txt'])

# Create full path list
path=[]
for i in range(0, len(file_names)):
    path.append(os.path.join(dir_names[i],file_names[i]))


# Set number of intervals as defined in COPASI
intervals=15000
# Set interval_size as defined in COPASI
interval_size=0.1

for i in range(0,np.size(path)):
    locs=simulate_locs.copasi2locs(path[i],interval_size,intervals)


