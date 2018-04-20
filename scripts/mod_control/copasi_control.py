
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
dir_names=['/fs/fs01/lv03/home/b_schwille/stehr/programs/copasi/reactions']
# Define names of locs_picked.hdf5 file
file_names=['taub1s&taub5s_kon2e6_n60&12_c10nM.txt']
#file_names.extend(['taub2s_kon2e6_n24_c10nM.txt'])
#file_names.extend(['taub2s_kon2e6_n24_c20nM.txt'])
#file_names.extend(['taub2s_kon2e6_n24_c50nM.txt'])
#file_names.extend(['taub2s_kon2e6_n24_c100nM.txt'])

# Create full path list
path=[]
for i in range(0, len(file_names)):
    path.append(os.path.join(dir_names[i],file_names[i]))


# Set number of intervals as defined in COPASI
intervals=100000
# Set interval_size as defined in COPASI
interval_size=0.05

for i in range(0,np.size(path)):
#    locs=simulate_locs.copasi2locs(path[i],interval_size,intervals)
    locs=simulate_locs.copasi2locs_double(path[i],interval_size,intervals)
