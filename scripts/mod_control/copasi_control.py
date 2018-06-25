
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
dir_names=['/fs/pool/pool-schwille-paint/Data/Simulation/18-06-06_copasi_bi_0-5s-5s/50k_exp100/']
# Define names of locs_picked.hdf5 file
file_names=['taub0-5s-5s_kon2-3e6_c10nM_exp100_f50k_N12-3.txt']
#file_names.extend(['taub2s_kon2e6_n24_c10nM.txt'])
#file_names.extend(['taub2s_kon2e6_n24_c20nM.txt'])
#file_names.extend(['taub2s_kon2e6_n24_c50nM.txt'])
#file_names.extend(['taub2s_kon2e6_n24_c100nM.txt'])

# Create full path list
path=[]
for i in range(0, len(file_names)):
    path.append(os.path.join(dir_names[i],file_names[i]))


# Set number of intervals as defined in COPASI
intervals=50000
# Set interval_size as defined in COPASI
interval_size=0.1

for i in range(0,np.size(path)):
#    locs=simulate_locs.copasi2locs(path[i],interval_size,intervals)
    locs=simulate_locs.copasi2locs_double(path[i],interval_size,intervals)
