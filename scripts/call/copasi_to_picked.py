# Load modules
import os #platform independent paths
import numpy as np
import importlib
# Load custom modules
import copasi_convert 
importlib.reload(copasi_convert)
#%%
# Define folder of locs_picked.hdf5 file
dir_names=[]
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/z.simulations/19-06-19_Pm2_2B07/N3']*4)
# Define names of locs_picked.hdf5 file
file_names=[]
file_names.extend(['N3_2-5nM.txt'])
file_names.extend(['N3_5nM.txt'])
file_names.extend(['N3_10nM.txt'])
file_names.extend(['N3_20nM.txt'])

# Create full path list
path=[]
for i in range(0, len(file_names)):
    path.append(os.path.join(dir_names[i],file_names[i]))


# Set number of intervals as defined in COPASI
intervals=[18000,18000,9000,9000]
# Set interval_size as defined in COPASI
interval_size=0.2

for i in range(0,np.size(path)):
    locs=copasi_convert.copasi2locs(path[i],interval_size,intervals[i])
#    locs=copasi_convert.copasi2locs_double(path[i],interval_size,intervals)
