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
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/z.simulations/19-01-22_copasi_Pm2-8nt_error-meas-time/N4/36k']*3)
# Define names of locs_picked.hdf5 file
file_names=[]
file_names.extend(['N4_30nM_1.txt'])
file_names.extend(['N4_30nM_2.txt'])
file_names.extend(['N4_30nM_3.txt'])

# Create full path list
path=[]
for i in range(0, len(file_names)):
    path.append(os.path.join(dir_names[i],file_names[i]))


# Set number of intervals as defined in COPASI
intervals=36000
# Set interval_size as defined in COPASI
interval_size=0.2

for i in range(0,np.size(path)):
    locs=copasi_convert.copasi2locs(path[i],interval_size,intervals)
#    locs=copasi_convert.copasi2locs_double(path[i],interval_size,intervals)
