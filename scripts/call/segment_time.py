# Import modules
import os
import importlib


#### Define path to file
dir_names=[]
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/z.simulations/19-01-22_copasi_Pm2-8nt_error-meas-time/N1/36k/30nM']*3)

file_names=[]
file_names.extend(['N1_30nM_1_locs_picked.hdf5'])
file_names.extend(['N1_30nM_2_locs_picked.hdf5'])
file_names.extend(['N1_30nM_3_locs_picked.hdf5'])

#### Create list of paths
path=[os.path.join(dir_names[i],file_names[i]) for i in range(0,len(file_names))]

#%%
#### Segment
import pickprops_calls as props_call
# Reload modules
importlib.reload(props_call)

noFrames_seg=9000
for p in path:
    print(p)
    props_call.segment_time(p,9000)