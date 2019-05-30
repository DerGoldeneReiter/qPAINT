# Import modules
import os
import importlib


#### Define path to file
dir_names=[]
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/z.simulations/19-01-22_copasi_Pm2-8nt_error-meas-time/N12/36k/2-5nM']*3)
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/z.simulations/19-01-22_copasi_Pm2-8nt_error-meas-time/N12/36k/5nM']*3)
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/z.simulations/19-01-22_copasi_Pm2-8nt_error-meas-time/N12/36k/10nM']*3)
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/z.simulations/19-01-22_copasi_Pm2-8nt_error-meas-time/N12/36k/20nM']*3)
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/z.simulations/19-01-22_copasi_Pm2-8nt_error-meas-time/N12/36k/30nM']*3)

file_names=[]
file_names.extend(['N12_2-5nM_1_locs_picked.hdf5'])
file_names.extend(['N12_2-5nM_2_locs_picked.hdf5'])
file_names.extend(['N12_2-5nM_3_locs_picked.hdf5'])
file_names.extend(['N12_5nM_1_locs_picked.hdf5'])
file_names.extend(['N12_5nM_2_locs_picked.hdf5'])
file_names.extend(['N12_5nM_3_locs_picked.hdf5'])
file_names.extend(['N12_10nM_1_locs_picked.hdf5'])
file_names.extend(['N12_10nM_2_locs_picked.hdf5'])
file_names.extend(['N12_10nM_3_locs_picked.hdf5'])
file_names.extend(['N12_20nM_1_locs_picked.hdf5'])
file_names.extend(['N12_20nM_2_locs_picked.hdf5'])
file_names.extend(['N12_20nM_3_locs_picked.hdf5'])
file_names.extend(['N12_30nM_1_locs_picked.hdf5'])
file_names.extend(['N12_30nM_2_locs_picked.hdf5'])
file_names.extend(['N12_30nM_3_locs_picked.hdf5'])

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
    props_call.segment_time(p,noFrames_seg)