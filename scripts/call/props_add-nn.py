# Import modules
import os
import importlib

#### Define path to file
dir_names=[]
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/19-01-22_c-series_pm2_1-12DS/DS1+12-Pm2-8nt_c20nM_p35uW_1/19-01-22_JS/FS_picked'])

file_names=[]
file_names.extend(['all-Pm2-8nt_c20nM_p35uW_Pos0_locs_picked_props_ig1.hdf5'])

#### Create list of paths
path=[os.path.join(dir_names[i],file_names[i]) for i in range(0,len(file_names))]

#%%
#### Segment
#import pickprops_calls as props_call
import pickprops_calls as props_call
# Reload modules
importlib.reload(props_call)

for p in path:
    print(p)
    locs_props_nn=props_call.props_add_nn(p)
