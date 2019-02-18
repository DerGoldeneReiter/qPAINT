#################################################### Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os #platform independent paths
import importlib
# Load user defined functions
import pickprops_calls as props_call
importlib.reload(props_call)

############################################################## Define number of groups that are combined into single one
N=48
q=3
compress=True
############################################################## Define data paths
dir_names=[]
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/19-02-04_c-series_SDS_Pm2_NoPEG_POC/SDS-Pm2-8nt-NoPEG_c5nM_p35uW_1/19-02-04_JS']*6)
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/19-02-04_c-series_SDS_Pm2_NoPEG_POC/SDS-Pm2-8nt-NoPEG_c10nM_p35uW_1/19-02-04_JS']*6)
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/19-02-04_c-series_SDS_Pm2_NoPEG_POC/SDS-Pm2-8nt-NoPEG_c20nM_p35uW_1/19-02-04_JS']*6)

# Define names of locs_picked.hdf5 file
file_names=[]
file_names.extend(['SDS-Pm2-8nt-NoPEG_c5nM_p35uW_1_MMStack_Pos0.ome_locs_picked.hdf5'])
file_names.extend(['SDS-Pm2-8nt-NoPEG_c5nM_p35uW_1_MMStack_Pos0.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-NoPEG_c5nM_p35uW_2_MMStack_Pos1.ome_locs_picked.hdf5'])
file_names.extend(['SDS-Pm2-8nt-NoPEG_c5nM_p35uW_2_MMStack_Pos1.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-NoPEG_c5nM_p35uW_2_MMStack_Pos2.ome_locs_picked.hdf5'])
file_names.extend(['SDS-Pm2-8nt-NoPEG_c5nM_p35uW_2_MMStack_Pos2.ome_locs_picked_props_ig1.hdf5'])

file_names.extend(['SDS-Pm2-8nt-NoPEG_c10nM_p35uW_1_MMStack_Pos0.ome_locs_picked.hdf5'])
file_names.extend(['SDS-Pm2-8nt-NoPEG_c10nM_p35uW_1_MMStack_Pos0.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-NoPEG_c10nM_p35uW_1_MMStack_Pos1.ome_locs_picked.hdf5'])
file_names.extend(['SDS-Pm2-8nt-NoPEG_c10nM_p35uW_1_MMStack_Pos1.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-NoPEG_c10nM_p35uW_1_MMStack_Pos2.ome_locs_picked.hdf5'])
file_names.extend(['SDS-Pm2-8nt-NoPEG_c10nM_p35uW_1_MMStack_Pos2.ome_locs_picked_props_ig1.hdf5'])

file_names.extend(['SDS-Pm2-8nt-NoPEG_c20nM_p35uW_1_MMStack_Pos0.ome_locs_picked.hdf5'])
file_names.extend(['SDS-Pm2-8nt-NoPEG_c20nM_p35uW_1_MMStack_Pos0.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-NoPEG_c20nM_p35uW_1_MMStack_Pos1.ome_locs_picked.hdf5'])
file_names.extend(['SDS-Pm2-8nt-NoPEG_c20nM_p35uW_1_MMStack_Pos1.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-NoPEG_c20nM_p35uW_1_MMStack_Pos2.ome_locs_picked.hdf5'])
file_names.extend(['SDS-Pm2-8nt-NoPEG_c20nM_p35uW_1_MMStack_Pos2.ome_locs_picked_props_ig1.hdf5'])

############################################################## Define
#### Create list of paths
path_locs=[os.path.join(dir_names[i],file_names[i]) for i in np.arange(0,len(file_names),2)]
path_props=[os.path.join(dir_names[i],file_names[i]) for i in np.arange(1,len(file_names),2)]

for i in range(0,len(path_locs)):
    locs_combine,groups_map=props_call.combine_picks(path_locs[i],path_props[i],N,q,compress)


