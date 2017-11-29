# Load packages
import h5py as h5py #hdf5 handling
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os.path as ospath #platform independent paths
import importlib
from tqdm import tqdm
import sys

# Add function and script path
sys.path.append('modules')
sys.path.append('scripts')

# Load user defined functions
import locs_groupprops as l2grp

# Define path
dir_name='E:\\Projects\\qPAINT\\data\\20171110_Max_2x2Grid_varN\\20171011_Flo_LS'
file_name='P1_Cy3b_PPT_10mW_2_5nM_4C_Origami_1_MMStack_Pos0.ome_locs_picked_singlespots.hdf5'
path=ospath.join(dir_name,file_name)

# Number of frames in measurement
NoFrames=60000

# Open hdf5 file
locs_file=h5py.File(path,'r')
# Load dataset 'locs' into np.array 'locs'
locs=locs_file['locs'][...]
locs_file.close()
# Sort locs after groups and frames
locs.sort(order=['group','frame'],axis=0)

#%%
importlib.reload(l2grp)

NoAllGroups=np.amax(locs['group'])
g_list=np.arange(0,NoAllGroups,1)
NoGroups=np.size(g_list)

max_n_events=500

groupprops=np.empty(NoGroups,dtype=[('group','i4',1),('mean_x','f4',1),
                                    ('mean_y','f4',1),('mean_frame','f4',1),
                                    ('std_frame','f4',1),('mean_photons','f4',1),('G0','f4',1),
                                    ('n_locs','i4',1),('n_events','i4',1),
                                    ('tau_d','f4',(1,max_n_events)),('tau_b','f4',(1,max_n_events))
                                    ])

groupprops['group'][:]=g_list


for i in tqdm(range(0,NoGroups)):
    g=g_list[i]
    locs_g=locs[:][locs['group']==g]
    
    groupprops['mean_x'][i]=l2grp.get_mean_x(locs_g)
    groupprops['mean_y'][i]=l2grp.get_mean_y(locs_g)
    groupprops['mean_frame'][i]=l2grp.get_mean_frame(locs_g)
    groupprops['std_frame'][i]=l2grp.get_std_frame(locs_g)
    groupprops['mean_photons'][i]=l2grp.get_mean_photons(locs_g)
    groupprops['G0'][i]=l2grp.get_G0(locs_g,NoFrames)
    groupprops['n_locs'][i]=l2grp.get_n_locs(locs_g)
    groupprops['n_events'][i]=l2grp.get_n_events(locs_g,1)
    
    [tau_d,tau_b]=l2grp.get_tau(locs_g,1,max_n_events)
    groupprops['tau_d'][i]=tau_d
    groupprops['tau_b'][i]=tau_b

#%% Save groupprops in hdf5 file
groupprops_file = h5py.File(path.replace('.hdf5','_groupprops.hdf5'), "w")
groupprops_file.create_dataset("groupprops", np.shape(groupprops), dtype='f')
groupprops_file.close()

#%% Plot unfiltered group property distributions
plt.figure(1)
plt.clf()
plt.subplot(3,1,1)
plt.hist(groupprops['std_frame'],50)
plt.subplot(3,1,2)
plt.hist(groupprops['mean_frame'],50)
plt.subplot(3,1,3)
plt.hist(groupprops['n_events'],50)

#%% Filter groups
filter_mean_x=[150.,350.]
filter_mean_y=[150.,350.]
filter_mean_frame=[26000,36000]
filter_std_frame=[14000,20000]

# Spatially filter (center)
istrue_filter_mean_x=(groupprops['mean_x']>filter_mean_x[0])&(groupprops['mean_x']<filter_mean_x[1])
istrue_filter_mean_y=(groupprops['mean_y']>filter_mean_y[0])&(groupprops['mean_y']<filter_mean_y[1])
istrue_filter_spatial=istrue_filter_mean_x&istrue_filter_mean_y
istrue_filter_spatial_inverse=~np.array(istrue_filter_spatial)

# Filter according to signal stability
istrue_filter_mean_frame=(groupprops['mean_frame']>filter_mean_frame[0])&(groupprops['mean_frame']<filter_mean_frame[1])
istrue_filter_std_frame=(groupprops['std_frame']>filter_std_frame[0])&(groupprops['std_frame']<filter_std_frame[1])
istrue_filter_stable=istrue_filter_mean_frame&istrue_filter_std_frame

#istrue_all=istrue_filter_spatial&istrue_filter_stable
istrue_all=istrue_filter_spatial_inverse&istrue_filter_stable

groupprops_filter=groupprops[:][istrue_all]

#%% Plot unfiltered group property distributions
plt.figure(2)
plt.clf()
plt.subplot(3,1,1)
plt.plot(groupprops_filter['mean_x'],groupprops_filter['mean_y'],'.')
ax=plt.gca()
ax.set_xlim([0,500])
ax.set_ylim([0,500])
plt.subplot(3,1,2)
test_var=np.divide(groupprops_filter['G0'],np.square(groupprops_filter['n_events']))
test_var=np.divide(test_var,np.square(groupprops_filter['mean_photons']))
plt.hist(test_var,50)
plt.subplot(3,1,3)
plt.hist(groupprops_filter['n_events'],50)



