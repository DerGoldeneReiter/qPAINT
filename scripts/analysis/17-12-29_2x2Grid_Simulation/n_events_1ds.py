# Load packages
import h5py as h5py #hdf5 handling
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os #platform independent paths
import importlib
import sys
from tqdm import tqdm
import matplotlib
from matplotlib.ticker import FormatStrFormatter

# Load user defined functions
import locs_groupprops as l2grp
import file_formats as fifo

# Define folder of locs.hdf5 file
dir_name=['E:/Projects/qPAINT/data/17-12-12_Simulation2x2Grid/FS_18-01-03_1ds_1ignore_lin']*8

# Define names of locs.hdf5 file
file_names=['01_locs_picked_groupprops.hdf5']
file_names.extend(['02_locs_picked_groupprops.hdf5'])
file_names.extend(['03_locs_picked_groupprops.hdf5'])
file_names.extend(['04_locs_picked_groupprops.hdf5'])
file_names.extend(['05_locs_picked_groupprops.hdf5'])
file_names.extend(['06_locs_picked_groupprops.hdf5'])
file_names.extend(['07_locs_picked_groupprops.hdf5'])
file_names.extend(['08_locs_picked_groupprops.hdf5'])

# Define labels
labels=['5nM']
labels.extend(['20nM'])
labels.extend(['10nM'])
labels.extend(['1nM'])
labels.extend(['5nM'])
labels.extend(['20nM'])
labels.extend(['10nM'])
labels.extend(['1nM'])

# Create list of paths
path=list()
for i in range(0,len(file_names)):
    path.append(os.path.join(dir_name[i],file_names[i]))
    
# Delete variables not used
del file_names
del dir_name

# Set general font size
matplotlib.rcParams.update({'font.size': 15})

#%%
# Single docking sites, tau_b set to 5s=100frames
# tau_b and tau_c for various concentrations
figure_size=np.array([11.69,8.27])*1.5
f1=plt.figure(1,figsize=figure_size)
f1.clf()

# Define axarr for plotting, shared x and y axis
axarr = f1.subplots(2,4, sharex=True, sharey=True)
# Control hspace 
f1.subplots_adjust(hspace=0.05, left=0.1,right=0.87)

# Order in which plots are shown
order100=[3,0,2,1]
order12=[7,4,6,5]

# Axes limits
#axarr[0,1].set_xlim([0,300])
# Labels
for i in range(0,4):
    axarr[1,i].set_xlabel(r'$N_{events}$')
for i in range(0,2):
    axarr[i,0].set_ylabel('Counts')
# Titles
f1.suptitle('Picasso Simulate: Number of binding events\n 1 docking site per group')
for i in range(0,4):
    axarr[0,i].set_title(labels[order100[i]],Fontsize=17)

# Texts
x_text=250
y_text=85

for i in range(0,4):
    # Mean & std of n_events
    n_events_mean_100=np.mean(fifo.read_grouprops(path[order100[i]])['n_events'])
    n_events_std_100=np.std(fifo.read_grouprops(path[order100[i]])['n_events'])
    n_events_mean_12=np.mean(fifo.read_grouprops(path[order12[i]])['n_events'])
    n_events_std_12=np.std(fifo.read_grouprops(path[order12[i]])['n_events'])
    
    # Text annotations
    axarr[0,i].text(x_text,y_text,'%d (%d)'%(n_events_mean_100,n_events_std_100))
    axarr[1,i].text(x_text,y_text,'%d (%d)'%(n_events_mean_12,n_events_std_12))

    
    # ac_tau
    axarr[0,i].hist(fifo.read_grouprops(path[order100[i]])['n_events'],bins='fd')
    # tau_b
    axarr[1,i].hist(fifo.read_grouprops(path[order12[i]])['n_events'],bins='fd',color='olive')

    
#Legend
axarr[0,3].legend(['100 Frames'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
axarr[1,3].legend(['12 Frames'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)