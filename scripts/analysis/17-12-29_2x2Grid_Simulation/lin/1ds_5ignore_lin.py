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
dir_name='E:/Projects/qPAINT/data/17-12-12_Simulation2x2Grid/FS_18-01-03_1ds_5ignore_lin'

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

# List of concentrations as floats
concs=[float(label.replace('nM',''))*1e-9 for label in labels]

# Overall simulation paramters
CycleTime=50e-3 # Frame cycle time [s]
k_on=2e6 # On-rate [(Ms)^-1]
tau_bs=np.array([5.]*4+[0.6]*4) # Bright times [s]
tau_ds=np.array([(k_on*conc)**-1 for conc in concs]) # Dark times [s]
tau_cs=np.divide(1,np.power(tau_bs,-1)+np.power(tau_ds,-1)) # Autocorr tau [s]
# Convert times to frames
tau_bs=tau_bs/CycleTime 
tau_ds=tau_ds/CycleTime
tau_cs=tau_cs/CycleTime


# Create list of paths
path=list()
path_pick=list()
for file_name in file_names:
    path.append(os.path.join(dir_name,file_name))
    path_pick.append(os.path.join(dir_name,file_name.replace('picked_groupprops','pickprops')))
    
# Delete variables not used
del file_name
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
axarr = f1.subplots(4,5, sharex=True, sharey=True)
# Control hspace 
f1.subplots_adjust(hspace=0.05)
# Axes limits
axarr[0,1].set_xlim([0,300])
# Labels
axarr[3,0].set_xlabel(r'$\tau_c$'+' [Frames]')
axarr[3,1].set_xlabel(r'$\tau_b$'+' [Frames]')
axarr[3,2].set_xlabel(r'$\tau_{b,ignore}$'+' [Frames]')
axarr[3,3].set_xlabel(r'$length$'+' [Frames]')
axarr[3,4].set_xlabel(r'$length_{cdf}$'+' [Frames]')
for i in range(0,4):
    axarr[i,0].set_ylabel('Counts')
# Titles
f1.suptitle('Picasso Simulate: Bright times\n'+r'$\tau_b$'+'= 5 s (100 frames)\n'+
            '1 docking site per pick, ignore = 5, linearized')
#axarr[0,0].set_title(r'$\tau_c$',Fontsize=17)
#axarr[0,1].set_title(r'$\tau_b$')
#axarr[0,2].set_title(r'$\tau_{b,ignore}$')
#axarr[0,3].set_title(r'$length$',Fontsize=17)
#axarr[0,4].set_title(r'$length_{cdf}$',Fontsize=17)
# Texts
x_text=150
y_text=40


# Order in which plots are shown
order=[3,0,2,1]
for i in range(0,4):
    # Mean & std of ac_tau
    ac_tau_mean=np.mean(fifo.read_grouprops(path[order[i]])['ac_tau'])
    ac_tau_std=np.std(fifo.read_grouprops(path[order[i]])['ac_tau'])
    # Mean & std of tau_b
    tau_b_mean=np.mean(fifo.read_grouprops(path[order[i]])['tau_b'])
    tau_b_std=np.std(fifo.read_grouprops(path[order[i]])['tau_b'])
    # Mean & std of tau_b_ignore
    tau_b_ignore_mean=np.mean(fifo.read_grouprops(path[order[i]])['tau_b_ignore'])
    tau_b_ignore_std=np.std(fifo.read_grouprops(path[order[i]])['tau_b_ignore'])
    # Mean & std of len_mean
    length_mean=np.mean(fifo.read_pickprops(path_pick[order[i]])['len_mean'])
    length_std=np.std(fifo.read_pickprops(path_pick[order[i]])['len_mean'])
    # Mean & std of length_cdf
    length_cdf_mean=np.mean(fifo.read_pickprops(path_pick[order[i]])['length_cdf'])
    length_cdf_std=np.std(fifo.read_pickprops(path_pick[order[i]])['length_cdf'])
    
    # Text annotations
    axarr[i,0].text(x_text,y_text,'%d (%d)'%(ac_tau_mean,ac_tau_std))
    axarr[i,0].text(x_text,y_text+10,'%d'%(tau_cs[order[i]]),color='red')
    axarr[i,1].text(x_text,y_text,'%d (%d)'%(tau_b_mean,tau_b_std))
    axarr[i,2].text(x_text,y_text,'%d (%d)'%(tau_b_ignore_mean,tau_b_ignore_std))
    axarr[i,3].text(x_text,y_text,'%d (%d)'%(length_mean,length_std))
    axarr[i,4].text(x_text,y_text,'%d (%d)'%(length_cdf_mean,length_cdf_std))
    
    # ac_tau
    axarr[i,0].hist(fifo.read_grouprops(path[order[i]])['ac_tau'],bins='fd',label=labels[order[i]])
    # tau_b
    axarr[i,1].hist(fifo.read_grouprops(path[order[i]])['tau_b'],bins='fd',label=labels[order[i]])
    # tau_b_ignore
    axarr[i,2].hist(fifo.read_grouprops(path[order[i]])['tau_b_ignore'],bins='fd',label=labels[order[i]])
    # len_mean
    axarr[i,3].hist(fifo.read_pickprops(path_pick[order[i]])['len_mean'],bins='fd',label=labels[order[i]])
    # len_mean
    axarr[i,4].hist(fifo.read_pickprops(path_pick[order[i]])['length_cdf'],bins='fd',label=labels[order[i]])
    
    #Legend
    axarr[i,4].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

#%%
# Single docking sites, k_on=2e6(Ms)^-1 i.e. diff. tau_d for diff. cs
# tau_d for various concentrations
figure_size=np.array([11.69,8.27])*1.5
f1=plt.figure(1,figsize=figure_size)
f1.clf()

# Define axarr for plotting, shared x and y axis
axarr = f1.subplots(4,4, sharex=True, sharey=True)
# Control hspace 
f1.subplots_adjust(hspace=0.1)
# Axes limits
axarr[0,1].set_xlim([10,5e4])
# Axes scale
for ax in axarr.flat:
    ax.set_xscale('log')
# Labels
axarr[3,0].set_xlabel(r'$\tau_d$'+' [Frames]')
axarr[3,1].set_xlabel(r'$\tau_{d,cdf}$'+' [Frames]')
axarr[3,2].set_xlabel(r'$dark_{mean}$'+' [Frames]')
axarr[3,3].set_xlabel(r'$dark_{cdf}$'+' [Frames]')
for i in range(0,4):
    axarr[i,0].set_ylabel('Counts')
# Titles
f1.suptitle('Picasso Simulate: Dark time for '+r'$\tau_b$'+' = 5s'+'\n'+r'$\tau_d$'+' for '+r'$k_{on}=2e6(Ms)^{-1}$'+'\n'+
            '1 docking site per pick, ignore = 5, linearized')
# Texts
x_text=15
y_text=70


# Order in which plots are shown
order=[3,0,2,1]
for i in range(0,4):
    # Mean & std of tau_b
    tau_d_mean=np.mean(fifo.read_grouprops(path[order[i]])['tau_d'])
    tau_d_std=np.std(fifo.read_grouprops(path[order[i]])['tau_d'])
    # Mean & std of tau_b
    tau_d_cdf_mean=np.mean(fifo.read_grouprops(path[order[i]])['tau_d_cdf'])
    tau_d_cdf_std=np.std(fifo.read_grouprops(path[order[i]])['tau_d_cdf'])
    # Mean & std of dark_mean
    length_mean=np.mean(fifo.read_pickprops(path_pick[order[i]])['dark_mean'])
    length_std=np.std(fifo.read_pickprops(path_pick[order[i]])['dark_mean'])
    # Mean & std of dark_cdf
    length_cdf_mean=np.mean(fifo.read_pickprops(path_pick[order[i]])['dark_cdf'])
    length_cdf_std=np.std(fifo.read_pickprops(path_pick[order[i]])['dark_cdf'])
    
    # Text annotations
    axarr[i,0].text(x_text+150,y_text,'%d (%d)'%(tau_d_mean,tau_d_std))
    axarr[i,0].text(x_text,y_text,'%d, '%(tau_ds[order[i]]),color='red')
    axarr[i,1].text(x_text,y_text,'%d (%d)'%(tau_d_cdf_mean,tau_d_cdf_std))
    axarr[i,2].text(x_text,y_text,'%d (%d)'%(length_mean,length_std))
    axarr[i,3].text(x_text,y_text,'%d (%d)'%(length_cdf_mean,length_cdf_std))
    
    # tau_d
    axarr[i,0].hist(fifo.read_grouprops(path[order[i]])['tau_d'],bins='fd',label=labels[order[i]])
    # tau_d_cdf
    axarr[i,1].hist(fifo.read_grouprops(path[order[i]])['tau_d_cdf'],bins='fd',label=labels[order[i]])
    # dark_mean
    axarr[i,2].hist(fifo.read_pickprops(path_pick[order[i]])['dark_mean'],bins='fd',label=labels[order[i]])
    # dark_cdf
    axarr[i,3].hist(fifo.read_pickprops(path_pick[order[i]])['dark_cdf'],bins='fd',label=labels[order[i]])
    
     #Legend
    axarr[i,3].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

#%%
# Single docking sites, tau_b set to 600ms=12frames
# tau_b and tau_c for various concentrations
figure_size=np.array([11.69,8.27])*1.5
f1=plt.figure(1,figsize=figure_size)
f1.clf()

# Define axarr for plotting, shared x and y axis
axarr = f1.subplots(4,5, sharex=True, sharey=True)
# Control hspace 
f1.subplots_adjust(hspace=0.05)
# Axes limits
axarr[0,1].set_xlim([0,50])
# Labels
axarr[3,0].set_xlabel(r'$\tau_c$'+' [Frames]')
axarr[3,1].set_xlabel(r'$\tau_b$'+' [Frames]')
axarr[3,2].set_xlabel(r'$\tau_{b,ignore}$'+' [Frames]')
axarr[3,3].set_xlabel(r'$length$'+' [Frames]')
axarr[3,4].set_xlabel(r'$length_{cdf}$'+' [Frames]')
for i in range(0,4):
    axarr[i,0].set_ylabel('Counts')
# Titles
f1.suptitle('Picasso Simulate: Bright times\n'+r'$\tau_b$'+'= 0.6 s (12 frames)\n'+
            '1 docking site per pick, ignore = 5, linearized')
#axarr[0,0].set_title(r'$\tau_c$',Fontsize=17)
#axarr[0,1].set_title(r'$\tau_b$')
#axarr[0,2].set_title(r'$\tau_{b,ignore}$')
#axarr[0,3].set_title(r'$length$',Fontsize=17)
#axarr[0,4].set_title(r'$length_{cdf}$',Fontsize=17)
# Texts
x_text=25
y_text=30
# Order in which plots are shown
order=[7,4,6,5]
for i in range(0,4):
    # Mean & std of ac_tau
    ac_tau_mean=np.mean(fifo.read_grouprops(path[order[i]])['ac_tau'])
    ac_tau_std=np.std(fifo.read_grouprops(path[order[i]])['ac_tau'])
    # Mean & std of tau_b
    tau_b_mean=np.mean(fifo.read_grouprops(path[order[i]])['tau_b'])
    tau_b_std=np.std(fifo.read_grouprops(path[order[i]])['tau_b'])
    # Mean & std of tau_b_ignore
    tau_b_ignore_mean=np.mean(fifo.read_grouprops(path[order[i]])['tau_b_ignore'])
    tau_b_ignore_std=np.std(fifo.read_grouprops(path[order[i]])['tau_b_ignore'])
    # Mean & std of len_mean
    length_mean=np.mean(fifo.read_pickprops(path_pick[order[i]])['len_mean'])
    length_std=np.std(fifo.read_pickprops(path_pick[order[i]])['len_mean'])
    # Mean & std of length_cdf
    length_cdf_mean=np.mean(fifo.read_pickprops(path_pick[order[i]])['length_cdf'])
    length_cdf_std=np.std(fifo.read_pickprops(path_pick[order[i]])['length_cdf'])
    
    # Text annotations
    axarr[i,0].text(x_text,y_text,'%1.0f (%2.1f)'%(ac_tau_mean,ac_tau_std))
    axarr[i,0].text(x_text,y_text+10,'%1.0f'%(tau_cs[order[i]]),color='red')
    axarr[i,1].text(x_text,y_text,'%1.0f (%2.1f)'%(tau_b_mean,tau_b_std))
    axarr[i,2].text(x_text,y_text,'%1.0f (%2.1f)'%(tau_b_ignore_mean,tau_b_ignore_std))
    axarr[i,3].text(x_text,y_text,'%1.0f (%2.1f)'%(length_mean,length_std))
    axarr[i,4].text(x_text,y_text,'%1.0f (%2.1f)'%(length_cdf_mean,length_cdf_std))
    
    # ac_tau
    axarr[i,0].hist(fifo.read_grouprops(path[order[i]])['ac_tau'],bins='fd',label=labels[order[i]])
    # tau_b
    axarr[i,1].hist(fifo.read_grouprops(path[order[i]])['tau_b'],bins='fd',label=labels[order[i]])
    # tau_b_ignore
    axarr[i,2].hist(fifo.read_grouprops(path[order[i]])['tau_b_ignore'],bins='fd',label=labels[order[i]])
    # len_mean
    axarr[i,3].hist(fifo.read_pickprops(path_pick[order[i]])['len_mean'],bins='fd',label=labels[order[i]])
    # len_mean
    axarr[i,4].hist(fifo.read_pickprops(path_pick[order[i]])['length_cdf'],bins='fd',label=labels[order[i]])
    
    #Legend
    axarr[i,4].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


#%%
# Single docking sites, k_on=2e6(Ms)^-1 i.e. diff. tau_d for diff. cs
# tau_d for various concentrations
figure_size=np.array([11.69,8.27])*1.5
f1=plt.figure(1,figsize=figure_size)
f1.clf()

# Define axarr for plotting, shared x and y axis
axarr = f1.subplots(4,4, sharex=True, sharey=True)
# Control hspace 
f1.subplots_adjust(hspace=0.1)
# Axes limits
axarr[0,1].set_xlim([100,1e5])
# Axes scale
for ax in axarr.flat:
    ax.set_xscale('log')
# Labels
axarr[3,0].set_xlabel(r'$\tau_d$'+' [Frames]')
axarr[3,1].set_xlabel(r'$\tau_{d,cdf}$'+' [Frames]')
axarr[3,2].set_xlabel(r'$dark_{mean}$'+' [Frames]')
axarr[3,3].set_xlabel(r'$dark_{cdf}$'+' [Frames]')
for i in range(0,4):
    axarr[i,0].set_ylabel('Counts')
# Titles
f1.suptitle('Picasso Simulate: Dark time for '+r'$\tau_b$'+' = 0,6s'+'\n'+r'$\tau_d$'+' for '+r'$k_{on}=2e6(Ms)^{-1}$'+'\n'+
            '1 docking site per pick, ignore = 5, linearized')
# Texts
x_text=120
y_text=70


# Order in which plots are shown
order=[7,4,6,5]
for i in range(0,4):
    # Mean & std of tau_b
    tau_d_mean=np.mean(fifo.read_grouprops(path[order[i]])['tau_d'])
    tau_d_std=np.std(fifo.read_grouprops(path[order[i]])['tau_d'])
    # Mean & std of tau_b
    tau_d_cdf_mean=np.mean(fifo.read_grouprops(path[order[i]])['tau_d_cdf'])
    tau_d_cdf_std=np.std(fifo.read_grouprops(path[order[i]])['tau_d_cdf'])
    # Mean & std of dark_mean
    length_mean=np.mean(fifo.read_pickprops(path_pick[order[i]])['dark_mean'])
    length_std=np.std(fifo.read_pickprops(path_pick[order[i]])['dark_mean'])
    # Mean & std of dark_cdf
    length_cdf_mean=np.mean(fifo.read_pickprops(path_pick[order[i]])['dark_cdf'])
    length_cdf_std=np.std(fifo.read_pickprops(path_pick[order[i]])['dark_cdf'])
    
    # Text annotations
    axarr[i,0].text(x_text+700,y_text,'%d (%d)'%(tau_d_mean,tau_d_std))
    axarr[i,0].text(x_text,y_text,'%d, '%(tau_ds[order[i]]),color='red')
    axarr[i,1].text(x_text,y_text,'%d (%d)'%(tau_d_cdf_mean,tau_d_cdf_std))
    axarr[i,2].text(x_text,y_text,'%d (%d)'%(length_mean,length_std))
    axarr[i,3].text(x_text,y_text,'%d (%d)'%(length_cdf_mean,length_cdf_std))
    
    # tau_d
    axarr[i,0].hist(fifo.read_grouprops(path[order[i]])['tau_d'],bins='fd',label=labels[order[i]])
    # tau_d_cdf
    axarr[i,1].hist(fifo.read_grouprops(path[order[i]])['tau_d_cdf'],bins='fd',label=labels[order[i]])
    # dark_mean
    axarr[i,2].hist(fifo.read_pickprops(path_pick[order[i]])['dark_mean'],bins='fd',label=labels[order[i]])
    # dark_cdf
    axarr[i,3].hist(fifo.read_pickprops(path_pick[order[i]])['dark_cdf'],bins='fd',label=labels[order[i]])    
    
    #Legend
    axarr[i,3].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    