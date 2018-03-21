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
dir_name='E:/Projects/qPAINT/data/18-01-02_OwnSimulation_4ds/0ignore_lin'

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
labels=['1nM']
labels.extend(['5nM'])
labels.extend(['10nM'])
labels.extend(['20nM'])
labels.extend(['1nM'])
labels.extend(['5nM'])
labels.extend(['10nM'])
labels.extend(['20nM'])

# List of concentrations as floats
concs=[float(label.replace('nM',''))*1e-9 for label in labels]

# Overall simulation parameters
tau_bs=np.array([100.]*4+[12.]*4) # Bright time [Frames]
tau_ds=np.array([10000.,2000.,1000.,500.,10000.,2000.,1000.,500.]) # Dark time [Frames]
NoFrames=[100000]*8
NoDocks=[4]*8
NoGroups=[400]*8
# Calculate expected tau_c
tau_cs=np.divide(1,np.power(tau_bs,-1)+np.power(tau_ds,-1)) # Autocorr tau [s]


# Create list of paths
path=list()
for file_name in file_names:
    path.append(os.path.join(dir_name,file_name))
    
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
axarr = f1.subplots(4,3, sharex=True, sharey=True)
# Control hspace 
f1.subplots_adjust(hspace=0.05)
# Axes limits
axarr[0,1].set_xlim([0,300])
# Labels
axarr[3,0].set_xlabel(r'$\tau_c$'+' [Frames]')
axarr[3,1].set_xlabel(r'$\tau_b$'+' [Frames]')
axarr[3,2].set_xlabel(r'$\tau_{b,ignore}$'+' [Frames]')
for i in range(0,4):
    axarr[i,0].set_ylabel('Counts')
# Titles
f1.suptitle('Own Simulation: Bright times\n'+r'$\tau_b$'+'= 100 frames\n'+
            '4 docking sites per pick, ignore = 0, linearized')
# Texts
x_text=150
y_text=15


# Order in which plots are shown
order=[0,1,2,3]
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
    
    # Text annotations
    axarr[i,0].text(x_text,y_text,'%d (%d)'%(ac_tau_mean,ac_tau_std))
    axarr[i,0].text(x_text,y_text+4,'%d'%(tau_cs[order[i]]),color='red')
    axarr[i,1].text(x_text,y_text,'%d (%d)'%(tau_b_mean,tau_b_std))
    axarr[i,2].text(x_text,y_text,'%d (%d)'%(tau_b_ignore_mean,tau_b_ignore_std))
    
    # ac_tau
    axarr[i,0].hist(fifo.read_grouprops(path[order[i]])['ac_tau'],bins='fd',label=labels[order[i]])
    # tau_b
    axarr[i,1].hist(fifo.read_grouprops(path[order[i]])['tau_b'],bins='fd',label=labels[order[i]])
    # tau_b_ignore
    axarr[i,2].hist(fifo.read_grouprops(path[order[i]])['tau_b_ignore'],bins='fd',label=labels[order[i]])
    
    #Legend
    axarr[i,2].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

#%%
# Single docking sites, k_on=2e6(Ms)^-1 i.e. diff. tau_d for diff. cs
# tau_d for various concentrations
figure_size=np.array([11.69,8.27])*1.5
f1=plt.figure(1,figsize=figure_size)
f1.clf()

# Define axarr for plotting, shared x and y axis
axarr = f1.subplots(4,2, sharex=True, sharey=True)
# Control hspace 
f1.subplots_adjust(hspace=0.1)
# Axes limits
axarr[0,1].set_ylim([0,28])
# Axes scale
for ax in axarr.flat:
    ax.set_xscale('log')
# Labels
axarr[3,0].set_xlabel(r'$\tau_d$'+' [Frames]')
axarr[3,1].set_xlabel(r'$\tau_{d,cdf}$'+' [Frames]')
for i in range(0,4):
    axarr[i,0].set_ylabel('Counts')
# Titles
f1.suptitle('Own Simulation: Dark time for '+r'$\tau_b$'+' = 100 frames'+'\n'+r'$\tau_d$'+' for '+r'$k_{on}=2e6(Ms)^{-1}$'+'\n'+
            '4 docking sites per pick, ignore = 0, linearized')
# Texts
x_text=100
y_text=25


# Order in which plots are shown
order=[0,1,2,3]
for i in range(0,4):
    # Mean & std of tau_b
    tau_d_mean=np.mean(fifo.read_grouprops(path[order[i]])['tau_d'])
    tau_d_std=np.std(fifo.read_grouprops(path[order[i]])['tau_d'])
    # Mean & std of tau_b
    tau_d_cdf_mean=np.mean(fifo.read_grouprops(path[order[i]])['tau_d_cdf'])
    tau_d_cdf_std=np.std(fifo.read_grouprops(path[order[i]])['tau_d_cdf'])
    
    # Text annotations
    axarr[i,0].text(x_text+70,y_text,'%d (%d)'%(tau_d_mean,tau_d_std))
    axarr[i,0].text(x_text,y_text,'%d, '%(tau_ds[order[i]]/4),color='red')
    axarr[i,1].text(x_text,y_text,'%d (%d)'%(tau_d_cdf_mean,tau_d_cdf_std))
    
    # tau_d
    axarr[i,0].hist(fifo.read_grouprops(path[order[i]])['tau_d'],bins='fd',label=labels[order[i]])
    # tau_d_cdf
    axarr[i,1].hist(fifo.read_grouprops(path[order[i]])['tau_d_cdf'],bins='fd',label=labels[order[i]])
    
     #Legend
    axarr[i,1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

#%%
# Single docking sites, tau_b set to 600ms=12frames
# tau_b and tau_c for various concentrations
figure_size=np.array([11.69,8.27])*1.5
f1=plt.figure(1,figsize=figure_size)
f1.clf()

# Define axarr for plotting, shared x and y axis
axarr = f1.subplots(4,3, sharex=True, sharey=True)
# Control hspace 
f1.subplots_adjust(hspace=0.05)
# Axes limits
#axarr[0,1].set_xlim([0,50])
# Labels
axarr[3,0].set_xlabel(r'$\tau_c$'+' [Frames]')
axarr[3,1].set_xlabel(r'$\tau_b$'+' [Frames]')
axarr[3,2].set_xlabel(r'$\tau_{b,ignore}$'+' [Frames]')
for i in range(0,4):
    axarr[i,0].set_ylabel('Counts')
# Titles
f1.suptitle('Own Simulation: Bright times\n'+r'$\tau_b$'+'= 12 frames\n'+
            '4 docking sites per pick, ignore = 0, linearized')
# Texts
x_text=15
y_text=20
# Order in which plots are shown
order=[4,5,6,7]
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
    
    # Text annotations
    axarr[i,0].text(x_text,y_text,'%1.0f (%2.1f)'%(ac_tau_mean,ac_tau_std))
    axarr[i,0].text(x_text,y_text+4,'%1.0f'%(tau_cs[order[i]]),color='red')
    axarr[i,1].text(x_text,y_text,'%1.0f (%2.1f)'%(tau_b_mean,tau_b_std))
    axarr[i,2].text(x_text,y_text,'%1.0f (%2.1f)'%(tau_b_ignore_mean,tau_b_ignore_std))
    
    # ac_tau
    axarr[i,0].hist(fifo.read_grouprops(path[order[i]])['ac_tau'],bins='fd',label=labels[order[i]])
    # tau_b
    axarr[i,1].hist(fifo.read_grouprops(path[order[i]])['tau_b'],bins='fd',label=labels[order[i]])
    # tau_b_ignore
    axarr[i,2].hist(fifo.read_grouprops(path[order[i]])['tau_b_ignore'],bins='fd',label=labels[order[i]])
    
    #Legend
    axarr[i,2].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


#%%
# Single docking sites, k_on=2e6(Ms)^-1 i.e. diff. tau_d for diff. cs
# tau_d for various concentrations
figure_size=np.array([11.69,8.27])*1.5
f1=plt.figure(1,figsize=figure_size)
f1.clf()

# Define axarr for plotting, shared x and y axis
axarr = f1.subplots(4,2, sharex=True, sharey=True)
# Control hspace 
f1.subplots_adjust(hspace=0.1)
# Axes limits
axarr[0,1].set_ylim([0,31])
# Axes scale
for ax in axarr.flat:
    ax.set_xscale('log')
# Labels
axarr[3,0].set_xlabel(r'$\tau_d$'+' [Frames]')
axarr[3,1].set_xlabel(r'$\tau_{d,cdf}$'+' [Frames]')
for i in range(0,4):
    axarr[i,0].set_ylabel('Counts')
# Titles
f1.suptitle('Own Simulation: Dark time for '+r'$\tau_b$'+' = 12 frames'+'\n'+r'$\tau_d$'+' for '+r'$k_{on}=2e6(Ms)^{-1}$'+'\n'+
            '4 docking sites per pick, ignore = 0, linearized')
# Texts
x_text=120
y_text=27


# Order in which plots are shown
order=[4,5,6,7]
for i in range(0,4):
    # Mean & std of tau_b
    tau_d_mean=np.mean(fifo.read_grouprops(path[order[i]])['tau_d'])
    tau_d_std=np.std(fifo.read_grouprops(path[order[i]])['tau_d'])
    # Mean & std of tau_b
    tau_d_cdf_mean=np.mean(fifo.read_grouprops(path[order[i]])['tau_d_cdf'])
    tau_d_cdf_std=np.std(fifo.read_grouprops(path[order[i]])['tau_d_cdf'])
    
    # Text annotations
    axarr[i,0].text(x_text+80,y_text,'%d (%d)'%(tau_d_mean,tau_d_std))
    axarr[i,0].text(x_text,y_text,'%d, '%(tau_ds[order[i]]/4),color='red')
    axarr[i,1].text(x_text,y_text,'%d (%d)'%(tau_d_cdf_mean,tau_d_cdf_std))
    
    # tau_d
    axarr[i,0].hist(fifo.read_grouprops(path[order[i]])['tau_d'],bins='fd',label=labels[order[i]])
    # tau_d_cdf
    axarr[i,1].hist(fifo.read_grouprops(path[order[i]])['tau_d_cdf'],bins='fd',label=labels[order[i]])
    
    #Legend
    axarr[i,1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    