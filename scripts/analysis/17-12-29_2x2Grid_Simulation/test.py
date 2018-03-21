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
dir_name='E:/Projects/qPAINT/data/17-12-12_Simulation2x2Grid/FS_18-01-03_1ds_1ignore_lin'

# Define names of locs.hdf5 file
file_names=['01_locs_picked_groupprops.hdf5']
file_names.extend(['02_locs_picked_groupprops.hdf5'])
file_names.extend(['03_locs_picked_groupprops.hdf5'])
file_names.extend(['04_locs_picked_groupprops.hdf5'])

# Define labels
labels=['5nM']
labels.extend(['20nM'])
labels.extend(['10nM'])
labels.extend(['1nM'])

# Overall simulation parameters
concs=np.array([float(label.replace('nM',''))*1e-9 for label in labels]) # List of concentrations as floats
tau_bs=np.array([100.]*4) # Bright time [Frames]
tau_ds=np.array([10000.,2000.,1000.,500.]) # Dark time [Frames]
NoFrames=[100000]*4
NoDocks=[1]*4
NoGroups=[400]*4
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

# Reorder list according to concentrations
order=[3,0,2,1]
path=[path[i] for i in order]
labels=[labels[i] for i in order]

# Set general font size
matplotlib.rcParams.update({'font.size': 15})
#%%
# Initialize mean&std of important parameters
G0=np.empty([len(path),2])
Ginf=np.empty([len(path),2])
tau_c=np.empty([len(path),2])
tau_d=np.empty([len(path),2])

# Initialize mean&std of important parameters
for i in range(0,len(path)):
    G0[i,0]=np.mean(fifo.read_grouprops(path[i])['ac_A'])
    G0[i,1]=np.std(fifo.read_grouprops(path[i])['ac_A'])
    Ginf[i,0]=np.mean(fifo.read_grouprops(path[i])['ac_offset'])
    Ginf[i,1]=np.std(fifo.read_grouprops(path[i])['ac_offset'])
    tau_c[i,0]=np.mean(fifo.read_grouprops(path[i])['ac_tau'])
    tau_c[i,1]=np.std(fifo.read_grouprops(path[i])['ac_tau'])
    tau_d[i,0]=np.mean(fifo.read_grouprops(path[i])['tau_d'])
    tau_d[i,1]=np.std(fifo.read_grouprops(path[i])['tau_d'])
    
#%%
f1=plt.figure(1)
f1.clf()
ax=f1.subplots(1,1)
f1.subplots_adjust(left=0.15,bottom=0.15)
# Limits
#ax.set_xlim([-0.02,0.43])
#ax.set_ylim([-0.02,0.43])
# Labels
ax.set_xlabel(r'$\frac{\tau_{b,in}}{\tau_{d,in}}$',fontsize=20)
ax.set_ylabel(r'$\frac{G_{\inf,out}}{G_{0,out}}$',fontsize=20)

x=np.divide(tau_bs,tau_ds)
y=np.divide(Ginf[:,0],G0[:,0])
yerr=np.sqrt(np.power(np.divide(Ginf[:,1],Ginf[:,0]),2)+np.power(np.divide(G0[:,1],G0[:,0]),2))
yerr=np.multiply(y,yerr)

ax.errorbar(x,y,yerr=yerr,fmt='o',markersize=6,label='Own Simulation')
ax.plot(x,x,'-',color='black')

ax.legend(loc=2)

