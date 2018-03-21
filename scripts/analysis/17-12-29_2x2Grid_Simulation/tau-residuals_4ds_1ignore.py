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
dir_name=['E:/Projects/qPAINT/data/17-12-12_Simulation2x2Grid/FS_18-01-03_4ds_1ignore_lin']*8


# Define names of locs.hdf5 file
file_names=['01_locs_picked_groupprops.hdf5']
file_names.extend(['02_locs_picked_groupprops.hdf5'])
file_names.extend(['03_locs_picked_groupprops.hdf5'])
file_names.extend(['04_locs_picked_groupprops.hdf5'])
file_names.extend(['05_locs_picked_groupprops.hdf5'])
file_names.extend(['06_locs_picked_groupprops.hdf5'])
file_names.extend(['07_locs_picked_groupprops.hdf5'])
file_names.extend(['08_locs_picked_groupprops.hdf5'])
file_names=file_names

# Define labels
labels=['5nM']
labels.extend(['20nM'])
labels.extend(['10nM'])
labels.extend(['1nM'])
labels=labels+labels

# Create list of paths
path=list()
for i in range(0,len(file_names)):
    path.append(os.path.join(dir_name[i],file_names[i]))
    
# Delete variables not used
del file_names
del dir_name

# Reorder list according to concentrations
order=[7,4,6,5,3,0,2,1]
path=[path[i] for i in order]
labels=[labels[i] for i in order]

# Overall simulation parameters
concs=np.array([float(label.replace('nM',''))*1e-9 for label in labels]) # List of concentrations as floats
tau_bs=np.array([12.]*4+[100.]*4) # Bright time [Frames]
tau_ds=[10000.,2000.,1000.,500.]
tau_ds=np.array(tau_ds+tau_ds)
 # Dark time [Frames]
NoFrames=[100000]*8
NoDocks=[4]*8
# Calculate expected tau_c
tau_cs=np.divide(1,np.power(tau_bs,-1)+np.power(tau_ds,-1)) # Autocorr tau [s]
p_infs=np.divide(tau_bs,tau_ds)*NoDocks[0]

# Set general font size
matplotlib.rcParams.update({'font.size': 15})

# Initialize mean&std of important parameters
G0=np.empty([len(path),2])
Ginf=np.empty([len(path),2])
tau_c=np.empty([len(path),2])

tau_b=np.empty([len(path),2])
tau_b_ignore=np.empty([len(path),2])
tau_d=np.empty([len(path),2])

tau_b_lin=np.empty([len(path),2])
tau_b_lin_ignore=np.empty([len(path),2])
tau_d_lin=np.empty([len(path),2])

p_inf=np.empty([len(path),2])
ac_tau_b=np.empty([len(path),2])
ac_tau_d=np.empty([len(path),2])

# Initialize mean&std of important parameters
for i in range(0,len(path)):
    G0[i,0]=np.mean(fifo.read_grouprops(path[i])['ac_A'])
    G0[i,1]=np.std(fifo.read_grouprops(path[i])['ac_A'])
    Ginf[i,0]=np.mean(fifo.read_grouprops(path[i])['ac_offset'])
    Ginf[i,1]=np.std(fifo.read_grouprops(path[i])['ac_offset'])
    tau_c[i,0]=np.mean(fifo.read_grouprops(path[i])['ac_tau'])
    tau_c[i,1]=np.std(fifo.read_grouprops(path[i])['ac_tau'])
    
    tau_b[i,0]=np.mean(fifo.read_grouprops(path[i])['tau_b'])
    tau_b[i,1]=np.std(fifo.read_grouprops(path[i])['tau_b'])
    tau_b_ignore[i,0]=np.mean(fifo.read_grouprops(path[i])['tau_b_ignore'])
    tau_b_ignore[i,1]=np.std(fifo.read_grouprops(path[i])['tau_b_ignore'])
    tau_d[i,0]=np.mean(fifo.read_grouprops(path[i])['tau_d'])
    tau_d[i,1]=np.std(fifo.read_grouprops(path[i])['tau_d'])
    
    tau_b_lin[i,0]=np.mean(fifo.read_grouprops(path[i])['tau_b_lin'])
    tau_b_lin[i,1]=np.std(fifo.read_grouprops(path[i])['tau_b_lin'])
    tau_b_lin_ignore[i,0]=np.mean(fifo.read_grouprops(path[i])['tau_b_lin_ignore'])
    tau_b_lin_ignore[i,1]=np.std(fifo.read_grouprops(path[i])['tau_b_lin_ignore'])
    tau_d_lin[i,0]=np.mean(fifo.read_grouprops(path[i])['tau_d_lin'])
    tau_d_lin[i,1]=np.std(fifo.read_grouprops(path[i])['tau_d_lin'])
    
    p_inf[i,0]=np.mean(fifo.read_grouprops(path[i])['ac_p_inf'])
    p_inf[i,1]=np.std(fifo.read_grouprops(path[i])['ac_p_inf'])
    ac_tau_b[i,0]=np.mean(fifo.read_grouprops(path[i])['ac_tau_b'])
    ac_tau_b[i,1]=np.std(fifo.read_grouprops(path[i])['ac_tau_b'])
    ac_tau_d[i,0]=np.mean(fifo.read_grouprops(path[i])['ac_tau_d'][fifo.read_grouprops(path[i])['ac_tau_d']<NoFrames[i]]) # Disregard all dark times longer than total number of frames
    ac_tau_d[i,1]=np.std(fifo.read_grouprops(path[i])['ac_tau_d'][fifo.read_grouprops(path[i])['ac_tau_d']<NoFrames[i]])
    
#%% 
# Is G quotient equal to p_inf??
figsize=[8,6]
f1=plt.figure(1, figsize=figsize)
f1.clf()
ax=f1.subplots(1,1)
f1.subplots_adjust(left=0.15,bottom=0.15)
# Limits
#ax.set_xlim([-0.02,0.43])
#ax.set_ylim([-0.02,0.43])
# Labels
ax.set_xlabel(r'$\frac{\tau_{b,in}}{\tau_{d,in}}$',fontsize=20)
ax.set_ylabel(r'$\frac{G_{\inf,out}}{G_{0,out}}$',fontsize=20)
f1.suptitle(r'$\frac{G_{\inf}}{G_{0}}$'+' : Own simulation, 4 docking sites')

pltrange=np.arange(0,4)
ax.errorbar(concs[pltrange]*1e9,p_inf[pltrange,0]-p_infs[pltrange]*0.8,
            yerr=p_inf[pltrange,1],
            fmt='o',markersize=6,label='12frames',capsize=5)
pltrange=np.arange(4,8)
ax.errorbar(concs[pltrange]*1e9,p_inf[pltrange,0]-p_infs[pltrange]*0.69,
            yerr=p_inf[pltrange,1],
            fmt='o',markersize=6,label='100frames',capsize=5)
ax.axhline(0,dashes=(1,1),color='black')

ax.legend(loc=3)

#%%
# tau_c ??
figsize=[8,6]
f1=plt.figure(1, figsize=figsize)
f1.clf()
ax=f1.subplots(2,1,sharex=True)
f1.subplots_adjust(left=0.11,bottom=0.1,right=0.8)

# Labels
ax[1].set_xlabel(r'$\tau_{c,in}$',fontsize=20)
ax[0].set_ylabel(r'$\tau_{c,out}$',fontsize=20)
ax[1].set_ylabel(r'$\tau_{c,out}$',fontsize=20)
# Title
f1.suptitle(r'$\tau_c$'+' : Picasso simulation, 4 docking sites')

pltrange=np.arange(0,4)
ax[0].errorbar(concs[pltrange]*1e9,tau_c[pltrange,0]-tau_cs[pltrange],
            yerr=tau_c[pltrange,1],
            fmt='o',markersize=6,label='12frames',capsize=5)
ax[0].plot(concs[pltrange]*1e9,tau_cs[pltrange]*0.1,
  '--',color='black',label='< 10%')
ax[0].plot(concs[pltrange]*1e9,-tau_cs[pltrange]*0.1,
  '--',color='black')
ax[0].axhline(0,color='black')
ax[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

pltrange=np.arange(4,8)
ax[1].errorbar(concs[pltrange]*1e9,tau_c[pltrange,0]-tau_cs[pltrange],
            yerr=tau_c[pltrange,1],
            fmt='o',markersize=6,label='100frames',capsize=5)
ax[1].plot(concs[pltrange]*1e9,tau_cs[pltrange]*0.1,
  '--',color='black',label='< 10%')
ax[1].plot(concs[pltrange]*1e9,-tau_cs[pltrange]*0.1,
  '--',color='black')
ax[1].axhline(0,color='black')
ax[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#%%
# tau_d absolute error
figsize=np.array([11.69,8.27])*1.5
f1=plt.figure(1, figsize=figsize)
f1.clf()
axarr=f1.subplots(2,1,sharex=True)
f1.subplots_adjust(left=0.1,bottom=0.1,right=0.8)

# Scales
axarr[0].set_xscale('log')
axarr[0].set_yscale('symlog')
axarr[1].set_xscale('log')
axarr[1].set_yscale('symlog')
# Titles
f1.suptitle('Dark time, '+r'$\tau_d$'+' : Own simulation, 4 docking sites')
axarr[0].set_title(r'$\tau_b$'+' = 20 frames', loc='left')
axarr[1].set_title(r'$\tau_b$'+' = 100 frames',loc='left')
# Labels
axarr[1].set_xlabel('Concentration [nM]',fontsize=15)
axarr[0].set_ylabel('Residual',fontsize=15)
axarr[1].set_ylabel('Residual',fontsize=15)

# 20 frames
# ac 
axarr[0].errorbar(concs[0:8]*1e9,ac_tau_d[0:8,0]-tau_ds[0:8],yerr=ac_tau_d[0:8,1],fmt='o',markersize=6,label='acPAINT',color='black',capsize=5)
# qPAINT old
axarr[0].errorbar(concs[0:8]*1e9,tau_d[0:8,0]-tau_ds[0:8]/4,yerr=tau_d[0:8,1],fmt='o',markersize=6,label='qPAINT, old',color='green',capsize=5)
# qPAINT new
axarr[0].errorbar(concs[0:8]*1e9,tau_d_lin[0:8,0]-tau_ds[0:8]/4,yerr=tau_d_lin[0:8,1],fmt='o',markersize=6,label='qPAINT, new',color='blue',capsize=5)
# Draw zero line
axarr[0].axhline(0,color='black',dashes=(1,1))
axarr[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

# 20 frames
# ac 
axarr[1].errorbar(concs[8:16]*1e9,ac_tau_d[8:16,0]-tau_ds[8:16],yerr=ac_tau_d[8:16,1],fmt='o',markersize=6,label='acPAINT',color='black',capsize=5)
# qPAINT old
axarr[1].errorbar(concs[8:16]*1e9,tau_d[8:16,0]-tau_ds[8:16]/4,yerr=tau_d[8:16,1],fmt='o',markersize=6,label='qPAINT, old',color='green',capsize=5)
# qPAINT new
axarr[1].errorbar(concs[8:16]*1e9,tau_d_lin[8:16,0]-tau_ds[8:16]/4,yerr=tau_d_lin[8:16,1],fmt='o',markersize=6,label='qPAINT, new',color='blue',capsize=5)
# Draw zero line
axarr[1].axhline(0,color='black',dashes=(1,1))
axarr[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

#%%
# tau_b
figsize=np.array([11.69,8.27])*1.5
f1=plt.figure(1, figsize=figsize)
f1.clf()
axarr=f1.subplots(2,1,sharex=True)
f1.subplots_adjust(left=0.1,bottom=0.1,right=0.8)

# Scales
axarr[0].set_xscale('log')
axarr[0].set_yscale('linear')
axarr[1].set_xscale('log')
axarr[1].set_yscale('linear')
# Titles
f1.suptitle('Dark time, '+r'$\tau_d$'+' : Own simulation, 4 docking sites')
axarr[0].set_title(r'$\tau_b$'+' = 20 frames', loc='left')
axarr[1].set_title(r'$\tau_b$'+' = 100 frames',loc='left')
# Labels
axarr[1].set_xlabel('Concentration [nM]',fontsize=15)
axarr[0].set_ylabel('Relative residual',fontsize=15)
axarr[1].set_ylabel('Relative residual',fontsize=15)

# 20 frames
# ac 
axarr[0].errorbar(concs[0:8]*1e9,np.divide(ac_tau_d[0:8,0]-tau_ds[0:8],tau_ds[0:8]),
     yerr=np.divide(ac_tau_d[0:8,1],tau_ds[0:8]),
     fmt='o',markersize=6,label='acPAINT',color='black',capsize=5)
# qPAINT old
axarr[0].errorbar(concs[0:8]*1e9,np.divide(tau_d[0:8,0]-tau_ds[0:8]/4,tau_ds[0:8]/4),
     yerr=np.divide(tau_d[0:8,1],tau_ds[0:8]/4),
     fmt='o',markersize=6,label='qPAINT, old',color='green',capsize=5)
# qPAINT new
axarr[0].errorbar(concs[0:8]*1e9,np.divide(tau_d_lin[0:8,0]-tau_ds[0:8]/4,tau_ds[0:8]/4),
     yerr=np.divide(tau_d_lin[0:8,1],tau_ds[0:8]/4),
     fmt='o',markersize=6,label='qPAINT, new',color='blue',capsize=5)

# Draw zero line
axarr[0].axhline(0,color='black')
axarr[0].axhline(0.1,color='black',dashes=(1,1),label='< 10%')
axarr[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

# 100 frames
# ac 
axarr[1].errorbar(concs[8:16]*1e9,np.divide(ac_tau_d[8:16,0]-tau_ds[8:16],tau_ds[8:16]),
     yerr=np.divide(ac_tau_d[8:16,1],tau_ds[8:16]),
     fmt='o',markersize=6,label='acPAINT',color='black',capsize=5)
# qPAINT old
axarr[1].errorbar(concs[8:16]*1e9,np.divide(tau_d[8:16,0]-tau_ds[8:16]/4,tau_ds[8:16]/4),
     yerr=np.divide(tau_d[8:16,1],tau_ds[8:16]/4),
     fmt='o',markersize=6,label='qPAINT, old',color='green',capsize=5)
# qPAINT new
axarr[1].errorbar(concs[8:16]*1e9,np.divide(tau_d_lin[8:16,0]-tau_ds[8:16]/4,tau_ds[8:16]/4),
     yerr=np.divide(tau_d_lin[8:16,1],tau_ds[8:16]/4),
     fmt='o',markersize=6,label='qPAINT, new',color='blue',capsize=5)

# Draw zero line
axarr[1].axhline(0,color='black')
axarr[1].axhline(0.1,color='black',dashes=(1,1),label='< 10%')
axarr[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#%%
# tau_b
figsize=np.array([11.69,8.27])*1.5
f1=plt.figure(1, figsize=figsize)
f1.clf()
axarr=f1.subplots(2,1,sharex=True)
f1.subplots_adjust(left=0.1,bottom=0.1,right=0.8)

# Scales
axarr[0].set_xscale('log')
axarr[0].set_yscale('linear')
axarr[1].set_xscale('log')
axarr[1].set_yscale('linear')
# Titles
f1.suptitle('Bright time, '+r'$\tau_b$'+' : Own simulation, 4 docking sites')
axarr[0].set_title(r'$\tau_b$'+' = 20 frames', loc='left')
axarr[1].set_title(r'$\tau_b$'+' = 100 frames',loc='left')
# Labels
axarr[1].set_xlabel('Concentration [nM]',fontsize=15)
axarr[0].set_ylabel(r'$\tau_b$',fontsize=20)
axarr[1].set_ylabel(r'$\tau_b$',fontsize=20)

pltrange=np.arange(0,4)
# 12 frames
# ac 
axarr[0].errorbar(concs[pltrange]*1e9,ac_tau_b[pltrange,0]-tau_bs[0],
     yerr=ac_tau_b[pltrange,1],
    fmt='o',markersize=6,label='ac',color='black',capsize=5)
# qPAINT old
axarr[0].errorbar(concs[pltrange]*1e9,tau_b[pltrange,0]-tau_bs[0],
     yerr=tau_b[pltrange,1],fmt='o',
     markersize=6,label='exp',color='green',capsize=5)
# qPAINT old
axarr[0].errorbar(concs[pltrange]*1e9,tau_b_lin[pltrange,0]-tau_bs[0],
     yerr=tau_b_lin[pltrange,1],
     fmt='o',markersize=6,label='lin',color='blue',capsize=5)
# Draw zero line
axarr[0].axhline(0,color='black')
axarr[0].axhline(2,color='black',dashes=(1,1),label='< 10%')
axarr[0].axhline(-2,color='black',dashes=(1,1))
axarr[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

pltrange=np.arange(4,8)
# 100 frames
# ac 
axarr[1].errorbar(concs[pltrange]*1e9,ac_tau_b[pltrange,0]-tau_bs[4],
     yerr=ac_tau_b[pltrange,1],
    fmt='o',markersize=6,label='ac',color='black',capsize=5)
# qPAINT old
axarr[1].errorbar(concs[pltrange]*1e9,tau_b[pltrange,0]-tau_bs[4],
     yerr=tau_b[pltrange,1],fmt='o',
     markersize=6,label='exp',color='green',capsize=5)
# qPAINT old
axarr[1].errorbar(concs[pltrange]*1e9,tau_b_lin[pltrange,0]-tau_bs[4],
     yerr=tau_b_lin[pltrange,1],
     fmt='o',markersize=6,label='lin',color='blue',capsize=5)
# Draw zero line
axarr[1].axhline(0,color='black')
axarr[1].axhline(10,color='black',dashes=(1,1),label='< 10%')
axarr[1].axhline(-10,color='black',dashes=(1,1))
axarr[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
