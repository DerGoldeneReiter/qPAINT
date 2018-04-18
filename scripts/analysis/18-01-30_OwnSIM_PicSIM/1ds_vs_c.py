# Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os #platform independent paths
import matplotlib as mpl

# Load user defined functions
import file_formats as fifo

##################################################################################################### Set paths to files and labels
# Define folder of locs.hdf5 file
dir_name=[r'E:\Projects\qPAINT\data\18-01-30_OwnSimulation_1ds\12frames']*4+\
[r'E:\Projects\qPAINT\data\18-01-30_OwnSimulation_1ds\100frames']*4+\
[r'E:\Projects\qPAINT\data\17-12-12_Simulation2x2Grid\12frames\ng500_box7_1ds']*4+\
[r'E:\Projects\qPAINT\data\17-12-12_Simulation2x2Grid\100frames\ng500_box7_1ds']*4
# Define names of locs.hdf5 file
file_names=['01_locs_picked_groupprops.hdf5']
file_names.extend(['02_locs_picked_groupprops.hdf5'])
file_names.extend(['03_locs_picked_groupprops.hdf5'])
file_names.extend(['04_locs_picked_groupprops.hdf5'])
file_names=2*file_names+2*[file_names[i] for i in [3,0,2,1]]
# Define labels
labels=['1nM']
labels.extend(['5nM'])
labels.extend(['10nM'])
labels.extend(['20nM'])
labels=4*labels
# Create list of paths
path=[]
for i in range(0,len(file_names)):
    path.extend([os.path.join(dir_name[i],file_names[i])])
    
#################################################################################################### Simulation input
concs=np.array([float(label.replace('nM',''))*1e-9 for label in labels]) # List of concentrations as floats
tau_bs=np.array([12.]*4+[100.]*4+[12.]*4+[100.]*4) # Bright time [Frames]
tau_ds=[10000.,2000.,1000.,500.]
tau_ds=np.array(4*tau_ds)
NoFrames=[100000]*16
NoDocks=np.array([1]*16)
tau_cs=np.divide(1,np.power(tau_bs,-1)+np.power(tau_ds,-1)) # Autocorr tau [s]
p_infs=NoDocks*np.divide(tau_bs,tau_ds)

#################################################################################################### Read, Filter, Means&Stds
groupprops_list=[]
groupprops_list_filter=[]
means=np.zeros([len(path),],dtype=fifo.read_groupprops(path[i]).dtype)
stds=np.zeros([len(path),],dtype=fifo.read_groupprops(path[i]).dtype)

for i in range(0,len(path)):  
    ################# Unfiltered groupprops to list
    groupprops=fifo.read_groupprops(path[i])# Single groupprops
    groupprops_list.extend(groupprops)# List containing all groupprops arrays (unfiltered)

    ################# Filter
    groupprops=groupprops[:][groupprops['ac_tau_d']<NoFrames[i]] # Remove all entries for which ac_tau_d>NoFrames
    
    ################# Filtered groupprops to list
    groupprops_list_filter.extend(groupprops)# List containing all groupprops arrays (unfiltered)
    
    ################# Mean and std of all fields in groupprops   
    for name in groupprops.dtype.names: 
        means[name][i]=np.mean(groupprops[name])
        stds[name][i]=np.std(groupprops[name])

#################################################################################################### Set plot style
fontsize=15
# Figure
mpl.rcParams['figure.figsize']=[7,7]
# Linewidths
mpl.rcParams['axes.linewidth']=1
mpl.rcParams['lines.linewidth']=3
mpl.rcParams['lines.markersize']=8
mpl.rcParams['errorbar.capsize']=6
mpl.rcParams['xtick.major.width']=1
# Font
mpl.rcParams['xtick.labelsize']=fontsize
mpl.rcParams['ytick.labelsize']=fontsize
mpl.rcParams['figure.titlesize']=fontsize
mpl.rcParams['legend.fontsize']=fontsize
mpl.rcParams['axes.titlesize']=fontsize
mpl.rcParams['axes.labelsize']=fontsize
mpl.rcParams['font.sans-serif']='Times New Roman'
#%% 
################################################################################################# p_inf
plt.clf()
f,axarr=plt.subplots(2,1,num=1)
f.subplots_adjust(left=0.13,bottom=0.1,top=0.93,right=0.9,hspace=0.25)

f.suptitle(r'$p_{\infty}$'+' : OUT vs. IN')
axarr[1].set_xlabel(r'$p_{\infty,in}$')
axarr[0].set_ylabel(r'$p_{\infty,out}$')
axarr[1].set_ylabel(r'$p_{\infty,out}$')
# 12 frames
pltrange=range(0,4)
axarr[0].errorbar(p_infs[pltrange],means['ac_p_inf'][pltrange],yerr=stds['ac_p_inf'][pltrange],fmt='o',label='OwnSIM',color='blue')
pltrange=range(8,12)
axarr[0].errorbar(p_infs[pltrange],means['ac_p_inf'][pltrange],yerr=stds['ac_p_inf'][pltrange],fmt='o',label='PicSIM',color='red')
axarr[0].plot(p_infs[pltrange],p_infs[pltrange],color='black')
axarr[0].legend()
# 100 frames
pltrange=range(4,8)
axarr[1].errorbar(p_infs[pltrange],means['ac_p_inf'][pltrange],yerr=stds['ac_p_inf'][pltrange],fmt='o',label='OwnSIM',color='blue')
pltrange=range(12,16)
axarr[1].errorbar(p_infs[pltrange],means['ac_p_inf'][pltrange],yerr=stds['ac_p_inf'][pltrange],fmt='o',label='PicSIM',color='red')
axarr[1].plot(p_infs[pltrange],p_infs[pltrange],color='black')
axarr[1].legend()

#%% 
################################################################################################# tau_c
plt.clf()
f,axarr=plt.subplots(2,1,num=1)
f.subplots_adjust(left=0.13,bottom=0.1,top=0.93,right=0.9,hspace=0.25)

f.suptitle(r'$\tau_{c}$'+' : OUT vs. IN')
axarr[1].set_xlabel(r'$\tau_{c,in}$')
axarr[0].set_ylabel(r'$\tau_{c,out}$')
axarr[1].set_ylabel(r'$\tau_{c,out}$')
# 12 frames
pltrange=range(0,4)
axarr[0].errorbar(tau_cs[pltrange],means['ac_tau'][pltrange],yerr=stds['ac_tau'][pltrange],fmt='o',label='OwnSIM',color='blue')
pltrange=range(8,12)
axarr[0].errorbar(tau_cs[pltrange],means['ac_tau'][pltrange],yerr=stds['ac_tau'][pltrange],fmt='o',label='PicSIM',color='red')
axarr[0].plot(tau_cs[pltrange],tau_cs[pltrange],color='black')
axarr[0].legend()
# 100 frames
pltrange=range(4,8)
axarr[1].errorbar(tau_cs[pltrange],means['ac_tau'][pltrange],yerr=stds['ac_tau'][pltrange],fmt='o',label='OwnSIM',color='blue')
pltrange=range(12,16)
axarr[1].errorbar(tau_cs[pltrange],means['ac_tau'][pltrange],yerr=stds['ac_tau'][pltrange],fmt='o',label='PicSIM',color='red')
axarr[1].plot(tau_cs[pltrange],tau_cs[pltrange],color='black')
axarr[1].legend()

#%% 
################################################################################################# tau_b
plt.clf()
f,axarr=plt.subplots(2,1,num=1)
f.subplots_adjust(left=0.13,bottom=0.1,top=0.93,right=0.9,hspace=0.25)

f.suptitle(r'$\tau_{b}$'+' : OUT vs. IN')
axarr[1].set_xlabel('Concentration [nM]')
axarr[0].set_ylabel(r'$\tau_{b,out}$')
axarr[1].set_ylabel(r'$\tau_{b,out}$')
# 12 frames
pltrange=range(0,4)
axarr[0].errorbar(concs[pltrange]*1e9,means['tau_b_lin'][pltrange],yerr=stds['tau_b_lin'][pltrange],fmt='o',label='OwnSIM',color='blue')
pltrange=range(8,12)
axarr[0].errorbar(concs[pltrange]*1e9,means['tau_b_lin'][pltrange],yerr=stds['tau_b_lin'][pltrange],fmt='o',label='PicSIM',color='red')
axarr[0].axhline(12,color='black')
axarr[0].legend()
# 100 frames
pltrange=range(4,8)
axarr[1].errorbar(concs[pltrange]*1e9,means['tau_b_lin'][pltrange],yerr=stds['tau_b_lin'][pltrange],fmt='o',label='OwnSIM',color='blue')
pltrange=range(12,16)
axarr[1].errorbar(concs[pltrange]*1e9,means['tau_b_lin'][pltrange],yerr=stds['tau_b_lin'][pltrange],fmt='o',label='PicSIM',color='red')
axarr[1].axhline(100,color='black')
axarr[1].legend()

#%% 
################################################################################################# tau_d
plt.clf()
f,axarr=plt.subplots(2,1,num=1)
f.subplots_adjust(left=0.13,bottom=0.1,top=0.93,right=0.9,hspace=0.25)

f.suptitle(r'$\tau_{d}$'+' : OUT vs. IN')
axarr[1].set_xlabel(r'$\tau_{d,in}$')
axarr[0].set_ylabel(r'$\tau_{d,out}$')
axarr[1].set_ylabel(r'$\tau_{d,out}$')
# 12 frames
pltrange=range(0,4)
axarr[0].errorbar(tau_ds[pltrange],means['tau_d_lin'][pltrange],yerr=stds['tau_d_lin'][pltrange],fmt='o',label='OwnSIM',color='blue')
pltrange=range(8,12)
axarr[0].errorbar(tau_ds[pltrange],means['tau_d_lin'][pltrange],yerr=stds['tau_d_lin'][pltrange],fmt='o',label='PicSIM',color='red')
axarr[0].plot(tau_ds[pltrange],tau_ds[pltrange],color='black')
axarr[0].legend()
# 100 frames
pltrange=range(4,8)
axarr[1].errorbar(tau_ds[pltrange],means['tau_d_lin'][pltrange],yerr=stds['tau_d_lin'][pltrange],fmt='o',label='OwnSIM',color='blue')
pltrange=range(12,16)
axarr[1].errorbar(tau_ds[pltrange],means['tau_d_lin'][pltrange],yerr=stds['tau_d_lin'][pltrange],fmt='o',label='PicSIM',color='red')
axarr[1].plot(tau_ds[pltrange],tau_ds[pltrange],color='black')
axarr[1].legend()

#%% 
################################################################################################# n_events
plt.clf()
f,axarr=plt.subplots(1,2,num=1,sharey=True)
f.subplots_adjust(left=0.13,bottom=0.1,top=0.85,right=0.9,hspace=0.25)

f.suptitle(r'$N_{events}$'+' : OUT vs. IN')
axarr[0].set_title('Bright time = 12 frames')
axarr[1].set_title('Bright time = 100 frames')
axarr[0].set_xlabel('Concentration [nM]')
axarr[1].set_xlabel('Concentration [nM]')
axarr[0].set_ylabel(r'$N_{events}$')

# 12 frames
pltrange=range(0,4)
axarr[0].errorbar(concs[pltrange]*1e9,means['n_events'][pltrange],yerr=stds['n_events'][pltrange],fmt='o',label='OwnSIM',color='blue')
pltrange=range(8,12)
axarr[0].errorbar(concs[pltrange]*1e9,means['n_events'][pltrange],yerr=stds['n_events'][pltrange],fmt='o',label='PicSIM',color='red')
axarr[0].plot(concs[pltrange]*1e9,np.divide(NoFrames[0],tau_bs[pltrange]+tau_ds[pltrange]),color='black')
axarr[0].legend()
# 100 frames
pltrange=range(4,8)
axarr[1].errorbar(concs[pltrange]*1e9,means['n_events'][pltrange],yerr=stds['n_events'][pltrange],fmt='o',label='OwnSIM',color='blue')
pltrange=range(12,16)
axarr[1].errorbar(concs[pltrange]*1e9,means['n_events'][pltrange],yerr=stds['n_events'][pltrange],fmt='o',label='PicSIM',color='red')
axarr[1].plot(concs[pltrange]*1e9,np.divide(NoFrames[0],tau_bs[pltrange]+tau_ds[pltrange]),color='black')
axarr[1].legend()