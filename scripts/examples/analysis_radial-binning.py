
#################################################### Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting

import os #platform independent paths
import importlib

# Load user defined functions
import locs_groupprops as l2grp
import file_formats as fifo
import fitfunc
# Reload modules
importlib.reload(l2grp)
importlib.reload(fitfunc)
# Set plot style
plt.style.use('classic')

###################################################Experimental settings
CycleTime=[0.1]*4 # Aquisition cycle time [s]
# Define labels
labels=['test']*3


##################################################################################################### File read in
# Define folder of locs.hdf5 file
dir_names=[] 
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/D134/18-06-14/s02_10nM_IS+10xtele_P15_exp100_fr25k_noautofoc_1/18-06-15_JS/']) ##TRX/PCA/PCD
#dir_names.extend(['/fs/pool/pool-schwille-paint/Data/D134/18-06-14/s01_10nM_OxCat+10xtele_P15_exp100_fr25k_noautofoc_1/18-06-15_JS/']) #CatOx autofoc off

# Define names of locs_picked.hdf5 file
file_names=[] 
file_names.extend(['s02_10nM_IS+10xtele_P15_exp100_fr25k_noautofoc_1_MMStack_Pos0.ome_locs_picked_groupprops.hdf5']) #TRX/PCA/PCD
#file_names.extend(['s01_10nM_OxCat+10xtele_P15_exp100_fr25k_noautofoc_1_MMStack_Pos0.ome_locs_picked_groupprops.hdf5']) #CatOx autofoc off

# Create list of paths
path=[]
for i in range(0,len(file_names)):
    path.extend([os.path.join(dir_names[i],file_names[i])])
    
#################################################################################################### Read, Filter, Means&Stds
groupprops_list=[]
groupprops_list_filter=[]
groupprops_list_filter_out=[]
means=np.zeros([len(path),],dtype=fifo.read_groupprops(path[i]).dtype)
stds=np.zeros([len(path),],dtype=fifo.read_groupprops(path[i]).dtype)


for p in range(0,len(path)):  
    ################# Unfiltered groupprops to list
    groupprops=fifo.read_groupprops(path[p])# Single groupprops
    groupprops_list.append(groupprops)# List containing all groupprops arrays (unfiltered)

    ################# Set filters
    # Set single filters, always one single filter has to be active!
    istrue_single=[]
    istrue_single.append(groupprops_list[p]['std_frame']>6800)
    istrue_single.append(groupprops_list[p]['mean_frame']>1000)
    istrue_single.append(groupprops_list[p]['mean_frame']<100000)
#    istrue_single.append(groupprops_list[p]['mono_tau']<150) # Uncomment for usage
    # Combine single filters to overall filter (and not)
    istrue=np.all(istrue_single,axis=0)
    istrue_not=~np.array(istrue)
    
    ################# Filtered groupprops to list
    groupprops_list_filter.append(groupprops_list[p][:][istrue]) # Create list of filtered groupprops
    groupprops_list_filter_out.append(groupprops_list[p][:][istrue_not]) # Create list of out-filtered groupprops
    
#    ################# Mean and std of all fields in filtered groupprops 
    for name in groupprops.dtype.names: 
        means[name][p]=np.mean(groupprops_list_filter[p][name])
        stds[name][p]=np.std(groupprops_list_filter[p][name])#/np.sqrt(len(groupprops_list_filter[p])) #Uncomment to switch between std and standard error of the mean

################################################### Show filtering   
p=0
#field='std_frame'
#
#f=plt.figure(num=10)
#f.clear()
#ax=f.add_subplot(1,1,1)
#ax.hist(groupprops_list_filter[p][field][:],bins='fd',label='p=%.0i, '%(p)+field)
#ax.legend(loc=1)

################################################### Radial& Intensity binning
# enter laser power measured behind objective
power=1.9e-3 # 20% 561 D134 [W]
# enter pixel size 
px_size=13e-6 # [m]
# Number of bins
NoBins=12

# parameters from Gaussian fit to laser intensity profile
x_0=239.34 # profile center in x
y_0=236.70 # profile center in y
theta=43.76 # rotation angle (clockwise)
sig_x=219.12/np.sqrt(2) # effective beam diameter in x
sig_y=253.15/np.sqrt(2) # effective beam diameter in x

popt_I=[1,x_0,y_0,sig_x,sig_y,theta,0]
popt_I[0]=(2*power)/(np.pi*(np.mean(popt_I[3:4])*px_size*1e-2)**2)*1e-7 #[kW/cm^2]


Sample_I_bin=np.empty([len(groupprops_list_filter),NoBins,1])
I_bin=np.empty([len(groupprops_list_filter),NoBins,2])
ac_tau_I_bin=np.empty([len(groupprops_list_filter),NoBins,2])
ac_A_I_bin=np.empty([len(groupprops_list_filter),NoBins,2])
tau_b_lin_I_bin=np.empty([len(groupprops_list_filter),NoBins,2])
tau_d_lin_I_bin=np.empty([len(groupprops_list_filter),NoBins,2])


for p in range(0,len(groupprops_list)):
    # Get intensity for every group
    I=fitfunc.gaussian_2D((groupprops_list_filter[p]['mean_x'],groupprops_list_filter[p]['mean_y']),*popt_I)
    # Bin Intensity in NoBins
    bins_I=np.linspace(np.min(I),np.max(I),NoBins+1)
    # Get groups in I according to binning
    digitized_I=np.digitize(I,bins_I)
    
    Sample_I_bin[p,:,0] = np.array([len(digitized_I[digitized_I == i]) for i in range(1,NoBins+1)])
    I_bin[p,:,0] = np.array([I[digitized_I == i].mean() for i in range(1, len(bins_I))])
    I_bin[p,:,1] = np.array([I[digitized_I == i].std() for i in range(1, len(bins_I))])
    
    ac_tau_I_bin[p,:,0] = np.array([groupprops_list_filter[p]['ac_tau'][digitized_I == i].mean() for i in range(1, len(bins_I))])
    ac_tau_I_bin[p,:,1] = np.array([groupprops_list_filter[p]['ac_tau'][digitized_I == i].std() for i in range(1, len(bins_I))])
    ac_A_I_bin[p,:,0] = np.array([groupprops_list_filter[p]['ac_A'][digitized_I == i].mean() for i in range(1, len(bins_I))])
    ac_A_I_bin[p,:,1] = np.array([groupprops_list_filter[p]['ac_A'][digitized_I == i].std() for i in range(1, len(bins_I))])    
    tau_b_lin_I_bin[p,:,0] = np.array([groupprops_list_filter[p]['tau_b_lin'][digitized_I == i].mean() for i in range(1, len(bins_I))])
    tau_b_lin_I_bin[p,:,1] = np.array([groupprops_list_filter[p]['tau_b_lin'][digitized_I == i].std() for i in range(1, len(bins_I))])
    tau_d_lin_I_bin[p,:,0] = np.array([groupprops_list_filter[p]['tau_d_lin'][digitized_I == i].mean() for i in range(1, len(bins_I))])
    tau_d_lin_I_bin[p,:,1] = np.array([groupprops_list_filter[p]['tau_d_lin'][digitized_I == i].std() for i in range(1, len(bins_I))])        
    

################################################### Plotting
#Define title for plots    
title='T=23C, P=20%/1.9mW, 3x tele, 10nM P1modACy3B'
#%%    
################################################### ac_tau
f=plt.figure(num=1)
f.clear()
ax=f.add_subplot(1,1,1)


for p in range(0, len(groupprops_list)):
    ax.errorbar(I_bin[p,:,0],ac_tau_I_bin[p,:,0]*CycleTime[p],yerr=np.divide(ac_tau_I_bin[p,:,1],np.sqrt(Sample_I_bin[p,:,0]))*CycleTime[p],label=labels[p])

ax.set_xscale('log')
#ax.set_ylim([0,7])

ax.set_xlabel(r'Intensity $[kW/cm^2]$')
ax.set_ylabel(r'ac_tau [s]')
ax.set_title(title,fontsize=8)
ax.legend()


#%%    
################################################### ac_A
f=plt.figure(num=2)
f.clear()
ax=f.add_subplot(1,1,1)


for p in range(0, len(groupprops_list)):
    ax.errorbar(I_bin[p,:,0],ac_A_I_bin[p,:,0],yerr=np.divide(ac_A_I_bin[p,:,1],np.sqrt(Sample_I_bin[p,:,0])),label=labels[p])

ax.set_xscale('log')
#ax.set_ylim(ylim)

ax.set_xlabel(r'Intensity $[kW/cm^2]$')
ax.set_ylabel(r'ac_A')
ax.set_title(title,fontsize=8)
ax.legend(loc=2)
#%%    
################################################### tau_b_lin
f=plt.figure(num=3)
f.clear()
ax=f.add_subplot(1,1,1)


for p in range(0, len(groupprops_list)):
    ax.errorbar(I_bin[p,:,0],tau_b_lin_I_bin[p,:,0]*CycleTime[p],yerr=np.divide(tau_b_lin_I_bin[p,:,1],np.sqrt(Sample_I_bin[p,:,0]))*CycleTime[p],label=labels[p])

ax.set_xscale('log')
#ax.set_ylim([0,7])

ax.set_xlabel(r'Intensity $[kW/cm^2]$')
ax.set_ylabel(r'tau_b_lin [s]')
ax.set_title(title,fontsize=8)
ax.legend()


#%%    
################################################### tau_d_lin
f=plt.figure(num=4)
f.clear()
ax=f.add_subplot(1,1,1)


for p in range(0, len(groupprops_list)):
    ax.errorbar(I_bin[p,:,0],tau_d_lin_I_bin[p,:,0]*CycleTime[p],yerr=np.divide(tau_d_lin_I_bin[p,:,1]*CycleTime[p],np.sqrt(Sample_I_bin[p,:,0])),label=labels[p])

ax.set_xscale('log')
#ax.set_ylim(ylim)

ax.set_xlabel(r'Intensity $[kW/cm^2]$')
ax.set_ylabel(r'tau_d_lin [s]')
ax.set_title(title,fontsize=8)
ax.legend(loc=2)

