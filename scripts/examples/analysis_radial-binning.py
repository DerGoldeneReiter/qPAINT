
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
# Load style sheet for plots
plt.style.use('~/qPAINT/styles/FoM.mplstyle')

###################################################Experimental settings
CycleTime=0.2 # Aquisition cycle time [s]
# Define labels
# Define labels
labels=[]
labels.extend(['flat @ 2.7mW'])
labels.extend(['gauss @ 2.7mW'])
labels.extend(['flat @ 5.5mW'])
labels.extend(['gauss @ 5.5mW'])
labels.extend(['flat @ 15.5mW'])
labels.extend(['gauss @ 15.5mW'])
labels.extend(['flat @ 35mW'])
labels.extend(['gauss @ 35mW'])
labels.extend(['flat @ 75.5mW'])
labels.extend(['gauss @ 75.5mW'])
##################################################################################################### File read in
# Define folder of locs.hdf5 file
dir_names=[] 
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/D042/18-08-10_FlatGaussPower/id15-80_P1-Cy3b-4nM_p100mW-9deg_flat_1/18-08-10_FS'])
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/D042/18-08-10_FlatGaussPower/id15-80_P1-Cy3b-4nM_p100mW-9deg_gauss_1/18-08-10_FS'])
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/D042/18-08-10_FlatGaussPower/id15-80_P1-Cy3b-4nM_p100mW-11deg_flat_1/18-08-10_FS/'])
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/D042/18-08-10_FlatGaussPower/id15-80_P1-Cy3b-4nM_p100mW-11deg_gauss_1/18-08-10_FS/'])
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/D042/18-08-10_FlatGaussPower/id15-80_P1-Cy3b-4nM_p100mW-17deg_flat_1/18-08-10_FS/'])
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/D042/18-08-10_FlatGaussPower/id15-80_P1-Cy3b-4nM_p100mW-17deg_gauss_1/18-08-10_FS/'])
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/D042/18-08-10_FlatGaussPower/id15-80_P1-Cy3b-4nM_p100mW-25deg_flat_1/18-08-10_FS/'])
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/D042/18-08-10_FlatGaussPower/id15-80_P1-Cy3b-4nM_p100mW-25deg_gauss_1/18-08-10_FS/'])
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/D042/18-08-10_FlatGaussPower/id15-80_P1-Cy3b-4nM_p100mW-50deg_flat_1/18-08-10_FS/'])
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/D042/18-08-10_FlatGaussPower/id15-80_P1-Cy3b-4nM_p100mW-50deg_gauss_1/18-08-10_FS/'])


# Define names of locs_picked.hdf5 file
file_names=[] 
file_names.extend(['id15-80_P1-Cy3b-4nM_p100mW-9deg_flat_1_MMStack.ome_locs_picked_groupprops.hdf5']) 
file_names.extend(['id15-80_P1-Cy3b-4nM_p100mW-9deg_gauss_1_MMStack.ome_locs_picked_groupprops.hdf5']) 
file_names.extend(['id15-80_P1-Cy3b-4nM_p100mW-11deg_flat_1_MMStack.ome_locs_picked_groupprops.hdf5']) 
file_names.extend(['id15-80_P1-Cy3b-4nM_p100mW-11deg_gauss_1_MMStack.ome_locs_picked_groupprops.hdf5']) 
file_names.extend(['id15-80_P1-Cy3b-4nM_p100mW-17deg_flat_1_MMStack.ome_locs_picked_groupprops.hdf5']) 
file_names.extend(['id15-80_P1-Cy3b-4nM_p100mW-17deg_gauss_1_MMStack.ome_locs_picked_groupprops.hdf5']) 
file_names.extend(['id15-80_P1-Cy3b-4nM_p100mW-25deg_flat_1_MMStack.ome_locs_picked_groupprops.hdf5']) 
file_names.extend(['id15-80_P1-Cy3b-4nM_p100mW-25deg_gauss_1_MMStack.ome_locs_picked_groupprops.hdf5'])
file_names.extend(['id15-80_P1-Cy3b-4nM_p100mW-50deg_flat_1_MMStack.ome_locs_picked_groupprops.hdf5']) 
file_names.extend(['id15-80_P1-Cy3b-4nM_p100mW-50deg_gauss_1_MMStack.ome_locs_picked_groupprops.hdf5'])


# Create list of paths
path=[]
for i in range(0,len(file_names)):
    path.extend([os.path.join(dir_names[i],file_names[i])])
    
#################################################################################################### Read, Filter, Means&Stds
groupprops_list=[]
groupprops_list_filter=[]
groupprops_list_filter_out=[]
means=np.zeros([len(path),],dtype=fifo.read_locs(path[i]).dtype)
stds=np.zeros([len(path),],dtype=fifo.read_locs(path[i]).dtype)

low_std_frame=[3100,3200,3200,3200,3200,3200,3200,3200,3000,2800]
low_mean_frame=[5300,5300,5300,5300,5200,5200,5200,4700,4700,4000]
up_mean_frame=[8000,7600,7600,7800,7700,7700,7700,7500,7800,8000]
up_mono_chi=[0.035,0.035,0.035,0.035,0.035,0.035,0.035,0.035,0.035,0.035]
up_mono_tau=[28,28,28,28,28,28,23,23,16,20]

for p in range(0,len(path)):  
    ################# Unfiltered groupprops to list
    groupprops=fifo.read_locs(path[p])# Single groupprops
    groupprops_list.append(groupprops)# List containing all groupprops arrays (unfiltered)

    ################# Set filters
    # Set single filters, always one single filter has to be active!
    istrue_single=[]
    istrue_single.append(groupprops_list[p]['std_frame']>low_std_frame[p])
    istrue_single.append(groupprops_list[p]['mean_frame']>low_mean_frame[p])
    istrue_single.append(groupprops_list[p]['mean_frame']<up_mean_frame[p])
    istrue_single.append(groupprops_list[p]['mono_chi']<up_mono_chi[p])
    istrue_single.append(groupprops_list[p]['mono_tau']<up_mono_tau[p])
    # Combine single filters to overall filter (and not)
    istrue=np.all(istrue_single,axis=0)
    istrue_not=~np.array(istrue)
    
    ################# Filtered groupprops to list
    groupprops_list_filter.append(groupprops_list[p][:][istrue]) # Create list of filtered groupprops
    groupprops_list_filter_out.append(groupprops_list[p][:][istrue_not]) # Create list of out-filtered groupprops
    
    ################# Mean and std of all fields in filtered groupprops 
    for name in groupprops.dtype.names: 
        means[name][p]=np.nanmean(groupprops_list_filter[p][name])
        stds[name][p]=np.nanstd(groupprops_list_filter[p][name])#/np.sqrt(len(groupprops_list_filter[p])) #Uncomment to switch between std and standard error of the mean

################################################### Show filtering   
p=9
field='mono_tau'

#f=plt.figure(num=10)
#f.clear()
#ax=f.add_subplot(1,1,1)
#ax.hist(groupprops_list_filter[p][field][:],bins='fd',label='p=%.0i, '%(p)+field)
#ax.legend(loc=1)

################################################### Radial binning
# Number of bins
NoBins=5
r0=[512,512]

n_groups_r_bin=np.empty([len(groupprops_list_filter),NoBins,1])
mono_tau_r_bin=np.empty([len(groupprops_list_filter),NoBins,2])
tau_b_r_bin=np.empty([len(groupprops_list_filter),NoBins,2])
tau_d_r_bin=np.empty([len(groupprops_list_filter),NoBins,2])


for p in range(0,len(groupprops_list)):
    # Get distance from center of illumination
#    r=np.sqrt(np.power(groupprops_list_filter[p]['mean_x']-r0[0],2)+np.power(groupprops_list_filter[p]['mean_y']-r0[1],2))
    r=(np.power(groupprops_list_filter[p]['mean_x']-r0[0],2)+np.power(groupprops_list_filter[p]['mean_y']-r0[1],2))
    # Bin Intensity in NoBins
    bins_r=np.linspace(0,470**2,NoBins+1)
#    bins_r=np.array([100,150,200,300,400,500])
    # Get groups in r according to binning
    digitized_r=np.digitize(r,bins_r)
    # Get number of groups in each bin 
    n_groups_r_bin[p,:,0] = np.array([len(digitized_r[digitized_r == i]) for i in range(1, len(bins_r))])
    # Get mean&std of mono_tau in each r-bin
    mono_tau_r_bin[p,:,0] = np.array([groupprops_list_filter[p]['mono_tau'][digitized_r == i].mean() for i in range(1, len(bins_r))])
    mono_tau_r_bin[p,:,1] = np.array([groupprops_list_filter[p]['mono_tau'][digitized_r == i].std() for i in range(1, len(bins_r))])
    # Get mean&std of tau_b_lin_ignore in each r-bin
    tau_b_r_bin[p,:,0] = np.array([groupprops_list_filter[p]['tau_b_ignore'][digitized_r == i].mean() for i in range(1, len(bins_r))])
    tau_b_r_bin[p,:,1] = np.array([groupprops_list_filter[p]['tau_b_ignore'][digitized_r == i].std() for i in range(1, len(bins_r))])
    # Get mean&std of tau_d_lin_ignore in each r-bin
    tau_d_r_bin[p,:,0] = np.array([groupprops_list_filter[p]['n_events'][digitized_r == i].mean() for i in range(1, len(bins_r))])
    tau_d_r_bin[p,:,1] = np.array([groupprops_list_filter[p]['n_events'][digitized_r == i].std() for i in range(1, len(bins_r))])      
    

################################################### Plotting
pltrange=[0,2,4,5,6,7,8,9]    
################################################### mono_tau
f=plt.figure(num=1)
f.clear()
f.subplots_adjust(bottom=0.1,left=0.13)
ax=f.add_subplot(1,1,1)


for p in pltrange:
    x=np.sqrt(bins_r[1:])
    y=mono_tau_r_bin[p,:,0]*CycleTime
    yerr=np.divide(mono_tau_r_bin[p,:,1],np.sqrt(n_groups_r_bin[p,:,0]))*CycleTime
    ax.errorbar(x,y,yerr=yerr,label=labels[p])

ax.set_xlabel(r'$r-r_{0}$ [px]')
ax.set_ylabel(r'$\tau_c$ [s]')
ax.legend(loc=4)

################################################### tau_b_lin_ignore
f=plt.figure(num=2)
f.clear()
f.subplots_adjust(bottom=0.1,left=0.13)
ax=f.add_subplot(1,1,1)


for p in pltrange:
    x=np.sqrt(bins_r[1:])
    y=tau_b_r_bin[p,:,0]*CycleTime
    yerr=np.divide(tau_b_r_bin[p,:,1],np.sqrt(n_groups_r_bin[p,:,0]))*CycleTime
    ax.errorbar(x,y,yerr=yerr,label=labels[p])

ax.set_xlabel(r'$r-r_{0}$ [px]')
ax.set_ylabel(r'$\tau_b$ [s]')
ax.legend(loc=1)

################################################### tau_b_lin_ignore
pltrange=[4,5]  
f=plt.figure(num=3)
f.clear()
f.subplots_adjust(bottom=0.1,left=0.13)
ax=f.add_subplot(1,1,1)


for p in pltrange:
    x=np.sqrt(bins_r[1:])
    y=tau_d_r_bin[p,:,0]#*CycleTime
    yerr=np.divide(tau_d_r_bin[p,:,1],np.sqrt(n_groups_r_bin[p,:,0]))#*CycleTime
    ax.errorbar(x,y,yerr=yerr,label=labels[p])

ax.set_xlabel(r'$r-r_{0}$ [px]')
ax.set_ylabel(r'$\tau_d$ [s]')
ax.legend(loc=1)