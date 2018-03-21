# Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os #platform independent paths
import importlib
# Load user defined functions
import locs_groupprops as l2grp
import file_formats as fifo
import fitfunc
importlib.reload(fifo)
importlib.reload(l2grp)
importlib.reload(fitfunc)
#Set plot style
plt.style.use(r'E:\Flo\repos\qPAINT\styles\FoM.mplstyle')

#T=22
dir_name=[r'E:\Projects\qPAINT\data\18-01-12\sample01_pseries10_convampf_3x_T22_1\JS_18-01-16\ng=200'] #P=10
dir_name.extend([r'E:\Projects\qPAINT\data\18-01-12\sample01_pseries20_2_convampf_3x_T22_1\JS_18-01-12\ng=350']) #P=20-2
dir_name.extend([r'E:\Projects\qPAINT\data\18-01-12\sample01_pseries30_convampf_3x_T22_1\JS_18-01-16\ng=400']) #P=30
dir_name.extend([r'E:\Projects\qPAINT\data\18-01-12\sample01_pseries40_convampf_3x_T22_1\JS_18-01-16\ng=400']) #P=40-2
dir_name.extend([r'E:\Projects\qPAINT\data\18-01-12\sample01_pseries100_convampf_3x_T22_1\JS_18-01-16\ng=400']) #P=100
#T=22
file_name=['sample01_pseries10_convampf_3x_T22_1_MMStack_Pos0.ome_locs_picked_groupprops.hdf5'] #P=10
file_name.extend(['sample01_pseries20_2_convampf_3x_T22_1_MMStack_Pos0.ome_locs_picked_groupprops.hdf5']) #P=20-2
file_name.extend(['sample01_pseries30_convampf_3x_T22_1_MMStack_Pos0.ome_locs_picked_groupprops.hdf5']) #P=30
file_name.extend(['sample01_pseries40_convampf_3x_T22_1_MMStack_Pos0.ome_locs_picked_groupprops.hdf5']) #P=40-2
file_name.extend(['sample01_pseries100_convampf_3x_T22_1_MMStack_Pos0.ome_locs_picked_groupprops.hdf5']) #P=100
#T=25
dir_name.extend([r'E:\Projects\qPAINT\data\18-01-12\sample01_p40_convampf_3x_T25_1\JS_18-01-16\ng=400']) # P=40
file_name.extend(['sample01_p40_convampf_3x_T25_1_MMStack_Pos0.ome_locs_picked_groupprops.hdf5'])
#T=32
dir_name.extend([r'E:\Projects\qPAINT\data\18-01-12\sample01_p40_convampf_3x_T32_1\JS_18-01-16\ng=400']) # P=40 
file_name.extend(['sample01_p40_convampf_3x_T32_1_MMStack_Pos0.ome_locs_picked_groupprops.hdf5'])

# Define labels
labels=['0.5mW']
labels.extend(['2mW'])
labels.extend(['4mW'])
labels.extend(['7mW'])
labels.extend(['19mW'])
labels.extend(['7mW'])
labels.extend(['7mW'])

CycleTime=0.1

# Create list of paths
path=[]
for i in range(0,len(file_name)):
    path.extend([os.path.join(dir_name[i],file_name[i])])

groupprops_list=[]
for i in range(0,len(path)):
    groupprops_list.append(fifo.read_groupprops(path[i]))
    

################################################### Filter:       
groupprops_list_filter=[]
groupprops_list_filter_out=[]

for p in range(0,len(groupprops_list)):
    # Single Filters
    istrue_ac_chi=groupprops_list[p]['ac_chisquare']<0.03
    istrue_ac_tau=groupprops_list[p]['ac_tau']<50
    # Overall filter
    istrue=istrue_ac_chi&istrue_ac_tau
    istrue_not=~np.array(istrue)
    
    groupprops_list_filter.append(groupprops_list[p][:][istrue])
    groupprops_list_filter_out.append(groupprops_list[p][:][istrue_not])

################################################### Show filtering   
file=0
field='ac_tau'

#f=plt.figure(num=10)
#f.clear()
#f.suptitle(field)
#ax=f.add_subplot(1,2,1)
#ax.hist(groupprops_list[file][field][:],bins='fd')
#ax.set_title('Unfiltered')
#ax=f.add_subplot(1,2,2)
#ax.hist(groupprops_list_filter[file][field][:],bins='fd')
#ax.set_title('Filtered')
#plt.show()


################################################### Radial& Intensity binning
power=np.array([0.45,1.9,4.2,7.2,19,7.2,7.2])*1e-3 # Total measured power after objective for every measurement [W] 
popt_I=[1,268.935,230.337,192.413,167.164,45.73,0] # Gaussian beam parameters (same for every measuremnt)

NoBins=13
Sample_I_bin=np.empty([len(groupprops_list),NoBins,1])
I_bin=np.empty([len(groupprops_list),NoBins,2])
tau_c_I_bin=np.empty([len(groupprops_list),NoBins,2])
p_inf_I_bin=np.empty([len(groupprops_list),NoBins,2])
ac_tau_b_I_bin=np.empty([len(groupprops_list),NoBins,2])
ac_tau_d_I_bin=np.empty([len(groupprops_list),NoBins,2])

tau_b_I_bin=np.empty([len(groupprops_list),NoBins,2])
tau_d_I_bin=np.empty([len(groupprops_list),NoBins,2])

for p in range(0,len(groupprops_list)):
    # Update peak intensity for every measurement
    popt_I[0]=(2*power[p])/(np.pi*(np.mean(popt_I[3:4])*16e-6*1e-2)**2)*1e-7 #[kW/cm^2]
    # Get intensity for every group
    I=fitfunc.gaussian_2D((groupprops_list_filter[p]['mean_x'],groupprops_list_filter[p]['mean_y']),*popt_I)
    # Bin Intensity in NoBins
    bins_I=np.linspace(np.min(I),np.max(I),NoBins+1)
    # Get groups in I according to binning
    digitized_I=np.digitize(I,bins_I)
    
    Sample_I_bin[p,:,0] = np.array([len(digitized_I[digitized_I == i]) for i in range(1,NoBins+1)])
    I_bin[p,:,0] = np.array([I[digitized_I == i].mean() for i in range(1, len(bins_I))])
    I_bin[p,:,1] = np.array([I[digitized_I == i].std() for i in range(1, len(bins_I))])
    
    tau_c_I_bin[p,:,0] = np.array([groupprops_list_filter[p]['ac_tau'][digitized_I == i].mean() for i in range(1, len(bins_I))])
    tau_c_I_bin[p,:,1] = np.array([groupprops_list_filter[p]['ac_tau'][digitized_I == i].std() for i in range(1, len(bins_I))])
    p_inf_I_bin[p,:,0] = np.array([groupprops_list_filter[p]['ac_p_inf'][digitized_I == i].mean() for i in range(1, len(bins_I))])
    p_inf_I_bin[p,:,1] = np.array([groupprops_list_filter[p]['ac_p_inf'][digitized_I == i].std() for i in range(1, len(bins_I))])    
    ac_tau_b_I_bin[p,:,0] = np.array([groupprops_list_filter[p]['ac_tau_b'][digitized_I == i].mean() for i in range(1, len(bins_I))])
    ac_tau_b_I_bin[p,:,1] = np.array([groupprops_list_filter[p]['ac_tau_b'][digitized_I == i].std() for i in range(1, len(bins_I))])
    ac_tau_d_I_bin[p,:,0] = np.array([groupprops_list_filter[p]['ac_tau_d'][digitized_I == i].mean() for i in range(1, len(bins_I))])
    ac_tau_d_I_bin[p,:,1] = np.array([groupprops_list_filter[p]['ac_tau_d'][digitized_I == i].std() for i in range(1, len(bins_I))])        
    
    tau_b_I_bin[p,:,0] = np.array([groupprops_list_filter[p]['tau_b_lin_ignore'][digitized_I == i].mean() for i in range(1, len(bins_I))])
    tau_b_I_bin[p,:,1] = np.array([groupprops_list_filter[p]['tau_b_lin_ignore'][digitized_I == i].std() for i in range(1, len(bins_I))])
    tau_d_I_bin[p,:,0] = np.array([groupprops_list_filter[p]['tau_d_lin'][digitized_I == i].mean() for i in range(1, len(bins_I))])
    tau_d_I_bin[p,:,1] = np.array([groupprops_list_filter[p]['tau_d_lin'][digitized_I == i].std() for i in range(1, len(bins_I))])    

    
#%% ################################################
ylim=[0.6,2.28]
################################################### Autocorrelation: tau
f=plt.figure(num=1)
f.clear()
ax=f.add_subplot(1,1,1)
for p in range(1, len(groupprops_list)):
    ax.errorbar(I_bin[p,:,0],tau_c_I_bin[p,:,0]*CycleTime,yerr=np.divide(tau_c_I_bin[p,:,1],np.sqrt(Sample_I_bin[p,:,0]))*CycleTime,label=labels[p])

ax.set_xscale('log')
ax.set_ylim(ylim)

ax.set_xlabel(r'Intensity $[kW/cm^2]$')
ax.set_ylabel(r'Autocorrelation $\tau_c$ [s]')
ax.legend(loc=3)

################################################### Picasso: BRIGHT
f=plt.figure(num=2)
f.clear()
ax=f.add_subplot(1,1,1)
for p in range(1, len(groupprops_list)):
    ax.errorbar(I_bin[p,:,0],tau_b_I_bin[p,:,0]*CycleTime,yerr=np.divide(tau_b_I_bin[p,:,1],np.sqrt(Sample_I_bin[p,:,0]))*CycleTime,label=labels[p])

ax.set_xscale('log')
ax.set_ylim(ylim)

ax.set_xlabel(r'Intensity $[kW/cm^2]$')
ax.set_ylabel(r'Picasso $\langle\tau_b\rangle$ [s]')
ax.legend(loc=3)

################################################### Temperature dependence
f=plt.figure(num=3)
f.clear()
ax=f.add_subplot(1,1,1)
select=[3,5,6]
ax.errorbar([22,25,32],tau_c_I_bin[select,0,0]*CycleTime,yerr=np.divide(tau_c_I_bin[select,0,1],np.sqrt(Sample_I_bin[select,0,0]))*CycleTime,label=r'$\tau_c$')
ax.errorbar([22,25,32],tau_b_I_bin[select,0,0]*CycleTime,yerr=np.divide(tau_b_I_bin[select,0,1],np.sqrt(Sample_I_bin[select,0,0]))*CycleTime,label=r'$\langle\tau_b\rangle$')

ax.set_xlabel(r'Set Temperature $[C^{\circ}]$')
ax.set_ylabel(r'$\tau_c$ and $\langle\tau_b\rangle$ [s]')
ax.legend(loc=1)

################################################### Autocorrelation: tau
f=plt.figure(num=4)
f.clear()
ax=f.add_subplot(1,1,1)
for p in range(1, len(groupprops_list)):
    ax.errorbar(I_bin[p,:,0],p_inf_I_bin[p,:,0],yerr=np.divide(p_inf_I_bin[p,:,1],np.sqrt(Sample_I_bin[p,:,0])),label=labels[p])

ax.set_xscale('log')
ax.set_xlabel(r'Intensity $[kW/cm^2]$')
ax.set_ylabel(r'Autocorrelation $p_{\infty}$')
ax.legend(loc=3)