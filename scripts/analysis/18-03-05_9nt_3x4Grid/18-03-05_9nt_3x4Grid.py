# Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os #platform independent paths
import importlib
import matplotlib as mpl

# Load user defined functions
import locs_groupprops as l2grp
import file_formats as fifo
import fitfunc
importlib.reload(fifo)
importlib.reload(l2grp)
importlib.reload(fitfunc)

#Set plot style
#plt.close('all')
plt.style.use(r'E:\Flo\repos\qPAINT\styles\FoM.mplstyle')

##################################################################################################### Set paths to files and labels
# Define folder of locs.hdf5 file
dir_name=[r'E:\Projects\qPAINT\data\18-02-21\sample1_P4mW_exp100ms_1\JS_18-03-05'] # 5nM
dir_name.extend([r'E:\Projects\qPAINT\data\18-02-21\sample1_P4mW_exp200ms_10nM_1\JS_18-02-21']) # 10nM
dir_name.extend([r'E:\Projects\qPAINT\data\18-02-21\sample1_P4mW_exp100ms_20nM_1\JS_18-02-22']) #20nM
# Define names of locs.hdf5 file
file_names=['sample2_P4mW_exp200ms_1_MMStack_Pos0.ome_locs_picked_groupprops.hdf5'] # 5nM
file_names.extend(['sample1_P4mW_exp200ms_10nM_1_MMStack_Pos0.ome_locs_render_picked_groupprops.hdf5']) # 10nM
file_names.extend(['sample1_P4mW_exp100ms_20nM_1_MMStack_Pos0.ome_locs_render_picked_groupprops.hdf5']) # 20nM
# Define labels
labels=['5nM']
labels.extend(['10nM'])
labels.extend(['20nM'])
CycleTime=[0.1,0.2,0.1]

# Create list of paths
path=[]
for i in range(0,len(file_names)):
    path.extend([os.path.join(dir_name[i],file_names[i])])

groupprops_list=[]
for i in range(0,len(path)):
    groupprops_list.append(fifo.read_groupprops(path[i]))
    

################################################### Filter:       
groupprops_list_filter=[]
groupprops_list_filter_out=[]

for p in range(0,len(groupprops_list)):
    # Single Filters
    istrue_ac_chi=groupprops_list[p]['ac_chisquare']<0.03   
    # Overall filter
    istrue=istrue_ac_chi
    istrue_not=~np.array(istrue)
    
    groupprops_list_filter.append(groupprops_list[p][:][istrue])
    groupprops_list_filter_out.append(groupprops_list[p][:][istrue_not])
    
################################################### Show filtering   
file=0
field='ac_chisquare'

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
power=4e-3 #[W]
popt_I=[1,279.98,244.13,399.99/np.sqrt(2),449.25/np.sqrt(2),43.60,0]
popt_I[0]=(2*power)/(np.pi*(np.mean(popt_I[3:4])*13e-6*1e-2)**2)*1e-7 #[kW/cm^2]

NoBins=12
Sample_I_bin=np.empty([len(groupprops_list),NoBins,1])
I_bin=np.empty([len(groupprops_list),NoBins,2])
tau_c_I_bin=np.empty([len(groupprops_list),NoBins,2])
p_inf_I_bin=np.empty([len(groupprops_list),NoBins,2])
ac_tau_b_I_bin=np.empty([len(groupprops_list),NoBins,2])
ac_tau_d_I_bin=np.empty([len(groupprops_list),NoBins,2])

tau_b_I_bin=np.empty([len(groupprops_list),NoBins,2])
tau_d_I_bin=np.empty([len(groupprops_list),NoBins,2])

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

################################################### Plotting
xlim=[4e-2,0.23]
################################################### Autocorrelation: tau
f=plt.figure(num=1)
f.clear()
ax=f.add_subplot(1,1,1)
for p in range(0, len(groupprops_list)):
    ax.errorbar(I_bin[p,:,0],tau_c_I_bin[p,:,0]*CycleTime[p],yerr=np.divide(tau_c_I_bin[p,:,1],np.sqrt(Sample_I_bin[p,:,0]))*CycleTime[p],label=labels[p])

ax.set_xscale('log')
ax.set_xlim(xlim)

ax.set_xlabel(r'Intensity $[kW/cm^2]$')
ax.set_ylabel(r'Autocorrelation $\tau_c$ [s]')
ax.legend(bbox_to_anchor=(1, 1))

################################################### Autocorrelation: p_inf
f=plt.figure(num=2)
f.clear()
ax=f.add_subplot(1,1,1)
for p in range(0, len(groupprops_list)):
    ax.errorbar(I_bin[p,:,0],p_inf_I_bin[p,:,0],yerr=np.divide(p_inf_I_bin[p,:,1],np.sqrt(Sample_I_bin[p,:,0])),label=labels[p])

ax.set_xscale('log')
ax.set_xlim(xlim)

ax.set_xlabel(r'Intensity $[kW/cm^2]$')
ax.set_ylabel(r'Autocorrelation $p_{\infty}$')
ax.legend(bbox_to_anchor=(1, 1))
#################################################### Picasso: BRIGHT
#f=plt.figure(num=3)
#f.clear()
#ax=f.add_subplot(1,1,1)
#for p in range(0, len(groupprops_list)):
#    ax.errorbar(I_bin[p,:,0],tau_b_I_bin[p,:,0]*CycleTime[p],yerr=np.divide(tau_b_I_bin[p,:,1],np.sqrt(Sample_I_bin[p,:,0]))*CycleTime[p],label=labels[p])
#
#ax.set_xlabel(r'Intensity $[kW/cm^2]$')
#ax.set_ylabel(r'Picasso $\tau_b$ [s]')
#ax.legend(bbox_to_anchor=(1,1))
#
#################################################### Picasso: DARK
#f=plt.figure(num=4)
#f.clear()
#ax=f.add_subplot(1,1,1)
#for p in range(0, len(groupprops_list)):
#    ax.errorbar(I_bin[p,:,0],tau_d_I_bin[p,:,0]*CycleTime[p],yerr=np.divide(tau_d_I_bin[p,:,1],np.sqrt(Sample_I_bin[p,:,0]))*CycleTime[p],label=labels[p])
#
#ax.set_xlabel(r'Intensity $[kW/cm^2]$')
#ax.set_ylabel(r'Picasso $\tau_d$ [s]')
#ax.legend(bbox_to_anchor=(1,1))
#
#################################################### Autocorrelation: BRIGHT
#f=plt.figure(num=5)
#f.clear()
#
#ax=f.add_subplot(1,1,1)
#for p in range(0, len(groupprops_list)):
#    ax.errorbar(I_bin[p,:,0],ac_tau_b_I_bin[p,:,0]*CycleTime[p],yerr=np.divide(ac_tau_b_I_bin[p,:,1],np.sqrt(Sample_I_bin[p,:,0]))*CycleTime[p],label=labels[p])
#
#ax.set_xlabel(r'Intensity $[kW/cm^2]$')
#ax.set_ylabel(r'Autocorrelation $\tau_b$ [s]')
#ax.legend(bbox_to_anchor=(1,1))
#
#################################################### Autocorrelation: DARK
#f=plt.figure(num=6)
#f.clear()
#ax=f.add_subplot(1,1,1)
#for p in range(0, len(groupprops_list)):
#    ax.errorbar(I_bin[p,:,0],ac_tau_d_I_bin[p,:,0]*CycleTime[p],yerr=np.divide(ac_tau_d_I_bin[p,:,1],np.sqrt(Sample_I_bin[p,:,0]))*CycleTime[p],label=labels[p])  
#ax.set_xlabel(r'Intensity $[kW/cm^2]$')
#ax.set_ylabel(r'Autocorrelation $\tau_d$ [s]')
#ax.legend(bbox_to_anchor=(1,1))





































