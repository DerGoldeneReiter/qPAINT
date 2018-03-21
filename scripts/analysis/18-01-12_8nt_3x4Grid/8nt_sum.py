# Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os #platform independent paths
import matplotlib as mpl
import importlib
# Load user defined functions
import locs_groupprops as l2grp
import file_formats as fifo
import fitfunc
importlib.reload(fifo)
importlib.reload(l2grp)
importlib.reload(fitfunc)

#T=22
dir_name=[r'E:\Projects\qPAINT\data\18-01-12\sample01_pseries10_convampf_3x_T22_1\JS_18-01-16\crop4sum_ng=300'] #P=10
dir_name.extend([r'E:\Projects\qPAINT\data\18-01-12\sample01_pseries20_2_convampf_3x_T22_1\JS_18-01-16\crop4sum_ng=1200']) #P=20-2
dir_name.extend([r'E:\Projects\qPAINT\data\18-01-12\sample01_pseries30_convampf_3x_T22_1\JS_18-01-16\crop4sum_ng=750']) #P=30
#dir_name.extend([r'E:\Projects\qPAINT\data\18-01-12\sample01_p40_convampf_3x_1\JS_18-01-16\crop4sum_ng1000']) #P=40 (1st)
dir_name.extend([r'E:\Projects\qPAINT\data\18-01-12\sample01_pseries40_convampf_3x_T22_1\JS_18-01-16\crop4sum_ng1200']) #P=40 (2nd)
dir_name.extend([r'E:\Projects\qPAINT\data\18-01-12\sample01_pseries100_convampf_3x_T22_1\JS_18-01-16\crop4sum_ng=1000']) #P=100
#T=22
file_name=['sample01_pseries10_convampf_3x_T22_1_MMStack_Pos0.ome_locs_sum_groupprops.hdf5'] #P=10
file_name.extend(['sample01_pseries20_2_convampf_3x_T22_1_MMStack_Pos0.ome_locs_sum_groupprops.hdf5']) #P=20-2 
file_name.extend(['sample01_pseries30_convampf_3x_T22_1_MMStack_Pos0.ome_locs_sum_groupprops.hdf5']) #P=30
#file_name.extend(['sample01_p40_convampf_3x_1_MMStack_Pos0.ome_locs_sum_groupprops.hdf5']) #P=40 (1st)
file_name.extend(['sample01_pseries40_convampf_3x_T22_1_MMStack_Pos0.ome_locs_sum_groupprops.hdf5']) #P=40 (2nd)
file_name.extend(['sample01_pseries100_convampf_3x_T22_1_MMStack_Pos0.ome_locs_sum_groupprops.hdf5']) #P=100
##T=25
dir_name.extend([r'E:\Projects\qPAINT\data\18-01-12\sample01_p40_convampf_3x_T25_1\JS_18-01-16\crop4sum_ng=1500']) # P=40
file_name.extend(['sample01_p40_convampf_3x_T25_1_MMStack_Pos0.ome_locs_sum_groupprops.hdf5'])
##T=32
dir_name.extend([r'E:\Projects\qPAINT\data\18-01-12\sample01_p40_convampf_3x_T32_1\JS_18-01-16\crop4sum_ng=1000']) # P=40 
file_name.extend(['sample01_p40_convampf_3x_T32_1_MMStack_Pos0.ome_locs_sum_groupprops.hdf5'])


power=np.array([0.45,1.9,4.2,7.2,19.0,7.2,7.2])*1e-3 # Total power in W
peak_I=(2*power)/(np.pi*(120*16e-6*1e-2)**2)*1e-7

label=[r'$I_p$'+' = %.2f'%((2*p)/(np.pi*(120*16e-6*1e-2)**2)*1e-7)+r' $\frac{kW}{cm^2}$' for p in power]
#label[3]=label[3]+' (1st)'
#label[4]=label[4]+' (2nd)'
label[5]=label[5]+' (25)'
label[6]=label[6]+' (32)'

path=[]
for i in range(0,len(file_name)):
    path.append(os.path.join(dir_name[i],file_name[i]))
    
groupprops_list=[]    
for i in range(0,len(path)):
    groupprops_list.append(fifo.read_grouprops(path[i]))


################################################### Filter:
# 1) Frames    
# 2) tau_c value        
groupprops_list_filter=[]


for p in range(0,len(groupprops_list)):
    istrue_tau_out=groupprops_list[p]['ac_tau']<40
    
    # Overall filter
    istrue=istrue_tau_out
    
    groupprops_list_filter.append(groupprops_list[p][:][istrue])

############################################# Tau vs power calculation
x0=268.75
y0=230.65    
sigma_x=118.53
sigma_y=136.23
theta=-6.09

NoBins=15

I_bin=np.empty([len(label),NoBins,2])
tau_c_I_bin=np.empty([len(label),NoBins,2])

for p in range(0,len(path)):
    # Get tau_c for each point
    tau_c=groupprops_list_filter[p]['ac_tau']
    # Get intensity for each (x,y) value according to measurement
    I=fitfunc.gaussian_2D((groupprops_list_filter[p]['mean_x'],groupprops_list_filter[p]['mean_y']), peak_I[p], x0, y0, sigma_x, sigma_y, theta,0)
    
    # Bin I 
    bins_I=np.linspace(np.min(I),np.max(I),NoBins+1)
    digitized_I=np.digitize(I,bins_I)
    
    # Get binned r & tau
    I_bin[p,:,0] = np.array([I[digitized_I == i].mean() for i in range(1, len(bins_I))])
    I_bin[p,:,1] = np.array([I[digitized_I == i].std() for i in range(1, len(bins_I))])
    tau_c_I_bin[p,:,0] = np.array([tau_c[digitized_I == i].mean() for i in range(1, len(bins_I))])
    tau_c_I_bin[p,:,1] = np.array([tau_c[digitized_I == i].std() for i in range(1, len(bins_I))]) 

#%% 
# Show filtered radial tau_cs for T=22
plt.clf()    
figsize=[np.sqrt(2)*7,7]    
f,ax=plt.subplots(1,1,num=1,figsize=figsize)
f.subplots_adjust(left=0.08,bottom=0.12,right=0.72)
ax.set_xscale('log')
ax.set_xlabel('Intensity '+r'$[\frac{kW}{cm^2}]$')
ax.set_ylabel(r'$\tau_c$'+' [s]')
ax.set_xlim([1e-2,4])
ax.set_ylim([0.1,2.42])

pltrange=[range(0,NoBins)]*8
#pltrange[1]=[range(1,NoBins)]
for p in range(0,len(path)):
    ax.errorbar(I_bin[p,pltrange[p],0],tau_c_I_bin[p,pltrange[p],0]*0.1,label=label[p])
    
ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)

#%%
# Give weighted arithmetic means for temperatures in certain intensity range
# Temperatures
Ts=[22,25,32]
# Define I ranges for all temperatures
I_range=[[7.5e-2,2e-1],[1e-1,4e-1],[1e-1,2]]
# Define range in I_bin corresponding to temperature
T_idx=[[range(0,6)],[6],[7]]


tau_c_T=np.empty([3,2])
for t in range(0,3):
    istrue_T_range=(I_bin[T_idx[t],:,0]>=I_range[t][0])&(I_bin[T_idx[t],:,0]<=I_range[t][1])
    tau_c_T[t,:]=fitfunc.w_mean(tau_c_I_bin[T_idx[t],:,0][istrue_T_range],tau_c_I_bin[T_idx[t],:,1][istrue_T_range])

figsize=[np.sqrt(2)*7,7]    
f,ax=plt.subplots(1,1,num=2,figsize=figsize)
f.subplots_adjust(left=0.08,bottom=0.12,right=0.72)

ax.errorbar(Ts,tau_c_T[:,0],yerr=tau_c_T[:,1])