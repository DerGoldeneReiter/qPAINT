# Load packages
import h5py as h5py #hdf5 handling
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import matplotlib as mpl
import os #platform independent paths
import importlib
# Load user defined functions
import locs_groupprops as l2grp
import file_formats as fifo
import fitfunc
importlib.reload(fifo)
importlib.reload(l2grp)
importlib.reload(fitfunc)

#single docks
dir_name=[r'E:\Projects\qPAINT\data\18-01-12\sample01_p40_convampf_3x_T25_1\JS_18-01-16\ng=400_single'] 
file_name=['sample01_p40_convampf_3x_T25_1_MMStack_Pos0.ome_locs_render_picked_groupprops.hdf5']

#N docks
dir_name.extend([r'E:\Projects\qPAINT\data\18-01-12\sample01_p40_convampf_3x_T25_1\JS_18-01-16\ng=400_N'])
file_name.extend(['sample01_p40_convampf_3x_T25_1_MMStack_Pos0.ome_locs_picked_groupprops.hdf5'])


power=np.array([7.2,7.2])*1e-3 # Total power [W]
peak_I=(2*power)/(np.pi*(120*16e-6*1e-2)**2)*1e-7 # peak Intensity [kW/cm^2]

label=[r'$I_p$'+' = %.2f'%((2*p)/(np.pi*(120*16e-6*1e-2)**2)*1e-7)+r' $\frac{kW}{cm^2}$' for p in power]

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
meanframe_mean_list=[]
meanframe_std_list=[]
stdframe_mean_list=[]
stdframe_std_list=[]

for p in range(0,len(groupprops_list)):
    # Mean of mean_frame
    meanframe_mean=np.mean(groupprops_list[p]['mean_frame'])
    meanframe_mean_list.append(meanframe_mean)
    # Std of mean_frame
    meanframe_std=np.std(groupprops_list[p]['mean_frame'])
    meanframe_std_list.append(meanframe_std)
    
    # Mean of std_frame
    stdframe_mean=np.mean(groupprops_list[p]['std_frame'])
    stdframe_mean_list.append(stdframe_mean)
    # Std of std_frame
    stdframe_std=np.std(groupprops_list[p]['std_frame'])
    stdframe_std_list.append(stdframe_std)
    
    # Define filters
    istrue_meanframe=np.absolute(groupprops_list[p]['mean_frame']-meanframe_mean)<meanframe_std
    istrue_stdframe=np.absolute(groupprops_list[p]['std_frame']-stdframe_mean)<stdframe_std
    istrue_tau_out=groupprops_list[p]['ac_tau']<100
    
    # Overall filter
    istrue=istrue_meanframe&istrue_stdframe&istrue_tau_out
    
    groupprops_list_filter.append(groupprops_list[p][:][istrue])

############################################# Tau vs power calculation
x0=268.75
y0=230.65    
sigma_x=118.53
sigma_y=136.23
theta=-6.09

NoBins=20
I_bin=np.empty([len(label),NoBins,2])
tau_c_I_bin=np.empty([len(label),NoBins,2])
p_inf_I_bin=np.empty([len(label),NoBins,2])

r_bin=np.empty([len(label),NoBins,2])
tau_c_r_bin=np.empty([len(label),NoBins,2])
p_inf_r_bin=np.empty([len(label),NoBins,2])
N_r_bin=np.empty([len(label),NoBins,2])

for p in range(0,len(path)):
    # Get tau_c for each point
    tau_c=groupprops_list_filter[p]['ac_tau']
    # Get p_inf for each point
    p_inf=groupprops_list_filter[p]['ac_p_inf']
    # Get p_inf for each point
    N=groupprops_list_filter[p]['NoDocks']
    # Radial distance from center of illumination
    r=np.sqrt(np.power(groupprops_list_filter[p]['mean_x']-x0,2)+np.power(groupprops_list_filter[p]['mean_y']-y0,2))
    # Get intensity for each (x,y) value according to measurement
    I=fitfunc.gaussian_2D((groupprops_list_filter[p]['mean_x'],groupprops_list_filter[p]['mean_y']), peak_I[p], x0, y0, sigma_x, sigma_y, theta,0)
    
    # Bin I 
    bins_I=np.linspace(np.min(I),np.max(I),NoBins+1)
    digitized_I=np.digitize(I,bins_I)
    bins_r=np.linspace(0,np.max(r),NoBins+1)
    digitized_r=np.digitize(r,bins_r)
    
    # Get binned r & tau
    # Autocorr
    I_bin[p,:,0] = np.array([I[digitized_I == i].mean() for i in range(1, len(bins_I))])
    I_bin[p,:,1] = np.array([I[digitized_I == i].std() for i in range(1, len(bins_I))])
    r_bin[p,:,0] = np.array([r[digitized_r == i].mean() for i in range(1, len(bins_r))])
    r_bin[p,:,1] = np.array([r[digitized_r == i].std() for i in range(1, len(bins_r))])
    
    tau_c_I_bin[p,:,0] = np.array([tau_c[digitized_I == i].mean() for i in range(1, len(bins_I))])
    tau_c_I_bin[p,:,1] = np.array([tau_c[digitized_I == i].std() for i in range(1, len(bins_I))])
    p_inf_I_bin[p,:,0] = np.array([p_inf[digitized_I == i].mean() for i in range(1, len(bins_I))])
    p_inf_I_bin[p,:,1] = np.array([p_inf[digitized_I == i].std() for i in range(1, len(bins_I))])
    
    tau_c_r_bin[p,:,0] = np.array([tau_c[digitized_r == i].mean() for i in range(1, len(bins_r))])
    tau_c_r_bin[p,:,1] = np.array([tau_c[digitized_r == i].std() for i in range(1, len(bins_r))])
    p_inf_r_bin[p,:,0] = np.array([p_inf[digitized_r == i].mean() for i in range(1, len(bins_r))])
    p_inf_r_bin[p,:,1] = np.array([p_inf[digitized_r == i].std() for i in range(1, len(bins_r))])
    N_r_bin[p,:,0] = np.array([N[digitized_r == i].mean() for i in range(1, len(bins_r))])
    N_r_bin[p,:,1] = np.array([N[digitized_r == i].std() for i in range(1, len(bins_r))])
#%% Show frame distribution
plt.clf()
f,ax=plt.subplots(2,2,num=1)
f.suptitle('8nt vs T: Unfiltered frames')
for p in range(0,len(groupprops_list)):
    ax[0,p].hist(groupprops_list[p]['mean_frame'][:],bins='fd')
    ax[1,p].hist(groupprops_list[p]['std_frame'][:],bins='fd')
    
#%% Show filtered frame distribution    
plt.clf()    
f,ax=plt.subplots(2,2,num=1)
f.suptitle('8nt vs T: Filtered frames')
for p in range(0,len(groupprops_list)):
    ax[0,p].hist(groupprops_list_filter[p]['mean_frame'][:],bins='fd')
    ax[1,p].hist(groupprops_list_filter[p]['std_frame'][:],bins='fd') 
    
#%% ################################################
#################################################### Bright times    
# Show filtered radial tau_c
plt.clf()    
figsize=[np.sqrt(2)*7,7]    
f,ax=plt.subplots(1,1,num=1,figsize=figsize)
f.subplots_adjust(left=0.08,bottom=0.12,right=0.72)
ax.set_xscale('log')
ax.set_xlabel('Intensity '+r'$[\frac{kW}{cm^2}]$')
ax.set_ylabel(r'$\tau_c$'+' [s]')
#ax.set_xlim([1e-2,4])
#ax.set_ylim([0.1,2.42])

pltrange=[range(1,NoBins)]*len(path)
for p in range(0,len(path)):
    ax.errorbar(I_bin[p,pltrange[p],0],tau_c_I_bin[p,pltrange[p],0],yerr=tau_c_I_bin[p,pltrange[p],1],label=label[p])
    
ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)

#%% 
# Show filtered radial p_inf
plt.clf()
figsize=[np.sqrt(2)*7,7]    
f,ax=plt.subplots(1,1,num=1,figsize=figsize)
f.subplots_adjust(left=0.08,bottom=0.12,right=0.72)
#ax.set_xscale('log')
ax.set_xlabel('Intensity '+r'$[\frac{kW}{cm^2}]$')
#ax.set_ylabel(r'$\tau_b$'+' [s]')
#ax.set_xlim([1e-2,4])
#ax.set_ylim([0.1,2.42])

pltrange=[range(0,NoBins)]*len(path)
for p in range(0,len(path)):
    ax.errorbar(I_bin[p,pltrange[p],0],p_inf_I_bin[p,pltrange[p],0],yerr=p_inf_I_bin[p,pltrange[p],1],label=label[p])
    
ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)

#%% 
# Show filtered p_inf vs r
plt.clf()
figsize=[np.sqrt(2)*7,7]    
f,ax=plt.subplots(1,1,num=1,figsize=figsize)
f.subplots_adjust(left=0.08,bottom=0.12,right=0.72)
#ax.set_xscale('log')
ax.set_xlabel('r [px]')
#ax.set_ylabel(r'$\tau_b$'+' [s]')
#ax.set_xlim([1e-2,4])
#ax.set_ylim([0.1,2.42])

pltrange=[range(0,NoBins)]*len(path)
for p in range(0,len(path)):
    ax.errorbar(r_bin[p,pltrange[p],0],p_inf_r_bin[p,pltrange[p],0],yerr=p_inf_r_bin[p,pltrange[p],1],label=label[p])
    
ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)

#%% 
# Show filtered p_inf vs I
plt.clf()
figsize=[np.sqrt(2)*7,7]    
f,ax=plt.subplots(1,1,num=1,figsize=figsize)
f.subplots_adjust(left=0.08,bottom=0.12,right=0.72)

ax.errorbar(I_bin[0,:,0],np.divide(p_inf_I_bin[1,:,0],p_inf_I_bin[0,:,0]),label=label[p])

#%% 
# Show filtered p_inf_N/p_inf_1
plt.clf()
figsize=[np.sqrt(2)*7,7]    
f,ax=plt.subplots(1,1,num=1,figsize=figsize)
f.subplots_adjust(left=0.08,bottom=0.12,right=0.72)

ax.errorbar(r_bin[0,:,0],np.divide(p_inf_r_bin[1,:,0],p_inf_r_bin[0,:,0]),label=label[p])

#%% 
# Show filtered NoDocks
plt.clf()
figsize=[np.sqrt(2)*7,7]    
f,ax=plt.subplots(1,1,num=1,figsize=figsize)
f.subplots_adjust(left=0.08,bottom=0.12,right=0.72)

ax.errorbar(r_bin[1,:,0],N_r_bin[1,:,0],yerr=N_r_bin[1,:,1],label=label[p])

#%%
# Get all groups with r<r_lim
r_lim=80
istrue_center=np.sqrt(np.power(groupprops_list_filter[p]['mean_x']-x0,2)+np.power(groupprops_list_filter[p]['mean_y']-y0,2))<r_lim
groupprops=groupprops_list_filter[1][:][istrue_center]

groupprops_filter=np.zeros(np.size(groupprops),dtype=[('group','i4',1),
                                ('r','f4'),# Location
                                ('NoDocks','f4',1), # ac dynamics
                                ('mean_photons','f4',1) # Photons
                                ])
groupprops_filter['group']=groupprops['group']
groupprops_filter['r']=np.sqrt(np.power(groupprops['mean_x']-x0,2)+np.power(groupprops['mean_y']-y0,2))
groupprops_filter['NoDocks']=groupprops['NoDocks']
groupprops_filter['mean_photons']=groupprops['mean_photons']

#%%
# Show filtered NoDocks
fontsize=20
mpl.rcParams['figure.figsize']=[16,8]
mpl.rcParams['axes.linewidth']=1
mpl.rcParams['lines.linewidth']=3
mpl.rcParams['lines.markersize']=8
mpl.rcParams['errorbar.capsize']=6
mpl.rcParams['xtick.major.width']=1


plt.clf()
f,ax=plt.subplots(1,2,num=1,gridspec_kw = {'width_ratios':[2, 4]})  
#f.set_size_inches(np.sqrt(2)*7, 7, forward=True)
f.subplots_adjust(left=0.05,bottom=0.09,right=0.98,wspace=0.3)


ax[0].hist(groupprops_filter['NoDocks'],bins='fd')
ax[0].set_title('Center of illumination [r < 80 px]\nNumber of docking sites',fontsize=fontsize)
ax[0].set_ylabel('Counts',fontsize=fontsize)
ax[0].set_xlabel('N',fontsize=fontsize)
ax[0].set_xlim([5,16])
ax[0].set_xticks(np.arange(5,16,2))
ax[0].text(11,34,'N = %.1f \n      (%.1f)'%(np.mean(groupprops_filter['NoDocks']),np.std(groupprops_filter['NoDocks'])),fontsize=fontsize)

label=[r'$p_{\infty,1}$',r'$p_{\infty,Origami}$']
for p in range(0,len(path)):
    ax[1].errorbar(r_bin[p,pltrange[p],0],p_inf_r_bin[p,pltrange[p],0],yerr=p_inf_r_bin[p,pltrange[p],1],label=label[p])
ax[1].set_title('Radial dependence of '+r'$p_{\infty}$'+' from center of illumination',fontsize=fontsize)
ax[1].set_ylabel(r'$p_{\infty}$',fontsize=fontsize+5)
ax[1].set_xlabel('Distance r from center of illumination [px]',fontsize=fontsize)
ax[1].legend(loc=0,fontsize=fontsize-2)

#%%
groupprops_file = h5py.File(path[1].replace('.hdf5','_filter.hdf5'), "w")
dset=groupprops_file.create_dataset("locs", np.shape(groupprops_filter), dtype=groupprops_filter.dtype)
dset[...]=groupprops_filter
groupprops_file.close() 