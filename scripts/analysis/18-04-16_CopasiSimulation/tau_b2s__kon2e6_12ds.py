# Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os #platform independent paths
import matplotlib as mpl

# Load user defined functions
import file_formats as fifo

# Style sheet
plt.style.use('/fs/fs01/lv03/home/b_schwille/stehr/qPAINT/styles/FoM.mplstyle')

##################################################################################################### Set paths to files and labels
# Define folder of locs.hdf5 file
dir_names=['/fs/pool/pool-schwille-paint/Data/Simulation/18-04-16_copasi/']*5
# Define names of locs_picked.hdf5 file
file_names=['taub2s_kon2e6_n12_c5nM_locs_picked_groupprops.hdf5']
file_names.extend(['taub2s_kon2e6_n12_c10nM_locs_picked_groupprops.hdf5'])
file_names.extend(['taub2s_kon2e6_n12_c20nM_locs_picked_groupprops.hdf5'])
file_names.extend(['taub2s_kon2e6_n12_c50nM_locs_picked_groupprops.hdf5'])
file_names.extend(['taub2s_kon2e6_n12_c100nM_locs_picked_groupprops.hdf5'])
# Define labels
labels=['5nM']
labels.extend(['10nM'])
labels.extend(['20nM'])
labels.extend(['50nM'])
labels.extend(['100nM'])
# Create list of paths
path=[]
for i in range(0,len(file_names)):
    path.extend([os.path.join(dir_names[i],file_names[i])])
    
#################################################################################################### Simulation input
interval_size=0.1    
concs=np.array([float(label.replace('nM',''))*1e-9 for label in labels]) # List of concentrations as floats
tau_bs=np.array([20.]*5) # Bright time [Frames]
tau_ds=np.divide(1,2e6*concs)/interval_size
NoFrames=[15001]*5
NoDocks=np.array([12]*5)
tau_cs=np.divide(1,np.power(tau_bs,-1)+np.power(tau_ds,-1)) # Autocorr tau [s]

#################################################################################################### Read, Filter, Means&Stds
groupprops_list=[]
groupprops_list_filter=[]
means=np.zeros([len(path),],dtype=fifo.read_groupprops(path[i]).dtype)
stds=np.zeros([len(path),],dtype=fifo.read_groupprops(path[i]).dtype)

for i in range(0,len(path)):  
    ################# Unfiltered groupprops to list
    groupprops=fifo.read_groupprops(path[i])# Single groupprops
    groupprops_list.append(groupprops)# List containing all groupprops arrays (unfiltered)

    ################# Filter
#    groupprops=groupprops[:][groupprops['tau_d_misc']>0] # Remove all entries for which ac_tau_d>NoFrames
    
    ################# Filtered groupprops to list
    groupprops_list_filter.append(groupprops)# List containing all groupprops arrays (unfiltered)
    
    ################# Mean and std of all fields in groupprops   
    for name in groupprops.dtype.names: 
        means[name][i]=np.mean(groupprops[name])
        stds[name][i]=np.std(groupprops[name])#/np.sqrt(len(groupprops))

################################################### Show filtering   
file=2
field='tau_d_misc'

f=plt.figure(num=10)
f.clear()
f.suptitle(field)
ax=f.add_subplot(1,2,1)
ax.hist(groupprops_list[file][field][:],bins='fd')
ax.set_title('Unfiltered')
ax=f.add_subplot(1,2,2)
ax.hist(groupprops_list_filter[file][field][:],bins='fd')
ax.set_title('Filtered')
plt.show()

#%%
# tau_b
f=plt.figure(num=1)
f.clear()
ax=f.add_subplot(1,1,1)

pltrange=np.arange(0,5)
field='ac_tau_b'
label='ACF'
color='blue'
ax.errorbar(concs[pltrange]*1e9,np.divide(means[field][pltrange]-tau_bs[pltrange],tau_bs[pltrange])*100,yerr=np.divide(stds[field][pltrange],tau_bs[pltrange])*100,fmt='o',label=label,color=color)
field='tau_b'
label='ECDF'
color='red'
ax.errorbar(concs[pltrange]*1e9,np.divide(means[field][pltrange]-tau_bs[pltrange],tau_bs[pltrange])*100,yerr=np.divide(stds[field][pltrange],tau_bs[pltrange])*100,fmt='o',label=label,color=color)
ax.axhline(0,color='black')
# Legend
ax.legend(loc=2)
# Scales
ax.set_xscale('log')
#ax.set_yscale('symlog')
# Labels
ax.set_xlabel('Concentration [nM]')
ax.set_ylabel(r'$\langle\tau_b\rangle$ relative error [%]')
#%%
# tau_d
f=plt.figure(num=2)
f.clear()
ax=f.add_subplot(1,1,1)

pltrange=np.arange(0,5)
field='ac_tau_d'
label='ACF'
color='blue'
ax.errorbar(concs[pltrange]*1e9,np.divide(means[field][pltrange]-tau_ds[pltrange],tau_ds[pltrange])*100,yerr=np.divide(stds[field][pltrange],tau_ds[pltrange])*100,fmt='o',label=label,color=color)
field='tau_d'
label='ECDF'
color='red'
ax.errorbar(concs[pltrange]*1e9,np.divide(means[field][pltrange]-tau_ds[pltrange]/NoDocks[0],tau_ds[pltrange]/NoDocks[0])*100,yerr=np.divide(stds[field][pltrange],tau_ds[pltrange]/NoDocks[0])*100,fmt='o',label=label,color=color)
ax.axhline(0,color='black')
# Legend
ax.legend(loc=2)
# Scales
ax.set_xscale('log')
# Labels
ax.set_xlabel('Concentration [nM]')
ax.set_ylabel(r'$\langle\tau_d\rangle$ relative error [%]')

#%%
# n_events
f=plt.figure(num=3)
f.clear()
ax=f.add_subplot(1,1,1)

pltrange=np.arange(0,5)
field='n_events'
label=r'$N_{docks}=12$'
color='black'
ax.errorbar(concs[pltrange]*1e9,means[field][pltrange],yerr=stds[field][pltrange],fmt='o',label=label,color=color)
# Legend
ax.legend(loc=0)
# Scales
ax.set_xscale('log')
# Labels
ax.set_xlabel('Concentration [nM]')
ax.set_ylabel(r'$n_{events}$')

#%%
# NoDocks
f=plt.figure(num=4)
f.clear()
ax=f.add_subplot(1,1,1)

pltrange=np.arange(0,5)
field='NoDocks_misc'
label=r'$N_{docks}=12$'
color='green'
ax.errorbar(concs[pltrange]*1e9,means[field][pltrange],yerr=stds[field][pltrange],fmt='o',label=label,color=color)
ax.axhline(NoDocks[0],color='black')
# Legend
ax.legend(loc=1)
# Scales
ax.set_xscale('log')
# Labels
ax.set_xlabel('Concentration [nM]')
ax.set_ylabel(r'$N_{docks,mix}$')

#%%
# G_0
f=plt.figure(num=5)
f.clear()
ax=f.add_subplot(1,1,1)

pltrange=np.arange(0,5)
field='ac_A'
label=r'$G_{0}$'
color='blue'
ax.errorbar(concs[pltrange]*1e9,means[field][pltrange],yerr=stds[field][pltrange],fmt='o',label=label,color=color)
ax.plot(concs[pltrange]*1e9,np.divide(tau_ds[pltrange],tau_bs[pltrange]*NoDocks[0]),label='Theory',color='black')
# Legend
ax.legend(loc=0)
# Scales
ax.set_xscale('log')
# Labels
ax.set_xlabel('Concentration [nM]')
ax.set_ylabel(r'Autocorrelation Amplitude $G_{0}$')