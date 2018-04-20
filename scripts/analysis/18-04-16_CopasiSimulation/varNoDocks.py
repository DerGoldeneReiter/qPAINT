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
dir_names=['/fs/pool/pool-schwille-paint/Data/Simulation/18-04-16_copasi/']*6
# Define names of locs_picked.hdf5 file
file_names=['taub2s_kon2e6_n1_c10nM_locs_picked_groupprops.hdf5']
file_names.extend(['taub2s_kon2e6_n4_c10nM_locs_picked_groupprops.hdf5'])
file_names.extend(['taub2s_kon2e6_n12_c10nM_locs_picked_groupprops.hdf5'])
file_names.extend(['taub2s_kon2e6_n24_c10nM_locs_picked_groupprops.hdf5'])
file_names.extend(['taub2s_kon2e6_n36_c10nM_locs_picked_groupprops.hdf5'])
file_names.extend(['taub2s_kon2e6_n48_c10nM_locs_picked_groupprops.hdf5'])
# Create list of paths
path=[]
for i in range(0,len(file_names)):
    path.extend([os.path.join(dir_names[i],file_names[i])])
    
#################################################################################################### Simulation input
interval_size=0.1
c=10e-9
kon=2e6
tau_b=20
NoDocks=np.array([1,4,12,24,36,48])
tau_ds=np.divide((1/(kon*c))/interval_size,NoDocks)
NoFrames=[15001]*6
tau_c=np.divide(tau_b*tau_ds[0],tau_b+tau_ds[0])

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

#%%
# tau_c & tau_b vs NoDocks
f=plt.figure(num=1)
f.clear()
ax=f.add_subplot(1,1,1)

pltrange=np.arange(0,6)
field='ac_tau'
label=r'$\tau_c$'
color='blue'
ax.errorbar(NoDocks[pltrange],np.divide(means[field][pltrange]-tau_c,tau_c)*100,yerr=np.divide(stds[field][pltrange],tau_c)*100,fmt='o',label=label,color=color)
field='tau_b'
label=r'$\tau_b$'
color='red'
ax.errorbar(NoDocks[pltrange],np.divide(means[field][pltrange]-tau_b,tau_b)*100,yerr=np.divide(stds[field][pltrange],tau_b)*100,fmt='o',label=label,color=color)
ax.axhline(0,color='black')
# Legend
ax.legend(loc=2)
# Labels
ax.set_xlabel(r'$N_{docks}$')
ax.set_ylabel(r'$\langle\tau\rangle$ relative error [%]')
#%%
# G_0 & tau_d vs NoDocks
f=plt.figure(num=2)
f.clear()
ax=f.add_subplot(1,1,1)

pltrange=np.arange(0,6)
field='ac_A'
label=r'$N_{docks}(G_0)$'
color='blue'
y=np.divide(means[field][0],means[field][pltrange])
yerr=y*np.divide(stds[field][pltrange],means[field][pltrange])
ax.errorbar(NoDocks[pltrange],y,yerr=yerr,fmt='o',label=label,color=color)
field='tau_d'
label=r'$N_{docks}(\tau_d)$'
color='red'
y=np.divide(means[field][0],means[field][pltrange])
yerr=y*np.divide(stds[field][pltrange],means[field][pltrange])
ax.errorbar(NoDocks[pltrange],y,yerr=yerr,fmt='o',label=label,color=color)


# Legend
ax.legend(loc=2)
# Labels
ax.set_xlabel(r'Input $N_{docks}$')
ax.set_ylabel(r'Output $N_{docks}$')