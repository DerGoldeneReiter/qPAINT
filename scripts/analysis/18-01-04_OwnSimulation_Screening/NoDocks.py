# Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os #platform independent paths
import matplotlib as mpl
# Load user defined functions
import file_formats as fifo
# Load style sheet
plt.style.use(r'E:\Flo\repos\qPAINT\styles\FoM.mplstyle')

# Define folder of locs.hdf5 file
dir_name=[r'E:\Projects\qPAINT\data\18-01-03_OwnSimulation_4ds\20frames']*8+\
[r'E:\Projects\qPAINT\data\18-01-03_OwnSimulation_12ds\20frames']*8


# Define names of locs.hdf5 file
file_names=['01_locs_picked_groupprops.hdf5']
file_names.extend(['02_locs_picked_groupprops.hdf5'])
file_names.extend(['03_locs_picked_groupprops.hdf5'])
file_names.extend(['04_locs_picked_groupprops.hdf5'])
file_names.extend(['05_locs_picked_groupprops.hdf5'])
file_names.extend(['06_locs_picked_groupprops.hdf5'])
file_names.extend(['07_locs_picked_groupprops.hdf5'])
file_names.extend(['08_locs_picked_groupprops.hdf5'])
file_names=file_names+file_names
# Define labels
labels=['0.5nM']
labels.extend(['2.5nM'])
labels.extend(['5nM'])
labels.extend(['10nM'])
labels.extend(['25nM'])
labels.extend(['50nM'])
labels.extend(['75nM'])
labels.extend(['100nM'])
labels=labels+labels
# Create list of paths
path=list()
for i in range(0,len(file_names)):
    path.append(os.path.join(dir_name[i],file_names[i]))
    
# Overall simulation parameters
concs=np.array([float(label.replace('nM',''))*1e-9 for label in labels]) # List of concentrations as floats
tau_bs=np.array([20.]*20) # Bright time [Frames]
tau_ds=[10000.,2000.,1000.,500.,200,100,(2/3)*100,50]
tau_ds=np.array(tau_ds+tau_ds+tau_ds[0:4])
 # Dark time [Frames]
NoFrames=[100000]*20
NoDocks=np.array([4]*8+[12]*8+[40]*4)
# Calculate expected tau_c
tau_cs=np.divide(1,np.power(tau_bs,-1)+np.power(tau_ds,-1)) # Autocorr tau [s]
# Calculate expected p_inf
p_infs=NoDocks*np.divide(tau_bs,tau_ds)


# Initialize mean&std of important parameters
ac_NoDocks=np.empty([len(path),2])
for i in range(0,len(path)):
    ac_NoDocks[i,0]=np.mean(fifo.read_groupprops(path[i])['NoDocks'])
    ac_NoDocks[i,1]=np.std(fifo.read_groupprops(path[i])['NoDocks'])/10

#%%
f=plt.figure(num=1)
f.clear()
f.subplots_adjust(hspace=0.3)
# N=4
ax=f.add_subplot(2,1,1)
pltrange=np.arange(1,8) 
ax.errorbar(concs[pltrange]*1e9,ac_NoDocks[pltrange,0],yerr=ac_NoDocks[pltrange,1],fmt='o',label=r'$N_{in}=4$',color='blue')
ax.axhline(NoDocks[pltrange[0]],color='blue')
ax.set_xscale('log')
ax.set_ylabel(r'$N_{out}$')
ax.legend(loc=1)

# N=12
ax=f.add_subplot(2,1,2)
pltrange=np.arange(9,16) 
ax.errorbar(concs[pltrange]*1e9,ac_NoDocks[pltrange,0],yerr=ac_NoDocks[pltrange,1],fmt='o',label=r'$N_{in}=12$',color='red')
ax.axhline(NoDocks[pltrange[0]],color='red')
ax.set_xscale('log')
ax.set_xlabel('Concentration [nM]')
ax.set_ylabel(r'$N_{out}$')
ax.legend(loc=1)


