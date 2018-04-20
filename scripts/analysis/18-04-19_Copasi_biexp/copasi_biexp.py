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
dir_names=['/fs/pool/pool-schwille-paint/Data/Simulation/18-04-19_copasi_biexp/']*2
# Define names of locs_picked.hdf5 file
file_names=['taub1s&taub5s_kon2e6_n12_c10nM_locs_picked_groupprops.hdf5']
file_names.extend(['taub1s&taub5s_kon2e6_n60&12_c10nM_locs_picked_groupprops.hdf5'])
#file_names.extend(['taub2s_kon2e6_n1_c20nM_locs_picked_groupprops.hdf5'])
#file_names.extend(['taub2s_kon2e6_n1_c50nM_locs_picked_groupprops.hdf5'])
#file_names.extend(['taub2s_kon2e6_n1_c100nM_locs_picked_groupprops.hdf5'])

# Create list of paths
path=[]
for i in range(0,len(file_names)):
    path.extend([os.path.join(dir_names[i],file_names[i])])
    
#################################################################################################### Simulation input
# read meta data
[TIFmeta,LOCmeta]=fifo.read_meta(path[0])
# Extract NoFrames from meta-data
NoFrames=TIFmeta['Frames'] # Number of frames in tif stack
# Exposure time
interval_size=0.05

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
#    groupprops=groupprops[:][groupprops['tau2']<300] 
#    groupprops=groupprops[:][groupprops['A1']<0.3]
#    groupprops=groupprops[:][groupprops['tau1']<40]
    
    ################# Filtered groupprops to list
    groupprops_list_filter.append(groupprops)# List containing all groupprops arrays (unfiltered)
    
    ################# Mean and std of all fields in groupprops   
    for name in groupprops.dtype.names: 
        means[name][i]=np.mean(groupprops[name])
        stds[name][i]=np.std(groupprops[name])#/np.sqrt(len(groupprops))

################################################### Show filtering   
file=1
field='tau1'
f=plt.figure(num=1)
f.clear()
ax=f.add_subplot(1,1,1)
ax.hist(groupprops_list_filter[file][field][:],bins='fd',label=field)
plt.legend(loc=1)
plt.show()

field='tau2'
f=plt.figure(num=2)
f.clear()
ax=f.add_subplot(1,1,1)
ax.hist(groupprops_list_filter[file][field][:],bins='fd',label=field)
plt.legend(loc=1)
plt.show()

field='A1'
f=plt.figure(num=3)
f.clear()
ax=f.add_subplot(1,1,1)
ax.hist(groupprops_list_filter[file][field][:],bins='fd',label=field)
plt.legend(loc=1)
plt.show()

field='A2'
f=plt.figure(num=4)
f.clear()
ax=f.add_subplot(1,1,1)
ax.hist(groupprops_list_filter[file][field][:],bins='fd',label=field)
plt.legend(loc=1)
plt.show()

#f=plt.figure(num=5)
#f.clear()
#ax=f.add_subplot(1,1,1)
#ax.hist(np.divide(groupprops_list_filter[file]['A2'][:],groupprops_list_filter[file]['A1'][:]),bins='fd',label='A2/A1')
#plt.legend(loc=1)
#plt.show()