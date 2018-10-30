# Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os #platform independent paths
import matplotlib as mpl 

# Load user defined functions
import file_formats as fifo
# Load style sheet for plots
plt.style.use('default')

###################################################Experimental settings
CycleTime=0.075 # Aquisition cycle time [s]

# Define labels
labels=(['6-3'])
labels.extend(['3-3']) 


##################################################################################################### File read in
# Define folder of locs.hdf5 file
dir_names=[]
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/D134/18-06-24/N3-3_10nM-P1modA_p30_T23_10MHz-g300_field1_1/18-06-24_FS/'])


# Define names of locs_picked.hdf5 file
file_names=[]
file_names.extend(['N3-3_10nM-P1modA_p30_T23_10MHz-g300_field1_1_MMStack_Pos0.ome_locs_picked_groupprops.hdf5'])


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

for p in range(0,len(path)):  
    ################# Unfiltered groupprops to list
    groupprops=fifo.read_locs(path[p])# Single groupprops
    groupprops_list.append(groupprops)# List containing all groupprops arrays (unfiltered)

    ################# Set filters
    # Set single filters, always one single filter has to be active!
    istrue_single=[]
    istrue_single.append(groupprops_list[p]['std_frame']>12000)
    istrue_single.append(groupprops_list[p]['mean_frame']>1000)
    istrue_single.append(groupprops_list[p]['mean_frame']<100000)
#    istrue_single.append(groupprops_list[p]['ac_tau']<25) # Uncomment for usage
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
p=0
field='std_frame'

f=plt.figure(num=10,figsize=[4,4])
f.clear()
f.suptitle('Filter result')
ax=f.add_subplot(1,1,1)
ax.hist(groupprops_list_filter[p][field][:],bins='fd',label=labels[p]+', '+field)
ax.legend(loc=1)

################################################### Plotting 

################################################### Histogram plot for one file in path

f=plt.figure(num=1,figsize=[5,4]) # Open and clear figure 1,add 1 subplot
f.clear()

ax=f.add_subplot(2,1,1) # Add subplot one to figure (row,column,index) assign to variable ax
p=0 #Set path
field='ac_tau' #set field
ax.hist(groupprops_list_filter[p][field],bins='fd',label=labels[p]) # Plot histogram on ax

ax.legend(loc=1) # Set legend and location
ax.set_xlim([0,500]) # Set x-axis limits
ax.set_xlabel('x')  #Set x-label
ax.set_ylabel('y')  #Set y-label

################################################### Trend of one variable over different files

f=plt.figure(num=2,figsize=[5,4]) # Open and clear figure 2,add 1 subplot
f.clear()

ax=f.add_subplot(2,1,1) # Add subplot one to figure (row,column,index) assign to variable ax
field='A1' #set field

x=np.arange(0,len(means))
ax.errorbar(x,means[field],yerr=stds[field],label='trend') # Plot histogram on ax

ax.legend(loc=0) # Set legend and location
ax.set_xlim([-.1,1.1]) # Set x-axis limits
ax.set_xlabel('x')  #Set x-label
ax.set_ylabel('y')  #Set y-label