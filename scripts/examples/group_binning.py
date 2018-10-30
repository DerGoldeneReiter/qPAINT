# Load packages
import numpy as np #numpy data formats and operators
import os #platform independent paths
import importlib
import copy
from tqdm import tqdm
# Load user defined functions
import file_formats as fifo
import fitfunc
import locs_groupprops as l2grp
importlib.reload(fifo)
importlib.reload(fitfunc)
importlib.reload(l2grp)
#%%
#####################################################################################################
bin=20
##################################################################################################### Set paths to files and labels
# Define folder of locs.hdf5 file
dir_names=['/fs/pool/pool-schwille-paint/Data/Simulation/18-06-04_copasi_bi_2s-6s/100k']*2
# Define names of locs_picked.hdf5 file
file_names=['taub2-6s_kon2-3e6_c10nM_exp100_f100k_N12-3_locs_picked_groupprops.hdf5'] # _groupprops file
file_names.extend(['taub2-6s_kon2-3e6_c10nM_exp100_f100k_N12-3_locs_picked.hdf5']) # _picked file

# Create list of paths
path=[]
for i in range(0,len(file_names)):
    path.extend([os.path.join(dir_names[i],file_names[i])])  
    

#################################################################################################### Read, Filter, Means&Stds
groupprops_list=[]
locs_list=[]
groupprops_list_filter=[]
groupprops_list_filter_out=[]
NoFrames=[]
for p in range(0,1):
    ################# Unfiltered groupprops to list
    groupprops=fifo.read_groupprops(path[p])# Single groupprops
    groupprops_list.append(groupprops)# List containing all groupprops arrays (unfiltered)
    locs=np.copy(fifo.read_locs(path[p+1])) # Single picked
    locs_list.append(np.copy(locs))# List containing all picked arrays (unfiltered)
    # Get TIF and Localize meta-data
    [TIFmeta,LOCmeta]=fifo.read_meta(path[p+1])
    # Number of frames in .tif stack used by Picasso Localize
    NoFrames.append(TIFmeta['Frames'])
    ################# Set filters
    istrue_single=[]
    # Set single filters, always one single filter has to be active!
    istrue_single.append(groupprops_list[p]['std_frame']>1000)
        
    # Combine single filters to overall filter (and not)
    istrue=np.all(istrue_single,axis=0)
    istrue_not=~np.array(istrue)
    
    ################# Filtered groupprops to list
    groupprops_list_filter.append(groupprops_list[p][:][istrue]) # Create list of filtered groupprops
    groupprops_list_filter_out.append(groupprops_list[p][:][istrue_not]) # Create list of out-filtered groupprops
################################################### Show filtering   
p=2
field='ac_chisquare'

#f=plt.figure(num=10,figsize=[7,7])
#f.clear()
#ax=f.add_subplot(1,1,1)
#ax.hist(groupprops_list_filter[p][field][:],bins='fd',label=labels[p]+', '+field)
#ax.legend(loc=1)
#################################################################################################### Binning of groups

locs_list_bin=[]

for p in range(0,1):
    ######################### Create look-up for groups to binned groups
    group_list=np.unique(groupprops_list_filter[p]['group']) # List of groups in filtered groupprops
    N_groups=len(group_list) # Total number of groups in group_list
    N_groups_bin=np.floor(N_groups/bin) # Integer number of binned groups 
    group2group_bin=np.zeros([len(group_list),2]) # Look-up array for groups to binned groups
    group2group_bin[:,0]=group_list # Assign original groups to first column
    # Assign binned group to second column
    new_id=1 # Binned groups will start with id=1, leftover will be id=0
    for i in np.int64(np.arange(0,bin*N_groups_bin,bin)):
        group2group_bin[i:i+bin,1]=new_id
        new_id=new_id+1
        
    # List of binned groups
    group_bin_list=np.unique(group2group_bin[:,1])
    # Binning of origamis, sum traces
    locs_bin=np.array([],dtype=locs.dtype)
    
    for gbin in tqdm(group_bin_list): # Loop through binned groups
        # Trace with bin channels
        trace_bin=np.zeros([NoFrames[p],bin])
        # Frames and photons for assignment to sum of binned traces
        framestrace=np.zeros([NoFrames[p],2])
        framestrace[:,0]=np.arange(0,NoFrames[p],1)
        for g in group2group_bin[group2group_bin[:,1]==gbin,0]: # Get original groups belonging to binned groups
            # Get locs of group
            locs_g=locs_list[p][:][locs_list[p]['group']==g]
            # Assign trace with bin channels
            trace_bin[:,np.mod(int(g),bin)]=l2grp.get_trace(locs_g,NoFrames[p])
        
        # Assign to sum of channels in trace_bin to framestrace
        framestrace[:,1]=np.sum(trace_bin,axis=1)
        # Remove all zero entries from framestrace
        framestrace=framestrace[framestrace[:,1]>0,:]
        # Create locs_g_bin
        locs_g_bin=np.zeros([len(framestrace)],dtype=locs_g.dtype)
        # Assign frames and photons and binned group number
        locs_g_bin['frame']=framestrace[:,0]
        locs_g_bin['photons']=framestrace[:,1]
        locs_g_bin['group']=np.ones(len(framestrace))*gbin
        # Stack to locs_bin
        locs_bin=np.hstack((locs_bin,locs_g_bin))
    
    locs_list_bin.append(locs_bin)
    
    # Save binned picked file
    savename=path[p+1].replace('picked','picked_bin%i'%(bin))
    fifo.save_locs(locs_list_bin[p],savename)
    
    # Create dictionary
    ADDdict={'Generated by':'group_binning','bin':bin}
    fifo.create_meta_locs2groupprops(path[p+1],ADDdict,extension='_bin%i.yaml'%(bin))

    # Save group2group_bin
    np.savetxt(savename.replace('.hdf5','.csv'),group2group_bin)
