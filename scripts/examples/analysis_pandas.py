#################################################### Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os #platform independent paths
import importlib
import pandas as pd 
# Load user defined functions
import var_io as io
import multitau
#%%
############################################################## Parameters & labels
labels=[]
labels.extend(['exp01'])

############################################################## Define data paths
#### Define folder of .hdf5 file
dir_names=[] 
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/z.PAINT-checks/18-11-02_SDS-P1-9nt/SDS9nt-80_P1-20nM_p100mW-8deg_flat_1/18-11-02_FS/'])
#### Define names of .hdf5 file
file_names=[] 
file_names.extend(['SDS9nt-80_P1-20nM_p100mW-8deg_flat_1_MMStack.ome_locs_picked_d1-0_props-1.hdf5'])

############################################################## Read in data
#### Create list of paths
path=[os.path.join(dir_names[i],file_names[i]) for i in range(0,len(file_names))]
#### Read in locs of path
locs_props=pd.concat([io.read_locs(p)[0] for p in path],keys=labels,names=['expID'])
#### Reset index for better viewing in spyder
locs_props=locs_props.reset_index()

#%%
############################################################# Filter
#### Copy locs props for filtering
locs_props_filter=locs_props.copy()

#### Set filters
locs_props_filter=locs_props_filter[locs_props_filter['std_frame']>5000]
locs_props_filter=locs_props_filter[locs_props_filter['mean_frame']>2000]

############################################################# Plotting
f=plt.figure(num=1)
f.clear()
ax=f.add_subplot(111)

#### Histogram
#ax.hist(locs_props_filter['mono_tau'],bins=np.arange(0,50,1))

#### Scatter
#ax.scatter(locs_props_filter['mono_tau'],locs_props_filter['mono_A'],alpha=0.2)

#### Plot of trace of group=g
#g=4
#trace=locs_props_filter['trace'].loc[locs_props_filter['group']==g].values[0]
#ax.plot(trace)




############################################################# Advanced

def ownfunc(df):
    #### Get bright time distribution
    tau_b_dist=df['tau_b_dist'].values[0]
    #### Define cirtical value
    tcrit=30
    #### Get number of bright times grater than critical value
    N_greater_tcrit=len(tau_b_dist[tau_b_dist>tcrit])
    #### Assign to pandas.Series
    s_out=pd.Series({'N_greater_tcrit':N_greater_tcrit})
    
    return s_out

props_of_props=locs_props_filter.groupby('group').apply(ownfunc)

