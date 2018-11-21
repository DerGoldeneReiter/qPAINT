#Test script for hdf5 handling in pandas package
#
#
#
#################################################### Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os #platform independent paths
import importlib
import pandas as pd
from tqdm import tqdm
# Load user defined functions
import pickprops as props
# Reload modules
importlib.reload(props)

############################################################## Define data
# Define folder of locs.hdf5 file
dir_names=[] 
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/z.PAINT-checks/18-11-07_SDS_P1-9nt_noPur/Schueder_SDS_9nt_npPur_P1-20nM_p100mW-8deg_flat_2/18-11-07_FS/'])
# Define names of locs_picked.hdf5 file
file_names=[] 
file_names.extend(['Schueder_SDS_9nt_npPur_P1-20nM_p100mW-8deg_flat_2_MMStack.ome_locs_picked.hdf5'])
# Create list of paths
path=[]
for i in range(0,len(file_names)):
    path.extend([os.path.join(dir_names[i],file_names[i])])

############################################################# HDF5 to pd.DataFrame
store=pd.HDFStore(path[0]) # Read hdf5 container
picked=store['locs'] # Load locs as pandas.DataFrame
picked=picked.abs() # Take absolute value of every variable
picked[['group','frame']]=picked[['group','frame']].astype(int) # Convert fields frame and group to integers as they can be used for indexing


############################################################# Apply functions to grouped data
NoFrames=18000

tqdm.pandas()
picked_props=picked.groupby(['group']).progress_apply(lambda df: props.get_ac(df,NoFrames)).reset_index()


#%%
############################################################# Test application to single pick df
g=1
picked_g=picked[picked['group']==g]
s_out=props.get_ac(picked_g,NoFrames)


f=plt.figure(num=1)
f.clear()
ax=f.add_subplot(2,1,1)
ax.plot(np.arange(0,NoFrames,1),s_out['trace'])
ax=f.add_subplot(2,1,2)
ax.plot(s_out['tau'],s_out['g'])
ax.set_xscale('log')

