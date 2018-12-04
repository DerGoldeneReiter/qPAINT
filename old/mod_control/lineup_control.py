# Load packages
import numpy as np #numpy data formats and operators
import os #platform independent paths
import matplotlib.pyplot as plt #plotting
import importlib
# Load user defined functions
import file_formats as fifo
import locs_groupprops as l2grp
importlib.reload(fifo)
importlib.reload(l2grp)
#%%
##################################################################################################### Set paths to files and labels
# Define folder of locs.hdf5 file
dir_names=[]
file_names=[]

dir_names.extend(['/fs/pool/pool-schwille-paint/Data/D042/18-08-14_FlatGauss_Dense_20nm/id15-80_P1-Cy3b-20nM_p200mW-38deg_flat_1/18-08-14_FS/']*2)

file_names.extend(['id15-80_P1-Cy3b-20nM_p200mW-38deg_flat_1_MMStack.ome_locs_render_picked_filter-in.hdf5'])
file_names.extend(['id15-80_P1-Cy3b-20nM_p200mW-38deg_flat_1_MMStack.ome_locs_render_filter-phot10k-35k_picked_filter-out.hdf5']) 


# Create list of paths
path=[]
for i in range(0,len(file_names)):
    path.extend([os.path.join(dir_names[i],file_names[i])])
    
locs_list_new=fifo.lineup_locs(path)




#%%
######################################### Select group
g=14
px2dist=130

######################################### Plotting
f=plt.figure(num=10,figsize=[5,5])
f.clear()
################################################################################ x vs y
ax=f.add_subplot(1,1,1)
ax.set_facecolor('k')

locs_g=locs_list_new[0][:][locs_list_new[0][:]['group']==g]
x=locs_g['x']*px2dist
y=locs_g['y']*px2dist
ax.scatter(x,y,s=100,marker='o',edgecolor='none',alpha=0.6,color='cyan')


