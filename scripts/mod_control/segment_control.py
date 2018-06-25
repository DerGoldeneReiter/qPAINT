#   Template for segmenting data sets
#
#   read in locs_picked.hdf5 files to segment data sets for analyzation
#   read in render.hdf5 files to segment data sets you want to visualize in render
#%%
# Load modules
import os #platform independent paths
# Load&Reload own modules
import importlib
import file_formats as fifo
importlib.reload(fifo)

#%%
# Define folder of locs_picked.hdf5 file
dir_names=['/fs/pool/pool-schwille-paint/Data/D042/18-06-22_1stOrigamiZyla/gauss_20mW_fr10k_exp200_1/18-06-22_FS/']

# Define names of locs_picked.hdf5 file
file_names=['gauss_20mW_fr10k_exp200_1_MMStack.ome_locs.hdf5']

# Create full path list
path=[]
for i in range(0, len(file_names)):
    path.append(os.path.join(dir_names[i],file_names[i]))

# Select first string of path
path=path[0]

# Input parameters for segmentation
NoFrames=fifo.read_meta(path)[0]['Frames'] # Number of frames in original file
NoSegments = 10         # In how many parts do you want to segment your file
NoFrames_new = int(NoFrames / NoSegments)  #how many frames will one segment contain

#%%    
# Segment a locs_picked-file into segements of NoFrames_new starting with frame=0 
locs_list = fifo.segment_locs(path, NoFrames_new, NoSegments) 
# Save segments in .hdf5 containers
fifo.save_locs_segments(locs_list,path)           
# Create a .yaml file for each segment
fifo.create_meta_segments(path,NoSegments,NoFrames_new)
