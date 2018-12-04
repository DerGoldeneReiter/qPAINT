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
dir_names=['/fs/pool/pool-schwille-paint/Data/p02.ds-bleaching/D042/18-11-13_p1Linker_OxCat/id15-80_P1-Linker-10nM_OxCat-p100mW-30d5deg_gauss_1/18-11-13_JS/18-11-15_FS/']

# Define names of locs_picked.hdf5 file
file_names=['id15-80_P1-Linker-10nM_OxCat-p100mW-30d5deg_gauss_1_MMStack.ome_locs_render_s1_picked_avg3_lineup.hdf5']

# Create full path list
path=[]
for i in range(0, len(file_names)):
    path.append(os.path.join(dir_names[i],file_names[i]))

# Select first string of path
path=path[0]
# Read yaml
dict_list=fifo.read_yaml(path)

# Input parameters for segmentation
NoFrames=dict_list[0]['Frames'] # Number of frames in original file
NoSegments = 5       # In how many parts do you want to segment your file
NoFrames_new = int(NoFrames / NoSegments)  #how many frames will one segment contain

# Change entry of yaml to new no. of frames
dict_list[0]['Frames']=NoFrames_new
# New dict for yaml
ADDdict=[{'Created by':'segment_control','Segments':NoSegments,'Frames per segment':NoFrames_new}]
#%%    
# Segment a locs_picked-file into segements of NoFrames_new starting with frame=0 
locs_list = fifo.segment_locs(path, NoFrames_new, NoSegments) 

for i in range(0, len(locs_list)):
    fifo.save_locs(locs_list[i],path.replace('.hdf5','_'+str(i+1) + '_of_' + str(int(NoSegments))+'.hdf5'))
    fifo.create_yaml(dict_list+ADDdict,path.replace('.hdf5','_'+str(i+1) + '_of_' + str(int(NoSegments))+'.hdf5'))
    
