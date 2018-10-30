# Load packages
import h5py as h5py #hdf5 handling
import yaml
import numpy as np
from tqdm import tqdm
import importlib
import copy

import locs_groupprops as l2grp
importlib.reload(l2grp)
############################################################### Read-in of hdf5-files
#%%
def read_locs(path):
    # Open locs_picked.hdf5 file
    locs_file=h5py.File(path,'r')
    # Load dataset 'locs' into np.array 'locs'
    locs=locs_file['locs'][...]
    locs_file.close()
    
    return locs

#%%
def read_pickprops(path):
    
    pickprops_file=h5py.File(path,'r')
    for name in pickprops_file: print(name)
    pickprops=pickprops_file['groups'][...]
    pickprops_file.close()
    
    return pickprops
############################################################### Save in hdf5-files
#%%
def save_locs(locs,path):
    # Save groupprops in hdf5 file
    locs_file = h5py.File(path,"w")
    dset=locs_file.create_dataset("locs", np.shape(locs), dtype=locs.dtype)
    dset[...]=locs
    locs_file.close() 
    return
#%%           
def save_locs_segments(locs_list,path):
    size = len(locs_list)
   
    for i in range(0, size):
        # Create savepath
        save_path = path.replace('.hdf5','_'+str(i+1) + '_of_' + str(size)+'.hdf5')
        # Create hdf5 container
        hf = h5py.File(save_path,'w')
        # Dump locs in container
        hf.create_dataset('locs', data = locs_list[i])       
        hf.close()
    return
############################################################### locs file manipulation (combine groups, segment by time,...)
#%%
def combine_locs(path,path_save):
   # Open first locs file 
   locs_file=h5py.File(path[0],'r')
   # Load dataset into locs
   locs=locs_file['locs'][...]
   
   for p in range(1,len(path)):
        # Open locs_picked.hdf5 file
        locs_file=h5py.File(path[p],'r')
        # Load dataset 'locs' into np.array 'locs'
        locs=np.append(locs,locs_file['locs'][...])
        locs_file.close()

   # Save combined locs
   locs_file_save = h5py.File(path_save, "w")
   dset=locs_file_save.create_dataset("locs", np.shape(locs), dtype=locs.dtype)
   dset[...]=locs
   locs_file_save.close()
   
   return locs
#%%
def segment_locs(path, NoFrames_new, NoSegments):
    # Read in locs-file
    locs=read_locs(path)
    # Sort after frames
    locs.sort(order=['frame'],axis=0)
    locs_list = []
    start = 0
    end = 0
    upper_boarder = NoFrames_new
    
    # run over all segments
    for i in range(0, int(NoSegments)): 
        # run over all entries with frame number within the selected segment
        for j in range(start, len(locs)):     
            # increment end for each entry with frame number within boarders
            if locs[j][0] < upper_boarder:
                end += 1
            # break, if the frame number exceeds the segment
            else:
                break
        # append a segment to the list by truncating the original array
        locs_list.append(locs[start:end])
        # adapt the frame number to start from 0 for each segment
        locs_list[i]['frame'][:] = locs_list[i]['frame'][:] % NoFrames_new
        # make the end point the new starting point for the next segment
        start = end
        # increase the upper boarder for the frame number for the next segment 
        upper_boarder += NoFrames_new
        
    return locs_list
#%%
def remove_field(locs,field):
    names = list(locs.dtype.names)
    if field in names:
        names.remove(field)
    locs_new = locs[names]
    return locs_new

#%%
def rotate_locs(locs,angle=0,origin=np.zeros(2)):
    locs=np.copy(locs)
    # Substract origin from locs x&y coordinates
    locs['x']=locs['x']-origin[0]
    locs['y']=locs['y']-origin[1]
    # Define 2dim rotation matrix for clockwise rotation in x-y plane
    angle = np.radians(angle)
    c, s = np.cos(angle), np.sin(angle)
    R = np.array(((c,-s), (s, c)))
    # Rotate localizations
    x_rot=R[0,0]*locs['x']+R[0,1]*locs['y']
    y_rot=R[1,0]*locs['x']+R[1,1]*locs['y']
    # Assign new x-y coordinates to locs
    locs['x']=x_rot
    locs['y']=y_rot
    # Shift to original position forsame field of view, i.e. add origin to x&y
    locs['x']=locs['x']+origin[0]
    locs['y']=locs['y']+origin[1]
    return locs
#%%
def lineup_locs(path,x_incr=3):
    # Create list of locs lead by path list
    locs_list=[]
    for p in range(0,len(path)):
        locs_list.append(read_locs(path[p]))# Single groupprops
    # Deepcopy locs_list
    locs_list_lineup=copy.deepcopy(locs_list)
    # List of groups in p=0
    group_list=np.unique(locs_list[0]['group'])
    
    x_add=3 # x-add-value that will be added to x to line up picks
    y_add=read_yaml(path[p])[0]['Height']/2 #y-add-value that will be added to y to put line in center
    
    # Loop over all groups
    for i in tqdm(range(0,np.size(group_list))):
        g=group_list[i]
        
        # Get locs for group=0 in locs_list[p=0]
        locs_g=locs_list[0][:][locs_list[0]['group']==g]
        x_subtract=l2grp.get_mean_x(locs_g) # Mean of x in group in p=0
        y_subtract=l2grp.get_mean_y(locs_g) # Mean of y in group in p=0
        
        # Lopp over paths
        for p in range(0,len(path)):
            # Boolean index for group
            istrue=locs_list[p]['group']==g
            # Align x to zero
            locs_list_lineup[p]['x'][istrue]=locs_list[p]['x'][istrue]-x_subtract
            # Align y to zero
            locs_list_lineup[p]['y'][istrue]=locs_list[p]['y'][istrue]-y_subtract
            # Add increment to x
            locs_list_lineup[p]['x'][istrue]=locs_list_lineup[p]['x'][istrue]+x_add
            # Add increment to y
            locs_list_lineup[p]['y'][istrue]=locs_list_lineup[p]['y'][istrue]+y_add
        # Increase x_add by increment
        x_add=x_add+x_incr
        
    # Remove group IDs, i.e. from locs_picked to locs
    locs_list_lineup_remove=[]
    for p in range (0,len(path)):
        locs_list_lineup_remove.append(remove_field(locs_list_lineup[p],'group'))
    
    ADDdict={'Generated by':'file_formats.lineup_locs'}
    # Save loop
    for p in range (0,len(path)):
        # Create dictionary
        dict_list=read_yaml(path[p]) # Get yaml dicts
        create_yaml(dict_list+ADDdict,path[p].replace('.hdf5','_lineup.yaml'))
        #Save locs
        save_locs(locs_list_lineup_remove[p],path[p].replace('.hdf5','_lineup.hdf5'))
            
    return locs_list_lineup
    
############################################################### Read and create meta data in yaml format
#%%
def read_yaml(readpath):
    # Read in .yaml file generated by picasso
    stream=open(str.replace(readpath,'.hdf5','.yaml'))
    doc=yaml.load_all(stream)
    
    dict_list=[]
    i=0
    for data in doc:
        dict_single=data
        dict_list.append(dict_single)
        i=i+1
    
    stream.close()
    
    return dict_list
#%%
def create_yaml(dict_list,savepath):    
    # open .yaml file
    with open(str.replace(savepath,'.hdf5','.yaml'),'w') as yaml_file:
        # Dump first dict without explicit start
        yaml.dump(dict_list[0],yaml_file,default_flow_style=False)
        # Dump all other dicts with explicit start
        for i in range(1,len(dict_list)): 
            yaml.dump(dict_list[i],yaml_file,default_flow_style=False,explicit_start=True)
        
    return

