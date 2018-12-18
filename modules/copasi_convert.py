import h5py as h5py #hdf5 handling
import numpy as np #numpy data formats and operators
from tqdm import tqdm
import pandas as pd

import var_io as io         
#%%
def copasi2locs(path,interval_size,intervals):
  
    # Get number of lines in .csv
    n_l=len(open(path).readlines())*0.5
    # Calculate number of iterations of experiments in csv file
    n_iter=int(n_l/(intervals+2))
    
    # Read in .txt as pandas TextFileReader object
    reader=pd.read_csv(path,
                       sep='\t',
                       header=None,
                       iterator=True
                       )
    
    # Define array for saving individual groups
    arr_all=np.zeros([0,6])
    # Iterate over all experiments
    for g in tqdm(range(0,n_iter)):
        # Save data to np.array
        arr=reader.get_chunk(2*(intervals+1)).values
        # Remove double entries
        arr=arr[np.arange(0,len(arr),2),:]
        # Remove zeros
        arr=arr[arr[:,1]!=0,:]
        # Convert time to frames
        arr[:,0]=arr[:,0]*int(1/interval_size)
        # Add column for x, y and bg
        arr=np.c_[arr,np.zeros(len(arr)),np.zeros(len(arr)),np.zeros(len(arr))]
        # Add column for group
        arr=np.c_[arr,np.ones(len(arr))*g]
        # Add arr to array of all arrays arr_all
        arr_all=np.r_[arr_all,arr]
    
     
    # Define dtype for locs
    dt=np.dtype([('frame', '<i4'), ('photons', 'f4'), ('x', 'f4'), ('y', 'f4'),('bg', 'f4'), ('group', 'i4')])
    # Create locs np.array
    locs=np.zeros(len(arr_all),dtype=dt)
    
    locs['frame']=arr_all[:,0]
    locs['photons']=arr_all[:,1]
    locs['x']=arr_all[:,2]
    locs['y']=arr_all[:,3]
    locs['bg']=arr_all[:,4]
    locs['group']=arr_all[:,5]
    
    savepath=path.replace('.txt','_locs_picked.hdf5')
    
    # Save locs in hdf5 file
    locs_file = h5py.File(savepath, "w")
    dset=locs_file.create_dataset("locs", np.shape(locs), dtype=locs.dtype)
    dset[...]=locs
    locs_file.close() 
    
    # Create .yaml file containing meta data
    TIFdict=[{'Camera': 'Simulation',
             'Frames': intervals+1,
             'Generated by': 'copasi',
             'NoGroups':n_iter}]
    LOCdict=[{'Generated by':'simulate_locs.copasi2locs'}]
    io.create_yaml(TIFdict+LOCdict,savepath)
    
           

    return locs

#%%
def copasi2locs_double(path,interval_size,intervals):
  
    # Get number of lines in .csv
    n_l=len(open(path).readlines())*0.5
    # Calculate number of iterations of experiments in csv file
    n_iter=int(n_l/(intervals+2))
    
    # Read in .txt as pandas TextFileReader object
    reader=pd.read_csv(path,
                       sep='\t',
                       header=None,
                       iterator=True
                       )
    
    # Define array for saving individual groups
    arr_all=np.zeros([0,5])
    # Iterate over all experiments
    for g in tqdm(range(0,n_iter)):
        # Save data to np.array
        arr=reader.get_chunk(2*(intervals+1)).values
        # Remove double entries, i.e. time1,time1,time2,time2,....
        arr=arr[np.arange(0,len(arr),2),:]
        # Sum of particle numbers of column 1&2 -> column 1
        arr[:,1]=arr[:,1]+arr[:,2]
        arr=arr[:,0:2]
        # Remove zeros, i.e. no localization
        arr=arr[arr[:,1]!=0,:]
        # Convert time to frames
        arr[:,0]=arr[:,0]*int(1/interval_size)
        # Add column for x and y
        arr=np.c_[arr,np.zeros(len(arr)),np.zeros(len(arr))]
        # Add column for group
        arr=np.c_[arr,np.ones(len(arr))*g]
        # Add arr to array of all arrays arr_all
        arr_all=np.r_[arr_all,arr]
    
     
    # Define dtype for locs
    dt=np.dtype([('frame', '<i4'), ('photons', 'f4'), ('x', 'f4'), ('y', 'f4'), ('group', 'i4')])
    # Create locs np.array
    locs=np.zeros(len(arr_all),dtype=dt)
    
    locs['frame']=arr_all[:,0]
    locs['photons']=arr_all[:,1]
    locs['x']=arr_all[:,2]
    locs['y']=arr_all[:,3]
    locs['group']=arr_all[:,4]
    
    savepath=path.replace('.txt','_locs_picked.hdf5')
    
    # Save locs in hdf5 file
    locs_file = h5py.File(savepath, "w")
    dset=locs_file.create_dataset("locs", np.shape(locs), dtype=locs.dtype)
    dset[...]=locs
    locs_file.close() 
    
    # Create .yaml file containing meta data
    TIFdict={'Camera': 'Simulation',
             'Frames': intervals+1,
             'Generated by': 'copasi',
             'NoGroups':n_iter}
    LOCdict={'Generated by':'simulate_locs.copasi2locs'}
    io.create_meta_simulate(savepath,TIFdict,LOCdict)
    
           

    return locs  