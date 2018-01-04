# Import modules
import numpy as np
from tqdm import tqdm
import h5py as h5py #hdf5 handling
import scipy.optimize
import multipletau
import importlib

# Import own modules in qPAINT
import fitfunc
importlib.reload(fitfunc)
import file_formats as fifo
importlib.reload(fifo)

#%%
def get_mean_x(locs):
    x=np.mean(locs['x'],axis=0)
    return x
#%%
def get_std_x(locs):
    x=np.std(locs['x'],axis=0)
    return x
#%%
def get_mean_y(locs):
    y=np.mean(locs['y'],axis=0)
    return y
#%%
def get_std_y(locs):
    y=np.std(locs['y'],axis=0)
    return y
#%%
def get_mean_frame(locs):
    frame=np.mean(locs['frame'],axis=0)
    return frame
#%%
def get_std_frame(locs):
    frame=np.std(locs['frame'],axis=0)
    return frame
#%%
def get_mean_photons(locs):
    photons=np.mean(locs['photons'],axis=0)
    return photons
#%%
def get_std_photons(locs):
    photons=np.std(locs['photons'],axis=0)
    return photons

#%%
def get_n_locs(locs):
    n_locs=np.size(locs['frame'])
    return n_locs

#%%
def get_trace(locs,NoFrames):
    # Define trace
    trace=np.zeros([1,NoFrames])
    # Add photons to trace for each frame
    trace[0,locs['frame']]=locs['photons']
    # Get unique frames and counts
    frames,frames_count=np.unique(locs['frame'],return_counts=True)
    # Frames with mutiple loaclizations
    frames_multi=frames[frames_count>1]
    # Number of frames with multiple localizations in group
    NoMultiFrames=np.size(frames_multi)
    # Replace single localization frame phtons by sum of multiple localization frame photons
    for i in range(0,np.size(frames_multi)):
        trace[0,frames_multi[i]]=np.sum(locs['photons'][locs['frame']==frames_multi[i]])
    
    
    return trace

#%%
def get_trace_bridged(locs,bridge_length,NoFrames):
    # Get trace, frames and photons 
    trace=get_trace(locs,NoFrames)
    frame=locs['frame']
    photons=locs['photons']
    # Get gap_length
    gap_length=frame[1:]-frame[0:-1]
    gap_length=np.float64(gap_length)
    # Set gap_length==1 to zero    
    gap_length[gap_length==1]=0
    # Set tau_ds>gap_length to zero
    gap_length[gap_length>bridge_length]=0
    # Make gap_length coincide with gap_starts in frame by adding 0 at beginning
    gap_length=np.concatenate([gap_length,[0]],axis=0)
    # Define start frames of gaps
    gap_start=frame[gap_length>0]
    # Define gap length
    gap_length=gap_length[gap_length>0]
    gap_length=np.uint32(gap_length)
    
    # Loop over gaps in localization frames
    gap_idx=np.empty([0],dtype=np.uint32)
    gap_value=np.empty([0],dtype=np.uint32)
    for i in range(0,np.size(gap_start)):
        # Inidices of gap to be bridged
        gap_frames=range(gap_start[i]+1,gap_start[i]+gap_length[i])
        # Inidices of gap to be bridged relative to gap_start frame        
        gap_frames_relative=range(1,gap_length[i])
        
        # Define parameters of linear gap function
        gap_function_grad=photons[frame==gap_start[i]+gap_length[i]]-photons[frame==gap_start[i]]
        gap_function_grad=gap_function_grad/gap_length[i]
        gap_function_offset=photons[frame==gap_start[i]]
        
        # Assign gap values according to linear function
        gap_value=np.hstack([gap_value,gap_function_grad*gap_frames_relative+gap_function_offset])
        # Assign gap indices
        gap_idx=np.hstack([gap_idx,gap_frames])
        
    trace[0,gap_idx]=gap_value
    
    return trace
  
#%%
def get_ac(locs,NoFrames):
    # Get trace 
    trace=get_trace(locs,NoFrames)
    # Multiple tau autocorrelation
    ac=multipletau.autocorrelate(np.transpose(trace),m=64, deltat=1,
                                 normalize=False,copy=False, dtype=np.float64())
    return ac 

#%%
def get_ac_bridged(locs,bridge_length,NoFrames):
    # Get trace 
    trace=get_trace_bridged(locs,bridge_length,NoFrames)
    # Multiple tau autocorrelation
    ac=multipletau.autocorrelate(np.transpose(trace),m=32, deltat=1,
                                 normalize=False,copy=False, dtype=np.float64())
    return ac

#%%
def get_ac_fit(locs,NoFrames,lastframe):
    # Get multiple tau autocorrelation
    ac=get_ac(locs,NoFrames)
    # Cut end tail of autocorrelation
    ac=ac[ac[:,0]<=lastframe,:]
    # Define start parameters for fit
    p0=np.empty([3])
    p0[0]=ac[1,1] # Amplitude
    p0[2]=ac[-1,1] # offset
    halfvalue=p0[2]+(p0[0]-p0[2])/2 # Value of half decay of ac
    p0[1]=np.argmin(np.abs(ac[:,1]-halfvalue)) # tau
    # Fit with monoexponential function
    popt,pcov=scipy.optimize.curve_fit(fitfunc.ac_monoexp,ac[1:,0],ac[1:,1],p0,method='lm')
    
    return popt

#%%
def get_ac_bridged_fit(locs,bridge_length,NoFrames,lastframe):
    # Get multiple tau autocorrelation
    ac=get_ac_bridged(locs,bridge_length,NoFrames)
    # Cut end tail of autocorrelation
    ac=ac[ac[:,0]<=lastframe,:]
    # Define start parameters for fit
    p0=np.empty([3])
    p0[0]=ac[1,1] # Amplitude
    p0[2]=ac[-1,1] # offset
    halfvalue=p0[2]+(p0[0]-p0[2])/2 # Value of half decay of ac
    p0[1]=np.argmin(np.abs(ac[:,1]-halfvalue)) # tau
    # Fit with monoexponential function
    popt,pcov=scipy.optimize.curve_fit(fitfunc.ac_monoexp,ac[1:,0],ac[1:,1],p0,method='lm')
    
    return popt
    
#%%
def get_tau(locs,ignore,**kwargs):
    # Get tau_d distribution
    dframes=locs['frame'][1:]-locs['frame'][0:-1]
    dframes=np.float64(dframes)    
    tau_d_dist=dframes[dframes>(ignore+1)] # tau_d distribution
    tau_d_dist=np.sort(tau_d_dist) # sort tau_d_dist
    
    # Get tau_b_distribution
    dframes[dframes<=(ignore+1)]=0 # Set bright frames to 0
    dframes[dframes>1]=1 # Set dark frames to 1
    dframes[dframes<1]=np.nan # Set bright frames to NaN
    
    mask_end=np.concatenate([dframes,[1]],axis=0) # Mask for end of events
    frames_end=locs['frame']*mask_end # Apply mask to frames to get end frames of events
    frames_end=frames_end[~np.isnan(frames_end)] 
    
    mask_start=np.concatenate([[1],dframes],axis=0) # Mask for start of events
    frames_start=locs['frame']*mask_start # Apply mask to frames to get start frames events
    frames_start=frames_start[~np.isnan(frames_start)]
    
    tau_b_dist=frames_end-frames_start+1 #tau_b distribution
    tau_b_dist=np.sort(tau_b_dist) # sort tau_b distribution
    
    # Remove all bright_frames<=1 from tau_b dist if bright_ignore set to True
    if kwargs['bright_ignore']==True:
        tau_b_dist=tau_b_dist[tau_b_dist>1]
    
    # Number of binding events
    n_events=np.size(tau_b_dist)
    
    # Create ecdf (empirical cumulative distribution function) of tau_b distribution
    importlib.reload(fitfunc)
    tau_b_bins,tau_b_cdf=fitfunc.ecdf(tau_b_dist)
    # Create ecdf (empirical cumulative distribution function) of tau_d distribution
    tau_d_bins,tau_d_cdf=fitfunc.ecdf(tau_d_dist)

    # Fitting cdf of tau_b
    if np.size(tau_b_cdf)==1: # 1 binding event!
        tau_b_mean=tau_b_dist
    elif np.size(tau_b_cdf)==2: # 2 binding events!
        tau_b_mean=np.mean(tau_b_dist)        
    elif np.size(tau_b_cdf)>2:
        idx=np.abs(tau_b_cdf-(1-np.exp(-1))).argmin() # Find idx in cfd closest to 1-e^-1 value
        p0_b=tau_b_bins[idx] # Start value= bin at which cfd closest to 1-e^-1 value
        try:
            if kwargs['fit']=='exp':
                tau_b_mean,pcov_b=scipy.optimize.curve_fit(fitfunc.exp_cdf,tau_b_bins,tau_b_cdf,p0_b) # Fitting
            elif kwargs['fit']=='lin':
                tau_b_bins=tau_b_bins[0:-1] # Get rid of last value due to linearization with log
                tau_b_cdf=tau_b_cdf[0:-1]
                tau_b_mean,pcov_b=scipy.optimize.curve_fit(fitfunc.exp_cdf_linearized,tau_b_bins,-np.log(np.abs(1-tau_b_cdf)),p0_b) # Fitting
        except RuntimeError:
            tau_b_mean=p0_b
            
    # Fitting cdf of tau_d
    if np.size(tau_d_cdf)==0: # Single binding event no dark time-> set to zero!
        tau_d_mean=0 
        p0_d=0 
    elif np.size(tau_d_cdf)==1: # 1 dark time-> only dark time!
        tau_d_mean=tau_d_dist 
        p0_d=tau_d_mean
    elif np.size(tau_d_cdf)==2: # 2 dark times-> Mean of 2 dark times!
        tau_d_mean=np.mean(tau_d_dist)
        p0_d=np.mean(tau_d_dist)
    elif np.size(tau_d_cdf)>2: # At least 3 dark times
        idx=np.abs(tau_d_cdf-(1-np.exp(-1))).argmin() # Find idx in cfd closest to 1-e^-1 value
        p0_d=tau_d_bins[idx] # Start value= bin at which cfd closest to 1-e^-1 value
        
        if kwargs['fit']=='exp': # Fit to exponential function
            tau_d_mean,pcov_d=scipy.optimize.curve_fit(fitfunc.exp_cdf,tau_d_bins,tau_d_cdf,p0_d) # Fitting
        
        elif kwargs['fit']=='lin': # Linearized fit via -log(1-ecdf)
            tau_d_bins=tau_d_bins[0:-1] # Get rid of last value due to linearization with log
            tau_d_cdf=tau_d_cdf[0:-1]
            tau_d_mean,pcov_d=scipy.optimize.curve_fit(fitfunc.exp_cdf_linearized,tau_d_bins,-np.log(np.abs(1-tau_d_cdf)),p0_d) # Fitting


    # Output configuration
    try:
        if kwargs['out']=='param':
            return [n_events,float(tau_b_mean),float(tau_d_mean),float(p0_d)]
        elif kwargs['out']=='all':
            qTau=np.empty(1,dtype=[('n_events','i4',1),('tau_b_mean','f4',1),
                                   ('tau_b_dist','f4',np.shape(tau_b_dist)),('tau_b_bins','f4',np.shape(tau_b_bins)),
                                   ('tau_b_cdf','f4',np.shape(tau_b_cdf)),
                                   ('tau_d_mean','f4',1),('tau_d_dist','f4',np.shape(tau_d_dist)),
                                   ('tau_d_bins','f4',np.shape(tau_d_bins)),('tau_d_cdf','f4',np.shape(tau_d_cdf)),
                                   ])
            qTau['n_events']=n_events
            qTau['tau_b_mean']=tau_b_mean
            qTau['tau_b_dist']=tau_b_dist
            qTau['tau_b_bins']=tau_b_bins
            qTau['tau_b_cdf']=tau_b_cdf
            qTau['tau_d_mean']=tau_d_mean
            qTau['tau_d_dist']=tau_d_dist
            qTau['tau_d_bins']=tau_d_bins
            qTau['tau_d_cdf']=tau_d_cdf
            
            return qTau
    
    except KeyError:
        return [n_events,float(tau_b_mean),float(tau_d_mean),float(p0_d)]

#%%
def locs2groupprops(path,ac_lastframe,ignore_dark):
    # Read in locs_picked.hdf5 file generated by Picasso
    locs=fifo.read_locs(path)
    
    # Read in locs_picked.yaml file generated by Picasso
    # Get TIF and Localize meta-data
    [TIFmeta,LOCmeta]=fifo.read_meta(path)
    # Number of frames in .tif stack used by Picasso Localize
    NoFrames=TIFmeta['Frames']
    
    
    # List of groups in locs
    group_list=np.unique(locs['group'])
    
    # Initiate groupprops
    groupprops=np.zeros(np.size(group_list),dtype=[('group','i4',1),('mean_x','f4',1),('std_x','f4',1),
                                    ('mean_y','f4',1),('std_y','f4',1),('mean_frame','f4',1),
                                    ('std_frame','f4',1),('mean_photons','f4',1),
                                    ('std_photons','f4',1),('n_locs','i4',1),
                                    ('n_events','i4',1),('tau_d','f4',1),('tau_b','f4',1),('tau_d_cdf','f4',1),
                                    ('n_events_ignore','i4',1),('tau_b_ignore','f4',1),
                                    ('ac_A','f4',1),('ac_tau','f4',1),('ac_offset','f4',1)
                                    ])

    groupprops['group'][:]=group_list
    
    
    for i in tqdm(range(0,np.size(group_list))):
        g=group_list[i]
        locs_g=locs[:][locs['group']==g]
        
        groupprops['mean_x'][i]=get_mean_x(locs_g)
        groupprops['std_x'][i]=get_std_x(locs_g)
        groupprops['mean_y'][i]=get_mean_y(locs_g)
        groupprops['std_y'][i]=get_std_y(locs_g)
        groupprops['mean_frame'][i]=get_mean_frame(locs_g)
        groupprops['std_frame'][i]=get_std_frame(locs_g)
        groupprops['mean_photons'][i]=get_mean_photons(locs_g)
        groupprops['std_photons'][i]=get_std_photons(locs_g)
        groupprops['n_locs'][i]=get_n_locs(locs_g)
        
        [n_events,tau_b,tau_d,tau_d_cdf]=get_tau(locs_g,ignore_dark,out='param',bright_ignore=False,fit='lin')
        groupprops['n_events'][i]=n_events
        groupprops['tau_b'][i]=tau_b
        groupprops['tau_d'][i]=tau_d
        groupprops['tau_d_cdf'][i]=tau_d_cdf
        
        [n_events_ignore,tau_b_ignore,tau_d_ignore,tau_d_cdf]=get_tau(locs_g,ignore_dark,out='param',bright_ignore=True,fit='lin')
        groupprops['n_events_ignore'][i]=n_events_ignore
        groupprops['tau_b_ignore'][i]=tau_b_ignore
        
        [ac_A,ac_tau,ac_offset]=get_ac_fit(locs_g,NoFrames,ac_lastframe)
        groupprops['ac_A'][i]=ac_A
        groupprops['ac_tau'][i]=ac_tau
        groupprops['ac_offset'][i]=ac_offset
           
    # Save groupprops in hdf5 file
    groupprops_file = h5py.File(path.replace('.hdf5','_groupprops.hdf5'), "w")
    dset=groupprops_file.create_dataset("groupprops", np.shape(groupprops), dtype=groupprops.dtype)
    dset[...]=groupprops
    groupprops_file.close() 
    
    # Create meta_dictionary to add to yaml file
    ADDdict={'Generated by':'locs_groupprops.locs2groupprops','ac_lastframe':ac_lastframe,'ignore_dark':ignore_dark,'tau_fit':'lin'}
    fifo.create_meta_locs2groupprops(path,ADDdict)
       
    return groupprops