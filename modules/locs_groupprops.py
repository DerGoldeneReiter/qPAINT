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
import multitau
importlib.reload(multitau)

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
    
    # Correct for negative photon values
#    trace=np.abs(trace)
    return trace

#%%
def get_ac(locs,NoFrames,m=16):
    # Get trace 
    trace=get_trace(locs,NoFrames)
    # Multiple tau autocorrelation
    ac=multitau.autocorrelate(np.transpose(trace),m=m, deltat=1,
                                 normalize=False,copy=False, dtype=np.float64())
#    ac=ac[:-(m),:] # Cut off tail due to added zeros in non-normalized mode
    return ac 

#%%
def get_ac_fit(locs,NoFrames):
    # Get multiple tau autocorrelation
    ac=get_ac(locs,NoFrames)
    # Define start parameters for fit
    p0=np.empty([3])
    p0[0]=ac[1,1] # Amplitude
    p0[2]=ac[-1,1] # offset
    halfvalue=p0[2]+(p0[0]-p0[2])/2 # Value of half decay of ac
    p0[1]=np.argmin(np.abs(ac[:,1]-halfvalue)) # tau
    # Bounds for fit 
    lowbounds=np.array([0,0,0])
    upbounds=np.array([np.inf,np.inf,np.inf])
    # Fit with monoexponential function, disregard
    try:
        popt,pcov=scipy.optimize.curve_fit(fitfunc.ac_monoexp,ac[1:,0],ac[1:,1],p0,bounds=(lowbounds,upbounds),method='trf')
    except RuntimeError:
        popt=p0
        
    chisquare=np.sum(np.square(np.divide(fitfunc.ac_monoexp(ac[1:,0],*popt)-ac[1:,1],popt[0]+popt[2])))/len(ac)
    popt=np.append(popt,np.sqrt(chisquare))
    return popt

#%%
def get_ac_NoDocks(locs,NoFrames,p_inf_1):
    # Get ac fit parameters
    [G_0,tau_c,G_inf,chisquare]=get_ac_fit(locs,NoFrames)
    # p_inf_N
    p_inf_N=G_inf/G_0
    # NoDocks
    NoDocks=p_inf_N/p_inf_1[0]
    # err_NoDocks
    err_NoDocks=NoDocks*(p_inf_1[1]/p_inf_1[0])
    return [NoDocks,err_NoDocks]
   
    #%%
def get_misc_NoDocks(locs,NoFrames,ignore=1,**kwargs):
    # Get ac fit parameters
    [G_0,tau_c,G_inf,chisquare]=get_ac_fit(locs,NoFrames)
    p_inf=G_inf/G_0 # p_inf
    # Get Picasso fit parameters
    [n_events,tau_b_pi,tau_d_pi]=get_tau(locs,ignore,out='param',bright_ignore=True,fit='lin')
    
    # Get tau_b by using misc terms p_inf and tau_d_pi
    tau_b=p_inf*tau_d_pi
    tau_d=(tau_b*tau_c)/np.abs(tau_b-tau_c)
    NoDocks=tau_d/tau_d_pi
    return [tau_b,tau_d,NoDocks]

#%%
def get_ac_tau(locs,NoFrames,NoDocks):
    # Get ac fit parameters
    [G_0,tau_c,G_inf,chisquare]=get_ac_fit(locs,NoFrames)
    # p_inf
    p_inf=G_inf/G_0
    # tau_b extracetd by ac
    ac_tau_b=(1+p_inf/NoDocks)*tau_c
    # tau_d extracetd by ac
    ac_tau_d=(1+NoDocks/p_inf)*tau_c
    return [p_inf,ac_tau_b,ac_tau_d]
     
#%%
def get_tau(locs,ignore=1,**kwargs):
    ################################################################ Get tau_d distribution
    dframes=locs['frame'][1:]-locs['frame'][0:-1] # Differentiate frames
    dframes=np.float64(dframes)    
    tau_d_dist=dframes[dframes>2] # Remove all dark frame-distances >2 (i.e. 1 dark frame between two bright frames neglected)
    tau_d_dist=np.sort(tau_d_dist) # Sorted tau_d distribution
    
    ################################################################# Get tau_b_distribution
    dframes[dframes<=(ignore+1)]=0 # Set bright frames to 0, here ignore value from input argument is used!
    dframes[dframes>1]=1 # Set dark frames to 1
    dframes[dframes<1]=np.nan # Set bright frames to NaN
    
    mask_end=np.concatenate([dframes,[1]],axis=0) # Mask for end of events, add 1 at end
    frames_end=locs['frame']*mask_end # Apply mask to frames to get end frames of events
    frames_end=frames_end[~np.isnan(frames_end)] # get only non-NaN values, removal of bright frames
    
    mask_start=np.concatenate([[1],dframes],axis=0) # Mask for start of events, add one at start
    frames_start=locs['frame']*mask_start # Apply mask to frames to get start frames events
    frames_start=frames_start[~np.isnan(frames_start)] # get only non-NaN values, removal of bright frames
    
    tau_b_dist=frames_end-frames_start+1 # get tau_b distribution
    tau_b_dist=np.sort(tau_b_dist) # sort tau_b distribution
    
    # Remove all bright_frames<=1 from tau_b dist if bright_ignore set to True
    if kwargs['bright_ignore']==True:
        tau_b_dist=tau_b_dist[tau_b_dist>1]
    
    # Number of binding events
    n_events=np.size(tau_b_dist)
    
    ################################################################# Extract taus: Mean or ECDF-fit
    # Create ecdf (empirical cumulative distribution function) of tau_b distribution
    importlib.reload(fitfunc)
    tau_b_bins,tau_b_cdf=fitfunc.ecdf(tau_b_dist)
    # Create ecdf (empirical cumulative distribution function) of tau_d distribution
    tau_d_bins,tau_d_cdf=fitfunc.ecdf(tau_d_dist)

    ################ Extract tau_b_mean
    if np.size(tau_b_dist)<=5: # If size of tau_b distribution<=5 --> mean of distribution
        tau_b_mean=np.mean(tau_b_dist)      
    elif np.size(tau_b_dist)>5: # If size of tau_b distribution>5 --> Fit ECDF
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
            
    ################ Extract tau_d_mean
    if np.size(tau_d_dist)<=5: # If size of tau_d distribution<=5 --> mean of distribution
        tau_d_mean=np.mean(tau_d_dist)     
    elif np.size(tau_d_dist)>5: # If size of tau_d distribution>5 --> Fit ECDF
        idx=np.abs(tau_d_cdf-(1-np.exp(-1))).argmin() # Find idx in cfd closest to 1-e^-1 value
        p0_d=tau_d_bins[idx] # Start value= bin at which cfd closest to 1-e^-1 value
        try:
            if kwargs['fit']=='exp': # Fit to exponential function
                tau_d_mean,pcov_d=scipy.optimize.curve_fit(fitfunc.exp_cdf,tau_d_bins,tau_d_cdf,p0_d) # Fitting
            
            elif kwargs['fit']=='lin': # Linearized fit via -log(1-ecdf)
                tau_d_bins=tau_d_bins[0:-1] # Get rid of last value due to linearization with log
                tau_d_cdf=tau_d_cdf[0:-1]
                tau_d_mean,pcov_d=scipy.optimize.curve_fit(fitfunc.exp_cdf_linearized,tau_d_bins,-np.log(np.abs(1-tau_d_cdf)),p0_d) # Fitting
        except RuntimeError:
            tau_d_mean=p0_d

    # Output configuration
    try:
        if kwargs['out']=='param':
            return [n_events,float(tau_b_mean),float(tau_d_mean)]
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
        return [n_events,float(tau_b_mean),float(tau_d_mean)]


#%%
def locs2groupprops(path,ignore_dark=1,**kwargs):
    # Make sure that either NoDocks or p_inf_1 is in kwargs
    assert len(kwargs)==1, "Pass either NoDocks or p_inf_1 to function"
    ############################## Read in locs_picked.hdf5 file generated by Picasso
    locs=fifo.read_locs(path)
#    locs=locs[:][locs['frame']>999]
    # List of groups in locs
    group_list=np.unique(locs['group'])
    ############################## Read in locs_picked.yaml file generated by Picasso
    # Get TIF and Localize meta-data
    [TIFmeta,LOCmeta]=fifo.read_meta(path)
    # Number of frames in .tif stack used by Picasso Localize
    NoFrames=TIFmeta['Frames']
    ############################## Define ignore_dark for different cases
    if LOCmeta['Generated by']=='simulate_locs': # Check if locs is generated by simulate locs
        ignore_dark=0 # For simulations ignore_dark is set to zero
            
    ############################## Initiate groupprops
    groupprops=np.zeros(np.size(group_list),dtype=[('group','i4',1),
                                    ('mean_x','f4',1),('std_x','f4',1),('mean_y','f4',1),('std_y','f4',1), # Location
                                    ('mean_frame','f4',1),('std_frame','f4',1), # Frames
                                    ('mean_photons','f4',1),('std_photons','f4',1), # Photons
                                    ('n_locs','i4',1),('n_events','i4',1),('n_events_ignore','i4',1), # Statistics                                    
                                    ('tau_b','f4',1),('tau_b_ignore','f4',1),('tau_d','f4',1), # qPAINT dynamics
                                    ('tau_b_lin','f4',1),('tau_b_lin_ignore','f4',1),('tau_d_lin','f4',1), # qPAINT dynamics linearized
                                    ('ac_A','f4',1),('ac_tau','f4',1),('ac_offset','f4',1),('ac_chisquare','f4',1), # ac fit                                    
                                    ('ac_p_inf','f4',1),('ac_tau_b','f4',1),('ac_tau_d','f4',1), # ac dynamics
                                    ('NoDocks','f4',1),('err_NoDocks','f4',1)
                                    ])
    # Assign groups
    groupprops['group'][:]=group_list
    
    # Start loop over all groups
    for i in tqdm(range(0,np.size(group_list))):
        g=group_list[i]
        locs_g=locs[:][locs['group']==g]
        
        #Location
        groupprops['mean_x'][i]=get_mean_x(locs_g)
        groupprops['std_x'][i]=get_std_x(locs_g)
        groupprops['mean_y'][i]=get_mean_y(locs_g)
        groupprops['std_y'][i]=get_std_y(locs_g)
        #Frames
        groupprops['mean_frame'][i]=get_mean_frame(locs_g)
        groupprops['std_frame'][i]=get_std_frame(locs_g)
        # Photons
        groupprops['mean_photons'][i]=get_mean_photons(locs_g)
        groupprops['std_photons'][i]=get_std_photons(locs_g)
        
        # Old qPAINT dynamics
        # Binding times of 1 frame are not neglected in ECDF
        [n_events,tau_b,tau_d]=get_tau(locs_g,ignore_dark,out='param',bright_ignore=False,fit='exp')
        # Binding times of 1 frame are neglected in ECDF
        [n_events_ignore,tau_b_ignore,tau_d]=get_tau(locs_g,ignore_dark,out='param',bright_ignore=True,fit='exp')
        
        # Linearized qPAINT dynamics
        # Binding times of 1 frame are not neglected in ECDF
        [n_events,tau_b_lin,tau_d_lin]=get_tau(locs_g,ignore_dark,out='param',bright_ignore=False,fit='lin')
        # Binding times of 1 frame are neglected in ECDF
        [n_events_ignore,tau_b_lin_ignore,tau_d_lin]=get_tau(locs_g,ignore_dark,out='param',bright_ignore=True,fit='lin')
        
        # Statistics
        groupprops['n_locs'][i]=get_n_locs(locs_g)
        groupprops['n_events'][i]=n_events
        groupprops['n_events_ignore'][i]=n_events_ignore
        # qPAINT dynamics old
        groupprops['tau_b'][i]=tau_b
        groupprops['tau_b_ignore'][i]=tau_b_ignore
        groupprops['tau_d'][i]=tau_d
        # qPAINT dynamics new
        groupprops['tau_b_lin'][i]=tau_b_lin
        groupprops['tau_b_lin_ignore'][i]=tau_b_lin_ignore
        groupprops['tau_d_lin'][i]=tau_d_lin
        
        # ac fit
        [ac_A,ac_tau,ac_offset,ac_chisquare]=get_ac_fit(locs_g,NoFrames)
        groupprops['ac_A'][i]=ac_A
        groupprops['ac_tau'][i]=ac_tau
        groupprops['ac_offset'][i]=ac_offset
        groupprops['ac_chisquare'][i]=ac_chisquare
        # ac dynamics
        if 'p_inf_1' in kwargs:# Case: p_inf_1 is given 
            [NoDocks,err_NoDocks]=get_ac_NoDocks(locs_g,NoFrames,kwargs['p_inf_1'])
            groupprops['NoDocks'][i]=NoDocks
            groupprops['err_NoDocks'][i]=err_NoDocks
            mode='p_inf_1=[%e,%e]'%(kwargs['p_inf_1'][0],kwargs['p_inf_1'][1])
        elif 'NoDocks' in kwargs: # Case: NoDocks is given
            NoDocks=kwargs['NoDocks']
            err_NoDocks=0
            groupprops['NoDocks'][i]=NoDocks
            groupprops['err_NoDocks'][i]=err_NoDocks
            mode='NoDocks=%d'%(NoDocks)
           
        [ac_p_inf,ac_tau_b,ac_tau_d]=get_ac_tau(locs_g,NoFrames,NoDocks)
        groupprops['ac_p_inf'][i]=ac_p_inf
        groupprops['ac_tau_b'][i]=ac_tau_b
        groupprops['ac_tau_d'][i]=ac_tau_d
              
    # Save groupprops in hdf5 file
    groupprops_file = h5py.File(path.replace('.hdf5','_groupprops.hdf5'), "w")
    dset=groupprops_file.create_dataset("locs", np.shape(groupprops), dtype=groupprops.dtype)
    dset[...]=groupprops
    groupprops_file.close() 
    
    # Create meta_dictionary to add to yaml file
    ADDdict={'Generated by':'locs_groupprops.locs2groupprops','ignore_dark':ignore_dark,'mode':mode}
    fifo.create_meta_locs2groupprops(path,ADDdict)
       
    return groupprops

#%%
def sum2groupprops(path,NoDocks):
    # Read in locs_picked.hdf5 file generated by Picasso
    locs=fifo.read_sum(path)
    # Read in locs_picked.yaml file generated by Picasso
    # Get TIF and Localize meta-data
    [TIFmeta,LOCmeta]=fifo.read_meta(path)
    # Number of frames in .tif stack used by Picasso Localize
    NoFrames=TIFmeta['Frames']
    # List of groups in locs
    group_list=np.unique(locs['n_id'])
    # Initiate groupprops
    groupprops=np.zeros(np.size(group_list),dtype=[('group','i4',1),
                                    ('mean_x','f4',1),('std_x','f4',1),('mean_y','f4',1),('std_y','f4',1), # Location
                                    ('mean_frame','f4',1),('std_frame','f4',1), # Frames
                                    ('mean_photons','f4',1),('std_photons','f4',1), # Photons
                                    ('ac_A','f4',1),('ac_tau','f4',1),('ac_offset','f4',1), # ac fit
                                    ('ac_p_inf','f4',1),('ac_tau_b','f4',1),('ac_tau_d','f4',1),('NoDocks','f4',1) # ac dynamics
                                    ])
    # Assign groups
    groupprops['group'][:]=group_list
    
    # Start loop over all groups
    for i in tqdm(range(0,np.size(group_list))):
        g=group_list[i]
        locs_g=locs[:][locs['n_id']==g]
        
        #Location
        groupprops['mean_x'][i]=get_mean_x(locs_g)
        groupprops['std_x'][i]=get_std_x(locs_g)
        groupprops['mean_y'][i]=get_mean_y(locs_g)
        groupprops['std_y'][i]=get_std_y(locs_g)
        #Frames
        groupprops['mean_frame'][i]=get_mean_frame(locs_g)
        groupprops['std_frame'][i]=get_std_frame(locs_g)
        # Photons
        groupprops['mean_photons'][i]=get_mean_photons(locs_g)
        groupprops['std_photons'][i]=get_std_photons(locs_g)
        
        # ac fit
        [ac_A,ac_tau,ac_offset]=get_ac_fit(locs_g,NoFrames)
        groupprops['ac_A'][i]=ac_A
        groupprops['ac_tau'][i]=ac_tau
        groupprops['ac_offset'][i]=ac_offset
        
        # ac dynamics
        [ac_p_inf,ac_tau_b,ac_tau_d]=get_ac_tau(locs_g,NoFrames,NoDocks)
        groupprops['ac_p_inf'][i]=ac_p_inf
        groupprops['ac_tau_b'][i]=ac_tau_b
        groupprops['ac_tau_d'][i]=ac_tau_d
        groupprops['NoDocks'][i]=NoDocks
        
        
    # Save groupprops in hdf5 file
    groupprops_file = h5py.File(path.replace('.hdf5','_groupprops.hdf5'), "w")
    dset=groupprops_file.create_dataset("groupprops", np.shape(groupprops), dtype=groupprops.dtype)
    dset[...]=groupprops
    groupprops_file.close() 
    
    # Create meta_dictionary to add to yaml file
    ADDdict={'Generated by':'locs_groupprops.sum2groupprops'}
    fifo.create_meta_sum2groupprops(path,ADDdict)
       
    return groupprops