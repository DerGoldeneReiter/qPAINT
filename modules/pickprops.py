# Import modules
import numpy as np
import scipy.optimize
import importlib
import pandas as pd
from tqdm import tqdm
import dask.dataframe as dd
# Import own modules
import varfuncs
importlib.reload(varfuncs)
import multitau
importlib.reload(multitau)

#%%
def get_ac(df,NoFrames):
    """ 
    1) Convert localizations of single pick to trace
        - trace(frame)=0 if no localization in frame
        - else trace(frame)=photons
        - length of trace will be NoFrames
    2) Compute normalized multi-tau autocorrelation function of trace -> help(multitau.autocorrelate)
    3) Least square fit of function f(tau)=mono_A*exp(-tau/mono_tau)+1 to normalized autocorrelation function -> help(pickprops.fit_ac_single)
    
    Parameters
    ---------
    df : pandas.DataFrame
        Picked localizations for single group. Required fields are 'frame' and 'photons'
    NoFrames : int
        Number of frames of localized image stack
    
    Returns
    -------
    s_out : pandas.Series
        Length: 6
        Column:
            'group' : int
        Index: 
            'trace' : numpy.ndarray
                trace
            'tau' : numpy.ndarray
                tau of autocorr.
            'g' : numpy.ndarray
                g(tau) autocorrelation value
            'mono_A' : float64
                Amplitude of monoexponential fit function
            'mono_tau' : float64
                Correlation time of monoeponential fit function
            'mono_chi' : float64
                Square root of chisquare value of fit with sigma=1 for all data points
    """
    ###################################################### Prepare trace
    # Get absolute values of photons, since sometimes negative values can be found
    df['photons']=df['photons'].abs() 
    # Sum multiple localizations in single frame
    df_sum=df[['frame','photons']].groupby('frame').sum()
    # Define trace of length=NoFrames with zero entries
    trace=np.zeros(NoFrames)
    # Add (summed) photons to trace for each frame
    trace[df_sum.index.values]=df_sum['photons'].values
    
    ###################################################### Generate autocorrelation
    ac=multitau.autocorrelate(trace,m=16, deltat=1,
                                 normalize=True,copy=False, dtype=np.float64())
    
    ###################################################### Fit mono exponential decay to autocorrelation
    mono_A,mono_tau,mono_chi=fit_ac_mono(ac) # Get fit
    
    ###################################################### Assignment to series 
    s_out=pd.Series({'trace':trace,
                     'tau':ac[:,0],'g':ac[:,1], # Autocorrelation function
                     'mono_A':mono_A,'mono_tau':mono_tau,'mono_chi':mono_chi}) # mono exponential fit results
    
    return s_out

#%%
def fit_ac_mono(ac):
    """ 
    Least square fit of function f(tau)=mono_A*exp(-tau/mono_tau)+1 to normalized autocorrelation function.
    
    Parameters
    ---------
    ac : numpy.ndarray
        1st column should correspond to delay time tau of autocorrelation function.
        2nd column should correspond to value g(tau) of autocorrelation function
    
    Returns
    -------
    mono_A : float64
        Amplitude of monoexponential fit function
    mono_tau : float64
        Correlation time of monoeponential fit function
    mono_chi : float64
        Chisquare value of fit with sigma=1 for all data points   
    """
    ###################################################### Define start parameters
    p0=np.empty([2])
    p0[0]=ac[1,1] # Amplitude
    halfvalue=1.+(p0[0]-1.)/2 # Value of half decay of ac
    p0[1]=np.argmin(np.abs(ac[:,1]-halfvalue)) # tau
    ###################################################### Fit boundaries
    lowbounds=np.array([0,0])
    upbounds=np.array([np.inf,np.inf])
    ###################################################### Fit data
    try:
        popt,pcov=scipy.optimize.curve_fit(varfuncs.ac_monoexp,ac[:,0],ac[:,1],p0,bounds=(lowbounds,upbounds),method='trf')
    except RuntimeError:
        popt=p0
    except ValueError:
        popt=p0
    except TypeError:
        popt=p0
    
    ###################################################### Calculate chisquare    
    chisquare=np.sum(np.square(varfuncs.ac_monoexp(ac[:,0],*popt)-ac[:,1]))/(len(ac)-2)
    
    return popt[0],popt[1],np.sqrt(chisquare)

#%%
def get_tau(df,ignore=1):
    """ 
    1) Computes bright, dark-time distributions and number of events a la Picasso
    2) Least square fit of function f(tau)=1-exp(-t/tau) to experimental continuous distribution function (ECDF) 
        -> help(pickprops.fit_tau)
    
    Parameters
    ---------
    df : pandas.DataFrame
        Picked localizations for single group. Required columns are 'frame'.
    ignore : int
        Disrupted binding events by duration ignore will be treated as single event with bright time of 
        total duration of bridged events. Defaults to ignore=1.
        
    
    Returns
    -------
    s_out : pandas.Series
        Length: 9
        Column:
            'group' : int
        Index: 
            'tau_b_dist' : numpy.ndarray
                Distribution of bright times with ignore value taken into account
            'tau_b' : float64
                Fit result of exponential function to bright time ECDF
            'tau_b_lin' : float64
                Fit result of line with offset=0 to linearized bright time ECDF given by -ln(1-ECDF)
            'tau_d_dist' : numpy.ndarray
                Distribution of dark times with ignore value taken into account
            'tau_b' : float64
                Fit result of exponential function to dark time ECDF
            'tau_b_lin' : float64
                Fit result of line with offset=0 to linearized dark time ECDF given by -ln(1-ECDF)
            'n_events' : float64
                Number of binding events
            'mean_event_times' : float64
                Mean of event times
            'std_event_times' : float64
                Standard deviation of event times
    """
    
    frames=df['frame'].values # Get sorted frames as numpy.ndarray
    frames.sort()
    ################################################################ Get tau_d distribution
    dframes=frames[1:]-frames[0:-1] # Get frame distances i.e. dark times
    dframes=dframes.astype(float) # Convert to float values for later multiplications
    tau_d_dist=dframes[dframes>(ignore+1)] # Remove all dark times>ignore+1
    tau_d_dist=np.sort(tau_d_dist) # Sorted tau_d distribution
 
    ################################################################# Get tau_b_distribution
    dframes[dframes<=(ignore+1)]=0 # Set (bright) frames to 0 that have nnext neighbor distance <= ignore+1
    dframes[dframes>1]=1 # Set dark frames to 1
    dframes[dframes<1]=np.nan # Set bright frames to NaN
    
    mask_end=np.concatenate([dframes,[1]],axis=0) # Mask for end of events, add 1 at end
    frames_end=frames*mask_end # Apply mask to frames to get end frames of events
    frames_end=frames_end[~np.isnan(frames_end)] # get only non-NaN values, removal of bright frames
    
    mask_start=np.concatenate([[1],dframes],axis=0) # Mask for start of events, add one at start
    frames_start=frames*mask_start # Apply mask to frames to get start frames events
    frames_start=frames_start[~np.isnan(frames_start)] # get only non-NaN values, removal of bright frames
    
    tau_b_dist=frames_end-frames_start+1 # get tau_b distribution
    tau_b_dist=np.sort(tau_b_dist) # sort tau_b distribution
    
    ################################################################# Number of events and their timing
    n_events=float(np.size(tau_b_dist)) # Number of binding events
    event_times=(frames_start+frames_end)/2 # Events time distribution over trace
    mean_event_times=np.mean(event_times) # Get mean of event times
    std_event_times=np.std(event_times) # Get std fo event times

    ################ Extract tau's
    if n_events<=5: # If n_events <= 5 --> Set all parameters to mean of distribution
        # Bright time
        tau_b=np.mean(tau_b_dist)
        tau_b_lin=np.mean(tau_b_dist)
        # Dark time
        tau_d=np.mean(tau_d_dist)
        tau_d_lin=np.mean(tau_d_dist)
        
    elif n_events>5: # If n_events > 5 --> Fitting of ECDF with exponential function, calculate median
        # Bright time
        tau_b,tau_b_lin=fit_tau(tau_b_dist)
        # Dark time
        tau_d,tau_d_lin=fit_tau(tau_d_dist)

    ###################################################### Assignment to series 
    s_out=pd.Series({'tau_b_dist':tau_b_dist,'tau_b':tau_b,'tau_b_lin':tau_b_lin, # Bright times
                     'tau_d_dist':tau_d_dist,'tau_d':tau_d,'tau_d_lin':tau_d_lin, # Dark times
                     'n_events':n_events,'mean_event_times':mean_event_times,'std_event_times':std_event_times}) # Events and timing
    return s_out    
#%%     
def fit_tau(tau_dist):
    """ 
    Least square fit of function f(t)=1-exp(-t/tau) to experimental continuous distribution function (ECDF) of tau_dist
    -> help(varfuncs.get_ecdf) equivalent to:
        matplotlib.pyplot.hist(tau_dist,bins=numpy.unique(tau_dist),normed=True,cumulative=True)
    
    Parameters
    ---------
    tau_dist : numpy.ndarray
        1 dimensional array of bright or dark times t   
    Returns
    -------
    tau : float64
        tau as obtained by fitting f(t) to ECDF(t) with least squares. 
    tau_lin : float64
        tau as obtained by fitting f_lin(t)=t/tau to linearized data -ln(1-ECDF(t)) with least squares. 
        Last value of ECDF is omitted due to usage of logarithm.
    """
    ####################################################### Get ECDF
    tau_bins,tau_ecdf=varfuncs.get_ecdf(tau_dist)
    
    ####################################################### Define start paramter
    idx=np.abs(tau_ecdf-(1-np.exp(-1))).argmin() # Find idx in ecfd closest to 1-e^-1 value
    p0=tau_bins[idx] # Start value = bin at which cfd closest to 1-e^-1 value
    
    ####################################################### Fit ECDF
    try:
        # Fit exponential directly
        tau,pcov=scipy.optimize.curve_fit(varfuncs.ecdf_exp,tau_bins,tau_ecdf,p0) 
        tau=tau[0]
        # Fit linearized data, omitt last value of ECDF du to linearization with log function
        tau_lin,pcov=scipy.optimize.curve_fit(varfuncs.ecdf_exp_lin,tau_bins[0:-1],-np.log(np.abs(1-tau_ecdf[0:-1])),p0)
        tau_lin=tau_lin[0]
    except RuntimeError:
        tau=p0
        tau_lin=p0
    except ValueError:
        tau=p0
        tau_lin=p0
    except TypeError:
        tau=p0
        tau_lin=p0
        
    return tau,tau_lin

#%%
def get_other(df):
    """ 
    Get mean and std values for a single group.
    
    Parameters
    ---------
    df : pandas.DataFrame
        Picked localizations for single group. Required columns are 'frame'.

    Returns
    -------
    s_out : pandas.Series
        Length: 6
        Column:
            'group' : int
        Index: 
            'mean_frame' : flaot64
                Mean of frames for all localizations in group
            'x' : float64
                Mean x position
            'y' : float64
                Mean y position
            'mean_photons' : float64
                Mean of photons for all localizations in group
            'std_frame' : flaot64
                Standard deviation of frames for all localizations in group
            'std_photons' : flaot64
                Standar deviation of photons for all localizations in group
    """
    # Get mean values
    s_mean=df[['frame','x','y','photons']].mean()
    mean_idx={'frame':'mean_frame','x':'x','y':'y','photons':'mean_photons'}
    # Get std values
    s_std=df[['frame','photons']].std()
    std_idx={'frame':'std_frame','photons':'std_photons'}
    # Combine output
    s_out=pd.concat([s_mean.rename(mean_idx),s_std.rename(std_idx)])
    
    return s_out

#%%
def get_props(df,NoFrames,ignore):
    """ 
    Wrapper function to combine:
        - pickprops.get_ac(df,NoFrames)
        - pickprops.get_tau(df,ignore)
        - pickprops.get_other(df)
    
    Parameters
    ---------
    df : pandas.DataFrame
        'locs' of locs_picked.hdf5 as given by Picasso
    NoFrames: int
        Number of frames of image stack corresponding to locs_picked.hdf5
    ignore: int
        Ignore as defined in props.get_tau
    Returns
    -------
    s : pandas.DataFrame
        Columns as defined in individual functions get_ac, get_tau. Index corresponds to 'group'.
    """
    
    # Call individual functions   
    s_ac=get_ac(df,NoFrames)
    s_tau=get_tau(df,ignore)
    s_other=get_other(df)
    # Combine output
    s_out=pd.concat([s_ac,s_tau,s_other])
    
    return s_out

#%%
def apply_props(df,NoFrames,ignore): 
    """
    Applies pick_props.get_props(df,NoFrames,ignore) to each group in non-parallelized manner. Progressbar is shown under calculation.
    """
    tqdm.pandas() # For progressbar under apply
    df_props=df.groupby('group').progress_apply(lambda df: get_props(df,NoFrames,ignore))
    
    return df_props

#%%
def apply_props_dask(df,NoFrames,ignore,NoPartitions): 
    """
    Applies pick_props.get_props(df,NoFrames,ignore) to each group in parallelized manner using dask by splitting df into 
    various partitions. No progressbar is shown.
    """
    ########### Partionate df using dask for parallelized computation
    df=dd.from_pandas(df,npartitions=NoPartitions) 
    ########### Define apply_props for dask which will be applied to different partitions of df
    def apply_props_2part(df,NoFrames,ignore): return df.groupby('group').apply(lambda df: get_props(df,NoFrames,ignore))
    ########### Map apply_props_2part to every partition of df for parallelized computing
    df_props=df.map_partitions(apply_props_2part,NoFrames,ignore).compute(scheduler='processes')
    
    return df_props













