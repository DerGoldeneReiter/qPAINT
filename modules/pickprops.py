# Import modules
import numpy as np
import scipy.optimize
import importlib
import pandas as pd

# Import own modules in qPAINT
import fitfunc
importlib.reload(fitfunc)
import multitau
importlib.reload(multitau)

#%%
def get_ac(df,NoFrames):
    """ 
    1) Convert localizations of single pick to trace
        - trace(frame)=0 if no localization in frame
        - else trace(frame)=photons
        - length of trace will be NoFrames
    2) Compute normalized multi-tau autocorrelation function of trace (according to custom multitau function)
    
    Parameters
    ---------
    df : pandas.DataFrame
        Picked localizations for single group. Required fields are 'frame' and 'photons'
    NoFrames : int
        Number of frames of localized image stack
    
    Returns
    -------
    s : pandas.Series
        Length: 2
        Column:
            'group' : int
        Index: 
            'trace' : numpy.ndarray
                trace
            'tau' : numpy.ndarray
                tau of autocorr.
            'g' : numpy.ndarray
                g(tau) autocorrelation value
    """
    ###################################################### Prepare trace
    # Sum multiple localizations in single frame
    df_sum=df[['frame','photons']].groupby('frame').sum()
    # Define trace of length=NoFrames with zero entries
    trace=np.zeros(NoFrames)
    # Add (summed) photons to trace for each frame
    trace[df_sum.index.values]=df_sum['photons'].values
    
    ###################################################### Generate autocorrelation
    ac=multitau.autocorrelate(trace,m=16, deltat=1,
                                 normalize=True,copy=False, dtype=np.float64())
    
    ###################################################### Fit single exponential decay to autocorrelation
    mono_A,mono_tau=fit_ac_single(ac)
    
    ###################################################### Assignment to series 
    s_out=pd.Series({'trace':trace,
                     'tau':ac[:,0],'g':ac[:,1], # Autocorrelation function
                     'mono_A':mono_A,'mono_tau':mono_tau}) # mono exponential fit results
    
    return s_out

#%%
def fit_ac_single(ac):
    
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
        popt,pcov=scipy.optimize.curve_fit(fitfunc.ac_monoexp,ac[:,0],ac[:,1],p0,bounds=(lowbounds,upbounds),method='trf')
    except RuntimeError:
        popt=p0
    except ValueError:
        popt=p0
    except TypeError:
        popt=p0
        
    chisquare=np.sum(np.square(np.divide(fitfunc.ac_monoexp(ac[:,0],*popt)-ac[:,1],popt[0]+1.)))/(len(ac)-len(popt))
    popt=np.append(popt,np.sqrt(chisquare))
    
    return popt[0],popt[1]