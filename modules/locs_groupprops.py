
import numpy as np

def get_mean_x(locs):
    x=np.mean(locs['x'],axis=0)
    return x
def get_mean_y(locs):
    y=np.mean(locs['y'],axis=0)
    return y
def get_mean_frame(locs):
    frame=np.mean(locs['frame'],axis=0)
    return frame
def get_mean_photons(locs):
    photons=np.mean(locs['photons'],axis=0)
    return photons

def get_std_frame(locs):
    frame=np.std(locs['frame'],axis=0)
    return frame

def get_G0(locs,NoFrames):
    G0=np.sum(np.square(locs['photons']),axis=0)*(1/NoFrames)
    return G0

def get_n_locs(locs):
    n_locs=np.size(locs['frame'])
    return n_locs

def get_n_events(locs,ignore):
    dframes=locs['frame'][1:-1]-locs['frame'][0:-2]
    tau_d=dframes[dframes>(ignore+1)]
    n_events=np.size(tau_d)+1
    return n_events

def get_tau(locs,ignore,max_n_events):
    tau_d_unisize=np.zeros([1,max_n_events]) # Create uniform sized tau distributons for storing
    tau_b_unisize=np.zeros([1,max_n_events]) #
    
    dframes=locs['frame'][1:]-locs['frame'][0:-1]
    dframes=np.float64(dframes)    
    tau_d=dframes[dframes>(ignore+1)] # tau_d distribution
    
    
    tau_d_unisize[0,:np.size(tau_d)]=tau_d # tau_d distribution with zeros padded to get to length no max_n_events
    
    # Get tau_b
    dframes[dframes<=(ignore+1)]=0 # Set bright frames to 0
    dframes[dframes>1]=1 # Set dark-frames to 1
    dframes[dframes<1]=np.nan # Set bright frames to 0
    
    mask_end=np.concatenate([dframes,[1]],axis=0) # Mask for end of events
    frames_end=locs['frame']*mask_end # Apply mask to frames to get end frames of events
    frames_end=frames_end[~np.isnan(frames_end)]
    
    mask_start=np.concatenate([[1],dframes],axis=0) # Mask for start of events
    frames_start=locs['frame']*mask_start # Apply mask to frames to get start frames events
    frames_start=frames_start[~np.isnan(frames_start)]
    
    
    tau_b=frames_end-frames_start+1 #tau_b distribution
    tau_b_unisize[0,:np.size(tau_b)]=tau_b # tau_b distribution with zeros padded to get to length no max_n_events
    
    return [tau_d_unisize,tau_b_unisize]