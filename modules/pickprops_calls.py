# Import modules
import numpy as np
import pandas as pd
import importlib

#%%
def group_groups(df,supsize,mode='concat'):
    
    import multitau
    import pickprops as props
    
    #### Function definitions 
    def concat_groups_kinetics(df):
        s_out=pd.Series({i:np.hstack(df[i].values) for i in df.columns})
        #### Drop 'sup_group'
        s_out=s_out.drop('sup_group')
        
        #### Apply AC functions
        ac=multitau.autocorrelate(s_out['trace'],m=16, deltat=1,
                                 normalize=True,copy=False, dtype=np.float64())
        mono_A,mono_tau,mono_chi=props.fit_ac_mono(ac) # Get fit
        
        #### Apply ECDF functions
        tau_b,tau_b_lin=props.fit_tau(s_out['tau_b_dist'])
        tau_d,tau_d_lin=props.fit_tau(s_out['tau_d_dist'])
        
        #### Assign new values
        s_out['tau']=ac[:,0]
        s_out['g']=ac[:,1]
        s_out['mono_A']=mono_A
        s_out['mono_tau']=mono_tau
        
        s_out['tau_b']=tau_b
        s_out['tau_b_lin']=tau_b_lin
        s_out['tau_d']=tau_d
        s_out['tau_d_lin']=tau_d_lin
        
        return s_out  
    
    #### Create sup-group index. Partionate groups into sup_groups of size=supsize. 
    #### Sup-group index will start with 1. Last group
    N=len(df)
    sup_groups=np.repeat(np.arange(1,int(N/supsize)+1),supsize)
    #### Assing group index of zero to remaining groups
    sup_groups.resize(N)
    
    #### Copy part of df and add sup-groups column, remove sup_groups samller tha required size
    df_sup=df[['group','trace','tau_b_dist','tau_d_dist']]
    df_sup['sup_group']=sup_groups
    df_sup=df_sup[df_sup['sup_group']>0]
    
    #### Apply functions
    df_out=df_sup.groupby('sup_group').apply(concat_groups_kinetics)
   
    return df_out

#%%
def segment_time(path,noFrames_seg):
    
    #### Load modules
    import var_io
    importlib.reload(var_io)
    
    #### Load locs and yaml
    locs,meta=var_io.read_locs(path)
    noFrames=meta[0]['Frames'] # Number of frames in locs
    #### Define time segments
    noSegments=np.floor(noFrames/noFrames_seg).astype(int) # no. of segments
    start_frame=[i*noFrames_seg for i in range(0,noSegments)] # start frames of new segments
    end_frame=[i*noFrames_seg for i in range(1,noSegments+1)] # end frames of new segment
    #### Assign noFrames_new to meta_seg for saving
    meta_seg=meta
    meta_seg[0]['Frames']=noFrames_seg
    
    #### Split locs file according to time segments
    for i in range(0,noSegments):
        print('_%i_of_%i'%(i+1,noSegments))
        #### Split locs
        locs_seg=locs.loc[(locs.frame>=start_frame[i])&(locs.frame<end_frame[i])]
        #### Substract start_frame
        locs_seg.loc[:,'frame']=locs_seg.loc[:,'frame'].values-start_frame[i]
        var_io.save_locs(locs_seg,meta_seg,path,savename_ext='_%i_of_%i'%(i+1,noSegments))
        
    return 

