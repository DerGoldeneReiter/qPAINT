# Import modules
import numpy as np
import pandas as pd

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




