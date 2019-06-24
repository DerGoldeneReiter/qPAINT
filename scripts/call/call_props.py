#Call script to module: pickprops
#
############################################################# Set parameters
NoPartitions=30

conc=[2.5,5,10,20] # Imager concentration [nM]
ignore=0
savename_ext='_props_ig%i'%(ignore)
omit_dist=True
kin_filter=False

###### Dictonary content for .yaml file
props_dict={'Generated by':'pickprops.get_props',
            'ignore':ignore,
            'omit_dist':omit_dist,
            'kin_filter': kin_filter}

############################################################## Define data
dir_names=[]
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/z.simulations/19-06-19_Pm2_2B07/N3']*4)

file_names=[]
file_names=[]
file_names.extend(['N3_2-5nM_locs_picked.hdf5'])
file_names.extend(['N3_5nM_locs_picked.hdf5'])
file_names.extend(['N3_10nM_locs_picked.hdf5'])
file_names.extend(['N3_20nM_locs_picked.hdf5'])

#################################################### Load packages
import os #platform independent paths
import importlib
# Load user defined functions
import pickprops as props
import var_io as io
import pickprops_calls as props_call
# Reload modules
importlib.reload(props)
importlib.reload(io)
importlib.reload(props_call)
#%%
############################################################# Read locs, apply props & save locs
######### Create list of paths
path=[os.path.join(dir_names[i],file_names[i]) for i in range(0,len(file_names))]

######### Read-Apply-Save loop
for i in range(0,len(path)):
    ######### File read in
    locs,locs_yaml=io.read_locs(path[i])
    
    ######### Get number of frames
    NoFrames=locs_yaml[0]['Frames']
    
    ######### Apply non-parallelized props
#    locs_props=props.apply_props(locs,conc[i],NoFrames,ignore)
    
    ######### Apply parallelized props
    locs_props=props.apply_props_dask(locs,conc[i],NoFrames,ignore,NoPartitions)
    
    ######### Drop objects for saving if omit=True
    if omit_dist:
        print('... removing distributionjs from output')
        locs_props=locs_props.drop(['trace','tau','g','tau_b_dist','tau_d_dist'],axis=1)
    
    if kin_filter:
        print('... applying kinetic filter')
        locs_props=props._kin_filter(locs_props)
    
    ######### Add nearest neigbour pick and distance
    print('... calculating nearest neighbour')
    locs_props=props_call.props_add_nn(locs_props)
    
    ######### Save .hdf5 and .yaml of locs_props    
    io.save_locs(locs_props,[locs_yaml,props_dict],path[i],savename_ext)
    

#%%
#import matplotlib.pyplot as plt
#import varfuncs
#import numpy as np
#
#f=plt.figure(num=11,figsize=[4,3])
#f.subplots_adjust(bottom=0.1,top=0.99,left=0.2,right=0.99)
#f.clear()
#
##### Autocorrelation
#ax=f.add_subplot(311)
#for g in [364]:
#    print(locs_props.loc[g,'mono_tau'],locs_props.loc[g,'mono_tau_lin'])
#    ax.plot(locs_props.loc[g,'tau'],
#            varfuncs.ac_monoexp(locs_props.loc[g,'tau'],locs_props.loc[g,'mono_A'],locs_props.loc[g,'mono_tau']),
#            '-',lw=2,c='r')
#    ax.plot(locs_props.loc[g,'tau'],
#            varfuncs.ac_monoexp(locs_props.loc[g,'tau'],locs_props.loc[g,'mono_A_lin'],locs_props.loc[g,'mono_tau_lin']),
#            '-',lw=2,c='b')
#    ax.plot(locs_props.loc[g,'tau'],locs_props.loc[g,'g'],'*')
#    ax.axhline(1,ls='--',lw=2,color='k')
#ax.set_xscale('symlog')
##### Trace
#ax=f.add_subplot(312)
#ax.plot(locs_props.loc[g,'trace'])
##### tau_d_dist
#ax=f.add_subplot(313)
#x=varfuncs.get_ecdf(locs_props.loc[g,'tau_d_dist'])[0]
#y=varfuncs.get_ecdf(locs_props.loc[g,'tau_d_dist'])[1]
#x_fit=np.arange(0,np.max(x),0.1)
#ax.plot(x,y)
#ax.plot(x_fit,varfuncs.ecdf_exp(x_fit,locs_props.loc[g,'tau_d']))

