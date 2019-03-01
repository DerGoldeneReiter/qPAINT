#################################################### Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os #platform independent paths
import pandas as pd 
import importlib
# Load user defined functions
import var_io as io
import kinetics
importlib.reload(kinetics)
plt.style.use('~/qPAINT/styles/FoM.mplstyle')

############################################################## Parameters & labels
#### Aquistion cycle time (s)
CycleTime=0.224 
#### Concentrations (nM) and repetitons
cs=[2.5]*3+[5]*3+[10]*3+[20]*3+[30]*3+[50]*3
rs=[1,2,3]*6
labels=['%4.1fnM-%i'%(c,r) for c,r in zip(cs,rs)]
#### Percentile for filtering (%)
q=5
#### Statistics used for fitting
use_tau_stat='mean'
use_A_stat='median'
#### Define outliers
outliers=[]
outliers=['50.0nM-1','50.0nM-2','50.0nM-3']
#### Saving
savedir='/fs/pool/pool-schwille-paint/Analysis/p04.lb-FCS/zz.Pm2-8nt/z.sum/data/'
savename='N01_Gel_B_01'
############################################################## Define data
dir_names=[]
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/19-02-14_c-series_SDS-Pm2-8nt/SDS-Pm2-8nt-GEL_c2-5nM_1/19-02-16_JS']*3)
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/19-02-14_c-series_SDS-Pm2-8nt/SDS-Pm2-8nt-GEL_c5nM_1/19-02-15_FS']*3)
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/19-02-14_c-series_SDS-Pm2-8nt/SDS-Pm2-8nt-GEL_c10nM_1/19-02-14_FS']*3)
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/19-02-14_c-series_SDS-Pm2-8nt/SDS-Pm2-8nt-GEL_c20nM_1/19-02-15_FS']*3)
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/19-02-14_c-series_SDS-Pm2-8nt/SDS-Pm2-8nt-GEL_c30nM_1/19-02-15_FS']*3)
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/p04.lb-FCS/19-02-14_c-series_SDS-Pm2-8nt/SDS-Pm2-8nt-GEL_c50nM_1/19-02-15_FS']*3)

file_names=[]
file_names.extend(['SDS-Pm2-8nt-GEL_c2-5nM_1_MMStack_Pos0.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c2-5nM_1_MMStack_Pos1.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c2-5nM_1_MMStack_Pos2.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c5nM_1_MMStack_Pos0.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c5nM_1_MMStack_Pos1.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c5nM_1_MMStack_Pos2.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c10nM_1_MMStack_Pos0.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c10nM_1_MMStack_Pos1.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c10nM_1_MMStack_Pos2.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c20nM_1_MMStack_Pos0.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c20nM_1_MMStack_Pos1.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c20nM_1_MMStack_Pos2.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c30nM_1_MMStack_Pos0.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c30nM_1_MMStack_Pos1.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c30nM_1_MMStack_Pos2.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c50nM_1_MMStack_Pos0.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c50nM_1_MMStack_Pos1.ome_locs_picked_props_ig1.hdf5'])
file_names.extend(['SDS-Pm2-8nt-GEL_c50nM_1_MMStack_Pos2.ome_locs_picked_props_ig1.hdf5'])

############################################################## Read in data
#### Create list of paths
path=[os.path.join(dir_names[i],file_names[i]) for i in range(0,len(file_names))]
#### Read in locs of path
locs_props=pd.concat([io.read_locs(p)[0] for p in path],keys=labels,names=['expID'])
locs_props.drop(outliers,inplace=True)
############################################################## Filter
X=kinetics._filter(locs_props,q)

############################################################## Statistics of observables
X_stats=kinetics._stats(X,CycleTime)

############################################################## Concentration series fitting
X,X_stats,X_fit=kinetics._fit_conc(X,X_stats,use_tau_stat,use_A_stat)

############################################################## Saving
X_stats.to_hdf(os.path.join(savedir,savename+'_stats.h5'),key='stats')
X_fit.to_hdf(os.path.join(savedir,savename+'_fit.h5'),key='fit')
############################################################## Plotting

#%%
#### Show histogram for expID
expID='30.0nM-1'
f=plt.figure(num=10,figsize=[7,8])
f.subplots_adjust(left=0.1,right=0.99,bottom=0.04,top=0.95)
ax=f.add_subplot(111)
f.clear()
X.loc[expID,['std_frame','mean_frame',
                      'mono_tau','mono_A',
                      'n_events','n_ac']].hist(bins='fd',ax=ax)

#### Kinetics plot
kinetics._plot(X_stats,X_fit,use_tau_stat,use_A_stat)
