#################################################### Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os #platform independent paths
import pandas as pd 
import scipy.stats as stats
# Load user defined functions
import var_io as io
# Load plt style
plt.style.use('~/qPAINT/styles/kinetics.mplstyle')

############################################################## Pre definitions
#### Create empty filter DataFrame
pre=pd.DataFrame(columns=['low_std_frame','low_mean_frame','up_mean_frame','up_mono_tau','up_mono_A','low_kde'])
#### Create filter DataFrame
#### Used fields for scatter plot and kde filtering
field1='mono_tau'
field2='mono_A'

############################################################## Experimental ID, filter settings and save paths
#### Working name for this evaluation
expID='YourID'

pre.loc[expID,'low_std_frame']=12e3
pre.loc[expID,'low_mean_frame']=20e3
pre.loc[expID,'up_mean_frame']=30e3
pre.loc[expID,'up_mono_tau']=50
pre.loc[expID,'up_mono_A']=8
pre.loc[expID,'low_kde']=0.2

savepath='/fs/pool/pool-schwille-paint/Analysis/z.PAINT-checks/18-12-04_kinSDS_test'

############################################################## Define data paths
dir_names=[]
dir_names.extend(['/fs/pool/pool-schwille-paint/Data/z.PAINT-checks/18-11-27_Pm_Pmr_check_FSc/SDS-Pm2-8nt-2xT_Im-8nt-20nM_p100mW-9deg!!!/18-11-29_FS'])
# Define names of locs_picked.hdf5 file
file_names=[]
file_names.extend(['SDS-Pm2-8nt-2xT_Im-8nt-20nM_p100mW-9deg_1_MMStack.ome_locs_render_picked_1985_d1_props-4.hdf5'])

############################################################## Read in data and prep
#### Create list of paths
path=os.path.join(dir_names[0],file_names[0])
#### Read in locs of path and corresponding yaml
locs_props,locs_props_yaml=io.read_locs(path)
#### Add expID to locs_props
locs_props['expID']=expID
############################################################# Apply pre-filters
#### Copy locs props to X for filtering
X=locs_props.copy()
#### Remove NaNs
X.dropna(axis=0,inplace=True)
#### Make pre-filters happen!
X.drop(X[(X.std_frame<pre.loc[expID,'low_std_frame'])&(X.expID==expID)].index,inplace=True)
X.drop(X[((X.mean_frame<pre.loc[expID,'low_mean_frame'])|(X.mean_frame>pre.loc[expID,'up_mean_frame']))&(X.expID==expID)].index,inplace=True)
X.drop(X[(X.mono_tau>pre.loc[expID,'up_mono_tau'])].index,inplace=True)
X.drop(X[(X.mono_A>pre.loc[expID,'up_mono_A'])&(X.expID==expID)].index,inplace=True)
############################################################# Apply kernel density filtering (kde)
xy=X.loc[X.expID==expID,[field1,field2]].values.transpose()
kde=stats.gaussian_kde(np.log(xy))(np.log(xy))
kde=kde/np.max(kde)
X.loc[X.expID==expID,'kde']=kde
del xy, kde
#### Make kde-filter happen!
X.drop(X[(X.kde<pre.loc[expID,'low_kde'])&(X.expID==expID)].index,inplace=True)
############################################################ Prepare results
#### Create results Series
results=pd.DataFrame(index=['mono_A','mono_tau','tau_b_ac','tau_d_ac','mean_photons','bg'])
#### Add column with tau_d_ac, tau_b_ac to X
X.loc[X.expID==expID,'tau_b_ac']=(X.mono_A+1)*(X.mono_tau/X.mono_A)
X.loc[X.expID==expID,'tau_d_ac']=(X.mono_A+1)*(X.mono_tau)
#### Save means in results DataFrame
for i in results.index: results.loc[i,'mean']=X.loc[X.expID==expID,i].mean()
for i in results.index: results.loc[i,'err_mean']=X.loc[X.expID==expID,i].std()/np.sqrt(len(X))
del i
results.to_csv(savepath+'/kinSDS_results.csv',index_label=expID)
############################################################# Plot distributions
#### Show effect of kde-filter in scatter plot
f=plt.figure(num=1,figsize=[4,3])
f.clear()
f.subplots_adjust(top=0.96,bottom=0.17,left=0.14)
ax=f.add_subplot(111)
plth=ax.scatter(X.loc[X.expID==expID,field1],X.loc[X.expID==expID,field2],c=X.loc[X.expID==expID,'kde'])
plt.colorbar(plth)
ax.set_xlabel(field1)
ax.set_ylabel(field2)
plt.savefig(savepath+'/kinSDS_scatter.png')
#### Show effect of kde-filter in histograms
f=plt.figure(num=2,figsize=[8,4])
f.subplots_adjust(top=0.93,bottom=0.06,left=0.055,right=0.995)
ax=f.add_subplot(111)
f.clear()
X.loc[X.expID==expID,['std_frame','mean_frame','mono_tau','mono_A','mean_photons','bg','tau_b_ac','tau_d_ac']].hist(bins='fd',ax=ax,layout=(2,4))
plt.savefig(savepath+'/kinSDS_hist.png')
