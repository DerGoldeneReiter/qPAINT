# Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting

import os #platform independent paths
import importlib
# Load user defined functionsf
import locs_groupprops as l2grp
import file_formats as fifo
import fitfunc
# Reload modules
importlib.reload(l2grp)
importlib.reload(fitfunc)
# Set plot style
plt.style.use('/fs/home/stehr/qPAINT/styles/FoM.mplstyle')
#%%
######################################### Read in data
# Define folder of locs_picked.hdf5 file
dir_names=['/fs/pool/pool-schwille-paint/Data/D042/18-06-27_biexp/N12-3_align/']
# Define names of locs_picked.hdf5 file
file_names=['12vs3_ac_picked.hdf5']

path=os.path.join(dir_names[0],file_names[0])
# read in locs files
locs=fifo.read_locs(path)
# Extract NoFrames from meta-data
NoFrames=fifo.read_meta(path)[0]['Frames']

#%%
######################################### Select group
g=115
locs_g=locs[:][locs['group']==g]

######################################### Parameters
ignore_dark=1
px2dist=160
bright_ignore=False
fit_mode='exp'
########################################## Apply locs2groupprops module
# Get tau dynamics according to Jungmann
qTau=l2grp.get_tau(locs_g,ignore_dark,out='all',bright_ignore=bright_ignore,fit=fit_mode)
# Get trace
trace=l2grp.get_trace(locs_g,NoFrames)
# Get mutiple tau autocorrelation
ac=l2grp.get_ac(locs_g,NoFrames,8)
# Fit autocorrelation with mono-exponential
popt=l2grp.get_ac_fit(locs_g,NoFrames)
# Fit autocorrelation with bi-exponential
popt_bi=l2grp.get_ac_fit_bi(locs_g,NoFrames)

######################################### Plotting
f=plt.figure(num=10,figsize=[8,10])
f.clear()
f.suptitle('Group = %i'%(g))
plt.subplots_adjust(top=0.85, left=0.1, right=0.95,bottom=0.08, wspace=0.4, hspace=0.45)
################################################################################ Autocorrelation
ax=f.add_subplot(3,2,1)
ax.set_title(r'$\tau_c$ = '+'%.0f, '%(popt[1])+
             r'$A$ = '+'%.2f, '%(popt[0])+'\n'
             r'$\tau_{c,1}$ = '+'%.0f, '%(popt_bi[1])+
             r'$A_1$ = '+'%.2f, '%(popt_bi[0])+'\n'
             r'$\tau_{c,2}$ = '+'%.0f, '%(popt_bi[3])+
             r'$A_2$ = '+'%.2f, '%(popt_bi[2])
             )
ax.plot(ac[1:,0],ac[1:,1],'bx')
ax.plot(ac[1:,0],fitfunc.ac_monoexp(ac[1:,0],*popt[:-1]),color='red')
# Test double exponential
ax.plot(ac[1:,0],fitfunc.ac_biexp(ac[1:,0],*popt_bi[:-1]),color='black',ls='-')
ax.plot(ac[1:,0],fitfunc.ac_monoexp(ac[1:,0],*popt_bi[0:2]),color='black',ls='--')
ax.plot(ac[1:,0],fitfunc.ac_monoexp(ac[1:,0],*popt_bi[2:4]),color='black',ls='--')

ax.set_xscale('log')
ax.set_ylabel('G [a.u]')
############################################################################### tau_d _ecdf linearized
ax=f.add_subplot(3,2,3)
ax.set_title(r'$\tau_d$ = '+'%.0f, '%(qTau['tau_d_mean'][0])+
             r'$n_{events}$ = '+'%.0i'%(qTau['n_events'][0]))
ax.plot(qTau['tau_d_bins'][0],-np.log(1-qTau['tau_d_cdf'][0]),'x')
ax.plot(qTau['tau_d_bins'][0],qTau['tau_d_bins'][0]/qTau['tau_d_mean'][0])
ax.set_ylabel('-log[1-ECDF(dark)]')
################################################################################ tau_b ecdf linearized
ax=f.add_subplot(3,2,5)
ax.set_title(r'$\tau_b$ = '+'%.2f'%(qTau['tau_b_mean'][0]))
ax.plot(qTau['tau_b_bins'][0],-np.log(1-qTau['tau_b_cdf'][0]),'x')
ax.plot(qTau['tau_b_bins'][0],qTau['tau_b_bins'][0]/qTau['tau_b_mean'][0])
ax.set_xlabel('Frame')
ax.set_ylabel('-log[1-ECDF(bright)]')
################################################################################ x vs y
ax=f.add_subplot(3,2,2)
x=locs_g['x']
x=(x-np.mean(x))*px2dist # Define distances relative to mean x&y
y=locs_g['y']
y=(y-np.mean(y))*px2dist
ax.scatter(x,y,s=100,marker='o',edgecolor='none',alpha=0.1,color='white')

ax.set_facecolor('black')
ax.set_title('Super-Resolution')
ax.set_xlabel('x [nm]')
ax.set_ylabel('y [nm]')
ax.set_xlim([-px2dist,px2dist])
ax.set_ylim([-px2dist,px2dist])
################################################################################ Photon statistics
ax=f.add_subplot(3,2,4)
ax.hist(locs_g['photons'][:],bins='fd')
ax.set_ylabel('Counts')
ax.set_xlabel('Photons')
################################################################################ Trace
ax=f.add_subplot(3,2,6)
ax.plot(np.arange(0,NoFrames),trace[0,:],'b-')
ax.set_xlim([0,float(NoFrames)])
ax.set_ylabel('Photons')
ax.set_xlabel('Frame')

plt.show()


f=plt.figure(num=1)
f.clear()
ax=f.add_subplot(1,1,1)
CycleTime=0.05
ax.plot(ac[1:,0]*CycleTime,ac[1:,1],'bx',label='data')
ax.plot(ac[1:,0]*CycleTime,fitfunc.ac_monoexp(ac[1:,0],*popt[:-1]),color='red',label=r'$G_{mono}(\tau)$')
# Test double exponential
ax.plot(ac[1:,0]*CycleTime,fitfunc.ac_biexp(ac[1:,0],*popt_bi[:-1]),color='black',ls='-',label=r'$G_{bi}(\tau)$')
ax.plot(ac[1:,0]*CycleTime,fitfunc.ac_monoexp(ac[1:,0],*popt_bi[0:2]),color='black',ls='--',label=r'$A_1 e^{\tau/\tau_1}+1$')
ax.plot(ac[1:,0]*CycleTime,fitfunc.ac_monoexp(ac[1:,0],*popt_bi[2:4]),color='black',ls='--',label=r'$A_2 e^{\tau/\tau_2}+1$')
ax.set_xscale('log')
ax.set_xlabel(r'$\tau$ [s]')
ax.set_ylabel(r'G ($\tau$) []')
ax.legend(loc=1)

#f=plt.figure(num=2)
#f.clear()
#ax=f.add_subplot(1,1,1)
#ax.plot(qTau['tau_d_bins'][0],-np.log(1-qTau['tau_d_cdf'][0]),'x',label='data')
#ax.plot(qTau['tau_d_bins'][0],qTau['tau_d_bins'][0]/qTau['tau_d_mean'][0],label='fit')
#ax.set_ylabel(r'-ln[1-ECDF($\tau_d$)]')
#ax.set_xlabel('Time [frames]')
#ax.legend(loc=2)
#
#f=plt.figure(num=3)
#f.clear()
#ax=f.add_subplot(1,1,1)
#ax.plot(qTau['tau_b_bins'][0],-np.log(1-qTau['tau_b_cdf'][0]),'x',label='data')
#ax.plot(qTau['tau_b_bins'][0],qTau['tau_b_bins'][0]/qTau['tau_b_mean'][0],label='fit')
#ax.set_ylabel(r'-ln[1-ECDF($\tau_b$)]')
#ax.set_xlabel('Time [frames]')
#ax.legend(loc=2)
#
#f=plt.figure(num=4)
#f.clear()
#plt.subplots_adjust(left=0.2)
#ax=f.add_subplot(1,1,1)
#ax.plot(np.arange(0,NoFrames),trace[0,:]+6000,'b-')
#ax.plot(np.arange(0,NoFrames),np.sign(np.float16(trace[0,:]))*5000,'r-')
#ax.set_ylabel('Photons')
#ax.set_xlabel('Time')
#ax.set_xlim([1680,1865])
#ax.set_ylim([-1000,19900])
