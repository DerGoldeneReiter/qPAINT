# Load packages
import h5py as h5py #hdf5 handling
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os #platform independent paths
import importlib
import sys
import scipy.optimize
import multipletau
import matplotlib.gridspec as gridspec
import matplotlib

# Load user defined functions
import locs_groupprops as l2grp
import file_formats
import fitfunc
import simulate_locs

#Set path
dir_name='E:/Projects/qPAINT/data/18-01-02_OwnSimulation_1ds'
file_name='01'
path=os.path.join(dir_name,file_name)

# Set parameters for simulation
NoDocks=4
NoFrames=100000
NoGroups=1

CycleTime=50e-3 # Exposure time [s]
tau_b=5 # Bright time [s]
c_I=20e-9    # Imager concentration [M]
k_on=2e6 # On-rate [(Ms)^-1]
tau_d=1/(k_on*c_I) # Dark time
tau_c=1/((1/tau_b)+k_on*c_I) # Expected autocorrelation tau_c  [s]
# Convert to frames
tau_b=tau_b/CycleTime
tau_d=tau_d/CycleTime
tau_c=tau_c/CycleTime

# Create simulation locs
importlib.reload(simulate_locs)
locs=simulate_locs.simulate_locs(path,tau_b,tau_d,NoFrames,NoDocks,NoGroups)

#%%
# Parameters for dynamics extraction
ignore=0 # ignore dark frames identical to Picasso
lastframe=6e3 # ac is just fitted up to lastframe to reduce artifacts

# Select group
g=0
locs_g=locs[locs['group']==g]

# Reload modules
importlib.reload(l2grp)
importlib.reload(fitfunc)

# Get tau dynamics according to Jungmann
qTau=l2grp.get_tau(locs_g,ignore,out='all',bright_ignore=False)
# Get trace
trace=l2grp.get_trace(locs_g,NoFrames)
# Get mutiple tau autocorrelation
ac=l2grp.get_ac(locs_g,NoFrames)
# Fit autocorrelation with monoexponential
popt=l2grp.get_ac_fit(locs_g,NoFrames,lastframe)
norm_factor=popt[0]+popt[2] # Normalization factor for normalization of amplitude to 1


# Plotting
plt.figure(1)
plt.clf()
#plt.suptitle(r'$\tau_b$'+' = 0.6 s (12 frames), 10 nM')
# Define gridspec for plotting
gs_up=gridspec.GridSpec(2,2)
gs_up.update(top=0.92, bottom=0.55, wspace=0.2, hspace=0.2)
ax_tau_b_dist=plt.subplot(gs_up[0, 0])
ax_tau_b_cfd=plt.subplot(gs_up[1, 0])
ax_tau_d_dist=plt.subplot(gs_up[0, 1])
ax_tau_d_cfd=plt.subplot(gs_up[1,1])
gs_low=gridspec.GridSpec(2,2)
gs_low.update(top=0.47, bottom=0.04,hspace=0.2)
ax_ac=plt.subplot(gs_low[0, :])
ax_trace=plt.subplot(gs_low[1, :])

# tau_b distribution
ax_tau_b_dist.hist(qTau['tau_b_dist'][0],bins='fd')
ax_tau_b_dist.set_title(r'$\tau_b$ = '+'%.2f'%(qTau['tau_b_mean'][0])+' frames ')
ax_tau_b_dist.set_ylabel('Counts')
#tau_b cfd
ax_tau_b_cfd.plot(qTau['tau_b_bins'][0],qTau['tau_b_cfd'][0],'x')
ax_tau_b_cfd.plot(qTau['tau_b_bins'][0],fitfunc.exp_cfd(qTau['tau_b_bins'][0],qTau['tau_b_mean'][0]))
ax_tau_b_cfd.set_ylabel('Probability')
ax_tau_b_cfd.set_ylim([0,1.1])
ax_tau_b_cfd.set_xscale('log')

# tau_d distribution
ax_tau_d_dist.hist(qTau['tau_d_dist'][0],bins='fd')
ax_tau_d_dist.set_title(r'$\tau_d$ = '+'%.2f'%(qTau['tau_d_mean'][0])+' frames ')
# tau_d _cfd
ax_tau_d_cfd.plot(qTau['tau_d_bins'][0],qTau['tau_d_cfd'][0],'x')
ax_tau_d_cfd.plot(qTau['tau_d_bins'][0],fitfunc.exp_cfd(qTau['tau_d_bins'][0],qTau['tau_d_mean'][0]))
ax_tau_d_cfd.set_ylim([0,1.1])
ax_tau_d_cfd.set_xscale('log')

# Autocorrelation
ax_ac.plot(ac[:,0],ac[:,1]/norm_factor,'bx')
ax_ac.plot(ac[1:,0],fitfunc.ac_monoexp(ac[1:,0],popt[0],popt[1],popt[2])/norm_factor,'r')
ax_ac.set_xscale('log')
ax_ac.set_xlim([0.9,lastframe])
#ax_ac.set_ylim([0,1.1])
ax_ac.set_ylabel('G [a.u]')
ax_ac.set_title(r'$\tau_c$ = '+'%.2f'%(popt[1]) +' frames ')

# Trace
ax_trace.plot(np.arange(0,NoFrames),trace[0,:],'b')
ax_trace.set_xlim([0,float(NoFrames)])
ax_trace.set_ylabel('Counts')
ax_trace.set_xlabel('Frames')
#ax_trace.set_xlim([56250,56300])



