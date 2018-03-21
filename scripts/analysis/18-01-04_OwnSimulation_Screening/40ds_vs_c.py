# Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os #platform independent paths
import matplotlib as mpl
# Load user defined functions
import file_formats as fifo
# Style sheet
plt.style.use(r'E:\Flo\repos\qPAINT\styles\FoM.mplstyle')

# Define folder of locs.hdf5 file
dir_name=[r'E:\Projects\qPAINT\data\18-01-03_OwnSimulation_12ds\20frames']*8
# Define names of locs.hdf5 file
file_names=['01_locs_picked_groupprops.hdf5']
file_names.extend(['02_locs_picked_groupprops.hdf5'])
file_names.extend(['03_locs_picked_groupprops.hdf5'])
file_names.extend(['04_locs_picked_groupprops.hdf5'])
# Define labels
labels=['0.5nM']
labels.extend(['2.5nM'])
labels.extend(['5nM'])
labels.extend(['10nM'])
labels.extend(['25nM'])
labels.extend(['50nM'])
labels.extend(['75nM'])
labels.extend(['100nM'])
labels=labels
# Create list of paths
path=list()
for i in range(0,len(file_names)):
    path.append(os.path.join(dir_name[i],file_names[i]))
    

# Overall simulation parameters
concs=np.array([float(label.replace('nM',''))*1e-9 for label in labels]) # List of concentrations as floats
tau_bs=np.array([20.]*8) # Bright time [Frames]
tau_ds=[10000.,2000.,1000.,500.,200,100,(2/3)*100,50]
tau_ds=np.array(tau_ds)
 # Dark time [Frames]
NoFrames=[100000]*8
NoDocks=np.array([40]*8)
# Calculate expected tau_c
tau_cs=np.divide(1,np.power(tau_bs,-1)+np.power(tau_ds,-1)) # Autocorr tau [s]
# Calculate expected p_inf
p_infs=NoDocks*np.divide(tau_bs,tau_ds)

# Initialize mean&std of important parameters
G0=np.empty([len(path),2])
Ginf=np.empty([len(path),2])
tau_c=np.empty([len(path),2])

tau_b=np.empty([len(path),2])
tau_b_ignore=np.empty([len(path),2])
tau_d=np.empty([len(path),2])

tau_b_lin=np.empty([len(path),2])
tau_b_lin_ignore=np.empty([len(path),2])
tau_d_lin=np.empty([len(path),2])

p_inf=np.empty([len(path),2])
ac_tau_b=np.empty([len(path),2])
ac_tau_d=np.empty([len(path),2])
ac_NoDocks=np.empty([len(path),2])

# Initialize mean&std of important parameters
for i in range(0,len(path)):
    tau_c[i,0]=np.mean(fifo.read_groupprops(path[i])['ac_tau'])
    tau_c[i,1]=np.std(fifo.read_groupprops(path[i])['ac_tau'])/10
    
    tau_b_lin[i,0]=np.mean(fifo.read_groupprops(path[i])['tau_b_lin'])
    tau_b_lin[i,1]=np.std(fifo.read_groupprops(path[i])['tau_b_lin'])/10
    tau_d_lin[i,0]=np.mean(fifo.read_groupprops(path[i])['tau_d_lin'])
    tau_d_lin[i,1]=np.std(fifo.read_groupprops(path[i])['tau_d_lin'])/10
    
    p_inf[i,0]=np.mean(fifo.read_groupprops(path[i])['ac_p_inf'])
    p_inf[i,1]=np.std(fifo.read_groupprops(path[i])['ac_p_inf'])/10
    ac_tau_b[i,0]=np.mean(fifo.read_groupprops(path[i])['ac_tau_b'])
    ac_tau_b[i,1]=np.std(fifo.read_groupprops(path[i])['ac_tau_b'])/10
    ac_tau_d[i,0]=np.mean(fifo.read_groupprops(path[i])['ac_tau_d'][fifo.read_groupprops(path[i])['ac_tau_d']<NoFrames[i]]) # Disregard all dark times longer than total number of frames
    ac_tau_d[i,1]=np.std(fifo.read_groupprops(path[i])['ac_tau_d'][fifo.read_groupprops(path[i])['ac_tau_d']<NoFrames[i]])/10

#%% 
## p_inf
#f=plt.figure(num=1)
#f.clear()
#f.subplots_adjust(left=0.2)
#ax=f.add_subplot(1,1,1)
#pltrange=np.arange(1,8)
## ac 
#ax.errorbar(concs[pltrange]*1e9,p_inf[pltrange,0]-p_infs[pltrange],yerr=p_inf[pltrange,1],fmt='o',label=r'$p_{\infty,out}-p_{\infty,in}$',color='blue')
## Draw 10% line
#ax.plot(concs[pltrange]*1e9,p_infs[pltrange]*0.1,'--',color='black',label='10% of '+r'$p_{\infty,in}$')
#ax.plot(concs[pltrange]*1e9,-p_infs[pltrange]*0.1,'--',color='black')
## Draw zero line
#ax.axhline(0,color='black')
#
## Legend
#ax.legend(loc=2)
## Scales
#ax.set_xscale('log')
#ax.set_yscale('linear')
## Labels
#ax.set_xlabel('Concentration [nM]')
#ax.set_ylabel('Residual')
##%%
## tau_c
#f=plt.figure(num=2)
#f.clear()
#f.subplots_adjust(left=0.2)
#ax=f.add_subplot(1,1,1)
#
#pltrange=np.arange(1,8)
## ac 
#ax.errorbar(concs[pltrange]*1e9,tau_c[pltrange,0]-tau_cs[pltrange],yerr=tau_c[pltrange,1],fmt='o',label=r'$\tau_{c,out}-\tau_{c,in}$',color='blue')
## Draw 10% line
#ax.plot(concs[pltrange]*1e9,tau_cs[pltrange]*0.1,'--',color='black',label='< 10% of '+r'$\tau_{c,in}$')
#ax.plot(concs[pltrange]*1e9,-tau_cs[pltrange]*0.1,'--',color='black')
## Draw zero line
#ax.axhline(0,color='black')
#
## Legend
#ax.legend(loc=1)
## Scales
#ax.set_xscale('log')
#ax.set_yscale('linear')
## Labels
#ax.set_xlabel('Concentration [nM]')
#ax.set_ylabel('Residual')


#%%
# tau_d
f=plt.figure(num=3)
f.clear()

ax=f.add_subplot(1,1,1)

pltrange=np.arange(0,4)
# ac 
ax.errorbar(concs[pltrange]*1e9,np.divide(ac_tau_d[pltrange,0]-tau_ds[pltrange],tau_ds[pltrange])*100,yerr=np.divide(ac_tau_d[pltrange,1],tau_ds[pltrange])*100,fmt='o',label=r'$\tau_{d,out}^{ac}-\tau_{d,in}$',color='blue')
# qPAINT
ax.errorbar(concs[pltrange]*1e9,np.divide(tau_d_lin[pltrange,0]-tau_ds[pltrange]/NoDocks[0],tau_ds[pltrange]/NoDocks[0])*100,yerr=np.divide(tau_d_lin[pltrange,1],tau_ds[pltrange]/NoDocks[0])*100,fmt='o',label=r'$\tau_{d,out}-\tau_{d,in}$',color='red')
ax.axhline(0,color='black')

# Legend
#ax.legend(loc=1)
# Scales
ax.set_xscale('log')
# Labels
ax.set_xlabel('Concentration [nM]')
ax.set_ylabel(r'$\langle\tau_d\rangle$ relative error [%]')

#%%
# tau_b
f=plt.figure(num=4)
f.clear()
ax=f.add_subplot(1,1,1)
pltrange=np.arange(0,4)
# ac 
ax.errorbar(concs[pltrange]*1e9,np.divide(ac_tau_b[pltrange,0]-tau_bs[0],tau_bs[0])*100,yerr=np.divide(ac_tau_b[pltrange,1],tau_bs[0])*100,fmt='o',label=r'$\tau_{b,out}^{ac}-\tau_{b,in}$',color='blue')
# qPAINT
ax.errorbar(concs[pltrange]*1e9,np.divide(tau_b_lin[pltrange,0]-tau_bs[0],tau_bs[0])*100,yerr=np.divide(tau_b_lin[pltrange,1],tau_bs[0])*100,fmt='o',label=r'$\tau_{b,out}-\tau_{b,in}$',color='red')
# Draw zero line
ax.axhline(0,color='black')


# Legend
#ax.legend(loc=1)
# Scales
ax.set_xscale('log')

# Labels
ax.set_xlabel('Concentration [nM]')
ax.set_ylabel(r'$\langle\tau_b\rangle$ relative error [%]')
