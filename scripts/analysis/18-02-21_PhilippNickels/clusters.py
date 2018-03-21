# Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import matplotlib as mpl
import os #platform independent paths
import importlib
from sklearn.cluster import KMeans 
from matplotlib.ticker import FormatStrFormatter
from cycler import cycler

# Load user defined functions
import locs_groupprops as l2grp
import file_formats as fifo
import fitfunc
importlib.reload(fifo)
importlib.reload(l2grp)
importlib.reload(fitfunc)

####################################################  Load data
dir_name=[r'E:\Projects\qPAINT\data\17-10-12_40xGrid_8nt10nt_PhilippNickel\all']
file_name=['170921_all_75pM-P3-Cy3b_2-5mW_30ms_1_MMStack_Pos0.ome_locs_picked_groupprops.hdf5']

path=[]
for i in range(0,len(file_name)):
    path.append(os.path.join(dir_name[i],file_name[i]))
    
groupprops_list=[]    
for i in range(0,len(path)):
    groupprops_list.append(fifo.read_groupprops(path[i]))
    

################################################### Filter:   
# 2) tau_c value        
groupprops_list_filter=[]
groupprops_list_filter_out=[]

for p in range(0,len(groupprops_list)):
    istrue_frame_mean=np.abs(groupprops_list[p]['mean_frame']-np.mean(groupprops_list[p]['mean_frame']))<15000
    istrue_frame_std=groupprops_list[p]['std_frame']>1000
    istrue_n_events=groupprops_list[p]['n_events']>5
    istrue_ac_chi=groupprops_list[p]['ac_chisquare']<0.051
    
    # Overall filter
    istrue=istrue_frame_mean&istrue_frame_std&istrue_n_events&istrue_ac_chi
    istrue_not=~np.array(istrue)
    
    groupprops_list_filter.append(groupprops_list[p][:][istrue])
    groupprops_list_filter_out.append(groupprops_list[p][:][istrue_not])

################################################### Show filtering  
# Show frame distribution
    
#f=plt.figure(1)
#f.clear()
#f.suptitle('Unfiltered frames')
#for p in range(0,len(groupprops_list)):
#    ax=f.add_subplot(2,len(groupprops_list),1)
#    ax.hist(groupprops_list[p]['mean_frame'][:],bins='fd')
#    ax=f.add_subplot(2,len(groupprops_list),2)
#    ax.hist(groupprops_list[p]['std_frame'][:],bins='fd')
#plt.show()
#
## Show filtered frame distribution
#f,ax=plt.subplots(2,len(groupprops_list),num=2)
#f.clear()
#f.suptitle('Filtered frames')
#for p in range(0,len(groupprops_list)):
#    ax=f.add_subplot(2,len(groupprops_list),1)
#    ax.hist(groupprops_list_filter[p]['mean_frame'][:],bins='fd')
#    ax=f.add_subplot(2,len(groupprops_list),2)
#    ax.hist(groupprops_list_filter[p]['std_frame'][:],bins='fd')
#plt.show()

################################################### Cluster ac
# Define 2 dimensional parameter space
X_ac=np.zeros([np.size(groupprops_list_filter),2])
X_ac[:,0]=np.log10(groupprops_list_filter[p]['ac_tau'][:])
X_ac[:,1]=np.log10(groupprops_list_filter[p]['n_events'][:])
kmeans_ac = KMeans(n_clusters=4).fit(X_ac)

################################################### Cluster pi
# Define 2 dimensional parameter space
X_ac=np.zeros([np.size(groupprops_list_filter),2])
X_ac[:,0]=np.log10(groupprops_list_filter[p]['tau_b_lin'][:])
X_ac[:,1]=np.log10(groupprops_list_filter[p]['n_events'][:])
kmeans_pi = KMeans(n_clusters=4).fit(X_ac)

#%%
################################################### Plot clustered pi
colors=['blue','red','orange','deeppink']
order=[0,1,2,3]

f=plt.figure(3,figsize=[7,7])
f.clear()
f.suptitle(r'$\tau_{b,lin}$ vs. $N_{events}$'+'\nCluster via k-means')
ax=f.add_subplot(1,1,1)
for c in  range(0,4):
    ax.scatter(groupprops_list_filter[0]['tau_b_lin'][kmeans_pi.labels_==c]*0.03,groupprops_list_filter[0]['n_events'][kmeans_pi.labels_==c],c=colors[order[c]])
ax.scatter(np.power(10,kmeans_pi.cluster_centers_[:,0])*0.03,np.power(10,kmeans_pi.cluster_centers_[:,1]),marker='x',color='black',s=100)
for i in range(0,len(groupprops_list_filter[p])):
    ax.text(groupprops_list_filter[0]['tau_b_lin'][i]*0.03, groupprops_list_filter[0]['n_events'][i],str(groupprops_list_filter[0]['group'][i]),color='black',fontsize=12)
ax.set_xscale('log')
ax.set_xlim([0.1,30])
ax.set_yscale('log')
ax.set_ylim([5,90])
ax.set_xlabel(r'$\tau_{b,lin}$'+' [s]')
ax.set_ylabel(r'$N_{events}$')
ax.yaxis.set_minor_formatter(FormatStrFormatter("%.f"))
ax.yaxis.set_major_formatter(FormatStrFormatter("%.f"))
plt.show()


order=[0,1,2,3]
f=plt.figure(4,figsize=[7,7])
f.clear()
f.suptitle(r'$\tau_c$ vs. $N_{events}$'+'\nCluster via k-means')
ax=f.add_subplot(1,1,1)
for c in  range(0,4):
    ax.scatter(groupprops_list_filter[0]['ac_tau'][kmeans_ac.labels_==c]*0.03,groupprops_list_filter[0]['n_events'][kmeans_ac.labels_==c],c=colors[order[c]])
ax.scatter(np.power(10,kmeans_ac.cluster_centers_[:,0])*0.03,np.power(10,kmeans_ac.cluster_centers_[:,1]),marker='x',color='black',s=100)
for i in range(0,len(groupprops_list_filter[p])):
    ax.text(groupprops_list_filter[0]['ac_tau'][i]*0.03, groupprops_list_filter[0]['n_events'][i],str(groupprops_list_filter[0]['group'][i]),color='black',fontsize=12)


ax.set_xscale('log')
ax.set_xlim([0.1,30])
ax.set_yscale('log')
ax.set_ylim([5,90])
ax.set_xlabel(r'$\tau_{c}$'+' [s]')
ax.set_ylabel(r'$N_{events}$')
ax.yaxis.set_minor_formatter(FormatStrFormatter("%.f"))
ax.yaxis.set_major_formatter(FormatStrFormatter("%.f"))
plt.show()


f=plt.figure(5,figsize=[7,7])
f.clear()
f.suptitle(r'$\tau_c$ vs. $p_{\infty}$'+'\nCluster via k-means')
ax=f.add_subplot(1,1,1)

for c in  range(0,4):
    ax.scatter(groupprops_list_filter[0]['ac_tau'][kmeans_ac.labels_==c]*0.03,groupprops_list_filter[0]['ac_p_inf'][kmeans_ac.labels_==c],c=colors[order[c]])
for i in range(0,len(groupprops_list_filter[p])):
    ax.text(groupprops_list_filter[0]['ac_tau'][i]*0.03, groupprops_list_filter[0]['ac_p_inf'][i],str(groupprops_list_filter[0]['group'][i]),color='black',fontsize=12)
    
ax.set_xscale('log')
ax.set_xlim([0.1,30])
ax.set_yscale('log')
ax.set_ylim([1e-26,3])
ax.set_xlabel(r'$\tau_{c}$'+' [s]')
ax.set_ylabel(r'$p_{\infty}$')

#order=[0,1,2,3]
#f=plt.figure(5,figsize=[8,8])
#f.clear()
#f.suptitle(r'$\tau_c$ vs. $N_{events}$'+'\nCluster via k-means')
#ax=f.add_subplot(1,1,1)
#for c in  range(0,4):
#    ax.scatter(groupprops_list_filter[0]['ac_tau'][kmeans_ac.labels_==c]*0.03,groupprops_list_filter[0]['n_events'][kmeans_ac.labels_==c],c=colors[order[c]])
#ax.scatter(np.power(10,kmeans_ac.cluster_centers_[:,0])*0.03,np.power(10,kmeans_ac.cluster_centers_[:,1]),marker='x',color='black',s=100)
#ax.scatter(groupprops_list_filter_out[0]['ac_tau']*0.03,groupprops_list_filter_out[0]['n_events'],marker='+',c='black')
##for i in range(0,len(groupprops_list_filter[p])):
##    ax.text(groupprops_list_filter[0]['ac_tau'][i]*0.03, groupprops_list_filter[0]['n_events'][i],str(groupprops_list_filter[0]['group'][i]),color='black',fontsize=12)
#for i in range(0,len(groupprops_list_filter_out[p])):
#    ax.text(groupprops_list_filter_out[0]['ac_tau'][i]*0.03, groupprops_list_filter_out[0]['n_events'][i],str(groupprops_list_filter_out[0]['group'][i]),color='black',fontsize=13)    
#
#ax.set_xscale('log')
#ax.set_xlim([0.04,50])
#ax.set_yscale('log')
#ax.set_ylim([1,90])
#ax.set_xlabel(r'$\tau_{c}$'+' [s]')
#ax.set_ylabel(r'$N_{events}$')
#ax.yaxis.set_minor_formatter(FormatStrFormatter("%.f"))
#ax.yaxis.set_major_formatter(FormatStrFormatter("%.f"))
#plt.show()