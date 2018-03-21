# Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import matplotlib as mpl
import os #platform independent paths
from mpl_toolkits.mplot3d import Axes3D

# Load user defined functions
import picasso_io
importlib.reload(picasso_io)
#%%
#file_name=['%icm_MMStack.ome.tif'%(i) for i in range(0,66,5)]
file_name=[r'Profile_MMStack.ome.tif']
dir_name=r'E:\Projects\SetupD042\data\18-02-28_piShaper_BeamCleaning_40-125Tele_150-400Mag'
path=[os.path.join(dir_name,name) for name in file_name]

# Load images of paths in stack
stack=[]
stack_meta=[]
for p in path:
    imag,info=picasso_io.load_tif(p)
    imag=imag[0,:,:]
    stack.extend([imag])
    stack_meta.extend([info])
        
labels=['%i cm'%(i) for i in range(0,66,5)]


mpl.rcParams['font.sans-serif']='Times New Roman'

f=plt.figure(num=11,figsize=[9,8])
f.clear()
f.suptitle('Lineplot along x through center')
ax=f.add_subplot(111)
ax.plot(stack[0][550,:])

f=plt.figure(num=10,figsize=[9,8])
f.clear()
for i in range(0,len(stack)):
    ax=f.add_subplot(1,1,i+1)
    stack[i][stack[i]>20000]=20000
    im=ax.imshow(stack[i],cmap='hot')
    
#    ax.set_title('%s'%(labels[i]))
    ax.set_xticks([])
    ax.set_yticks([])
#    ax.set_xlim([300,980])
#    ax.set_ylim([200,824])

f.subplots_adjust(top=0.9,bottom=0.01, left=0.05, right=0.85, hspace=0.01,wspace=0.0)
f.suptitle(r'$\pi$Shaper: Output magnified with -40mm/125mm 3x telescope'+'\n0.375x telecentric imaging (d = 15.4mm)')
cax = plt.axes([0.88, 0.035, 0.05, 0.85])
cbar=f.colorbar(im,cax=cax)

plt.show()

