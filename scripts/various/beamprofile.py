# Load packages
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os #platform independent paths
import scipy.optimize

# Load user defined functions
import fitfunc
import picasso_io

#%%
dir_name=r'E:\Projects\qPAINT\data\18-01-12\sample01_pseries100_convampf_3x_T22_1'
file_name=r'sample01_pseries100_convampf_3x_T22_1_MMStack_Pos0.ome.tif'
path=os.path.join(dir_name,file_name)
stack,info=picasso_io.load_tif(path)

#%%
stack_sum=np.sum(stack[0:2000,:,:],axis=0)/np.shape(stack[0:2000,:,:])[0] # Integrate over all frames, divide by frame number
#stack_sum=stack # Use if just one image

#%%
# Number of rows& columns in 
Nr=np.shape(stack_sum)[0]
Nc=np.shape(stack_sum)[1]
# Meshgrid for rows and columns
rc_grid=np.empty([Nr,Nc,2])
rc_grid[:,:,0],rc_grid[:,:,1]=np.meshgrid(np.arange(0,Nc,1),np.arange(0,Nr,1))

# Fit
# Inital guess
A=np.max(stack_sum[:])
x0=250#center[0]
y0=250#center[1]
sigma_x=100
sigma_y=100
theta=50
offset=np.min(stack_sum[:])
p0=[A,x0,y0,sigma_x,sigma_y,theta,offset]
lowbounds=[0,0,0,0,0,0,0]
upbounds=[np.inf,len(stack_sum),len(stack_sum),len(stack_sum),len(stack_sum),90,A]

# Fit
popt,pcov=scipy.optimize.curve_fit(fitfunc.gaussian_2D,(rc_grid[:,:,0].ravel(),rc_grid[:,:,1].ravel()),stack_sum.ravel(),p0=p0,bounds=(lowbounds,upbounds))
# Create fitted data
stack_fitted=fitfunc.gaussian_2D((rc_grid[:,:,0],rc_grid[:,:,1]),*popt).reshape(Nr,Nc)

#%%
plt.clf()
f,ax=plt.subplots(1,num=1)
plt.subplots_adjust(top=0.8)
f.suptitle('Beam Profile: Jung, 561nm, 2x2 binning, 1 px =  6.5 microns\n'+\
           r'$x_0$'+' = %.2f px, '%(popt[1])+\
           r'$y_0$'+' = %.2f px, '%(popt[2])+\
           r'$\theta$'+' = %.2f\n'%(popt[5])+\
           r'$\sigma_x$'+' = %.2f px '%(popt[3]*np.sqrt(2))+r'$[1/e^2]$ ,'+\
           r'$\sigma_y$'+' = %.2f px '%(popt[4]*np.sqrt(2))+r'$[1/e^2]$ ,'+'\n'\
           )

ax.imshow(stack_sum)
ax.contour(rc_grid[:,:,0],rc_grid[:,:,1],stack_fitted, 8, colors='w')

#%%
row=260

plt.clf()
f,ax=plt.subplots(1,num=1)
plt.subplots_adjust(top=0.83)

f.suptitle('Line plot at row %i'%(row))
ax.plot(stack_sum[row,:],'.')
ax.plot(stack_fitted[row,:])