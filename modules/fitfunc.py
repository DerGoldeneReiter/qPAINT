
import numpy as np

def exp_cdf(x,p):
    y=1-np.exp(-x/p)
    return y

def exp_pdf(x,p):
    y=1-(1/p)*np.exp(-x/p)
    return y

def ac_monoexp(x,A,tau):
    y=A*np.exp(-x/tau)+1.
    return y

def ac_biexp(x,As,taus,Al,taul):
    y=As*np.exp(-x/taus)+Al*np.exp(-x/taul)+1.
    return y

def ecdf(x):
    xs,ys=np.unique(x,return_counts=True) # Give unique values and counts in x
    ys=np.cumsum(ys) # Empirical cfd of dist without bias from binning
    ys = ys/len(x) # normalize so sum == 1
    return (xs,ys)

def exp_cdf_linearized(x,p):
    y=np.divide(x,p)
    return y

def poly_1(x,off,grad):
    y=off+grad*x
    return y

def gaussian_2D(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    (x, y) = xdata_tuple                                                        
    xo = float(xo)                                                              
    yo = float(yo)                                                              
    a = (np.cos(theta)**2)/(1*sigma_x**2) + (np.sin(theta)**2)/(1*sigma_y**2)   
    b = -(np.sin(2*theta))/(2*sigma_x**2) + (np.sin(2*theta))/(2*sigma_y**2)    
    c = (np.sin(theta)**2)/(1*sigma_x**2) + (np.cos(theta)**2)/(1*sigma_y**2)   
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)         
                        + c*((y-yo)**2)))                                   
    return g.ravel()

def w_mean(mean,std):
    # Set values of stds==0 to maximum std
    std[std==0]=np.max(std)
    # Define weigths
    w=np.divide(1,std)
    # Sum of weights
    w_sum=np.sum(w)
    # Weighted mean
    w_mean=np.sum(np.multiply(w,mean))/w_sum
    # Weighted std
    w_std=np.sqrt(np.sum(np.multiply(w,np.power(mean-w_mean,2)))/w_sum)
    return [w_mean,w_std]