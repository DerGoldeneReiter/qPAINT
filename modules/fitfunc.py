
import numpy as np

def exp_cdf(x,p):
    y=1-np.exp(-x/p)
    return y

def exp_pdf(x,p):
    y=1-(1/p)*np.exp(-x/p)
    return y

def ac_monoexp(x,A,tau,off):
    y=A*np.exp(-x/tau)+off
    return y

def ecdf(x):
    xs,ys=np.unique(x,return_counts=True) # Give unique values and counts in x
    ys=np.cumsum(ys) # Empirical cfd of dist without bias from binning
    ys = ys/len(x) # normalize so sum == 1
    return (xs,ys)

def exp_cdf_linearized(x,p):
    y=np.divide(x,p)
    return y
