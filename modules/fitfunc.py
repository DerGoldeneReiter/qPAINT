from scipy.stats import binom
import numpy as np
from scipy.interpolate import interp1d

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
    ys = ys/ys[-1] # normalize so sum == 1
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

def gaussian_1D(x, *p):
	a, b, c, d = p
	y = a*np.exp(-np.power((x - b), 2.)/(2. * c**2.)) + d

	return y

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

def get_e_ecdf(x):
    xs,ys=ecdf(x)
    value=1-1/np.exp(1)
    index=np.argmin(np.abs(ys-value))
    e_ecdf=xs[index]
    return e_ecdf

def tau_c_of_conc(conc,koff,kon):
    y=np.divide(1,koff+kon*conc)
    return y

def A_of_conc(conc,koff,kon,N):
    y=np.divide(koff,N*kon*conc)
    return y

def binomial_single(N,p):
    # Success number k
    k=np.arange(0,N+1,1)
    # Create empty array for success number k and probability
    k_prob=np.zeros([len(k),2])
    # Assign success number k to first row
    k_prob[:,0]=k
    # Assign probability of success to second row
    k_prob[:,1]=binom.pmf(k,N,p)
    return k_prob

def binomial_stoich(N1,N2,p,bins):
    # Get probability to detect k docking sites for each species
    k1_prob=binomial_single(N1,p)
    k2_prob=binomial_single(N2,p)
    # Init array for ratio an probabilities
    Nratio=np.zeros([len(k1_prob)*len(k2_prob),2])
    
    Nratio_index=0
    for i in range(0,len(k1_prob)):
        for j in range(0,len(k2_prob)):
            Nratio[Nratio_index,0]=k1_prob[i,0]/k2_prob[j,0]
            Nratio[Nratio_index,1]=k1_prob[i,1]*k2_prob[j,1]
            Nratio_index=Nratio_index+1
    # Remove infinity entries
    Nratio=Nratio[np.invert(np.isinf(Nratio[:,0]))]
    # Normalization after removing infs like in real experiment
    Nratio[:,1]=Nratio[:,1]/np.sum(Nratio[:,1])
    # Initialize Nratio_binnned
    Nratio_binned=np.zeros([len(bins)-1,2])
    # Loop thorugh bins and search for values
    for i in range(1,len(bins)):
        # Start of bin (included)
        Nratio_bin_start=bins[i-1]
        # End of bin (excluded)
        Nratio_bin_end=bins[i]
        # Boolean for half open interval, start included, end excluded
        istrue_interval=(Nratio[:,0]>=Nratio_bin_start)&((Nratio[:,0]<Nratio_bin_end))
        # Assign start value of interval 
        Nratio_binned[i-1,0]=Nratio_bin_start
        # Sum probabilities
        Nratio_binned[i-1,1]=np.sum(Nratio[istrue_interval,1])
        
    return Nratio_binned[:,1] 

def binomial_A_pdf(bins,N,taub,p,taud):
    # Get success of k under N ds with detetction probability p
    k_pdf=binomial_single(N,p)
    # Get inverse for calculation of A of autocorrelation
    k_pdf[:,0]=1/k_pdf[:,0]
    # Remove division by zero for k=0
    k_pdf=k_pdf[1:,:]
    # re-normalization after removing inf entries
    k_pdf[:,1]=k_pdf[:,1]/np.sum(k_pdf[:,1])
    # Flip k_prob to get ascending 1/k
    k_pdf=np.flip(k_pdf,axis=0)
    # Multiply k with taud/taub to get A
    k_pdf[:,0]=k_pdf[:,0]*(taud/taub)
    A_pdf=k_pdf.copy()
    # Interpolation function for discrete pdf
    fA_pdf=interp1d(A_pdf[:,0],A_pdf[:,1])
    # Apply interpolation to binning
    A_pdf_inter=np.zeros([len(bins),2])
    # Assign bins to A_pdf_inter
    A_pdf_inter[:,0]=bins
    # Search for bin values inside possible A_pdf values
    istrue_bins=(bins>=A_pdf[0,0])&(bins<A_pdf[-1,0])
    # Assign interpolated values within bins
    A_pdf_inter[istrue_bins,1]=fA_pdf(bins[istrue_bins])
    A_pdf_inter[:,1]=A_pdf_inter[:,1]/np.sum(A_pdf_inter[:,1])
    
    return A_pdf_inter[:,1]
