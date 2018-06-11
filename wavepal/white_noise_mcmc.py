import numpy as np
from scipy.stats import chi2 as chi2distr
from scipy.stats import gamma as gammadistr

def white_noise_mcmc(alpha_gamma,beta_gamma,smoothing_length,nmcmc):

	""" white_noise_mcmc generates chi-square samples, weigthed by a variance whose inverse follows a gamma distribution, for use with MCMC computations related to white noise in 'Wavepal'. More details in:
		'A General Theory on Spectral Analysis for Irregularly Sampled Time Series. I. Frequency Analysis', G. Lenoir and M. Crucifix
		Inputs:
		- alpha_gamma [float]: 'alpha' parameter of the gamma distribution
		- beta_gamma [float]: scale parameter of the gamma distribution
		- smoothing_length [int]: 2*smoothing_length is the length of each sample
		- nmcmc [int]: number of samples
		Outputs:
		- mywn [numpy array of floats - dim=(2*smoothing_length,nmcmc)]: the weighted chi-square samples.
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""

	# generate white noise processes for MCMC
	mywn=chi2distr.rvs(1,size=(smoothing_length*2,nmcmc))
	varwn_distr=1./gammadistr.rvs(alpha_gamma,scale=beta_gamma,size=nmcmc)
	for k in range(smoothing_length*2):
		mywn[k,:]=mywn[k,:]*varwn_distr
	return mywn