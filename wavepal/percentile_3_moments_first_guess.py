import numpy as np
from scipy.special import polygamma
from scipy.optimize import fminbound
import sys

def percentile_3_moments_first_guess(g,M):
	
	""" percentile_3_moments_first_guess returns a first guess for 'percentile_3_moments' function, on the basis of a 2 moments approximation (that is a weighted chi-square approx.). For more details, see:
		'A General Theory on Spectral Analysis for Irregularly Sampled Time Series. I. Frequency Analysis', G. Lenoir and M. Crucifix
		as well as:
		E. W. Stacy and G. A. Mihram. Parameter estimation for a generalized gamma distribution. Technometrics, 7(3):349-358, 1965.
		Inputs:
		- g [float]: multiplicative weight for the chi-square distribution
		- M [float]: dof of the chi-square distribution
		Outputs:
		- alpha [float], beta [float] and delta [float]: the first guess parameters of the generalized gamma distribution.
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""
	
	d1=np.log(2.0*g)+polygamma(0,M/2.0)
	d2=polygamma(1,M/2.0)
	d3=polygamma(2,M/2.0)
	myfun=lambda x: np.absolute(polygamma(2,x)/polygamma(1,x)**1.5+np.absolute(d3/d2**1.5)) # abs -> to use "fminbound" (we use it because the polygamma fct must take positive entries) - (don't forget we look for the zeros)
	alpha,_,ierr,_=fminbound(myfun,0.0,10000.0,disp=1,full_output=True)
	try:
		assert ierr==0
	except AssertionError:
		print "Error in percentile_3_moments_first_guess.py: convergence failed with fminbound"
		sys.exit(1)
	delta=np.sqrt(polygamma(1,alpha))/d2
	beta=np.exp(d1+polygamma(0,alpha)/delta)

	return alpha,beta,delta
