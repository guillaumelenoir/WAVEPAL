from scipy.stats import chi2 as chi2distr
import numpy as np

def percentile_2_moments(trA,trAsq,proba):
	
	""" percentile_2_moments returns the approximate percentiles of a normal form x'*A*x (where x is a multivariate standard normal distribution and A is a real symmetric matrix) by conserving its first 2 moments. The normal form is approximated with a weighted chi square with M dof (or equivalently with a 2 params gamma distribution).
		Inputs:
		- trA [1-dim numpy array of floats]: Trace of Matrix A. Several values are allowed (put them in a vector) in order to compute the percentiles of several normal forms.
		- trAsq [1-dim numpy array of floats - size=trA.size]: Trace of A^2. This is the squared frobenius norm of A (compute it that way, it's much quicker)
		- proba [1-dim numpy array of floats]: the percentage at which the percentile is computed (for ex. 0.95) -> 0<=proba<=1. Several values are allowed (put them in a vector) in order to compute the percentiles at different percentages.
		Outputs:
		- percentile [numpy array of floats - dim=(trA.size,proba.size)]: the percentiles.
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""
	
	nproba=proba.size
	ntrA=trA.size
	percentile=np.zeros((ntrA,nproba))
	for k in range(nproba):
		for l in range(ntrA):
			if (trA[l]==0. and trAsq[l]==0.):
				continue
			else:
				percentile[l,k]=trAsq[l]/trA[l]*chi2distr.ppf(proba[k],(trA[l]**2)/trAsq[l])

	return percentile
