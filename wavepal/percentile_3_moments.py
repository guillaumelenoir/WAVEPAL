import numpy as np
from scipy.optimize import least_squares
from scipy.special import gamma
from scipy.stats import gengamma
from percentile_3_moments_first_guess import percentile_3_moments_first_guess
from tqdm import trange
import sys

def percentile_3_moments(trA,trAsq,trAcub,proba,MaxFunEvals):
	
	""" percentile_3_moments returns the approximate percentiles of a normal form x'*A*x (where x is a multivariate standard normal distribution and A is a real symmetric matrix) by conserving its first 3 moments. The normal form is approximated with a generalized gamma distribution (which has 3 parameters).
		Inputs:
		- trA [1-dim numpy array of floats]: Trace of Matrix A. Several values are allowed (put them in a vector) in order to compute the percentiles of several normal forms.
		- trAsq [1-dim numpy array of floats - size=trA.size]: Trace of A^2. 
		- trAcub [1-dim numpy array of floats - size=trA.size]: Trace of A^3.
		N.B.: Remember that tr(A)=sum of the eigenvalues - Tr(A^2)=sum of squared eigenvalues - etc. It is usually quicker to compute trA, trAsq and trAcub from the eigenvalues of A.
		- proba [1-dim numpy array of floats]: the percentage at which the percentile is computed (for ex. 0.95) -> 0<=proba<=1. Several values are allowed (put them in a vector) in order to compute the percentiles at different percentages.
		- MaxFunEvals: "max_nfev" option for "least_squares" - see python help of "scipy.optimize.least_squares"
		Outputs:
		- alpha [numpy array of floats - size=trA.size] 
		- beta [numpy array of floats - size=trA.size]
		- delta [numpy array of floats - size=trA.size]
		alpha, beta and delta are the parameters of the generalized gamma distribution. More details in:
		'A General Theory on Spectral Analysis for Irregularly Sampled Time Series. I. Frequency Analysis', G. Lenoir and M. Crucifix
		- percentile [numpy array of floats - dim=(trA.size,proba.size)]: the percentiles.
		-----------------------------
		WARNING FOR EXTRA USES:
		If TrA, trAsq and trAcub are vectors, the parameters alpha, beta and delta of the generalized gamma distribution are determined for the first entry of those vectors. They are then used as a first guess for the next entry of trA, trAsq and trAcub to determine the new values of alpha, beta and delta. Etc. We thus implicitely guess that alpha, beta and delta are changing slowly when going through the values of trA, trAsq and trAcub, which is quite realistic in our case (the confidence levels slowly vary along a frequency (periodogram) or along the translation time (scalogram)).
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""
	
	m=trA.size
	nproba=proba.size
	percentile=np.zeros((m,nproba))
	l=0
	while (trA[l]==0. and trAsq[l]==0. and trAcub[l]==0.):
		l+=1
	g=trAsq[l]/trA[l]
	R=trA[l]**2/trAsq[l]
	alpha0,beta0,delta0=percentile_3_moments_first_guess(g,R)
	c0=np.zeros(3); c0[0]=alpha0; c0[1]=beta0; c0[2]=delta0  # first guess for the first frequency
	alpha=np.zeros(m)
	beta=np.zeros(m)
	delta=np.zeros(m)
	myfun1=lambda x: x[1]*gamma(x[0]+1.0/x[2])/gamma(x[0])
	myfun2=lambda x: x[1]**2*gamma(x[0]+2.0/x[2])/gamma(x[0])
	myfun3=lambda x: x[1]**3*gamma(x[0]+3.0/x[2])/gamma(x[0])
	print "Root-searching for the coefficients of the generalized gamma distribution:"
	for k in trange(m):
		if (trA[k]==0. and trAsq[k]==0. and trAcub[k]==0.):
			continue
		else:
			moment1=trA[k]
			moment2=2.0*trAsq[k]+trA[k]**2
			moment3=8.0*trAcub[k]+6.0*trA[k]*trAsq[k]+trA[k]**3
			# NOTE: I must use a method which does not give negative parameters as a solution
			# => "least_squares" is suitable for that, because It allows the user to provide bounds on the parameter values
			# !!! absolute() for use with least_squares (don't forget we look for the zeros)
			F=lambda x: [np.absolute(myfun1(x)-moment1),np.absolute(myfun2(x)-moment2),np.absolute(myfun3(x)-moment3)]
			answ=least_squares(F,c0,bounds=(0.0*c0,1000.0*c0),ftol=1.e-15,xtol=1.e-09,max_nfev=MaxFunEvals)
			try:
				assert answ.status>0
			except AssertionError:
				print "Error in percentile_3_moments.py with function least_squares"
				sys.exit(1)
			alphak=answ.x[0]
			betak=answ.x[1]
			deltak=answ.x[2]
			c0=answ.x	# first guess for the next frequency (whose solution should be close to the current one)
			percentile[k,:]=gengamma.ppf(proba,alphak,deltak,scale=betak)
			alpha[k]=alphak
			beta[k]=betak
			delta[k]=deltak

	return alpha,beta,delta,percentile
