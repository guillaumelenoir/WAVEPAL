import numpy as np
import numpy.linalg as la
import sys

def detrending(t,data,pol_order):
	
	""" detrending returns the least-squares polynomial trend of a time series.
		Inputs:
		- t [1-dim numpy array of floats]: the times of the time series.
		- data [1-dim numpy array of floats - size=t.size]: the data values of the time series.
		- pol_order [int]: the order of the polynomial trend. pol_order=-1 means no trend (trend is imposed to zero). pol_order=0 means trend=constant=data average. pol_order=1 means linear detrending (trend=a*t+b). etc.
		Outputs:
		- trend [1-dim numpy array of floats - size=t.size]: the least-squares polynomial trend.
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""
	
	if pol_order==-1:
		trend=np.zeros(t.size)
	else:
		tt=t[:]/t[-1]  # for numerical stability (because we take the powers of t !!!)
		z=np.ones((tt.size,pol_order+1))
		for k in range(1,pol_order+1):
			z[:,k]=tt[:]**k
		amat=np.dot(np.transpose(z),z)
		bmat=np.dot(np.transpose(z),data)
		beta=la.solve(amat,bmat)
		# Check that the solution is correct
		try:
			assert np.allclose(np.dot(amat,beta),bmat)==True
		except AssertionError:
			print "Error when computing detrending: Matrix is singular"
			sys.exit(1)
		trend=np.dot(z,beta)

	return trend