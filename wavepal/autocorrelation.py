import numpy as np

def autocorrelation(x,lag_min,lag_max):
	
	""" autocorrelation computes the autocorrelation of a time series, according to the formulas given in https://en.wikipedia.org/wiki/Autocorrelation
		Inputs:
		- x [1-dim numpy array of floats]: the time series (the time is supposed to be on a regular grid with a timestep of 1)
		- lag_min [int]: the minimum lag
		- lag_max [int]: the maximum lag
		Outputs:
		- r [1-dim numpy array of floats]: the autocorrelation, of size = lag_max-lag_min+1
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""
	
	x=x-np.mean(x)
	n=x.size
	if lag_min==0 and lag_max==0:
		r=np.zeros(1)
		r[0]=1.0
	elif lag_min==0:
		r=np.zeros(lag_max+1)
		r[0]=1.0
		for k in range(1,lag_max+1):
			r[k]=np.dot(x[:-k],x[k:])/float(n-k)
		r[1:]/=np.var(x)
	else:
		r=np.zeros(lag_max-lag_min+1)
		for k in range(lag_min,lag_max+1):
			r[k-lag_min]=np.dot(x[:-k],x[k:])/float(n-k)
		r/=np.var(x)
	return r