import numpy as np

def gen_car1(t,alpha,sig):
	
	""" gen_car1 generates red noise time series with uneven time grid.
		Inputs:
		- t [1-dim numpy array of floats]: the times of the time series.
		- alpha [1-dim numpy array of floats]: values of the coefficient of the AR-1 process, as defined in
			'A General Theory on Spectral Analysis for Irregularly Sampled Time Series. I. Frequency Analysis', G. Lenoir and M. Crucifix
		- sig [1-dim numpy array of floats - same size as alpha]: values of the standard deviation of the white noise process.
		Outputs:
		- rednoise_vec [numpy array of floats - dim=(t.size,alpha.size)]: red noise time series. There alpha.size time series, corresponding to the values of the couple (alpha,sig).
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""

	n=t.size
	m=alpha.size
	rednoise_vec=np.zeros((n,m))
	multfact=np.sqrt(1.0/2.0/alpha)*sig
	rednoise_vec[0,:]=multfact*np.random.standard_normal(m)
	for k in range(1,n):
		rho=np.exp(-alpha*(t[k]-t[k-1]))
		rednoise_vec[k,:]=rho*rednoise_vec[k-1,:]+multfact*np.sqrt(1.0-rho**2)*np.random.standard_normal(m)

	return rednoise_vec