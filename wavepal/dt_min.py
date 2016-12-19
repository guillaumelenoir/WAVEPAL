import numpy as np

def dt_min(time):
	
	""" dt_min computes the smallest time step of an irregularly sampled time series. 
		The time step is defined as dt[k]=time[k]-time[k-1]
		Inputs:
		- time [1-dim numpy array of floats]: the times of the time series. 
		Outputs:
		- dt_min [float]: the smallest time step
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""
	
	dt_min=time[-1]-time[0]
	n=time.size
	for k in range(1,n):
		dt_min=np.minimum(dt_min,(time[k]-time[k-1]))
	return dt_min