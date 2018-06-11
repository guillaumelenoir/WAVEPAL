import numpy as np
from autocorrelation import autocorrelation

def decorrelation_length(x,min_autocorrelation):
	
	""" decorrelation_length returns the first occurrence lag at which the autocorrelation becomes smaller than min_autocorrelation.
		Inputs:
		- x [1-dim numpy array of floats]: the time series (the time is supposed to be on a regular grid with a timestep of 1)
		- min_autocorrelation [float]: the value of the autocorrelation for which we want the corresponding lag.
		Outputs:
		- mylength [int]: The first occurrence lag at which the autocorrelation becomes smaller than 'min_autocorrelation'.
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""
	
	n=x.size
	mylength=np.nan
	lag_min=0
	lag_max=9
	mybreak=False
	while lag_max<n:
		r=autocorrelation(x,lag_min,lag_max)
		for k in range(10):
			if r[k]<min_autocorrelation:
				mylength=lag_min+k
				mybreak=True
				break
		if mybreak==True:
			break
		lag_min+=10
		lag_max+=10
	return mylength