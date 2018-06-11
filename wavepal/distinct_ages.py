import numpy as np

def distinct_ages(t,y):
	
	""" distinct_ages returns a copy of a time series after having removed the times/ages multiple occurrences. It selects the first occurrence of a time/age if there are many.
		Example: consider the time series: t=[1,2,2,4,5,5,6] and x=[1.,8.,9.,7.,7.,5.,0.]. distinct_ages will return t=[1,2,4,5,6] and x=[1.,8.,7.,7.,0.].
		Inputs: 
		- t [1-dim numpy array of floats]: the times of the time series.
		- y [1-dim numpy array of floats - size=t.size]: the data values of the time series corresponding to t.
		Outputs:
		- tt [1-dim numpy array of floats]: the times of the time series. They are all distinct.
		- yy [1-dim numpy array of floats - size=tt.size]: the data values of the time series corresponding to tt.
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""
	
	n=t.size
	tt=np.zeros(1)
	tt[0]=t[0]
	yy=np.zeros(1)
	yy[0]=y[0]
	for k in range(1,n):
		if t[k]==t[k-1]:
			print "same time/age at ", t[k]
		else:
			tt=np.append(tt,t[k])
			yy=np.append(yy,y[k])
	return tt, yy