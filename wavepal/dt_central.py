import numpy as np

def dt_central(t):
	
	""" dt_central returns the central time step of a time vector t. It is defined as
		dt[k]=(t[k+1]-t[k-1])/2. for all the k except k=0 and k=t.size-1
		dt[0]=t[1]-t[0]
		dt[t.size-1]=t[-1]-t[-2]
		Inputs:
		- t [1-dim numpy array of floats]: the times.
		Outputs:
		- dt [1-dim numpy array of floats - size=t.size]: the central time step.
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""

	N=t.size
	dt=np.zeros(N)
	dt[0]=t[1]-t[0]
	dt[-1]=t[-1]-t[-2]
	for k in range(1,N-1):
		dt[k]=(t[k+1]-t[k-1])/2.
	return dt
