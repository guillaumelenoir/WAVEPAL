import numpy as np
from tapering_window import tapering_window

# computes the average WEIGHTED time-step (the weight = the taper window)

def dt_normal(t,D,mywindow):
	
	""" dt_normal returns the average weighted time step on a given WOSA segment. The weight is the taper window and is computed at mid-times:
		tmid[k]=(t[k]+t[k+1])/2.0
		The time step is: dt[k]=t[k+1]-t[k]
		Inputs:
		- t [1-dim numpy array of floats]: the times.
		- D [float]: the temporal length of the WOSA segment.
		- mywindow [int]: window choice for the windowing of the WOSA segment. See tapering_window.py for more details.
		Outputs:
		- dt [float]: the average weigthed time step.
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""

	n=t.size
	tmid=np.zeros(n-1)
	for k in range(n-1):
		tmid[k]=(t[k]+t[k+1])/2.0
	G=tapering_window(tmid,D,mywindow)
	dt_weighted=(t[1]-t[0])*G[0]
	for k in range(1,n-1):
		dt_weighted=dt_weighted+(t[k+1]-t[k])*G[k]
	dt=dt_weighted/np.sum(G)

	return dt
