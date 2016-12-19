import numpy as np
from scipy.special import iv

def tapering_window(time,D,mywindow):
	
	""" tapering_window returns the window for tapering a WOSA segment.
		Inputs:
		- time [1-dim numpy array of floats]: times along the WOSA segment.
		- D [float]: Temporal length of the WOSA segment.
		- mywindow [int]: Choice of tapering window:
			-> 1: Square window
			-> 2: Triangular window
			-> 3: sin window
			-> 4: sin**2 (Hanning) window
			-> 5: sin**3 window
			-> 6: sin**4 window
			-> 7: Hamming window, defined as 0.54-0.46*np.cos(2.0*np.pi*time/D)
			-> 8: 4-term Blackman-Harris window, with a0=0.35875 and a1=0.48829 and a2=0.14128 and a3=0.01168
			-> 9: Kaiser-Bessel window, with parameter alpha=2.5
			-> 10: Gaussian window, with standard dev. sigma=D/6.0
		The terminology and formulas come from:
		F. Harris. On the use of windows for harmonic analysis with the discrete fourier transform. Proceedings of the IEEE, 66(1):51-83, January 1978.
		WARNING: Provide the vector 'time' such that for all k=0,...,time.size-1, we have time[k]>=0 and time[k]<=D
		Outputs:
		- tapering_window [1-dim numpy array of floats - size=time.size]: the tapering window.
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""
	
	T=time.size
	if mywindow==1:
		tapering_window=np.ones(T)
	elif mywindow==2:
		tapering_window=1.0-np.absolute(time-D/2.0)/(D/2.0)
	elif mywindow==3:
		tapering_window=np.sin(np.pi*time/D)
	elif mywindow==4:
		tapering_window=(np.sin(np.pi*time/D))**2
	elif mywindow==5:
		tapering_window=(np.sin(np.pi*time/D))**3
	elif mywindow==6:
		tapering_window=(np.sin(np.pi*time/D))**4
	elif mywindow==7:
		tapering_window=0.54-0.46*np.cos(2.0*np.pi*time/D)
	elif mywindow==8:
		a0=0.35875
		a1=0.48829
		a2=0.14128
		a3=0.01168
		tapering_window=a0-a1*np.cos(2.0*np.pi*time/D)+a2*np.cos(4.0*np.pi*time/D)-a3*np.cos(6.0*np.pi*time/D)
	elif mywindow==9:
		alpha=2.5
		tapering_window=iv(0,np.pi*alpha*np.sqrt(1.0-((time-D/2.0)/(D/2.0))**2))
	elif mywindow==10:
		sig=D/6.0
		tapering_window=np.exp(-(time-D/2.0)**2/2.0/sig**2)
	else:
		print "Error: The window number you entered is not valid. Check input variable 'mywindow'."
		return
	return tapering_window
