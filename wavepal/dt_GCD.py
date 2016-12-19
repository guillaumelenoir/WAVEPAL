import numpy as np

def dt_GCD(t):
	
	"""	Computes the greatest common divisor of the time steps of the times in t
		Inputs: 
		- t [1-dim numpy array of ints]: the times.
		WARNING: the times must all be distinct, of INTEGER type, and sorted in ascending order. In principle, the times come from a data file, such that they have a limited precision (in other words, they are rational numbers). Consequently, they can be transformed under integer form (note that doing that automatically can be quite challenging with python).
		Outputs:
		- dt_GCD [int]: the greatest common divisor of the time steps of the times in t
		Example:
		- See example_dt_GCD.py
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""
	
	# check inputs
	try:
		assert np.issubsctype(t,int)
	except:
		print "Error at input 't': must be a numpy array of 'int' type"
		return
	# Main code
	dt=np.zeros(t.size-1,dtype=long)
	for k in range(dt.size):
		dt[k]=t[k+1]-t[k]
	dt_GCD=dt[0]
	for k in range(1,dt.size):
		dt_GCD=np.core._internal._gcd(dt_GCD,dt[k])    # np.core._internal._gcd(n1,n2) computes the GCD of the integers n1 and n2 - NB:  GCD(n1,n2,n3)=GCD(GCD(n1,n2),n3)
	
	return dt_GCD