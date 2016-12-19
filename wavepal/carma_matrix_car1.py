import numpy as np

def carma_matrix_car1(t,alpha,sig):
	
	""" carma_matrix_car1 computes the matrix KK' for an zero-mean CAR-1 process.
		Triangular matrix K is defined in
		'A General Theory on Spectral Analysis for Irregularly Sampled Time Series. I. Frequency Analysis', G. Lenoir and M. Crucifix
		The product KK' (where K' is the transpose of K) takes a very simple form in the case of an CAR-1 process.
		Inputs:
		- t [1-dim numpy array of floats]: the times of the time series
		- alpha [float]: the coefficient of the CAR-1 process (as defined in the above reference)
		- sig [float]: the standard deviation of the white noise process
		Outputs:
		- Matrix KK' [numpy array of floats, dimension=(n,n), where n is the size of t]
		-------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""
	
	n=t.size
	carma_mat=np.zeros((n,n))
	for k in range(n):
		carma_mat[k,k]=1.
		for l in range(k):
			carma_mat[k,l]=np.exp(-alpha*(t[k]-t[l]))
			carma_mat[l,k]=carma_mat[k,l]
	return sig**2/2./alpha*carma_mat