import numpy as np
import numpy.linalg as la

def carma_matrix(t,p,q,alpha,beta,varwn):

	""" carma_matrix computes the matrix K for an zero-mean CARMA(p,q) process with p>q and p>=1.
		Matrix K is defined in
		'A General Theory on Spectral Analysis for Irregularly Sampled Time Series. I. Frequency Analysis', G. Lenoir and M. Crucifix
		Inputs:
		- t [1-dim numpy array of floats]: the times of the time series
		- p and q [int]: the orders of the CARMA process. It is necessary that p>q and p>=1.
		- alpha [1-dim numpy array of floats]: vector of size p+1. It contains the CAR-coefficients, and is defined as [1., alpha_(p-1), ..., alpha_0] (following the convention defined in the article cited above). Pay attention to the order of those coefficients in the vector.
		- beta [1-dim numpy array of floats]: vector of size q+1. It contains the CMA-coefficients, and is defined as [1., beta_1, beta_2, ..., beta_q] (following the convention defined in the article cited above).
		- varwn [float]: the variance of the white noise process
		Outputs:
		- Matrix K [numpy array of floats, dimension=(n,n*p), where n is the size of t]
		-------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""

	bvec=np.zeros((p,1))
	bvec[range(q+1),0]=beta
	Umat=np.zeros((p,p),dtype='complex128')
	r=np.roots(alpha)
	print "roots of the autoregressive polynomial [must have a negative real part]: ", r
	Vmat=np.zeros((p,p))
	Vmatfact=np.zeros(p,dtype='complex128')
	for k in range(p):
		Vmatfact[k]=2.0*np.real(r[k])*np.prod((r[:k]-r[k])*(np.conj(r[:k])+r[k]))*np.prod((r[k+1:]-r[k])*(np.conj(r[k+1:])+r[k]))
	for k in range(p):
		for l in range(p):
			Umat[l,k]=r[k]**l
			Vmat[k,l]=np.real(np.sum((r**k)*((-r)**l)/Vmatfact))    # np.real(...) to be sure it is real (because it can be numerically sligthly complex)
	Vmat=-varwn*Vmat
	# N.B: Vmat is a real symmetric positive semi-definite matrix - consequently, mat1 is real positive or 0, and Wmat is real
	mat1,Wmat=la.eigh(Vmat)
	for k in range(p):
		Wmat[:,k]=Wmat[:,k]*np.sqrt(mat1[k])
	Umatinv=la.inv(Umat)
	kappa=Umatinv[:,-1]
	Ubvecprim=np.dot(np.transpose(bvec),Umat)
	n=t.size
	carma_matrix=np.zeros((n,n*p),dtype='complex128')
	# case i=0 (first row of CARMA matrix)
	cvec=np.dot(np.transpose(bvec),Wmat)
	for j in range(p):
		carma_matrix[0,j]=cvec[0,j]
	# case i>0 (other rows)
	UmatinvWmat=np.dot(Umatinv,Wmat)
	Ymat=np.zeros((p,p),dtype='complex128')
	P_j_mat=np.zeros((p,p,n),dtype='complex128')
	sigmamat=np.zeros((p,p),dtype='complex128')
	for i in range(1,n):
		for k in range(p):
			for l in range(k+1):
				sigmamat[k,l]=kappa[k]*np.conj(kappa[l])/(r[k]+np.conj(r[l]))*(1.0-np.exp((r[k]+np.conj(r[l]))*(t[i]-t[i-1])))
				sigmamat[l,k]=np.conj(sigmamat[k,l])
			Ymat[k,k]=np.exp(r[k]*(t[i]-t[0]))
		sigmamat=-varwn*sigmamat
		# N.B: Umat*sigmamat*Umat' is a real symmetric positive semi-definite matrix - consequently, mat1 is real positive or 0, and mat2 is real
		mat1,mat2=la.eigh(np.dot(Umat,np.dot(sigmamat,np.transpose(np.conj(Umat)))))
		for k in range(p):
			mat2[:,k]=mat2[:,k]*np.sqrt(mat1[k])
		P_j_mat[:,:,i]=la.solve(Umat,mat2)
		# case j=0 (first p columns)
		Kmat=np.dot(Ymat,UmatinvWmat)
		cvec=np.dot(Ubvecprim,Kmat)
		for l in range(p):
			carma_matrix[i,l]=cvec[0,l]
		# case j>0 (other columns)
		for j in range(1,i+1):
			for k in range(p):
				Ymat[k,k]=np.exp(r[k]*(t[i]-t[j]))
			Kmat=np.dot(Ymat,P_j_mat[:,:,j])
			cvec=np.dot(Ubvecprim,Kmat)
			for l in range(p):
				carma_matrix[i,p*j+l]=cvec[0,l]
	# carma_matrix is theoretically real, but the numerical result may be (very sligthly) complex
	return np.real(carma_matrix)
