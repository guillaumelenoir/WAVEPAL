import numpy as np
import numpy.linalg as la
from scipy.special import gamma
from scipy.special import gammainc
from scipy.optimize import brentq
from percentile_2_moments import percentile_2_moments
from percentile_3_moments import percentile_3_moments
#from tqdm import trange
import sys

def percentile_n_moments(traces,proba,ind_full,algo,MaxFunEvals):
	
	""" percentile_n_moments returns the approximate percentiles of a normal form x'*A*x (where x is a multivariate standard normal distribution and A is a real symmetric matrix) by conserving its first n moments. The normal form is approximated with a (generalized) gamma-polynomial distribution. More details in:
		'A General Theory on Spectral Analysis for Irregularly Sampled Time Series. I. Frequency Analysis', G. Lenoir and M. Crucifix
		Inputs:
		- traces [numpy array of floats]: matrix of dimension (m,n), the element (i,j) being tr(A_i^(j+1)) (where tr() is the trace of a matrix and indices start at 0). Several values are allowed for matrix A (thus denoted A_i) in order to compute the percentiles of several normal forms. 'traces' must have at least 2 columns, i.e. n>=2 (the minimum approx. is a 2-moments one). Remember that tr(A)=sum of the eigenvalues - Tr(A^2)=sum of squared eigenvalues - etc. It is usually quicker to compute tr(A), tr(A^2), tr(A^3), ... from the eigenvalues of A. However, there is an exception in the case n=2. In such case, it is quicker to compute the traces without going through the eigenvalues, with tr(A) and tr(A^2)=the squared frobenius norm of A.
		- proba [1-dim numpy array of floats]: the percentage at which the percentile is computed (for ex. 0.95) -> 0<=proba<=1. Several values are allowed (put them in a vector) in order to compute the percentiles at different percentages.
		- ind_full [list of ints]: vector of indices (0<index<m-1) at which the d-moments-conserved-percentiles will be returned, for 2<=d<=n. This is returned in output variable 'percentile_full'.
		- algo='gamma-polynomial' or algo='generalized-gamma-polynomial'. The first choice is more stable and is thus recommended in most cases.
		- MaxFunEvals: "max_nfev" option for "least_squares" - see python help of "scipy.optimize.least_squares". Used if algo='generalized-gamma-polynomial'.
		Outputs:
		- percentile [numpy array of floats - dim=(traces.shape[0],proba.size)]: the percentiles.
		- percentile_full [numpy array of floats - dim=(len(ind_full),traces.shape[1]-1,proba.size)]: percentiles for the normal forms corresponding to the indices given in 'ind_full' input variable. Percentiles with d moments conserved, 2<=d<=n, are returned. That allows to check the convergence of the d-moments approximation.
		-----------------------------
		WARNING FOR EXTRA USES:
		If 'traces' contains more than one row and if n>2, the percentiles are determined for the first row of traces. They are then used as a first guess for compuation of the percentiles corresponding to the next row of traces. Etc. We thus implicitely guess that the percentiles vary slowly when going through the rows of traces, which is quite realistic in our case (the confidence levels slowly vary along a frequency (periodogram) or along the translation time (scalogram)).
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""
	
	m=traces.shape[0]
	n_moments=traces.shape[1]
	nproba=proba.size
	p=len(ind_full)
	percentile_full=np.zeros((p,n_moments-1,nproba))
	myindex=range(n_moments+1)
	myindex_float=np.arange(float(n_moments+1))
	#print "Computes the analytical approximate confidence levels for each frequency:"
	if n_moments==2:
		percentile=percentile_2_moments(traces[:,0],traces[:,1],proba)
		for k in range(p):
			percentile_full[k,0,:]=percentile[ind_full[k],:]   # index 0 for the second moment
	elif n_moments==3 and algo=='generalized-gamma-polynomial':
		_,_,_,percentile=percentile_3_moments(traces[:,0],traces[:,1],traces[:,2],proba,MaxFunEvals)
		for k in range(p):
			percentile_full[k,0,:]=percentile_2_moments(traces[ind_full[k],0],traces[ind_full[k],1],proba)
			percentile_full[k,1,:]=percentile[ind_full[k],:]
	else:
		cum=np.zeros(traces.shape)
		cum[:,0]=traces[:,0]
		for s in range(1,n_moments):
			cum[:,s]=2**s*gamma(s+2)*traces[:,s]/float(s+1)
		moment=np.zeros((m,n_moments+1))
		moment[:,0]=1.0							# 0th order moment
		for h in range(1,n_moments+1):
			mymoment=np.zeros(m)
			for i in range(h):
				mymoment=mymoment+cum[:,h-i-1]*moment[:,i]/gamma(h-i)/gamma(i+1)
			moment[:,h]=mymoment*gamma(h)
		percentile2=percentile_2_moments(traces[:,0],traces[:,1],proba)
		if algo=='generalized-gamma-polynomial':
			alphad,betad,deltad,percentile3=percentile_3_moments(traces[:,0],traces[:,1],traces[:,2],proba,MaxFunEvals)
			l=0
			while (np.sum(traces[l,:]==0.0)==n_moments):
				l+=1
			c0=percentile3[l,:] # first guess
		elif algo=='gamma-polynomial':
			alphad=moment[:,1]**2/(moment[:,2]-moment[:,1]**2)
			betad=moment[:,2]/moment[:,1]-moment[:,1]
			l=0
			while (np.sum(traces[l,:]==0.0)==n_moments):
				l+=1
			c0=percentile2[l,:] # first guess
		nu=np.zeros((m,2*n_moments+1))
		if algo=='generalized-gamma-polynomial':
			for h in range(2*n_moments+1):
				nu[:,h]=betad**h*gamma(alphad+h/deltad)/gamma(alphad)
		elif algo=='gamma-polynomial':
			for h in range(2*n_moments+1):
				nu[:,h]=betad**h*gamma(alphad+h)/gamma(alphad)
		numat=np.zeros((n_moments+1,n_moments+1))
		percentile=np.zeros((m,nproba))
		compt=0
		if algo=='generalized-gamma-polynomial':
			#for h in trange(m):
			for h in range(m):
				if np.sum(traces[h,:]==0.0)==n_moments:
					continue
				for k in range(n_moments+1):
					for i in range(n_moments+1):
						numat[i,k]=nu[h,i+k]
				myvec=np.transpose(moment[h,:])
				ksi=la.solve(numat,myvec)
				# Check that the solution is correct
				try:
					assert np.allclose(np.dot(numat,ksi),myvec)==True
				except AssertionError:
					print "Error when computing the n-moments approx for the linear combination of chi-squares: Matrix is singular"
					sys.exit(1)
				myfact=ksi[:]*betad[h]**myindex*gamma(myindex_float[:]/deltad[h]+alphad[h])
				cumdist=lambda c0var: np.sum(myfact[:]*gammainc(myindex_float[:]/deltad[h]+alphad[h],c0var**deltad[h]/betad[h]**deltad[h]))
				for k in range(nproba):
					cumdisth=lambda c0var: cumdist(c0var)/gamma(alphad[h])-proba[k]
					c0_try_low=c0_try_high=c0[k]
					while True:
						c0_try_low=c0_try_low/2.0
						c0_try_high=c0_try_high*2.0
						cumdisth_low=cumdisth(c0_try_low)
						cumdisth_high=cumdisth(c0_try_high)
						if (cumdisth_low<0 and cumdisth_high>0):
							break
					c0[k]=brentq(cumdisth,c0_try_low,c0_try_high)
					percentile[h,k]=c0[k]
				if compt<p and h==ind_full[compt]:
					percentile_full[compt,0,:]=percentile2[h,:]
					percentile_full[compt,1,:]=percentile3[h,:]
					for k in range(4,n_moments):
						numatk=numat[:(k+1),:(k+1)]
						myvec=np.transpose(moment[h,:(k+1)])
						ksi=la.solve(numatk,myvec)
						try:
							assert np.allclose(np.dot(numatk,ksi),myvec)==True
						except AssertionError:
							print "Error when computing the n-moments approx for the linear combination of chi-squares: Matrix is singular"
							sys.exit(1)
						myindexk=range(k+1)
						myindexk_float=np.arange(float(k+1))
						myfact=ksi[:(k+1)]*betad[h]**myindexk*gamma(myindexk_float/deltad[h]+alphad[h])
						cumdist=lambda c0var: np.sum(myfact[:]*gammainc(myindexk_float[:]/deltad[h]+alphad[h],c0var**deltad[h]/betad[h]**deltad[h]))
						for l in range(nproba):
							cumdistk=lambda c0var: cumdist(c0var)/gamma(alphad[h])-proba[l]
							c0_try_low=c0_try_high=c0[l]
							while True:
								c0_try_low=c0_try_low/2.0
								c0_try_high=c0_try_high*2.0
								cumdistk_low=cumdistk(c0_try_low)
								cumdistk_high=cumdistk(c0_try_high)
								if (cumdistk_low<0 and cumdistk_high>0):
									break
							c0[l]=brentq(cumdistk,c0_try_low,c0_try_high)
							percentile_full[compt,k-2,l]=c0[l]
					percentile_full[compt,n_moments-2,:]=percentile[h,:]
					compt=compt+1
		elif algo=='gamma-polynomial':
			#for h in trange(m):
			for h in range(m):
				if np.sum(traces[h,:]==0.0)==n_moments:
					continue
				for k in range(n_moments+1):
					for i in range(n_moments+1):
						numat[i,k]=nu[h,i+k]
				myvec=np.transpose(moment[h,:])
				ksi=la.solve(numat,myvec)
				# Check that the solution is correct
				try:
					assert np.allclose(np.dot(numat,ksi),myvec)==True
				except AssertionError:
					print "Error when computing the n-moments approx for the linear combination of chi-squares: Matrix is singular"
					sys.exit(1)
				myfact=ksi[:]*betad[h]**myindex*gamma(myindex_float[:]+alphad[h])
				cumdist=lambda c0var: np.sum(myfact[:]*gammainc(myindex_float[:]+alphad[h],c0var/betad[h]))
				for k in range(nproba):
					cumdisth=lambda c0var: cumdist(c0var)/gamma(alphad[h])-proba[k]
					c0_try_low=c0_try_high=c0[k]
					while True:
						c0_try_low=c0_try_low/2.0
						c0_try_high=c0_try_high*2.0
						cumdisth_low=cumdisth(c0_try_low)
						cumdisth_high=cumdisth(c0_try_high)
						if (cumdisth_low<0 and cumdisth_high>0):
							break
					c0[k]=brentq(cumdisth,c0_try_low,c0_try_high)
					percentile[h,k]=c0[k]
				if compt<p and h==ind_full[compt]:
					percentile_full[compt,0,:]=percentile2[h,:]
					for k in range(3,n_moments):
						numatk=numat[:(k+1),:(k+1)]
						myvec=np.transpose(moment[h,:(k+1)])
						ksi=la.solve(numatk,myvec)
						try:
							assert np.allclose(np.dot(numatk,ksi),myvec)==True
						except AssertionError:
							print "Error when computing the n-moments approx for the linear combination of chi-squares: Matrix is singular"
							sys.exit(1)
						myindexk=range(k+1)
						myindexk_float=np.arange(float(k+1))
						myfact=ksi[:(k+1)]*betad[h]**myindexk*gamma(myindexk_float+alphad[h])
						cumdist=lambda c0var: np.sum(myfact[:]*gammainc(myindexk_float[:]+alphad[h],c0var/betad[h]))
						for l in range(nproba):
							cumdistk=lambda c0var: cumdist(c0var)/gamma(alphad[h])-proba[l]
							c0_try_low=c0_try_high=c0[l]
							while True:
								c0_try_low=c0_try_low/2.0
								c0_try_high=c0_try_high*2.0
								cumdistk_low=cumdistk(c0_try_low)
								cumdistk_high=cumdistk(c0_try_high)
								if (cumdistk_low<0 and cumdistk_high>0):
									break
							c0[l]=brentq(cumdistk,c0_try_low,c0_try_high)
							percentile_full[compt,k-2,l]=c0[l]
					percentile_full[compt,n_moments-2,:]=percentile[h,:]
					compt=compt+1
	return percentile,percentile_full
