import numpy as np
import numpy.linalg as la
from tapering_window import tapering_window

def LS_WOSA(time,myfreq,freq_ind,myprojvec,D,tau,Q,Q_true,myind_time,myind_freq,mywindow,pol_order,weight_WOSA):

	""" LS_WOSA computes the basic bricks for all the computations related to the periodogram in 'Wavepal' class. It is called once for each frequency.
		Inputs:
		- time [1-dim numpy array of floats]: the times of the time series
		- myfreq [float]: the frequency at which the periodogram is to be computed
		- freq_ind [int]: the index of the frequency at which the periodogram is to be computed.
		- myprojvec [numpy array of floats - dimension=(time.size,pol_order+3)]: array with content related to the trend. This is an output from 'trend_vectors', in 'Wavepal' class.
		- D [float]: the temporal length of the WOSA segments, as output from freq_analysis_prelims.
		- tau [1-dim numpy array of floats]: values of the times at which start the WOSA segments, as output from freq_analysis_prelims (variable 'tau_considered')
		- Q [int]: Number of WOSA segments, as output from freq_analysis_prelims (variable 'Qtot')
		- Q_true [1-dim numpy array of ints]: Number of WOSA segments for each frequency, as output from freq_analysis_prelims (variable 'Qvec')
		- myind_time [numpy array of ints - dim=(tau.size,2)]: min. and max. temporal indices (of the vector 'time') for each WOSA segment, as output from freq_analysis_prelims
		- myind_freq [list of size=tau.size]: Each entry of the list contains an array with the frequency indices (of the output vector 'freq') which are taken into account on the WOSA segment, as output from freq_analysis_prelims (variable 'myind_freq_full')
		- mywindow [int]: window choice for the windowing of the WOSA segments. See tapering_window.py for more details.
		- pol_order [int]: order of the polynomial trend. pol_order=-1 means no trend.
		- weight_WOSA [1-dim numpy array of floats - size=tau.size]: the weights for the weighted periodogram.
		Outputs:
		- M2 [numpy array of floats - dimension=(time.size,2*Q_true)]: array containing the vectors on which we perform the orthogonal projection, in order to compute the periodogram. See:
		'A General Theory on Spectral Analysis for Irregularly Sampled Time Series. I. Frequency Analysis', G. Lenoir and M. Crucifix
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""
	
	N=time.size
	M2=np.zeros((N,2*Q_true))
	ll=-1
	for l in range(Q):
		if(freq_ind not in myind_freq[l]):
			continue   # next iteration in the loop "for l in range(Q)"
		else:
			myindl_0=myind_time[l,0]
			myindl_1=myind_time[l,1]
			mytime=time[myindl_0:myindl_1+1]-tau[l]
			gvec=np.zeros(N)
			gvec[myindl_0:myindl_1+1]=tapering_window(mytime,D,mywindow)
			domega=2.0*np.pi*myfreq
			mycos=np.cos(domega*time)
			mycos[:]=gvec*mycos				# weighted cosine
			for p in range(pol_order+1):     # Gram-Schmidt
				h=myprojvec[:,p]
				mycos[:]-=np.dot(h,mycos)*h
			mycos[:]/=la.norm(mycos)
			myprojvec[:,pol_order+1]=mycos[:]
			mysin=np.sin(domega*time)
			mysin[:]=gvec*mysin
			for p in range(pol_order+2):     # Gram-Schmidt
				h=myprojvec[:,p]
				mysin[:]-=np.dot(h,mysin)*h
			mysin[:]/=la.norm(mysin)
			myprojvec[:,pol_order+2]=mysin[:]
			ll+=1
			M2[:,2*ll]=mycos[:]*weight_WOSA[l]
			M2[:,2*ll+1]=mysin[:]*weight_WOSA[l]
			
	return M2
