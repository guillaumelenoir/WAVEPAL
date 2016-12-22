import numpy as np
from dt_normal import dt_normal
from dt_central import dt_central
from tqdm import trange
from tapering_window import tapering_window
import copy

def freq_analysis_prelims(time,freq,D,betafact,mywindow,coverage,freq_min_bound,freq_max_bound,pol_degree,WOSA_segments,weighted_WOSA):
	
	""" freq_analysis_prelims returns some variables for the frequency analysis in 'Wavepal' class.
		Inputs:
		- time [1-dim numpy array of floats]: the times of the time series, distinct and in ascending order.
		- freq [1-dim numpy array of floats]: the frequencies, distinct and in ascending order.
		- D [float]: the temporal length of the WOSA segments.
		- betafact [float - value in [0.,1.[]: overlapping factor for the WOSA segments.
		- mywindow [int]: window choice for the windowing of the WOSA segments. See tapering_window.py for more details.
		- coverage [float between 0. and 100.]: minimal coverage (in percent) of the data points along the segment length. Below this value, a WOSA segment is not considered. 
		- freq_min_bound [str - value = 'yes' or 'no']: limit, or not, the lower bound of the frequency range, for each WOSA segment. More details in
		'A General Theory on Spectral Analysis for Irregularly Sampled Time Series. I. Frequency Analysis', G. Lenoir and M. Crucifix
		- freq_max_bound [str - value = 'yes' or 'no']: limit, or not, the upper bound of the frequency range, for each WOSA segment. More details in the above cited article.
		- pol_degree [int]: Degree of the polynomial trend. pol_degree=-1 means no trend.
		- WOSA_segments [int or str]: Choose the minimal number of WOSA segments to be present at each frequency to take it into account, thus defining the frequency range for the analysis.
			-> WOSA_segments='all': No restrictions on the number of segments per frequency.
			-> WOSA_segments='max': Consider only the frequencies for which the number of WOSA segments is maximal. This is the most restrictive case.
			-> WOSA_segments=None: Consider only the frequencies for which the number of WOSA segments is at least 10, or maximal if there are less than 10 segments.
			-> WOSA_segments=n (n is an integer): Consider only the frequencies for which the number of WOSA segments is at least n.
		- weighted_WOSA [str - value = 'yes' or 'no']: 'yes' if weighted periodogram, or 'no' if not.
		Outputs:
		- freq [1-dim numpy array of floats]: frequencies on which the analysis is to be performed. 
		- tau_considered [1-dim numpy array of floats]: values of the times at which start the WOSA segments. Note that those times are not necessarily equal to values of the vector 'time'.
		- myind_time [numpy array of ints - dim=(tau_considered.size,2)]: min. and max. temporal indices (of the vector 'time') for each WOSA segment.
		- myind_freq_full [list of size=tau_considered.size]: Each entry of the list contains an array with the frequency indices (of the output vector 'freq') which are taken into account on the WOSA segment.
		- myind_Q [1-dim numpy array of ints - size=tau_considered.size]: indices of the segments taken into account, on the basis of a regular set of WOSA segments.
		- D [float]: the temporal length of the WOSA segments. Re-estimated value.
		- Qvec [1-dim numpy array of ints - size=freq.size]: Number of WOSA segments for each frequency. 
		- Qtot [int]: Number of WOSA segments. Note that Qtot is not necessarily equal to Qvec.max(), because each WOSA segment does not necessarily holds all the frequencies.
		- weight_WOSA [1-dim numpy array of floats - size=tau_considered.size]: the weights for the weighted periodogram. 
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""

	N=time.size
	J=freq.size
	coverage/=100.
	dt_centr=dt_central(time)  	# Central time step
	Q=int(np.floor(np.absolute(time[-1]-time[0]-D)/(1.0-betafact)/D))+1
	if (float(N)/float(Q))<10.:
		print "WARNING: Less than 10 data points per WOSA segment (in average)"
	D=(time[-1]-time[0])/(1.0+(1.0-betafact)*(Q-1))       # Estimate D from Q -> New value for D
	print "Re-estimated D factor (WOSA): ", D
	tau=np.zeros(Q)
	for k in range(Q):
		tau[k]=time[0]+(1.0-betafact)*D*float(k)
	myind_time=np.zeros((Q,2),dtype=int)  # min. and max. temporal indices for each WOSA segment.
	myind_freq=np.zeros((Q,2),dtype=int)  # min. and max. frequency indices for each WOSA segment.
	Qvec=np.zeros(J,dtype=int)
	tau_considered=np.zeros(Q)
	weight_WOSA=np.zeros(Q)
	myind_Q=np.zeros(Q,dtype=int)
	kin=0
	ll=-1
	print "Preliminary steps for the WOSA periodogram:"
	for l in trange(Q):
		count=0
		myindl_0=-1
		myindl_1=-1
		for k in range(kin,N):
			if ((time[k]-tau[l])>=0.0 and (time[k]-tau[l])<=D):
				if count==0:
					myindl_0=k
					myindl_1=k
					count=1
				else:
					myindl_1=k
		# If myindl_0==-1, there is a gap in the data that is > D
		# If less than pol_degree+3 data points, Gram-Schmidt is KO
		if (myindl_0>-1 and (myindl_1-myindl_0)>(pol_degree+1)):
			kin=myindl_0
			mytime=time[myindl_0:myindl_1+1]-tau[l]
			mydt1=dt_normal(mytime,D,mywindow)
			gvec=tapering_window(mytime,D,mywindow)
			mydt2=np.sum(gvec*dt_centr[myindl_0:myindl_1+1])/np.sum(gvec)
			per_max=mytime[-1]-mytime[0]
			per_min=2.0*max(mydt1,mydt2)
			# Check that the coverage in time on the WOSA segment is big enough
			if(per_max/D<coverage):
				continue
			# Check the frequency bounds for the WOSA segment
			if(freq_min_bound=='yes'):
				kmin=-1
				for k in range(J):
					if (1.0/freq[k]<per_max):
						kmin=k
						break
			elif(freq_min_bound=='no'):
				kmin=0
			if(freq_max_bound=='yes'):
				kmax=-1
				for k in range(J-1,-1,-1):
					if (1.0/freq[k]>per_min):
						kmax=k
						break
			elif(freq_max_bound=='no'):
				kmax=J-1
			if (kmin==-1 or kmax==-1):
				continue    # next iteration for the loop "for l in range(Q)"
			else:  # Go on if the frequency bounds exist
				ll+=1
				tau_considered[ll]=tau[l]
				myind_Q[ll]=l
				myind_time[ll,0]=myindl_0
				myind_time[ll,1]=myindl_1
				myind_freq[ll,0]=kmin
				myind_freq[ll,1]=kmax
				for k in range(kmin,kmax+1):
					Qvec[k]+=1
				if weighted_WOSA=="yes":
					weight_WOSA[ll]=np.sqrt(2.*np.sum(gvec**2))/np.sum(gvec)
				else:
					weight_WOSA[ll]=1.
	Qtot=ll+1         # Total number of WOSA segments on which the periodogram is going to be computed
	if type(WOSA_segments) is str:
		if WOSA_segments.lower()=="all":
			myQ=1
		elif WOSA_segments.lower()=="max":
			myQ=Qvec.max()
	elif WOSA_segments is None:
		myQ=min(Qvec.max(),10)
	else:
		myQ=WOSA_segments
	myind_freq_min=np.amin(myind_freq[0:Qtot,0])
	myind_freq_max=np.amax(myind_freq[0:Qtot,1])
	# NB: on range(myind_freq_min,myind_freq_max+1), we are sure that each frequency is in at least 1 WOSA segment
	for k in range(myind_freq_min,myind_freq_max+1):
		if Qvec[k]>=myQ:
			k0=k
			break
	myind_freq_full=[None]*Qtot
	mylargeset=set(range(myind_freq_min,myind_freq_max+1))
	for l in range(Qtot):
		mysmallset=set(range(myind_freq[l,0],myind_freq[l,1]+1))
		myind_freq_full[l]=np.asarray(list(mysmallset.intersection(mylargeset)))-k0
	myrange=Qvec>=myQ
	freq=copy.copy(freq[myrange])
	Qvec=copy.copy(Qvec[myrange])
	# check if a segment does not contain any frequency and makes the corresponding index reshaping
	# It would be for sure possible to do it in a better, more "python", way
	myind_freq_full_new=[None]*Qtot
	myind_time_new=np.zeros((Qtot,2),dtype=int)
	tau_considered_new=np.zeros(Qtot)
	myind_Q_new=np.zeros(Qtot,dtype=int)
	weight_WOSA_new=np.zeros(Qtot)
	kk=-1
	for k in range(Qtot):
		if myind_freq_full[k].size>0:
			kk+=1
			myind_freq_full_new[kk]=myind_freq_full[k]
			myind_time_new[kk,:]=myind_time[k,:]
			tau_considered_new[kk]=tau_considered[k]
			myind_Q_new[kk]=myind_Q[k]
			weight_WOSA_new[kk]=weight_WOSA[k]
	Qtot=kk+1
	del myind_freq_full
	myind_freq_full=myind_freq_full_new[0:Qtot]
	del myind_time
	myind_time=myind_time_new[0:Qtot,:]
	del myind_time_new
	del tau_considered
	tau_considered=tau_considered_new[0:Qtot]
	del tau_considered_new
	del myind_Q
	myind_Q=myind_Q_new[0:Qtot]
	del myind_Q_new
	del weight_WOSA
	weight_WOSA=weight_WOSA_new[0:Qtot]
	del weight_WOSA_new
	# Is there a gap in the frequency range?
	if freq.size>1:
		freq_step=freq[1]-freq[0]
		eps=0.001
		if (freq[-1]-freq[0])<((freq.size-1)*freq_step-eps*freq_step):
			warning_not_continuous=" with discontinuities => Be careful when plotting the figures in fct of the frequency"
		else:
			warning_not_continuous=""
	else:
		warning_not_continuous=""
	print "Re-estimated frequency range: from ", freq[0]," to ", freq[-1], warning_not_continuous

	return freq,tau_considered,myind_time,myind_freq_full,myind_Q,D,Qvec,Qtot,weight_WOSA
