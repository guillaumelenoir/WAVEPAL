import numpy as np
from tqdm import trange
from dt_central import dt_central
import copy

def timefreq_analysis_prelims(time,tau,scale,w0,gauss_spread,eps,dt_GCD,shannonnyquistexclusionzone,weighted_CWT,smoothing_coeff,smoothing_type):

	""" timefreq_analysis_prelims returns some variables for the time-frequency analysis in 'Wavepal' class.
		Inputs:
		- time [1-dim numpy array of floats]: the times of the time series, distinct and in ascending order.
		- tau [1-dim numpy array of floats]: the times at which the CWT is to be computed, distinct and in ascending order.
		- scale [1-dim numpy array of floats]: the scales, distinct and in ascending order.
		- w0 [float]: the usual parameter for the Morlet wavelet controlling the time-frequency resolution. Minimal allowed value is w0=5.5.
		- gauss_spread [float]: parameter for the spread of gaussian. 2*gauss_spread*std (where std is the standard dev. of the gaussian) is the approximate SUPPORT (i.e. where the function is not zero) of a gaussian. Typical values are gauss_spread=3.0 (conservative choice) or sqrt(2.0) (value taken in Torrence and Compo, 1998, and some subsequent papers). This is used for the computation of the cone of influence and for the max. allowed scale.
		- eps [float]: parameter controlling the flexibility on the value of the border of the Shannon-Nyquist exclusion zone. Typical value is eps=1.0e-06. 
		- dt_GCD [float]: the greatest common divisor of the time steps of the times in 'time'. 
		- shannonnyquistexclusionzone [str - value=True or False]: activate or not the Shannon-Nyquist exclusion zone.
		- weighted_CWT [str - value=True or False]: True if weighted scalogram, or False if not.
		- smoothing_coeff [float]: smoothing the CWT is performed this way: at a given (tau,scale), the CWT_smoothed is the average of the CWT along neighbouring values of tau (at the same scale). The interval of those neighbouring values of tau is delimited by: tau-smoothing_coeff*std and tau+smoothing_coeff*std, where std is the standard deviation of the gaussian. Note that:
			-> the std depends on the scale.
			-> Near the edges of the time series or close the the Shannon-Nyquist exclusion zone, the full interval over tau-smoothing_coeff*std to tau+smoothing_coeff*std cannot be considered. In such case, the CWT_smoothed at (tau,scale) is either ignored (if smoothing_type='fixed') or the interval around tau is shortened (if smoothing_type='variable').
		- smoothing_type [str - value='fixed' or 'variable']: See above the explanations for 'smoothing_coeff'.
		Outputs:
		- scale [1-dim numpy array of floats]: the scales on which the analysis is to be performed. 
		- coi1 [1-dim numpy array of floats - size=scale.size]: cone of influence - left part
		- coi2 [1-dim numpy array of floats - size=scale.size]: cone of influence - right part
		- coi1_smooth [1-dim numpy array of floats - size=scale.size]: border of the forbidden zone on the left side. Only used if smoothing_type='fixed'.
		- coi2_smooth [1-dim numpy array of floats - size=scale.size]: border of the forbidden zone on the right side. Only used if smoothing_type='fixed'.
		- coi_smooth_ind [1-dim numpy array of ints - size=tau.size]: indices, of the vector 'tau', corresponding to coi1_smooth and coi2_smooth.
		- weight_cwt [numpy array of floats - dim=(tau.size,scale.size)]: the weights for the weighted scalogram.
		- scalelim1_ind [1-dim numpy array of ints - size=tau.size]: scale indices at the border of the Shannon-Nyquist exclusion zone. For a given tau, this is the first scale index for which the CWT is to be computed.
		- scalelim1_smooth [1-dim numpy array of floats - size=tau.size]: scales at the border of the Shannon-Nyquist exclusion zone when there is smoothing. N.B.: if smoothing_type='fixed' and if there is smoothing (smoothing_coeff>0), the Shannon zone is refined.
		- scalelim1_ind_smooth [1-dim numpy array of ints - size=tau.size]: indices of the scales in scalelim1_smooth, plus 1. For a given tau, this is the first scale index for which the CWT_smoothed is drawn.
		- Qmax [int]: maximum length over which smoothing is to be performed, expressed in number of indices of tau.
		- n_outside_scalelim1 [1-dim numpy array of ints - size=scale.size]: number of tau's outside the Shannon-Nyquist exclusion zone, for each scale. 
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""

	N=time.size
	J=scale.size
	Q=tau.size
	
	# Cone of influence
	coi1=time[0]+gauss_spread*w0*scale
	coi2=time[-1]-gauss_spread*w0*scale

	# Central time step
	dt_centr=dt_central(time)

	# Normal time step and timemid
	dt=np.zeros(N-1)
	timemid=np.zeros(N-1)
	for k in range(N-1):
		dt[k]=time[k+1]-time[k]
		timemid[k]=(time[k]+time[k+1])/2.

	# Weights for the CWT squared norm and Shannon-Nyquist exclusion zone
	print "Weights for the CWT squared norm and Shannon-Nyquist exclusion zone:"
	dcon=1./2./w0**2
	scalelim1=np.ones(Q)*dt_GCD/np.pi*(1.+eps)
	scalelim1_ind=np.zeros(Q,dtype=int)
	weight_cwt=-np.ones((Q,J))
	if shannonnyquistexclusionzone is True and weighted_CWT is True:
		for k in trange(Q):
			tau_k=tau[k]
			time_ind0=0
			time_ind1=N-1
			for l in range(J-1,-1,-1):
				scale_l=scale[l]
				time_ind0=np.argmin(np.absolute(time[time_ind0:]-tau_k+3.*w0*scale_l))+time_ind0
				time_ind1=np.argmin(np.absolute(time[:(time_ind1+1)]-tau_k-3.*w0*scale_l))
				ind0=max(time_ind0-1,0)			# take the index minus 1
				ind1=min(time_ind1+1,N-1)+1		# take the index plus 1
				mydt=dt[ind0:ind1-1]
				mytime=timemid[ind0:ind1-1]-tau_k
				hvec=np.exp(-dcon/scale_l**2*mytime**2)
				dtmp1=sum(hvec*mydt)/sum(hvec)
				mytime=time[ind0:ind1]-tau_k
				mydt=dt_centr[ind0:ind1]
				gvec=np.exp(-dcon/scale_l**2*mytime**2)
				dtmp2=sum(gvec*mydt)/sum(gvec)
				dtmp=max(dtmp1,dtmp2)
				if scale_l<=dtmp/np.pi*(1.+eps):
					scalelim1[k]=scale[l]
					scalelim1_ind[k]=l+1
					break
				weight_cwt[k,l]=np.sqrt(2.*sum(gvec**2)/sum(gvec)**2)
	elif shannonnyquistexclusionzone is False and weighted_CWT is True:
		for k in trange(Q):
			tau_k=tau[k]
			time_ind0=0
			time_ind1=N-1
			Ipass=0
			for l in range(J-1,-1,-1):
				scale_l=scale[l]
				time_ind0=np.argmin(np.absolute(time[time_ind0:]-tau_k+3.*w0*scale_l))+time_ind0
				time_ind1=np.argmin(np.absolute(time[:(time_ind1+1)]-tau_k-3.*w0*scale_l))
				ind0=max(time_ind0-1,0)			# take the index minus 1
				ind1=min(time_ind1+1,N-1)+1		# take the index plus 1
				mydt=dt[ind0:ind1-1]
				mytime=timemid[ind0:ind1-1]-tau_k
				hvec=np.exp(-dcon/scale_l**2*mytime**2)
				dtmp1=sum(hvec*mydt)/sum(hvec)
				mytime=time[ind0:ind1]-tau_k
				mydt=dt_centr[ind0:ind1]
				gvec=np.exp(-dcon/scale_l**2*mytime**2)
				weight_cwt[k,l]=np.sqrt(2.*sum(gvec**2)/sum(gvec)**2)
				dtmp2=sum(gvec*mydt)/sum(gvec)
				dtmp=max(dtmp1,dtmp2)
				if scale_l<=dtmp/np.pi*(1.+eps) and Ipass==0:
					scalelim1[k]=scale[l]
					Ipass=1
	elif shannonnyquistexclusionzone is True and weighted_CWT is False:
		for k in trange(Q):
			tau_k=tau[k]
			time_ind0=0
			time_ind1=N-1
			for l in range(J-1,-1,-1):
				scale_l=scale[l]
				time_ind0=np.argmin(np.absolute(time[time_ind0:]-tau_k+3.*w0*scale_l))+time_ind0
				time_ind1=np.argmin(np.absolute(time[:(time_ind1+1)]-tau_k-3.*w0*scale_l))
				ind0=max(time_ind0-1,0)			# take the index minus 1
				ind1=min(time_ind1+1,N-1)+1		# take the index plus 1
				mydt=dt[ind0:ind1-1]
				mytime=timemid[ind0:ind1-1]-tau_k
				hvec=np.exp(-dcon/scale_l**2*mytime**2)
				dtmp1=sum(hvec*mydt)/sum(hvec)
				mytime=time[ind0:ind1]-tau_k
				mydt=dt_centr[ind0:ind1]
				gvec=np.exp(-dcon/scale_l**2*mytime**2)
				dtmp2=sum(gvec*mydt)/sum(gvec)
				dtmp=max(dtmp1,dtmp2)
				if scale_l<=dtmp/np.pi*(1.+eps):
					scalelim1[k]=scale[l]
					scalelim1_ind[k]=l+1
					break
				weight_cwt[k,l]=1.
	elif shannonnyquistexclusionzone is False and weighted_CWT is False:
		for k in trange(Q):
			tau_k=tau[k]
			time_ind0=0
			time_ind1=N-1
			Ipass=0
			for l in range(J-1,-1,-1):
				scale_l=scale[l]
				time_ind0=np.argmin(np.absolute(time[time_ind0:]-tau_k+3.*w0*scale_l))+time_ind0
				time_ind1=np.argmin(np.absolute(time[:(time_ind1+1)]-tau_k-3.*w0*scale_l))
				ind0=max(time_ind0-1,0)			# take the index minus 1
				ind1=min(time_ind1+1,N-1)+1		# take the index plus 1
				mydt=dt[ind0:ind1-1]
				mytime=timemid[ind0:ind1-1]-tau_k
				hvec=np.exp(-dcon/scale_l**2*mytime**2)
				dtmp1=sum(hvec*mydt)/sum(hvec)
				mytime=time[ind0:ind1]-tau_k
				mydt=dt_centr[ind0:ind1]
				gvec=np.exp(-dcon/scale_l**2*mytime**2)
				dtmp2=sum(gvec*mydt)/sum(gvec)
				dtmp=max(dtmp1,dtmp2)
				weight_cwt[k,l]=1.
				if scale_l<=dtmp/np.pi*(1.+eps) and Ipass==0:
					scalelim1[k]=scale[l]
					Ipass=1

	# redefine the scale and other variables because new min scale
	min_scalelim1_ind=min(scalelim1_ind)
	scale=copy.copy(scale[min_scalelim1_ind:])
	coi1=copy.copy(coi1[min_scalelim1_ind:])
	coi2=copy.copy(coi2[min_scalelim1_ind:])
	weight_cwt=copy.copy(weight_cwt[:,min_scalelim1_ind:])
	scalelim1_ind[:]=scalelim1_ind[:]-min_scalelim1_ind

	# Smoothing the CWT => parameters
	#print "Parameters for smoothing the CWT:"
	if shannonnyquistexclusionzone is True:
		if smoothing_type=="fixed":
			scalelim1_ind_max=np.amax(scalelim1_ind)-1
			scalelim1_smooth=copy.copy(scalelim1)
			scalelim1_ind_smooth=copy.copy(scalelim1_ind)
			for k in range(Q):
				tau_k=tau[k]
				for l in range(scalelim1_ind_max,-1,-1):
					scale_l=scale[l]
					ind_left=np.argmin(np.absolute(tau-(tau_k-smoothing_coeff*w0*scale_l)))
					ind_right=np.argmin(np.absolute(tau-(tau_k+smoothing_coeff*w0*scale_l)))
					if np.sum(scalelim1[ind_left:(ind_right+1)]>=scale_l)>0:
						scalelim1_smooth[k]=scale[l]
						scalelim1_ind_smooth[k]=l+1
						break
			# redefine the scale and other variables because new min scale
			min_scalelim1_ind_smooth=min(scalelim1_ind_smooth)
			scale=copy.copy(scale[min_scalelim1_ind_smooth:])
			weight_cwt=copy.copy(weight_cwt[:,min_scalelim1_ind_smooth:])
			scalelim1_ind[:]=scalelim1_ind[:]-min_scalelim1_ind_smooth
			scalelim1_ind_smooth[:]=scalelim1_ind_smooth[:]-min_scalelim1_ind_smooth
			# redefine the scale and other variables because new max scale
			scalemax=(time[-1]-time[0])/2./w0/(gauss_spread+smoothing_coeff)
			scalemax_ind=np.argmin(np.absolute(scale-scalemax))
			scale=copy.copy(scale[:(scalemax_ind+1)])
			weight_cwt=copy.copy(weight_cwt[:,:(scalemax_ind+1)])
			scalelim1_ind[scalelim1_ind>scalemax_ind]=scalemax_ind
			scalelim1_ind_smooth[scalelim1_ind_smooth>scalemax_ind]=scalemax_ind
			# Redefine the cone of influence
			coi1_smooth=tau[0]+smoothing_coeff*w0*scale
			coi2_smooth=tau[-1]-smoothing_coeff*w0*scale
			coi1=time[0]+(smoothing_coeff+gauss_spread)*w0*scale
			coi2=time[-1]-(smoothing_coeff+gauss_spread)*w0*scale
			# Indices for cwt coi1_smooth and coi2_smooth
			J=scale.size
			coi_smooth_ind=np.zeros(Q,dtype=int)
			for k in range(Q):
				tau_k=tau[k]
				for l in range(J-1,scalelim1_ind_smooth[k]-1,-1):
					if tau_k>coi1_smooth[l] and tau_k<coi2_smooth[l]:
						coi_smooth_ind[k]=l
						break
		elif smoothing_type=="variable":
			scalelim1_smooth=copy.copy(scalelim1)
			scalelim1_ind_smooth=copy.copy(scalelim1_ind)
			coi1_smooth=time[0]*np.ones(coi1.size)
			coi2_smooth=time[-1]*np.ones(coi2.size)
			J=scale.size
			coi_smooth_ind=np.ones(Q,dtype=int)*(J-1)
	elif shannonnyquistexclusionzone is False:
		scalelim1_smooth=copy.copy(scalelim1)
		scalelim1_ind_smooth=copy.copy(scalelim1_ind)
		if smoothing_type=="fixed":
			# redefine the scale and other variables because new max scale
			scalemax=(time[-1]-time[0])/2./w0/(gauss_spread+smoothing_coeff)
			scalemax_ind=np.argmin(np.absolute(scale-scalemax))
			scale=copy.copy(scale[:(scalemax_ind+1)])
			weight_cwt=copy.copy(weight_cwt[:,:(scalemax_ind+1)])
			scalelim1_ind[scalelim1_ind>scalemax_ind]=scalemax_ind
			scalelim1_ind_smooth[scalelim1_ind_smooth>scalemax_ind]=scalemax_ind
			# Redefine the cone of influence
			coi1_smooth=tau[0]+smoothing_coeff*w0*scale
			coi2_smooth=tau[-1]-smoothing_coeff*w0*scale
			coi1=time[0]+(smoothing_coeff+gauss_spread)*w0*scale
			coi2=time[-1]-(smoothing_coeff+gauss_spread)*w0*scale
			# Indices for cwt coi1_smooth and coi2_smooth
			J=scale.size
			coi_smooth_ind=np.zeros(Q,dtype=int)
			for k in range(Q):
				tau_k=tau[k]
				for l in range(J-1,scalelim1_ind_smooth[k]-1,-1):
					if tau_k>coi1_smooth[l] and tau_k<coi2_smooth[l]:
						coi_smooth_ind[k]=l
						break
		elif smoothing_type=="variable":
			coi1_smooth=time[0]*np.ones(coi1.size)
			coi2_smooth=time[-1]*np.ones(coi2.size)
			J=scale.size
			coi_smooth_ind=np.ones(Q,dtype=int)*(J-1)

	# Computes the maximum length over which smoothing is to be performed - That naturally occurs at the highest scale
	Qmax=0
	for k in range(Q):
		tau_k=tau[k]
		# Indices for cwt smoothing
		ind_left=np.argmin(np.absolute(tau-(tau_k-smoothing_coeff*w0*scale[-1])))
		ind_right=np.argmin(np.absolute(tau-(tau_k+smoothing_coeff*w0*scale[-1])))
		if smoothing_type=="fixed":
			# average the cwt over range(ind_left,ind_right+1)
			Qmax=np.maximum(Qmax,ind_right+1-ind_left)
		if smoothing_type=="variable":
			Qmax=np.maximum(Qmax,np.sum(scalelim1[ind_left:(ind_right+1)]<scale_l))

	# number of tau[k]'s on which the CWT is to be computed (for each scale), i.e. outside the Shannon-Nyquist exclusion zone.
	n_outside_scalelim1=np.ones(J,dtype=int)*Q
	J=scale.size
	for l in range(J):
		count=0
		for k in range(Q):
			if l>=scalelim1_ind[k]:
				count+=1
		n_outside_scalelim1[l]=count

	return scale,coi1,coi2,coi1_smooth,coi2_smooth,coi_smooth_ind,weight_cwt,scalelim1_ind,scalelim1_smooth,scalelim1_ind_smooth,Qmax,n_outside_scalelim1
