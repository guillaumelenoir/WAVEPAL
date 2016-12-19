import numpy as np

def mov_av(t,x,l,type_av="t"):
	
	""" mov_av returns the moving average of a time series. 
		If type_av="n":
		For each data point, it makes the average on l points, i.e. enter in the average: the (l-1)/2 data points before, the (l-1)/2 data points after, and the current data point. l must then be odd. On the borders of the time series, the formula is adapted: the average is taken on less data points (all the available ones); see the code for more details.
		If type_av="t":
		For each data point x[k] at time t[k], it makes the average on a fixed time interval l. Enter in the average: all the data points s.t. their time is between max(t[0],t[k]-l/2). and min(t[-1],t[k]+l/2). The time series being possibly irregularly sampled, the number of data points on which the average is performed may differ from one interval to another. See the code for more details.
		Required Inputs:
		- t [1-dim numpy array of floats]: the times of the time series.
		- x [1-dim numpy array of floats - size=t.size]: the data points of the time series
		- l [int or float]: 
		If type_av="n", this the number of points on which the average is performed. If it is not odd, 1 is added to l. l must be of 'int' type.
		If type_av="t", this is the time interval on which tge average is performed. l must be of 'float' type.
		Optional Inputs:
		- type_av="t": is "t" (default) or "n".
		Outputs:
		- mov_av_x [1-dim numpy array of floats - size=t.size]: the moving averaged time series.
		RECOMMENDATION:
		Before running mov_av, run 'check_data' module of the 'Wavepal' class in order to remove the data points having the same times.
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""
			
	# check inputs
	try:
		assert np.issubsctype(t,float)
	except:
		print "Error at input 't': must be a numpy array of 'float' type"
		return
	try:
		assert np.issubsctype(x,float)
	except:
		print "Error at input 'x': must be a numpy array of 'float' type"
		return
	try:
		assert (type(type_av) is str) and ((type_av.lower()=="n") or (type_av.lower()=="t"))
	except AssertionError:
		print "Error at input 'type_av': must be 'n' or 't'"
		return
	type_av=type_av.lower()
	if type_av=="n":
		try:
			assert (type(l) is int) and l>0
		except AssertionError:
			print "Error at input 'l': must be of 'int' type and >0"
			return
	elif type_av=="t":
		try:
			assert (type(l) is float) and l>0.
		except AssertionError:
			print "Error at input 'l': must be of 'float' type and >0."
			return
	# Main code
	n=x.size
	mov_av_x=np.zeros(n)
	if type_av=="n":
		if l%2==0:
			lbis=l+1      # if even number => +1
		else:
			lbis=l
		m=(lbis-1)/2
		# on the borders of the time series
		for k in range(m):
			mov_av_x[k]=np.sum(x[:(k+m+1)])/float(k+1+m)
			mov_av_x[n-1-k]=np.sum(x[(n-1-k-m):])/float(k+1+m)
		# everywhere except on the border
		for k in range(m,n-m):
			mov_av_x[k]=np.sum(x[(k-m):(k+m+1)])/lbis
	elif type_av=="t":
		t_start=t[0]
		t_end=t[-1]
		m=l/2.0
		ind_left=0
		ind_right=0
		for k in range(n):
			t_left=max(t_start,t[k]-m)
			ind_left=ind_left+np.argmin(np.absolute(t[ind_left:]-t_left))
			t_right=min(t_end,t[k]+m)
			ind_right=ind_right+np.argmin(np.absolute(t[ind_right:]-t_right))
			mov_av_x[k]=np.sum(x[ind_left:ind_right+1])/float(ind_right+1-ind_left)
	return mov_av_x




