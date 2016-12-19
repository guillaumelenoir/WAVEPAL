import numpy as np

def mov_av(x,l):
	
	""" mov_av returns the moving average of a time series. For each data point, it makes the average on l points, i.e. enter in the average: the (l-1)/2 data points before, the (l-1)/2 data points after, and the current data point. l must then be odd. On the borders of the time series, the formula is adapted: the average is taken on less data points (all the available ones); see the code for more details.
		Inputs:
		- x [1-dim numpy array of floats]: the data points of the time series
		- l [int]: the number of points on which the average is performed. If it is not odd, 1 is added to l.
		Outputs:
		- mov_av_x [1-dim numpy array of floats - size=x.size]: the moving averaged time series.
		-----------------------------
		This is part of WAVEPAL
		(C) 2016 G. Lenoir"""
	
	n=x.size
	mov_av_x=np.zeros(n)
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
	return mov_av_x
