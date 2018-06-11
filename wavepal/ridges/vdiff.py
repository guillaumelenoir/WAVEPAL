import numpy as np

def vdiff(x,dim):

	"""VDIFF	Length-preserving first central difference.

   DX=VDIFF(X,DIM) differentiates X along dimension DIM using the first
   central difference; DX is the same size as X.
   
   !!! Only works for 1-D and 2-D arrays
   _____________________________________________________________________

   It uses the first forwards / first backwards difference at the first and last point, respectively.
   _____________________________________________________________________

   Usage:  x=vdiff(x,dim)
   __________________________________________________________________
   This is part of JLAB
   (C) 2000--2011 J.M. Lilly
   Rewritten in python 2.X by G. Lenoir, October 2016"""
 
	y=np.zeros(x.shape)
	if x.ndim==1:
		if dim==1:
			if x.size<=1:
				print "Error in vdiff.py: length too small to perform numerical derivatives"
				return
			else:
				y[0]=x[1]-x[0]
				for k in range(1,x.size-1):
					y[k]=(x[k+1]-x[k-1])/2.
				y[-1]=x[-1]-x[-2]
		else:
			print "Error in vdiff.py"
			return
	elif x.ndim==2:
		if dim==1:
			if x.shape[0]<=1:
				print "Error in vdiff.py: length too small to perform numerical derivatives"
				return
			else:
				y[0,:]=x[1,:]-x[0,:]
				for k in range(1,x.shape[0]-1):
					y[k,:]=(x[k+1,:]-x[k-1,:])/2.
				y[-1,:]=x[-1,:]-x[-2,:]
		elif dim==2:
			if x.shape[1]<=1:
				print "Error in vdiff.py: length too small to perform numerical derivatives"
				return
			else:
				y[:,0]=x[:,1]-x[:,0]
				for k in range(1,x.shape[1]-1):
					y[:,k]=(x[:,k+1]-x[:,k-1])/2.
				y[:,-1]=x[:,-1]-x[:,-2]
		else:
			print "Error in vdiff.py"
			return
	else:
		print "Error in vdiff.py"
		return
	return y
