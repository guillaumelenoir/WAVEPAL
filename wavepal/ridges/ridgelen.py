import numpy as np

def ridgelen(id,ir,jr,fr):

	"""	RIDGELEN  Wavelet ridge length expressed as number of full cycles.

	LEN=RIDGELEN(ID,IR,JR,FR) determines the length of the ridges given by ID, IR, JR, and FR, with LEN expressed in number of cycles completed along the ridge.

	ID, IR and JR are arrays of ridge indices (type: int), and FR is the frequency of the wavelet transform value along the ridge (type: float), as output by RIDGEWALK.

	LEN is a row vector with the same size as the input arrays.
	__________________________________________________________________
	This is part of JLAB
	(C) 2009--2011 J.M. Lilly
	Rewritten in python 2.X by G. Lenoir, October 2016

	RIDGELEN is a low-level function called by RIDGEWALK."""


	lr=ridgelenloop(id,ir,jr,fr)

	return lr



def ridgelenloop(id,ir,jr,fr):
	
	import copy

	myid=copy.copy(id)
	if myid.size==0:
		myid=np.cumsum(ir==-999,0)
	index=(ir!=-999)
	lr=np.nan*ir
	if index.size>0:
		lr[index]=ridgelen1(myid[index],ir[index],jr[index],fr[index])
	return lr



def ridgelen1(id,ir,jr,fr):

	from blocknum import blocknum
	import copy

	myfr=copy.copy(fr)
	myfr=myfr/(2.*np.pi)      # Convert to cyclic frequency
	num,a,b=blocknum(id)
	myfr[np.isnan(myfr)]=0.
	ar=np.cumsum(myfr,0)     # ridge age
	lena=np.absolute(ar[b]-ar[a])			# 1-D array
	len1=np.zeros(id.shape)
	len1[a]=lena
	len1=copy.copy(np.cumsum(len1))
	len2=np.zeros(id.shape)
	myvec=np.zeros(lena.size)
	myvec[0]=0.
	myvec[1:]=lena[0:-1]
	len2[a]=copy.copy(myvec)
	len2=copy.copy(np.cumsum(len2))
	lr=len1-len2

	return lr
