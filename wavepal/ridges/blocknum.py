import numpy as np

def blocknum(x):

	"""BLOCKNUM  Numbers the contiguous blocks of an array.

   Suppose X is a column vector which contains blocks of identical
   values, e.g. X=[0 0 0 1 1 3 2 2 2]

   N, A, B=BLOCKNUM(X) gives each contiguous block of identical values a
   unique number, in order beginning with 0.  Each elements of N
   specifies the number of the block to which the corresponding
   element of X belongs.

   In the above example,  N=[0 0 0 1 1 2 3 3 3]   (under a numpy array)

   It also returns numpy arrays A and B which are indices
   into the first and last point, respectively, of each block.  In
   the above example, A=[0 3 5 6] and B=[2 4 5 8]

   Usage: num,a,b=blocknum(x)
   _________________________________________________________________
   This is part of JLAB
   (C) 2000, 2004 J.M. Lilly
   Rewritten in python 2.X by G. Lenoir, October 2016"""

	import copy

	# In PYTHON: only works for x which is a 1-D array
	
	dx=np.diff(x)
	index=np.where(np.absolute(dx)>0)
	num=np.zeros(x.size,dtype=int)
	if index[0].size>0:
		index[0][:]+=1
		num[index]=1
	num=copy.copy(np.cumsum(num))
	a=np.zeros(1,dtype=int)
	a[0]=0
	myvec=np.where(np.diff(num)!=0)
	myvec[0][:]+=1
	a=np.append(a,myvec)
	b=np.where(np.diff(num)!=0)
	b=np.append(b,num.size-1)

	return num, a, b
