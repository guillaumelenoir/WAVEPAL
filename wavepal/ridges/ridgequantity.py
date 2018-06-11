import numpy as np

def ridgequantity(w,fs):

	"""RIDGEQUANTITY  The ``ridge quantity'' associated with a wavelet transform.

   WARNING: Type of ridge is here imposed to be "amplitude" (see the matlab version for extra choices).

   RQ,FW = RIDGEQUANTITY(W,FS) returns the ridge quantity RQ associated with the analytic wavelet transform W computed at frequencies FS. RQ is the same size as W. It also returns the transform frequency matrix FW. FW has the radian instantaneous frequency at each scale and is the same size as W.
 
   For details see

       Lilly and Olhede (2010).  On the analytic wavelet transform.
           IEEE Trans. Info. Theory.
   _____________________________________________________________________

   RIDGEQUANTITY is a low-level function called by RIDGEWALK.

   See also RIDGEWALK

   Usage: rq, fw = ridgequantity(w,fs)
   __________________________________________________________________
   This is part of JLAB
   (C) 2009 J.M. Lilly
   Rewritten in python 2.X by G. Lenoir, October 2016"""
	
	from vdiff import vdiff
	from frac import frac
 
	om=np.zeros(w.shape)
	for k in range(w.shape[0]):
		om[k,:]=fs[:]

	#This is dw/ds
	rq=om**2*frac(vdiff(w,2),vdiff(om,2))*frac(np.ones(w.shape),w)
	
	return rq, om
