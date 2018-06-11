import numpy as np

def ridgewalk(W,fs,N=1.5,chi=0.,alpha=1./4.):

	"""RIDGEWALK  Extract wavelet transform ridges
	
   IR,JR,XR,FR = RIDGEWALK(W,FS) where W is a wavelet transform matrix at frequecies FS, returns the wavelet ridges of transform W.

   The columns of W correspond to different frequencies, specified by the
   frequency array FS, at which the wavelet transform was performed.  Note
   that FS assumes a unit sample rate.

   The frequencies FS are expected to be ordered from highest to lowest.

   RIDGEWALK returns the following quantities along ridges

       IR     Ridge indices into rows of W (time)
       JR     Ridge indices into columns of W (scale)
       XR     Estimated signal along the ridge
       FR     Transform frequency values in radian frequency

   All output variables are vectors of the same length.  A NAN appears in
   each value at the bottom of every ridge to separate different ridges.
   _______________________________________________________________________

   Options

   RIDGEWALK(...,N=1.5,CHI=0.,alpha=1./4.) specifies options for the ridge computation.

        N  -- Removes all ridges of less than N periods in length
      CHI  -- Removes all small amplitude ridge points having |W|<CHI
	  ALPHA  -- Controls agressiveness of chaining across scales. Increase ALPHA to chain ridges more agressively across scales, or decrease ALPHA to supress chaining across scales.
   _______________________________________________________________________
   
   This is part of JLAB
   (C) 2004--2011 J.M. Lilly
   Rewritten in python 2.X by G. Lenoir, October 2016
   ___________________________________________________________________

   ALPHA is defined as a normalized frequency difference

         ALPHA  =  DOMEGA / OMEGA

   where OMEGA is the transform frequency, and DOMEGA is the difference
   between the frequency predicted for the next point based on the
   transform at a "tail", and the actual frequency at prospective "heads".

   The chaining parameter is defined in such a way that it does not
   need to be changed as time sampling or frequency sampling changes.
   However, for strongly chirping signals or weakly chirping, noisy
   signals, better performance may be obtianed by adjusting it.
   ______________________________________________________________

   Spurious ridge points

   RIDGEWALK has a continigency for rejecting spurious ridge points.
   These tend to occur on the flanks of interesting signals, and
   reflect the wavelet structure rather than the signal structure.

   See ISRIDGEPOINT for details."""

	from colbreaks import colbreaks
	fs_sorted_descending=np.sort(fs)[::-1]
	try:
	   assert np.prod(fs==fs_sorted_descending)==1.
	except AssertionError:
		print "Error: The frequencies FS should be sorted from highest to lowest."
		return
	mybool,rq,om=isridgepoint(W,fs,chi)
	print "RIDGEWALK chaining ridges."
	id,ir,jr,wr,fr=ridgechains(fs,N,mybool,W,om,alpha)
	if id.size>0:
		_,ir,jr,wr,fr=colbreaks(id,ir,jr,wr,fr)
		# I skip "ridgeinterp", which should be called here

	print "RIDGEWALK finished."

	return ir,jr,wr,fr



def isridgepoint(w,fs,chi):

	"""ISRIDGEPOINT  Finds wavelet ridge points.

   BOOL, RQ, OM = ISRIDGEPOINT(W,FS,CHI) where W is a wavelet transform matrix at
   *cyclic* frequecies FS, finds all ridge points of W with amplitudes
   |W| exceeding the amplitude cutoff A.

   BOOL is a matrix of the same size as W, which is equal to one for
   those elements of W which are ridge points, and zero otherwise.

   ISRIDGEPOINT rejects spurious ridge points.
   These tend to occur on the flanks of interesting signals, and
   reflect the wavelet structure rather than the signal structure.

   A ridge point is considered spurious if either it is located at an
   amplitude minima, or if the frequency anomaly (transform frequency
   minus scale frequency) is a maximum.

   See also RIDGEQUANTITY, RIDGEWALK.
   __________________________________________________________________
   This is part of JLAB
   (C) 2006--2009 J.M. Lilly
   Modified by G. Lenoir, October 2016"""

	from ridgequantity import ridgequantity
	
	print "RIDGEWALK looking for ridge points..."

	rq,om=ridgequantity(w,fs)
	rqm=np.zeros(rq.shape)
	rqm[:,0]=rq[:,-1]
	rqm[:,1:]=rq[:,0:-1]
	rqp=np.zeros(rq.shape)
	rqp[:,-1]=rq[:,0]
	rqp[:,0:-1]=rq[:,1:]

	#This is d/ds < 0 since scale decreases in columns
	mybool=np.logical_or(np.logical_and(rqm<0.,rqp>=0.),np.logical_and(rqm<=0.,rqp>0.))
	
	err=np.absolute(rq)

	#Ensure maximum not minimum
	vshiftboolp1=np.zeros(mybool.shape,dtype=bool)
	vshiftboolp1[:,-1]=mybool[:,0]
	vshiftboolp1[:,0:-1]=mybool[:,1:]
	vshiftboolm1=np.zeros(mybool.shape,dtype=bool)
	vshiftboolm1[:,0]=mybool[:,-1]
	vshiftboolm1[:,1:]=mybool[:,0:-1]
	vshifterrp1=np.zeros(err.shape)
	vshifterrp1[:,-1]=err[:,0]
	vshifterrp1[:,0:-1]=err[:,1:]
	vshifterrm1=np.zeros(err.shape)
	vshifterrm1[:,0]=err[:,-1]
	vshifterrm1[:,1:]=err[:,0:-1]
	mybool[np.logical_and(np.logical_and(mybool,vshiftboolp1),err>vshifterrp1)]=False
	mybool[np.logical_and(np.logical_and(mybool,vshiftboolm1),err>vshifterrm1)]=False
	
	del vshiftboolp1, vshiftboolm1, vshifterrp1, vshifterrm1

	bool1=np.logical_not(np.isnan(w))  #Remove NANs
	bool2=np.logical_not((np.absolute(w)<chi))  #Remove those less than cutoff amplitude
	mybool=mybool*bool1*bool2
	mybool[:,0]=False
	mybool[:,-1]=False

	print "RIDGEWALK found "+str(np.where(mybool)[0].size)+" ridge points."

	return mybool, rq, om



def ridgechains(fs,N,mybool,x,f,alpha):

	"""RIDGECHAINS  Forms ridge curves by connecting transform ridge points.

   ID, II, JJ, XR, FR = RIDGECHAINS(fs, N, BOOL, x, f, alpha) forms chains of ridge points
   of wavelet transform W.

   Ridge points are all points of W for which BOOL, a matrix of the
   same size as x, equals one.  Only ridges of at least N periods in
   length are returned.
   
   ID is a unique ID number assigned to each ridge.  II and JJ are
   the time- and scale-indices along the ridges.  XR is the wavelet
   transform along the ridge.
 
   All output variables are the same size.
 
   See also RIDGEWALK.
   __________________________________________________________________
   This is part of JLAB
   (C) 2006--2007 J.M. Lilly
   Modified by G. Lenoir, October 2016"""
 
	import copy
	from frac import frac
	from ridgelen import ridgelen
	from vdiff import vdiff
 
	if np.where(mybool)[0].size==0:
		id=np.zeros(0)
		ii=np.zeros(0)
		jj=np.zeros(0)
		xr=np.zeros(0)
		fr=np.zeros(0)
		return
    
	dfdt=vdiff(f,1)

	# WARNING: In matlab, 'find' returns indices under LINEAR form (i.e. it's a vector). In python, np.where returns indices under the same shape as the array, e.g. as tuples (i,j) for a 2-D matrix. np.ravel_multi_index allows to transform those tuples to a LINEAR form, but, on the contrary to matlab, iteration is performed on ROWS (in matlab, 'find' iterates on COLUMNS).

	where_bool=np.where(mybool)     # to get indices going through ROWS first
	ii=where_bool[0]    # NB: on the contrary to matlab, sorting to ascending order is direct here because we iterate on ROWS with np.where
	jj=where_bool[1]

	#Using new algorithm as of November 2007, faster and also prevents ridge breaking
	xr=x[where_bool]			  #Transform value along ridge
	fr=f[where_bool]              #Frequency along ridge
	fsr=fs[jj]					  #Scale frequency along ridge
	fr_next=fr+dfdt[where_bool]   #Predicted frequency at next point
	fr_prev=fr-dfdt[where_bool]   #Predicted frequency at previous point

	cumbool=np.cumsum(mybool,1)
	J=np.amax(cumbool[:,-1])

	indexmat=-999*np.ones((f.shape[0],J),dtype=int)
	nextindexmat=-999*np.ones((f.shape[0],J),dtype=int)
	iimat=-999*np.ones((f.shape[0],J),dtype=int)
	jjmat=-999*np.ones((f.shape[0],J),dtype=int)
	fsmat=np.nan*np.ones((f.shape[0],J))
	frmat=np.nan*np.ones((f.shape[0],J))
	fr_nextmat=np.nan*np.ones((f.shape[0],J))
	fr_prevmat=np.nan*np.ones((f.shape[0],J))

	#Indices for this point
	for k in range(ii.size):
		indexmat[ii[k],cumbool[where_bool][k]-1]=k
	nonanindex_mat=np.where(np.transpose(indexmat)!=-999)[::-1]   # to get indices going throgh COLUMNS first

	iimat[nonanindex_mat]=copy.copy(ii[indexmat[nonanindex_mat]])
	jjmat[nonanindex_mat]=copy.copy(jj[indexmat[nonanindex_mat]])
	fsmat[nonanindex_mat]=copy.copy(fsr[indexmat[nonanindex_mat]])
	frmat[nonanindex_mat]=copy.copy(fr[indexmat[nonanindex_mat]])
	fr_nextmat[nonanindex_mat]=copy.copy(fr_next[indexmat[nonanindex_mat]])
	fr_prevmat[nonanindex_mat]=copy.copy(fr_prev[indexmat[nonanindex_mat]])

	#Time difference from points here to next points
	vshiftiimat=np.zeros(iimat.shape)
	vshiftiimat[-1,:]=iimat[0,:]
	vshiftiimat[0:-1,:]=iimat[1:,:]
	dii=np.sum(vshiftiimat-iimat,1)
	del iimat, vshiftiimat

	#Scale frequency difference from points here to next points
	fsmat3=np.zeros((fsmat.shape[0],fsmat.shape[1],J))
	frmat3=np.zeros((frmat.shape[0],frmat.shape[1],J))
	for k in range(J):
		fsmat3[:,:,k]=copy.copy(fsmat)
		frmat3[:,:,k]=copy.copy(frmat)

	#Predicted minus actual frequency at this point
	fr_nextmat3=np.zeros((fr_nextmat.shape[0],fr_nextmat.shape[1],J))
	for k in range(J):
		fr_nextmat3[:,:,k]=copy.copy(fr_nextmat)
	vshiftfrmat3=np.zeros(frmat3.shape)
	vshiftfrmat3[-1,:,:]=frmat3[0,:,:]
	vshiftfrmat3[0:-1,:,:]=frmat3[1:,:,:]
	df1=np.transpose(vshiftfrmat3,(0,2,1))-copy.copy(fr_nextmat3)
	df1=frac(df1,frmat3)

	del fr_nextmat3, vshiftfrmat3

	#Expected minus actual frequency at next point
	fr_prevmat3=np.zeros((fr_prevmat.shape[0],fr_prevmat.shape[1],J))
	for k in range(J):
		fr_prevmat3[:,:,k]=copy.copy(fr_prevmat)
	vshiftfr_prevmat3=np.zeros(fr_prevmat3.shape)
	vshiftfr_prevmat3[-1,:,:]=fr_prevmat3[0,:,:]
	vshiftfr_prevmat3[0:-1,:,:]=fr_prevmat3[1:,:,:]
	df2=np.transpose(vshiftfr_prevmat3,(0,2,1))-copy.copy(frmat3)
	df2=frac(df2,frmat3)

	df=(np.absolute(df1)+np.absolute(df2))/2.

	df[df>alpha]=np.nan

	del fr_prevmat3, frmat3, df1, df2, vshiftfr_prevmat3

	#Keep when they are the same
	#Set df to nan except when min along one direction
	mindf=np.nanmin(df,axis=1)
	jjmin=np.zeros((df.shape[0],df.shape[2]),dtype=int)
	for k in range(df.shape[0]):
		for j in range(df.shape[2]):
			try:
				jjmin[k,j]=np.nanargmin(df[k,:,j])
			except ValueError:
				print "WARNING: ValueError with numpy.nanargmin: All-NaN slice encountered"
				jjmin[k,j]=0
	iimin=np.zeros((df.shape[0],df.shape[1]),dtype=int)
	myvec=range(df.shape[0])
	for k in range(df.shape[1]):
		iimin[:,k]=copy.copy(myvec)
	kkmin=np.zeros((df.shape[0],df.shape[1]),dtype=int)
	myvec=range(df.shape[1])
	for k in range(df.shape[0]):
		kkmin[k,:]=copy.copy(myvec)
	df=np.nan*df
	df[iimin,jjmin,kkmin]=copy.copy(mindf)     # N.B.: No need to make "squeeze" with python
	del myvec, iimin, jjmin, kkmin

	mindf=np.nanmin(df,axis=2)
	jjmin=np.zeros((df.shape[0],df.shape[1]),dtype=int)
	for k in range(df.shape[0]):
		for j in range(df.shape[1]):
			try:
				jjmin[k,j]=np.nanargmin(df[k,j,:])
			except ValueError:
				print "WARNING: ValueError with numpy.nanargmin: All-NaN slice encountered"
				jjmin[k,j]=0
	index=np.where(np.logical_not(np.isnan(np.transpose(mindf))))[::-1]

	if index[0].size>0:
		# N.B.: mindf, ii2 and jj2 are not used afterwards
		index0_where=np.where(index[0]<x.shape[0]-1)     # 1-D array
		ii2=copy.copy(index[0][index0_where])
		jj2=copy.copy(index[1][index0_where])
		index=[None]*2
		index[0]=np.zeros(ii2.size,dtype=int)
		index[1]=np.zeros(ii2.size,dtype=int)
		for k in range(ii2.size):
			index[0][k]=ii2[k]
			index[1][k]=jj2[k]
		index=tuple(index)
		nextindexmat[index]=indexmat[ii2+1,jjmin[index]]

	id=copy.copy(ii)*-999

	for i in range(nextindexmat.shape[0]):
		id[indexmat[i,np.where(indexmat[i,:]!=-999)]]=indexmat[i,np.where(indexmat[i,:]!=-999)]

	for i in range(nextindexmat.shape[0]):
		id[nextindexmat[i,np.where(nextindexmat[i,:]!=-999)]]=id[indexmat[i,nextindexmat[i,:]!=-999]]
	sorter=np.argsort(id,kind='mergesort')        # !!! kind='mergesort' is ABSOLUTELY necessary, to make it stable. Ridgewalk fails with the default 'quicksort' algo.
	id=np.sort(id)
	ii=copy.copy(ii[sorter])
	jj=copy.copy(jj[sorter])
	xr=copy.copy(xr[sorter])
	fr=copy.copy(fr[sorter])

	#Remove ridge lines of length shorter than a specified length
	lr=ridgelen(id,ii,jj,fr)
	mysorter=np.where(lr>=N)
	id=copy.copy(id[mysorter])
	ii=copy.copy(ii[mysorter])
	jj=copy.copy(jj[mysorter])
	xr=copy.copy(xr[mysorter])
	fr=copy.copy(fr[mysorter])
					
	print "RIDGEWALK pruning to "+str(id.size)+" ridge points."

	return id,ii,jj,xr,fr
