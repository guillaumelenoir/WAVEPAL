import numpy as np

def colbreaks(I1,I2,I3,I4,I5):

	"""COLBREAKS  Insert NANs into discontinuties in a vector.

	Usage: O1, O2, O3, O4, O5=COLBREAKS(I1, I2, I3, I4, I5)
	
	I1, I2 and I3 are vectors with integers. I4 and I5 are vectors with floats. They are all vectors of the same length.
	
	It inserts -999 into I1 at its discontinuties (and puts one -999 at the end). Example:
	A numpy array [2,2,2,3,3] becomes [2,2,2,-999,3,3,-999]
	
	It then uses I1 as a reference for the other vectors and inserts NANs or -999 into all the vectors at locations where I1 is discontinous. If type=="float" => inserts NAN (I4 and I5) ; type=="int" => inserts -999 (I2 and I3).
	
	_________________________________________________________________
	This is part of JLAB
	(C) 2000--2008 J.M. Lilly
	Rewritten in python 2.X by G. Lenoir, October 2016"""

	mybreak=np.diff(I1)
	mybreak=np.insert(mybreak,0,0)
	myind=np.where(mybreak!=0)
	I1_out=np.append(np.insert(I1,myind[0],-999),-999)
	I2_out=np.append(np.insert(I2,myind[0],-999),-999)
	I3_out=np.append(np.insert(I3,myind[0],-999),-999)
	I4_out=np.append(np.insert(I4,myind[0],np.nan),np.nan)
	I5_out=np.append(np.insert(I5,myind[0],np.nan),np.nan)

	return I1_out, I2_out, I3_out, I4_out, I5_out
