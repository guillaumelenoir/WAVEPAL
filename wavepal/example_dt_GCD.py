# Example for the function dt_GCD.

import numpy as np
from distinct_ages import distinct_ages
from dt_GCD import dt_GCD
import sys
from os import path

here = path.abspath(path.dirname(__file__))
data=np.genfromtxt(here+"/data/Cheng_speleothems_nature2016/d18O.txt")
t=data[:,0]
mydata=data[:,1]
t, mydata=distinct_ages(t,mydata)
# There is, for each time in the file, maximum 3 numbers after the coma, so that I multiply by alpha=1000.
alpha=1000.0   # note that 10000., 100000., etc. works as well
tt=np.zeros(t.size,dtype=long)
for k in range(t.size):
	tt[k]=long(np.round(t[k]*alpha))
dt_GCD_int=dt_GCD(tt)
print "dt_GCD is ", float(dt_GCD_int)/alpha