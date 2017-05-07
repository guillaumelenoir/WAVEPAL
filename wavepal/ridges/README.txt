All the code of this package is adapted from jlab (https://github.com/jonathanlilly/jLab). jlab is available in Matlab, and all the code relative to « ridgewalk.m » and its dependencies has been transformed to Python2.X. 

N.B.: To call directly this package from Python: import wavepal.ridges

WARNING:
————————
The transformed code (Matlab -> Python) is much less general than the original Matlab code:
- Only the AMPLITUDE ridges are available. Indeed, PHASE ridges are not convenient for irregularly sampled time series. 
- There is no ridges interpolation, i.e. the call to « ridgeinterp » at the end of ridgewalk is skipped. 
- Ridgewalk and its dependencies were transformed in a « minimal » way, i.e. this only works with ridgewalk for the amplitude ridge analysis of the CWT-amplitude of a 1-D signal. Compared to the original Matlab version, the functions in this folder have been greatly simplified or skipped (and replaced by a simple bunch of code where they were previously called). Consequently, MOST OF THOSE Python FUNCTIONS MUST NOT BE USED IN ANOTHER CONTEXT. 
- Description of the functions was strongly modified to match the Python version. 