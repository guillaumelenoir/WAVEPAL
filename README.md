# WAVEPAL 3.3 python package

WAVEPAL is a package, written in python2.X, that performs frequency (WOSA Lomb-Scargle-based periodogram) and time-frequency (Smoothed Morlet wavelet scalogram, as an extension of the on the Lomb-Scargle periodogram) analyses of irregularly sampled time series with a polynomial trend, WITHOUT INTERPOLATING the data. Moreover, it proposes significance testing against a large choice of background processes, i.e. CARMA(p,q) processes. It also provides tools for filtering the signal (band and ridges filtering).

--------------

## References

The provisional reference for the package is
> **G. Lenoir and M. Crucifix.** *Wavepal: A python software for wavelet analysis of irregularly sampled time series.* Poster presentation at the 1st CVAS workshop: What do we know about multicentennial, multimillenial variability?, Guesthouse of the University of Hamburg, November 2016.

and is available [here](http://www.elic.ucl.ac.be/users/lenoir/mywebsite/docs/poster_CVAS_2016.pdf). Two reference papers are currently under writing. 

When `WAVEPAL` performs significance testing against a CARMA(p,q) background noise, it calls the python package [carma_pack](https://github.com/brandonckelly/carma_pack) (except when p=q=0: the case of a white noise background is entirely analysed within `WAVEPAL`), whose reference is
**B. C. Kelly, A. C. Becker, M. Sobolewska, A. Siemiginowska, and P. Uttley.** *Flexible and scalable methods for quantifying stochastic variability in the era of massive time-domain astronomical data sets.* The Astrophysical Journal, 788(1):33, 2014.

Ridges filtering is based on functions that I have rewritten from the matlab package [jlab](http://www.jmlilly.net/jmlsoft.html). Main reference is
**J. Lilly and S. Olhede.** *On the analytic wavelet transform.* IEEE Transactions on Information Theory, 56(8):4135–4156, aug. 2010.

---------------

## Installation 

Python 2.X and all the dependencies of WAVEPAL are installed. Feel free to modify the installation file. First, you need to download the package on your machine. Then, follow the steps explained below. Installation ends with the message "INSTALLATION COMPLETED SUCCESSFULLY".

### On Linux (Ubuntu)

Just enter
```
sh Linux_install.sh
```
and that is all!

------------

## Examples

See the the file `Quick_start.ipynb` in the `examples/` folder. 

---------------

## Dependencies

WAVEPAL depends on the following python modules:
* numpy      (for core functionality)
* scipy      (for core functionality)
* matplotlib (for generating plots)
* acor       (for calculating the autocorrelation time scale of MCMC samples; this is used only by carma_pack)
* [tqdm](https://pypi.python.org/pypi/tqdm) (make progress meters)
* [carma_pack](https://github.com/brandonckelly/carma_pack) (MCMC sampler for performing Bayesian inference on continuous-ARMA (CARMA) models) 

Moreover, carma_pack requires the instatllation of [BOOST](http://www.boost.org) (for linking python and C++) and [ARMADILLO](http://arma.sourceforge.net) (C++ linear algebra library).

--------------

## Source Code

All the source code of the WAVEPAL package is provided in the 'wavepal/' folder. The main code is in Wavepal.py, that contains the class 'Wavepal'. However, you should not work directly with those files since you just need to import the package once it is installed (see the examples provided). 
The folder 'carmcmc' contains 2 sligthly modified python files from carma_pack, where I only modified some bunch of code related to the figures. 
The folder 'wavepal/ridges/' contains functions rewritten from [jlab](http://www.jmlilly.net/jmlsoft.html). They were transformed from matlab to python and strongly modified. 

## Acknowledgements

Many thanks to [Pierre-Yves Barriat](https://be.linkedin.com/in/pybarriat) for having built the installation files. 


