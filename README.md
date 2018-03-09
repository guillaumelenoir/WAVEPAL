# WAVEPAL python package

WAVEPAL is a package, written in python 2.X, that performs frequency and time-frequency analyses of irregularly sampled time series with a polynomial trend, **WITHOUT INTERPOLATING** the data. Frequency analysis is based on the Lomb-Scargle periodogram and WOSA smoothing method. Time-frequency analysis is performed with the smoothed Morlet wavelet scalogram, defined as an extension of the Lomb-Scargle periodogram. Moreover, the package proposes significance testing against a large choice of background processes, i.e. CARMA(p,q) processes. It also provides tools for filtering the signal (band and ridges filtering).

--------------

## References

The two references are:
> **G. Lenoir and M. Crucifix.** [*A general theory on frequency and time–frequency analysis of irregularly sampled time series based on projection methods – part 1: Frequency analysis.*](http://dx.doi.org/10.5194/npg-25-145-2018) Nonlinear Processes in Geophysics, 25(1):145–173, 2018.

> **G. Lenoir and M. Crucifix.** [*A general theory on frequency and time–frequency analysis of irregularly sampled time series based on projection methods – part 2: Extension to time–frequency analysis.*](http://dx.doi.org/10.5194/npg-25-175-2018) Nonlinear Processes in Geophysics, 25(1):175–200, 2018.

When `WAVEPAL` performs significance testing against a CARMA(p,q) background noise, it calls the python package [carma_pack](https://github.com/brandonckelly/carma_pack) (except when p=q=0: the case of a white noise background is entirely analysed within `WAVEPAL`), whose reference is
> **B. C. Kelly, A. C. Becker, M. Sobolewska, A. Siemiginowska, and P. Uttley.** *Flexible and scalable methods for quantifying stochastic variability in the era of massive time-domain astronomical data sets.* The Astrophysical Journal, 788(1):33, 2014.

Ridges filtering is based on functions that I have rewritten from the matlab package [jlab](https://github.com/jonathanlilly/jLab). The main reference is
> **J. Lilly and S. Olhede.** *On the analytic wavelet transform.* IEEE Transactions on Information Theory, 56(8):4135–4156, aug. 2010.

---------------

## Installation 

First, you need to download the package on your machine (click on "clone or download", download the zip file, and unzip it). Then, follow the steps explained below. Installation ends with the message "INSTALLATION COMPLETED SUCCESSFULLY".

### On Linux (Ubuntu)

Just enter
```
sh Linux_install.sh
```
and that is all!

Python 2.X and all the dependencies of WAVEPAL are installed. Feel free to modify the installation file. 

### On Mac OS X (tested on Mac OS 10.9 Mavericks, 10.10 Yosemite, and 10.12 Sierra)

You first need to install MacPorts from [this website](https://www.macports.org/install.php). You also need xcodebuild, but this should be installed if you use the [terminal](https://en.wikipedia.org/wiki/Terminal_(macOS)) (if not, you will be asked to install it the first time you use the terminal). 

Open the [terminal](https://en.wikipedia.org/wiki/Terminal_(macOS)) and change directory (with command `cd`) to the folder containing the file `MacOSX_install.sh`.

Just enter 
```
sh MacOSX_install.sh
```
and that is all!

Python 2.7 and all the dependencies of WAVEPAL are installed. Type `python2.7` instead of `python` at the terminal when you call python. Feel free to modify the installation file.

### Tips
In the last part of the installation process, a test file is run to check that all is okay. It is indicated by 
```
****************
* Test on data *
****************
```
Just before the test, the paths towards the built-in packages are showed, and python may choose to take old versions, that you had already installed. If you encounter troubles with some packages during the test, I recommend to delete old versions (whose path is given), and reinstall WAVEPAL.

If you reinstall WAVEPAL from the same folder downloaded from github at the first installation, first delete the following files and folder:
* build
* dist
* WAVEPAL.egg-info

and rename: 
* 'carmcmc_pack' to 'carmcmc'
* 'wavepal_pack' to 'wavepal'

------------

## Examples

Open the file `Quick_start.ipynb` in the `examples/` folder, or click [here](https://github.com/guillaumelenoir/WAVEPAL/blob/master/examples/Quick_start.ipynb)

---------------

## Dependencies

WAVEPAL depends on the following python modules:
* numpy      (for core functionality)
* scipy      (for core functionality)
* matplotlib (for generating plots)
* acor       (for calculating the autocorrelation time scale of MCMC samples; this is used only by carma_pack)
* [tqdm](https://pypi.python.org/pypi/tqdm) (for progress meters)
* [carma_pack](https://github.com/brandonckelly/carma_pack) (MCMC sampler for performing Bayesian inference on continuous-ARMA (CARMA) models) 

Moreover, carma_pack requires the installation of [BOOST](http://www.boost.org) (for linking python and C++) and [ARMADILLO](http://arma.sourceforge.net) (C++ linear algebra library).

--------------

## Source Code

All the source code of the WAVEPAL package is provided in the `wavepal/` folder. The main code is in `Wavepal.py`, that contains the class `Wavepal`. However, you should not work directly with those files since you just need to import the package once it is installed (see the examples provided). 

The folder `carmcmc/` contains 2 sligthly modified python files from [carma_pack](https://github.com/brandonckelly/carma_pack), where I only modified some bunch of code related to the figures. 

The folder `wavepal/ridges/` contains functions rewritten from [jlab](http://www.jmlilly.net/jmlsoft.html). They were transformed from matlab to python and strongly modified. 

## Acknowledgements

Many thanks to [Pierre-Yves Barriat](https://be.linkedin.com/in/pybarriat) for having built the installation files. 


