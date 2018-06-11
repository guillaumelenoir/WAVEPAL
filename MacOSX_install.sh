#!/bin/bash

mypath=$(pwd)

### NEED MACPORTS
# https://www.macports.org/install.php

type port >/dev/null 2>&1 || { echo >&2 "I require port but it's not installed.  Aborting."; exit 1; }

### NEED ROOT PRIVILEGES ###
# MacOS 10.12 Sierra
sudo port install python27
type xcodebuild >/dev/null 2>&1 || { echo >&2 "I require xcodebuild but it's not installed.  Aborting."; exit 1; }
sudo port install py27-gnureadline
sudo port install py27-numpy
sudo port install py27-scipy
sudo port install py27-matplotlib
sudo port install py27-acor
sudo port install armadillo
sudo port install py27-tqdm
############################

### install an old version of boost (recent versions not compatible with carma pack)
cd /opt/local    # it seems it does not work if it is run in $mypath
sudo git clone --single-branch https://github.com/macports/macports-ports.git macports-ports_temp_wavepal
cd macports-ports_temp_wavepal
sudo git checkout fb1a595a492c02cc70476aa71db416d20397d71c   # go back to the commit of version 1.59
cd devel/boost
sudo port install +python27
cd ../../..
sudo rm -r macports-ports_temp_wavepal

cd $mypath
echo "Install carma_pack"
git clone https://github.com/brandonckelly/carma_pack.git
rsync -av --delete carmcmc carma_pack/src
cd carma_pack/src
BOOST_DIR="/opt/local" ARMADILLO_DIR="/opt/local" NUMPY_DIR="/opt/local" python2.7 setup.py install --user
cd ../..

echo "Install WAVEPAL"
python2.7 setup.py install --user

mv carmcmc carmcmc_pack
mv wavepal wavepal_pack
sudo rm -r carma_pack

echo "****************"
echo "* Test on data *"
echo "****************"
cd test
python2.7 test.py
cd figures
rm -f *

#Done
