#!/bin/bash

### NEED ROOT PRIVILEGES ###
# Ubuntu 16.04 Xenial
sudo apt-get install git bison make g++ devscripts g++-multilib cmake python-pip python-dev libboost-all-dev python-scipy python-matplotlib libopenblas-base libarpack++2-dev libarmadillo-dev
# CentOs 7
#sudo yum install git bison make gcc-c++ cmake python-devel python-pip scipy python-matplotlib boost-devel blas-devel lapack-devel 
############################

pip install --upgrade pip --user
pip install acor --user
pip install tqdm --user

echo "Install carma_pack"
git clone https://github.com/brandonckelly/carma_pack.git
rsync -av --delete carmcmc carma_pack/src
cd carma_pack/src
BOOST_DIR="/usr" ARMADILLO_DIR="/usr" NUMPY_DIR="/usr" python setup.py install --user
cd ../..

echo "Install WAVEPAL"
python setup.py install --user

echo "****************"
echo "* Test on data *"
echo "****************"
cd test
python test.py
cd figures
rm -f *

#Done
