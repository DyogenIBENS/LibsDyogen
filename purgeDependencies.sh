#!/bin/bash
# https://askubuntu.com/questions/761592/unable-to-apt-get-dist-upgrade-on-a-persistent-ubuntu-16-04-usb
sudo apt-get remove --purge git 
sudo apt-get remove --purge cython
sudo apt-get remove --purge python-numpy
sudo apt-get remove --purge python-numpy
sudo apt-get remove --purge python-scipy
sudo apt-get remove --purge python-matplotlib
sudo apt-get remove --purge cmake 
sudo apt-get remove --purge g++ 
sudo apt-get remove --purge libpng12-dev 
sudo apt-get remove --purge zlib1g-dev
sudo apt-get remove --purge mpi
