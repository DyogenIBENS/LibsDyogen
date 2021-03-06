#!/bin/bash

# save initial working dir
iwd=$(pwd)

# Install dependencies
# we assume that apt-get is the package manager
# enable the 'universe' deposit
sudo add-apt-repository universe
sudo apt-get update
# core dependencies: python2.7 git cython
sudo apt-get install python2.7 
sudo apt-get install git 
sudo apt-get install cython 
# marginal dependencies: matplotlib scipy numpy
# you may need to change the name of some package depending on your distribution
sudo apt-get install python-numpy
# Next packages needs to enable the deposit 'Universe'
sudo apt-get install python-matplotlib 
sudo apt-get install python-scipy 

PATH_PARENT_ALL="/home/${USER}/Libs"
echo "LibsDyogen and plugged softwares will be installed in ${PATH_PARENT_ALL}" >&2
mkdir -p ${PATH_PARENT_ALL}

# Install LibsDyogen
echo 'Install LibsDyogen' >&2
PATH_LIBSDYOGEN="${PATH_PARENT_ALL}/LibsDyogen"
git clone https://github.com/DyogenIBENS/LibsDyogen ${PATH_LIBSDYOGEN}
# Then cythonise *.pyx
bash ${PATH_LIBSDYOGEN}/cythonisePyxFiles.sh ${PATH_LIBSDYOGEN}


# Install homology teams
echo 'Install homology teams' >&2
cd ${PATH_PARENT_ALL}
wget http://euler.slu.edu/~goldwasser/homologyteams/homologyteams-1.1.zip
unzip homologyteams-1.1.zip
cd homologyteams-1.1/src
# compile with gcc
make
# To plug homologyteams to LibsDyogen
# sed -i "/PATH_HOMOLOGYTEAMS_BIN =/c\PATH_HOMOLOGYTEAMS_BIN = \"${PATH_PARENT_ALL}/homologyteams-1.1/src/homologyteams\"" ${PATH_LIBSDYOGEN}/utils/myGeneTeams.py

# Install i-ADHoRe 3.0
echo 'Install i-ADHoRe 3.0' >&2
cd ${PATH_PARENT_ALL}
wget http://bioinformatics.psb.ugent.be/downloads/psb/i-adhore/i-adhore-3.0.01.tar.gz
tar -zxvf i-adhore-3.0.01.tar.gz
rm i-adhore-3.0.01.tar.gz
# finish the installation following the INSTALL file provided by the ADHoRe team
# less i-adhore-3.0.01/INSTALL
cd i-adhore-3.0.01
mkdir build
cd build
# core dependencies: cmake and g++ 
sudo apt-get install cmake 
sudo apt-get install g++
# marginal dependencies: libpng and zlib
sudo apt-get install libpng-dev 
sudo apt-get install zlib1g-dev
# We need to enable the deposit 'Universe' to install mpi
sudo apt-get install mpi
cmake ..
make
# You do not need to install it, skip the make install
# To plug i-adhore to LibsDyogen
# sed -i "/PATH_ADHORE_BIN =/c\PATH_ADHORE_BIN = \"${PATH_PARENT_ALL}/i-adhore-3.0.01/build/src/i-adhore\"" ${PATH_LIBSDYOGEN}/utils/myADHoRe.py

# Install Cyntenator
echo 'Install Cyntenator' >&2
cd ${PATH_PARENT_ALL}
# download the cyntenator files (pointed by https://www.bioinformatics.org/cyntenator/wiki/Main/HomePage)
wget -r -np -nH --cut-dirs=3 -R index.html https://bbc.mdc-berlin.de/svn/bioinformatics/Software/cyntenator/
cd cyntenator
# core dependency g++
sudo apt-get install g++
# compile
g++ -Wno-deprecated cyntenator.cpp localign.cpp genome.cpp flow.cpp species_tree.cpp -o cyntenator
# To plug cyntenator to LibsDyogen
# sed -i "/PATH_CYNTENATOR_BIN =/c\PATH_CYNTENATOR_BIN = \"${PATH_PARENT_ALL}/cyntenator/cyntenator\"" ${PATH_LIBSDYOGEN}/utils/myCyntenator.py

# reinit working dir
cd ${iwd}

# Next 2 uncommented lines have been moved from the top to here since LibsDyogen/enum.py conflicts with python3. (at least while building the docker container)
#    Traceback (most recent call last):
#      File "/usr/bin/add-apt-repository", line 11, in <module>
#        from softwareproperties.SoftwareProperties import SoftwareProperties, shortcut_handler
#      File "/usr/lib/python3/dist-packages/softwareproperties/SoftwareProperties.py", line 49, in <module>
#        from xml.sax.saxutils import escape
#      File "/usr/lib/python3.5/xml/sax/saxutils.py", line 6, in <module>
#        import os, urllib.parse, urllib.request
#      File "/usr/lib/python3.5/urllib/request.py", line 88, in <module>
#        import http.client
#      File "/usr/lib/python3.5/http/__init__.py", line 1, in <module>
#        from enum import IntEnum
#      File "/home/Dev/LibsDyogen/enum.py", line 66
#        raise NotImplementedError, \
#                                 ^
#    SyntaxError: invalid syntax
# Add the LibsDyogen root folder to the PYTHONPATH environment variable
echo "export PYTHONPATH=\"\${PYTHONPATH}:${PATH_LIBSDYOGEN}\"" >> ~/.bashrc
export PYTHONPATH=${PYTHONPATH}:${PATH_LIBSDYOGEN}

# Post-processing information
green='\e[0;32m'
NC='\e[0m' # No Color
echo -e "${green}Installation finished${NC}"
# Need to run bash to update PYTHONPATH, it is not possible to export a variable to the environment from a bash script
# https://stackoverflow.com/questions/16618071/can-i-export-a-variable-to-the-environment-from-a-bash-script-without-sourcing-i
echo -e "${green}please run 'export PYTHONPATH=\${PYTHONPATH}:${PATH_LIBSDYOGEN}' to update the PYTHONPATH environment variable${NC}"
echo -e "${green}next time you open a bash session, you won't need this since your ~/.bashrc will be loaded${NC}"
echo -e "${green}make sure that 'echo \${PYTHONPATH}' shows the path to the folder LibsDyogen${NC}"

