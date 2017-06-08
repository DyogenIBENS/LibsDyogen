 _     _ _         ____
| |   (_) |__  ___|  _ \ _   _  ___   __ _  ___ _ __
| |   | | '_ \/ __| | | | | | |/ _ \ / _` |/ _ \ '_ \
| |___| | |_) \__ \ |_| | |_| | (_) | (_| |  __/ | | |
|_____|_|_.__/|___/____/ \__, |\___/ \__, |\___|_| |_|
                         |___/       |___/

LibsDyogen version 1.01 (02/06/2017)
python v2.7 at least is needed

This code may be freely distributed and modified under the terms of the GNU General Public License version 3 (GPL v3)
and the CeCILL licence version 2 of the CNRS. These licences are contained in the files:
    1) LICENSE-GPL.txt (or http://www.gnu.org/licenses/gpl-3.0-standalone.html)
    2) LICENCE-CeCILL.txt (or http://www.cecill.info/licences/Licence_CeCILL_V2-en.html)

Copyright for this code is held jointly by the Dyogen (DYnamic and Organisation of GENomes) team
of the Institut de Biologie de l'Ecole Normale Supérieure (IBENS) 46 rue d'Ulm Paris, and the individual authors.

Copyright © 2017 IBENS/Dyogen and authors
Authors' contributions:
---------------------------
Matthieu MUFFATO: started the library and created: myGenomes, myTools, myPsOutput, myProteinTre, myPhylTree, myKaryoDrawer, myFile, myGenomes, myGraph, myMultiprocess, myMaths
Joseph LUCAS: Updated the library and organised it, created: myDiags, myLightGenomes, myCondor, myGeneTeams, myGenomesDrawer, myIntervals, myMapping, myProbas, mySvgDrawer and updated: myTools
Nga THI THUY NGUYEN: Updated some functions of myDiags and myGraph
Lucas Tittmann: Updated some functions of myDiags
Hugues ROEST CROLLIUS: supervisor

mail : hrc@ens.fr or jlucas@ens.fr

Installation:
-------------
# The easiest way to install LibsDyogen is to launch the remote script INSTALL.sh hosted on github
# This script will clone the github deposit itself
# Dependencies will be installed with the package manager "apt-get" of debian distributions
# LibsDyogen and its plugged softwares will be installed into /home/${USER}/Libs
# the sudo password is required for:
# 1) the installation of dependencies 
# 2) the addition of the LibsDyogen folder into the PYTHONPATH, editing the ~/.bashrc

# Install curl, if you don't have it
sudo apt-get update
sudo apt-get install curl
# Use curl to execute the remote file INSTALL.sh hosted on github
bash <(curl -s https://raw.githubusercontent.com/DyogenIBENS/LibsDyogen/master/INSTALL.sh)

# If it does not work follow the next indications 

# Install dependencies:
# git
# python 2.7
# cython
# scipy
# numpy
# matplotlib
sudo apt-get update
sudo apt-get install git python2.7 cython python-matplotlib python-scipy python-numpy

# Choose a path for the installation of the library and plugins (here we chose /home/<user>/Libs)
PATH_PARENT_ALL="/home/${USER}/Libs"
# Create the main folder for the installation
mkdir -p ${PATH_PARENT_ALL}

# Install LibsDyogen
# Clone the LibsDyogen library
PATH_LIBSDYOGEN="${PATH_PARENT_ALL}/LibsDyogen"
git clone https://github.com/DyogenIBENS/LibsDyogen ${PATH_LIBSDYOGEN}
# Add the LibsDyogen root folder to the PYTHONPATH environment variable
echo "export PYTHONPATH=\"\${PYTHONPATH}:${PATH_LIBSDYOGEN}\"" >> ~/.bashrc
export PYTHONPATH=${PYTHONPATH}:${PATH_LIBSDYOGEN}
# Then cythonise *.pyx
bash ${PATH_LIBSDYOGEN}/cythonisePyxFiles.sh ${PATH_LIBSDYOGEN}

Install optional plugins:
------------------------------------
# Install homology teams

cd ${PATH_PARENT_ALL}
wget http://euler.slu.edu/~goldwasser/homologyteams/homologyteams-1.1.zip
unzip homologyteams-1.1.zip
cd homologyteams-1.1/src
make
# To plug homologyteams to LibsDyogen, update the PATH_HOMOLOGYTEAMS_BIN variable in ${PATH_LIBSDYOGEN}/utils/myGeneTeams.py
# PATH_HOMOLOGYTEAMS_BIN = "<PATH_PARENT_ALL>/homologyteams-1.1/src/homologyteams"
# with <PATH_PARENT_ALL> the appropriate path : in our case it's /home/<user>/Libs, with <user> your user name

# Install i-ADHoRe 3.0

cd ${PATH_PARENT_ALL}
wget http://bioinformatics.psb.ugent.be/downloads/psb/i-adhore/i-adhore-3.0.01.tar.gz
tar -zxvf i-adhore-3.0.01.tar.gz
rm i-adhore-3.0.01.tar.gz
# finish the installation following the INSTALL file provided by the ADHoRe team
# less i-adhore-3.0.01/INSTALL
cd i-adhore-3.0.01
mkdir build
cd build
# you need cmake for this step (type "sudo apt-get install cmake", if you don't already have it)
sudo apt-get install cmake
cmake ..
make
# You do not need to install it, skip the make install
# Verify that i-adhore is working
cd ../testset/datasetI
../../build/src/i-adhore datasetI.ini
# To plug i-adhore to LibsDyogen, update the PATH_ADHORE_BIN variable in ${PATH_LIBSDYOGEN}/utils/myADHoRe.py
# PATH_ADHORE_BIN = "<PATH_PARENT_ALL>/i-adhore-3.0.01/build/src/i-adhore"
# with <PATH_PARENT_ALL> the appropriate path : in our case it's /home/<user>/Libs, with <user> your user name

# Install Cyntenator

cd ${PATH_PARENT_ALL}
# download the cyntenator files (pointed by https://www.bioinformatics.org/cyntenator/wiki/Main/HomePage)
wget -r -np -nH --cut-dirs=3 -R index.html https://bbc.mdc-berlin.de/svn/bioinformatics/Software/cyntenator/
cd cyntenator
# compile
g++ -Wno-deprecated cyntenator.cpp localign.cpp genome.cpp flow.cpp species_tree.cpp -o cyntenator
# Read the INSTALL file to find tests to check the installation, for instance try:
./cyntenator -t "(HSX.txt MMX.txt)"     -h phylo HSCFMM.blast  "((HSX.txt:1.2 MMX.txt:1.3):0.5 CFX.txt:2.5):1" > human_mouse
# To plug cyntenator to LibsDyogen, update the PATH_CYNTENATOR_BIN variable in ${PATH_LIBSDYOGEN}/utils/myCyntenator.py
# PATH_CYNTENATOR_BIN = "<PATH_PARENT_ALL>/cyntenator/cyntenator"
# with <PATH_PARENT_ALL> the appropriate path : in our case it's /home/<user>/Libs, with <user> your user name

About Licences:
---------------
# Justify the use of snippets from stack overflow
http://programmers.stackexchange.com/questions/12171/how-does-fair-use-apply-to-code-snippets

##########################
# Details about scripts/ #
##########################

# This folder contains several scripts that edit, convert species tree, gene trees or gene families

# Convert nhx (or .nwk, newick) gene trees to our tabular format (phylTree):
scripts/nhxGeneTrees2phylTreeGeneTrees.py data/geneTrees.example.nhx > res/geneTrees.protTree

# Convert a newick species tree into a phylTree species tree:
scripts/newickSpeciesTree2phylTreeSpeciesTree.py data/speciesTree.nwk > res/speciesTree.phylTree

# Extract the ancestral gene content (ancGene) from the gene trees:
scripts/ancGenesFromGeneTrees.py res/speciesTree.phylTree res/geneTrees.protTree -out:ancGenes=res/ancGenes.example.%s.list.bz2 > res/geneTrees.afterExtractingAncGenes.protTree

These ancGenes files define gene families. An "ancGene" is an ancestral gene, a gene of the ancestor of interest,
usually the MRCA (Most Recent Common Ancestor). All descendant genes of the same ancestral gene belong to the same
family.

Usually when two species Sa and Sb are compared, gene families are defined by ancGenes.<MRCA(Sa,Sb)>
