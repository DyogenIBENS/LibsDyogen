# LibsDyogen
[![DOI](https://zenodo.org/badge/45687957.svg)](https://zenodo.org/badge/latestdoi/45687957)

Collaborative python library used in the [DYOGEN team](http://www.ibens.ens.fr/?rubrique43&lang=fr) for studying the evolution of gene order in vertebrates.

## Installation
The easiest way to install LibsDyogen is to execute the remote script INSTALL.sh hosted on github.
This script will clone the github deposit itself.
Dependencies will be installed with the package manager "apt-get" of debian distributions.
LibsDyogen and its plugged softwares will be installed into /home/${USER}/Libs.
The sudo password is required to
* Enable the universe deposit
* Instal dependencies with apt-get
* Add the LibsDyogen folder into the PYTHONPATH, editing the ~/.bashrc

Install *curl*, if you don't have it
```
sudo apt-get update
sudo apt-get install curl
```

Use *curl* to execute the remote file INSTALL.sh hosted on github
```
bash <(curl -s https://raw.githubusercontent.com/DyogenIBENS/LibsDyogen/master/INSTALL.sh)
```

If it does not work follow the next instructions.

### Dependencies
Core dependencies
1. python 2.7
2. cython

Marginal dependencies
1. numpy
2. scipy
3. matplotlib

You also need git for cloning the git deposit.

### Detailed installation guidelines
Enable the deposit 'Universe' on debian/ubuntu distributions
```
sudo add-apt-repository universe
```

Install dependencies
```
sudo apt-get update
sudo apt-get install git python2.7 cython python-matplotlib python-scipy python-numpy
```
Depending on your distribution and version, you may need to change the names of packages

If this did not work, look at indications of authors
* matplotlib : https://matplotlib.org/users/installing.html
* scipy : https://www.scipy.org/install.html

Choose a path for the installation of the library and plugins.
Here we choose /home/<user>/Libs.
```
PATH_PARENT_ALL="/home/${USER}/Libs"
```
Create the main folder for the installation
```
mkdir -p ${PATH_PARENT_ALL}
```

#### Install LibsDyogen
Clone the LibsDyogen library
```
PATH_LIBSDYOGEN="${PATH_PARENT_ALL}/LibsDyogen"
git clone https://github.com/DyogenIBENS/LibsDyogen ${PATH_LIBSDYOGEN}
```
Add the LibsDyogen root folder to the PYTHONPATH environment variable
```
echo "export PYTHONPATH=\"\${PYTHONPATH}:${PATH_LIBSDYOGEN}\"" >> ~/.bashrc
export PYTHONPATH=${PYTHONPATH}:${PATH_LIBSDYOGEN}
```
Then cythonise *.pyx files (=compile C-like files)
```
bash ${PATH_LIBSDYOGEN}/cythonisePyxFiles.sh ${PATH_LIBSDYOGEN}
```

#### Install optional plugins
##### Install Homology Teams
```
cd ${PATH_PARENT_ALL}
wget http://euler.slu.edu/~goldwasser/homologyteams/homologyteams-1.1.zip
unzip homologyteams-1.1.zip
cd homologyteams-1.1/src
```
Compile sources with gcc
```
make
```
homolgyteams should now be plugged automatically to LibsDyogen.
If it is not plugged, update the PATH_HOMOLOGYTEAMS_BIN variable in ${PATH_LIBSDYOGEN}/utils/myGeneTeams.py.
```
PATH_HOMOLOGYTEAMS_BIN = "<PATH_PARENT_ALL>/homologyteams-1.1/src/homologyteams"
```
with <PATH_PARENT_ALL> the appropriate path : in our case it's /home/<user>/Libs, with <user> your user name.

##### Install i-ADHoRe 3.0
```
cd ${PATH_PARENT_ALL}
wget http://bioinformatics.psb.ugent.be/downloads/psb/i-adhore/i-adhore-3.0.01.tar.gz
tar -zxvf i-adhore-3.0.01.tar.gz
rm i-adhore-3.0.01.tar.gz
```
Finish the installation following the INSTALL file provided by the ADHoRe team
```
cd i-adhore-3.0.01
mkdir build
cd build
```
*Cmake* and *g++* are needed for this step.
```
sudo apt-get install cmake
sudo apt-get install g++
```
If you don't already have them.

*libpng* and *zlib* are marginal dependencies of i-ADHoRe 3.0.
```
sudo apt-get install libpng-dev
sudo apt-get install zlib1g-dev
```
Depending on your distribution and version, you may need to change the names of packages libpng-dev and zlib1g-dev.
```
cmake ..
make
```
You do not need to install it, skip ```make install```.

Verify that i-adhore is working
```
cd ../testset/datasetI
../../build/src/i-adhore datasetI.ini
```
i-adhore should be plugged automatically to LibsDyogen.
If it is not plugged properly, update the PATH_ADHORE_BIN variable in ${PATH_LIBSDYOGEN}/utils/myADHoRe.py.
```
PATH_ADHORE_BIN = "<PATH_PARENT_ALL>/i-adhore-3.0.01/build/src/i-adhore"
```
with <PATH_PARENT_ALL> the appropriate path : in our case it's /home/<user>/Libs, with <user> your user name.

##### Install Cyntenator
```
cd ${PATH_PARENT_ALL}
```
Download cyntenator sourcefiles (located in the folder https://www.bioinformatics.org/cyntenator/wiki/Main/HomePage)
```
wget -r -np -nH --cut-dirs=3 -R index.html https://bbc.mdc-berlin.de/svn/bioinformatics/Software/cyntenator/
cd cyntenator
```
Compile
```
g++ -Wno-deprecated cyntenator.cpp localign.cpp genome.cpp flow.cpp species_tree.cpp -o cyntenator
```
Verify that Cyntenator is working.
Read the INSTALL file to find tests to check the installation, for instance try
```
./cyntenator -t "(HSX.txt MMX.txt)" -h phylo HSCFMM.blast  "((HSX.txt:1.2 MMX.txt:1.3):0.5 CFX.txt:2.5):1" > human_mouse
```
cyntenator should be plugged automatically to LibsDyogen.
If it is not plugged, update the PATH_CYNTENATOR_BIN variable in ${PATH_LIBSDYOGEN}/utils/myCyntenator.py
```
PATH_CYNTENATOR_BIN = "<PATH_PARENT_ALL>/cyntenator/cyntenator"
```
with <PATH_PARENT_ALL> the appropriate path : in our case it's /home/<user>/Libs, with <user> your user name

## Usage

See each module for usage.
Core classes are mainly
* utils.myLighGenomes.LightGenome
* utils.myLighGenomes.Families
* utils.myLighGenomes.Genome
* utils.myPhylTree.PhylogeneticTree
* utils.myProteinTree.ProteinTree

### Details about the scripts/ folder

This folder contains several scripts that edit, convert species tree, gene trees or gene families.

Convert nhx (or .nwk, newick) gene trees to our tabular format (phylTree)
```
scripts/nhxGeneTrees2phylTreeGeneTrees.py data/geneTrees.example.nhx > res/geneTrees.protTree
```

Convert a newick species tree into a phylTree species tree
```
scripts/newickSpeciesTree2phylTreeSpeciesTree.py data/speciesTree.nwk > res/speciesTree.phylTree
```

Extract the ancestral gene content (ancGene) from the gene trees
```
scripts/ancGenesFromGeneTrees.py res/speciesTree.phylTree res/geneTrees.protTree -out:ancGenes=res/ancGenes.example.%s.list.bz2 > res/geneTrees.afterExtractingAncGenes.protTree
```

These ancGenes files define gene families. An "ancGene" is an ancestral gene, a gene of the ancestor of interest,
usually the MRCA (Most Recent Common Ancestor). All descendant genes of the same ancestral gene belong to the same
family.

Usually when two species Sa and Sb are compared, gene families are defined by ancGenes.<MRCA(Sa,Sb)>


## Update
If you want to keep LibsDyogen up to date, execute
```
cd ${PATH_LIBSDYOGEN}
git pull
```
often.

This will upgrade your local git deposit to the last commit.

If you want a more stable version, after `git pull`, you can downgrade to the latest tagged version (=stable release), just execute

1. Get tags from the github deposit: `git fetch --tags`
2. Get the latest tag name ``latestTag=$(git describe --tags `git rev-list --tags --max-count=1`)``
3. Checkout the latest tag: `git checkout $latestTag`

Otherwise you can
1. List all tagged versions: `git tag -l`
2. Checkout the version you want: `git checkout <tagName>`

## Contributing
If you want to contribute to this deposit please
1. Fork it
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request

## Credits
* Matthieu Muffato: started the library and created: myGenomes, myTools, myPsOutput, myProteinTre, myPhylTree, myKaryoDrawer, myFile, myGenomes, myGraph, myMultiprocess, myMaths
* Joseph Lucas: Updated the library and organised it, created: myDiags, myLightGenomes, myCondor, myGeneTeams, myGenomesDrawer, myIntervals, myMapping, myProbas, mySvgDrawer and updated: myTools
* Nga thi thuy Nguyen: Updated some functions of myDiags and myGraph
* Lucas Tittmann: Updated some functions of myDiags
* Hugues Roest Crollius: supervisor

## License
This code may be freely distributed and modified under the terms of the GNU General Public License version 3 (GPL v3)
and the CeCILL licence version 2 of the CNRS. These licences are contained in the files:
* LICENSE-GPL.txt (or http://www.gnu.org/licenses/gpl-3.0-standalone.html)
* LICENCE-CeCILL.txt (or http://www.cecill.info/licences/Licence_CeCILL_V2-en.html)
Copyright for this code is held jointly by the Dyogen (DYnamic and Organisation of GENomes) team
of the Institut de Biologie de l'Ecole Normale Supérieure (IBENS) 46 rue d'Ulm Paris, and the individual authors.

## Contacts

* [Joseph Lucas](jlucas@ens.fr)
* [Alexandra Louis](alouis@ens.fr)
* [Hugues Roest Crollius](hrc@ens.fr)