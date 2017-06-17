#!/bin/bash
#Launch all the commands in the README file and stops on errors if any
set -e
# Any subsequent command which fails will cause the shell script to exit
# immediately
red='\e[0;31m'
green='\e[0;32m'
NC='\e[0m' # No Color

#############################################################
#	Check integrity of pre-processing scripts in LibsDyogen #
#############################################################
mkdir -p test/res

preProcessCommandLines=(
# convet a .nhx tree into a protTree (forest of gene trees)
"scripts/nhxGeneTrees2phylTreeGeneTrees.py test/data/geneTrees.example.nhx > test/res/geneTrees.protTree"
# convet a .nwk tree into a phylTree
"scripts/newickSpeciesTree2phylTreeSpeciesTree.py test/data/speciesTree.nwk > test/res/speciesTree.phylTree"
# extract ancGenes (family)  from the species tree and the forest of gene trees
"scripts/familiesFromGeneTrees.py test/res/speciesTree.phylTree test/res/geneTrees.protTree -out:ancGenes=test/res/%s.families.bz2 > test/res/geneTreesAfterExtractingAncGenes.protTree"
)
for line in "${preProcessCommandLines[@]}"
	do
		echo -e "${green}${line}${NC}"
		eval ${line}
done
