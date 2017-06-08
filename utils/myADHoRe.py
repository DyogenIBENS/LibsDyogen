#!/usr/bin/python
# -*- coding: utf-8 -*-

# LibsDyogen version 1.0 (6/11/2015)
# python v2.7 at least is needed
# Copyright © 2015 IBENS/Dyogen : Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : jlucas@ens.fr
# Licences GLP v3 and CeCILL v2

# Wrapper for 'ADHoRe'
import collections
import glob
import os
import getpass # https://stackoverflow.com/questions/842059/is-there-a-portable-way-to-get-the-current-username-in-python
import subprocess
import sys
import myFile, myLightGenomes, myDiags

# change this line for each local installation
PATH_ADHORE_BIN = '/home/' + getpass.getuser() + '/Libs/i-adhore-3.0.01/build/src/i-adhore"

# [("out:Chromosomes",str,), ("withScaffolds",bool,False), ("minChromLength",int,1), ('removeSpeciesSpecificGenes',bool,True)]
# genome1 = myLightGenomes.LightGenome(arguments["genome1"])
# genome2 = myLightGenomes.LightGenome(arguments["genome2"])
# families = utils.myLightGenomes.Families(arguments["ancGenes"])
from utils import myTools


def ourGenomeToADHoReGenomeAndFamily(genome, families,
                                     outADHoReChromosomes="../i-adhore/Homo.sapiens/Genome.Homo.sapiens.Chr%s.list",
                                     removeSpeciesSpecificGenes=True):
    ADHoReFamily = []
    assert isinstance(genome, myLightGenomes.LightGenome)
    genome.removeUnofficialChromosomes()
    genome.removeChrsStrictlySmallerThan(1)
    for chr, chrom in genome.iteritems():
        chromosomeRewritten = []
        for g in chrom:
            #print >> sys.stderr, "process gene %s of chromosome %s" % (g.names[0],chrom)
            fam = families.getFamilyByName(g.n, default=None)
            if fam is not None: # Genes that are inherited
                chromosomeRewritten.append((g.n, g.s, fam.fn))
            else: # Genes species specific
                if removeSpeciesSpecificGenes:
                    continue
                else:
                    chromosomeRewritten.append((g.n, g.s, None))
        with open(outADHoReChromosomes % chr, 'w') as f_chromosomes:
            for (gn, gs, fn) in chromosomeRewritten:
                print >> f_chromosomes, gn + ("+" if gs == 1 else "-")
                ADHoReFamily.append((gn, (fn if fn != None else gn)))
    return ADHoReFamily

def ourGenomesToADHoReGenomesAndFamilies(genome1, genome2, families,
                                         outADHoReChromosomes="res/Genome.%s.Chr%s.list",
                                         outAHoReFamilies="res/families.csv"):
        assert isinstance(genome1, myLightGenomes.LightGenome)
        outADHoReChromosomes1 = outADHoReChromosomes % ('G1', '%s')
        ADHoReFamily1 = ourGenomeToADHoReGenomeAndFamily(genome1, families,
                                                         outADHoReChromosomes=outADHoReChromosomes1)
        outADHoReChromosomes2 = outADHoReChromosomes % ('G2', '%s')
        assert isinstance(genome2, myLightGenomes.LightGenome)
        ADHoReFamily2 = ourGenomeToADHoReGenomeAndFamily(genome2, families,
                                                         outADHoReChromosomes=outADHoReChromosomes2)
        ADHoReFamily = ADHoReFamily1 + ADHoReFamily2
        ADHoReFamily.sort(key=lambda x: x[1])
        testNoDuplicate = set()
        with open(outAHoReFamilies, 'w') as f_families:
            for (gn, fn) in ADHoReFamily:
                assert gn not in testNoDuplicate, gn
                print >> f_families, '\t'.join([gn, fn])
                testNoDuplicate.add(gn)
        return ((outADHoReChromosomes1, outADHoReChromosomes2), outAHoReFamilies)

# Configure ADHoRe
def printConfigurationFileADHoRe(g1, outADHoReChromosomes1, g2, outADHoReChromosomes2,
                                 outAHoReConfigurationFile, resADHoRePath, outAHoReFamilies,
                                 tandemGapMax=5, gapMax=5, pThreshold=0.001, vizualize='false', minimalLengthForSbs=3):

    with open(outAHoReConfigurationFile, 'w') as f:
        for (Gi, (genome, outADHoReChromosomes)) in enumerate([(g1, outADHoReChromosomes1), (g2, outADHoReChromosomes2)]):
            Gi = 'G%s' % (Gi+1)
            print >> f, 'genome= %s' % Gi
            for (chr, _) in genome.iteritems():
                print >> f, '%s_%s ' % (Gi, chr) + outADHoReChromosomes % (chr)
            print >> f, '\n'

        print >> f, 'blast_table=%s\n' % outAHoReFamilies

        print >> f, 'output_path=%s/' % resADHoRePath
        print >> f, 'table_type=family'
        print >> f, 'tandem_gap=%s' % tandemGapMax
        print >> f, 'gap_size=%s' % gapMax
        print >> f, 'cluster_gap=%s' % gapMax
        print >> f, 'q_value=0.90'
        print >> f, 'prob_cutoff=%s' % pThreshold
        print >> f, 'anchor_points= %s' % minimalLengthForSbs
        print >> f, 'level_2_only=true'
        print >> f, 'cluster_type=colinear'
        print >> f, 'visualizeGHM=%s' % vizualize
        print >> f, 'visualizeAlignment=false'
        print >> f, 'write_stats=true'
        print >> f, 'number_of_threads=4'
    return outAHoReConfigurationFile

def dictGeneName2ancGeneName(outAHoReFamilies):
    geneName2ancGeneName = {}
    for i, (geneName, ancGeneName) in enumerate(myFile.myTSV.readTabular(outAHoReFamilies, [str, str])):
        geneName2ancGeneName[geneName] = ancGeneName
    return geneName2ancGeneName

def genomeFromAdhoreGenomes(outADHoReGenes):

    def update(genome1, gtb1, tbName2Pos1):
        genome1[chr].append(myLightGenomes.OGene(geneName, s))
        assert gIdx == len(genome1[chr]) - 1
        if isTandemRepresentative:
            gtb1[chr].append([myLightGenomes.OGene(genomeN, s)])
            tbName2Pos1[genomeN] = (chr, len(gtb1[chr]) - 1)
        else:
            assert isRemaped
            assert tbName2Pos1[tandemReprName] == tbIdx
            gtb1[chr][tbIdx].append(myLightGenomes.OGene(genomeN, s))

    genome1 = myLightGenomes.LightGenome()
    gtb1 = myLightGenomes.LightGenome()
    tbName2Pos1 = {}

    genome2 = myLightGenomes.LightGenome()
    gtb2 = myLightGenomes.LightGenome()
    tbName2Pos2 = {}

    foo= set()
    # id\tgenome\tlist\tcoordinate\torientation\tremapped_coordinate\tis_tandem\tis_tandem_representative\ttandem_representativeremapped
    for i, (geneName, genomeN, chr, gIdx, s, tbIdx, isTandem, isTandemRepresentative, tandemReprName, isRemaped) in enumerate(
            myFile.myTSV.readTabular(outADHoReGenes, [str, str, str, str, str, str, str, str, str, str])):
        if i == 0:
            continue  #  ignore the first line (header)
        assert s == '+' or s == '-'
        isTandemRepresentative = True if isTandemRepresentative == '-1' else False
        isRemaped = True if isRemaped == '-1' else False
        tbIdx = int(tbIdx)
        if geneName in foo:
            print >> sys.stderr, 'Warning: twice the same ancestral gene in adhore genes.txt'
        foo.add(geneName)

        s = (+1 if s == '+' else -1)
        if genomeN == 'G1':
            update(genome1, gtb1, tbName2Pos1, tbName2Pos1)
        else:
            assert genomeN == 'G2'
            update(genome2, gtb2, tbName2Pos2, tbName2Pos2)

        # ENSG00000186092 Homo.sapiens H_1 0 +	0 0 0 'None' 0
        # geneName2Orientation[geneName] = s
        # remappedGenesAtGeneName[tandemReprName].add(geneName)

        return ((genome1, gtb1), (genome2, gtb2))

def lightGenomeFrombaseClusterOfADHoRe(outAHoReFamilies, outADHoReGenes, outADHoReAnchorpoints):
    ancSbGenome = myLightGenomes.LightGenome()
    # Created usefull dictionnaries
    geneName2ancGeneName = {}
    for i,(geneName, ancGeneName) in enumerate(myFile.myTSV.readTabular(outAHoReFamilies, [str, str])):
        geneName2ancGeneName[geneName] = ancGeneName

    test=set([])
    geneName2Orientation = {}
    for i,(geneName,_,_,_,s,_,_,_,isTandemRepresentative,_) in enumerate(myFile.myTSV.readTabular(outADHoReGenes, [str,str,str,str,str,str,str,str,str,str])):
        if i == 0:
            continue # ignore the first line (header)
        #id\tgenome list\tcoordinate\torientation\tremapped_coordinate\tis_tandem\tis_tandem_representative\ttandem_representativeremapped
        #ENSG00000186092 Homo.sapiens H_1 0 +	0 0 0 'None' 0
        assert s == '+' or s == '-'
        if geneName in test:
            print >> sys.stderr, 'Warning: twice the same ancestral gene in adhore genes.txt'
        test.add(geneName)
        s = (+1 if s=='+' else -1)
        geneName2Orientation[geneName] = s

    #Read the file anchorpoints.txt
    #id\tmultiplicon\tbasecluster\tgene_x\tgene_y\tcoord_x(inTBs)\tcoord_y(inTBs)\is_real_anchorpoint
    old_ancGene=None
    old_SBId=None
    for i,(_,_,SBId,tb1Name,tb2Name,_,_,_) in enumerate(myFile.myTSV.readTabular(outADHoReAnchorpoints, [str, str, str, str, str, str,str,str])):
        if i == 0:
            continue # ignore the first line (header)
        # find the name of the ancestral gene
        ancGene1 = geneName2ancGeneName[tb1Name]
        ancGene2 = geneName2ancGeneName[tb2Name]
        assert ancGene1 == ancGene2
        ancGene = ancGene1
        #Build a genome with sbs as contigs and only conserve one gene of a tb
        if SBId == old_SBId:
            if ancGene == old_ancGene: # Happens when there are paralogous hps in the same diagonal
                #print >> sys.stderr, "paralogous hps in the same diagonal, (it is rare in ADHoRe)"
                continue # This may create diags of size 1 that will need to be removed
        # No ancestral orientation reported, take the orientation of the first genome

        # verify that genes are ordered properly!
        ancSbGenome[SBId].append((ancGene,geneName2Orientation[tb1Name],tb1Name,tb2Name))
        old_ancGene = ancGene
        old_SBId = SBId

    #Last loop for taking care of sbs of length <=1, or length <=2 if the option 'useSbsLength2' is set to false
    newAncSbGenome=myLightGenomes.LightGenome()
    cptSB = 0
    for (SBId,ancSb) in ancSbGenome.items():
        for (ancGene,s,tb1Name,tb2Name) in ancSb:
            newAncSbGenome[cptSB].append(myLightGenomes.OGene(ancGene, s))
        cptSB += 1
    (ancSbGenome,newAncSbGenome) = (newAncSbGenome,ancSbGenome)
    del newAncSbGenome
    return ancSbGenome

def launchADHoRe(genome1, genome2, families, gapMax=5, tandemGapMax=5, pThreshold=0.001, minimalLengthForSbs=3,
                 filterType=myDiags.FilterType.InBothGenomes,
                 outADHoReChromosomes="../res/%s/Genome.%s.Chr%s.list",
                 outAHoReFamilies="../res/families.csv",
                 outAHoReConfigurationFile="../res/dataset_G1_G2.ini",
                 resADHoRePath="../res/resADHoRe",
                 pathIADHOREbin=PATH_ADHORE_BIN):
    filesToRemove = glob.glob(os.path.dirname(outADHoReChromosomes) + '/Genome.*.list')
    filesToRemove += [outAHoReFamilies]
    filesToRemove += [outAHoReConfigurationFile]
    filesToRemove += [resADHoRePath + '/alignment.txt']
    filesToRemove += [resADHoRePath + '/anchorpoints.txt']
    filesToRemove += [resADHoRePath + '/baseclusters.txt']
    filesToRemove += [resADHoRePath + '/collinear_portions.txt']
    filesToRemove += [resADHoRePath + '/duplicated_portions.txt']
    filesToRemove += [resADHoRePath + '/genes.txt']
    filesToRemove += [resADHoRePath + '/list_elements.txt']
    filesToRemove += [resADHoRePath + '/multiplicon_pairs.txt']
    filesToRemove += [resADHoRePath + '/multiplicons.txt']
    filesToRemove += [resADHoRePath + '/segments.txt']
    for filename in filesToRemove:
        if os.path.isfile(filename):
            # os.remove will only remove one file at a time
            os.remove(filename)

    assert gapMax >= 1
    assert tandemGapMax >= 2
    assert minimalLengthForSbs >= 3

    ((g1_tb, tb2g1, (nCL1, nGL1)), (g2_tb, tb2g2, (nCL2, nGL2))) = myDiags.editGenomes(genome1, genome2, families,
                                                                                       filterType=filterType,
                                                                                       labelWith='FamName',
                                                                                       tandemGapMax=tandemGapMax,
                                                                                       minChromLength=1,
                                                                                       keepOriginal=True)
    genes1Removed = set([g.n for (chr, chrom) in genome1.iteritems() for (gIdx, g) in enumerate(chrom) if ((chr not in tb2g1) or (gIdx not in tb2g1[chr].old))])
    genes2Removed = set([g.n for (chr, chrom) in genome2.iteritems() for (gIdx, g) in enumerate(chrom) if ((chr not in tb2g2) or (gIdx not in tb2g2[chr].old))])
    # print >> sys.stderr, 'nb genes1 removed = %s' % len(genes1Removed)
    # print >> sys.stderr, 'nb nGL1 = %s' % nGL1
    # print >> sys.stderr, 'nb genes2 removed = %s' % len(genes2Removed)
    # print >> sys.stderr, 'nb nGL2 removed = %s' % nGL2
    assert len(genes1Removed) == nGL1
    assert len(genes2Removed) == nGL2
    genome1.removeGenes(genes1Removed)
    genome2.removeGenes(genes2Removed)
    genome1.getGeneNames(checkNoDuplicates=True)
    genome2.getGeneNames(checkNoDuplicates=True)
    assert isinstance(genome1, myLightGenomes.LightGenome)
    # create ADHoRE input Data
    ((outADHoReChromosomes1, outADHoReChromosomes2), outAHoReFamilies) = \
        ourGenomesToADHoReGenomesAndFamilies(genome1, genome2, families,
                                             outADHoReChromosomes=outADHoReChromosomes,
                                             outAHoReFamilies=outAHoReFamilies)
    # create configuration file
    configurationFilePath = printConfigurationFileADHoRe(genome1, outADHoReChromosomes1, genome2, outADHoReChromosomes2,
                                                         outAHoReConfigurationFile, resADHoRePath, outAHoReFamilies ,
                                                         tandemGapMax=tandemGapMax, gapMax=gapMax, pThreshold=pThreshold, vizualize='false',
                                                         minimalLengthForSbs=minimalLengthForSbs)
    # Launch ADHore
    options = ' '
    bashComand = pathIADHOREbin + options + configurationFilePath
    subProcStdin = None  #  subprocess.PIPE, if you want it
    subProcStderr = None  #  subprocess.PIPE, if you want it
    proc = subprocess.Popen(bashComand, shell=True, stdin=subProcStdin, stdout=subprocess.PIPE, stderr=subProcStderr)
    # returns stdout and stderr
    adhoreStdOut, subProcStderr = proc.communicate(subProcStdin)
    # proc.returncode, might be interesting too
    # stdout is a str
    # Use this line to output the output of GRIMM
    print >> sys.stderr, subProcStderr

    outADHoReGenes = resADHoRePath + '/' + 'genes.txt'
    outADHoReAnchorpoints = resADHoRePath + '/' + 'anchorpoints.txt'

    # convert synteny blocks of ADHoRe into a genome
    baseClustersAsLightGenome = lightGenomeFrombaseClusterOfADHoRe(outAHoReFamilies, outADHoReGenes, outADHoReAnchorpoints)

    return baseClustersAsLightGenome

if __name__ == '__main__':
    os.chdir('/home/' + getpass.getuser() + '/Libs/PhylDiag/data')
    genome1 = myLightGenomes.LightGenome('genesST.Homo.sapiens.list.bz2')
    genome2 = myLightGenomes.LightGenome('genesST.Mus.musculus.list.bz2')
    families = myLightGenomes.Families('ancGenes.Euarchontoglires.list.bz2')
    print launchADHoRe(genome1, genome2, families, gapMax=5, tandemGapMax=5, pThreshold=0.001, minimalLengthForSbs=3,
                       outADHoReChromosomes="../res/Genome.%s.Chr%s.list",
                       outAHoReFamilies="../res/families.csv",
                       outAHoReConfigurationFile="../res/dataset_G1_G2.ini",
                       resADHoRePath="../res/resADHoRe",
                       pathIADHOREbin=PATH_ADHORE_BIN)
