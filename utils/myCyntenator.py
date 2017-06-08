#!/usr/bin/python
# -*- coding: utf-8 -*-

# LibsDyogen version 1.0 (6/11/2015)
# python v2.7 at least is needed
# Copyright © 2015 IBENS/Dyogen : Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : jlucas@ens.fr
# Licences GLP v3 and CeCILL v2

# Wrapper for 'Cyntenator'
import collections
import os
import getpass # https://stackoverflow.com/questions/842059/is-there-a-portable-way-to-get-the-current-username-in-python
import subprocess
import sys
import myFile, myLightGenomes, myDiags

# change this line for each local installation
PATH_CYNTENATOR_BIN = '/home/' + getpass.getuser() + '/Libs/cyntenator/cyntenator"

# mismatch de 0 si on veut des segments conservés !
# Cyntenator génère parfois des alignements du genre
#
# Align 1
# A A
# ---
# B E
# C D
# D C
# E B
# ---
# F F
#
# Les paires de gènes ancestraux B C D E sont 4 mismatches.
#
# pour obtenir les segments conservés, il faut trouver le juste milieu entre un
# mismatch=-100000 (ce qui générera aucun mismatch) et un mismatch pas trop
# élevé pour permettre les insertions de duplications disperses (qui génère un
# unique mismatch).

# 1) Soit on empeche les mismatch -> on détecte les micro-inversions mais chaque
# insertion de gene dupliqué en tandem brise artificiellement les segments
# conservés
# 2) Soit on les autorise les mismatch -> on autorise des insertions de gènes
# dupliqués en tandem mais on ne détecte pas les micro-inversions

# [("out:Chromosomes",str,), ("withScaffolds",bool,False), ("minChromLength",int,1), ('removeSpeciesSpecificGenes',bool,True)]
# genome1 = myLightGenomes.LightGenome(arguments["genome1"])
# genome2 = myLightGenomes.LightGenome(arguments["genome2"])
# families = utils.myLightGenomes.Families(arguments["ancGenes"])

def ourGenomeToCyntenatorGenomeAndFamily(genome, families,
                                         outCyntenatorGenome="genome.%s.list.cyntenator",
                                         removeSpeciesSpecificGenes=True):
    assert isinstance(genome, myLightGenomes.LightGenome)
    genomeRewritten = myLightGenomes.LightGenome()
    for chr, chrom in genome.iteritems():
        genomeRewritten[chr] = []
        for g in chrom:
            #print >> sys.stderr, "process gene %s of chromosome %s" % (g.names[0],chrom)
            fam = families.getFamilyByName(g.n, default=None)
            if fam is not None: # Genes that are inherited
                genomeRewritten[chr].append((g.n, g.s, fam.fn))
            else: # Genes species specific
                if removeSpeciesSpecificGenes:
                    continue
                else:
                    genomeRewritten[chr].append((g.n, g.s, None))
    #print >> sys.stderr, '!!! print %s' % outCyntenatorGenome
    with open(outCyntenatorGenome, 'w') as f_genome:
        print >> f_genome, "#genome"
        for chr in genomeRewritten:
            for (idxG, (gn, gs, fn)) in enumerate(genomeRewritten[chr]):
                #print >> sys.stderr, "%s %s %s %s %s" % (fn, chr, str(idxG), str(idxG+1), ("+" if gs == 1 else "-"))
                print >> f_genome, "%s %s %s %s %s" % (fn, chr, str(idxG), str(idxG+1), ("+" if gs == 1 else "-"))

def ourGenomesAndFamiliesToCyntenatorGenomes(genome1, genome2, families,
                                             outCyntenatorGenomes="res/genome.%s.list.cyntenator"):
    assert isinstance(genome1, myLightGenomes.LightGenome)
    outCyntenatorGenome1 = outCyntenatorGenomes % 'G1'
    ourGenomeToCyntenatorGenomeAndFamily(genome1, families,
                                         outCyntenatorGenome=outCyntenatorGenome1)
    outCyntenatorGenome2 = outCyntenatorGenomes % 'G2'
    assert isinstance(genome2, myLightGenomes.LightGenome)
    ourGenomeToCyntenatorGenomeAndFamily(genome2, families,
                                         outCyntenatorGenome=outCyntenatorGenome2)
    return (outCyntenatorGenome1, outCyntenatorGenome2)

def lightGenomeFromPairwiseAlignementOfCyntenator(outCyntenatorAlignment, families, includeMicroRearrangedGenesInSbs=True):
    ancSbGenome = myLightGenomes.LightGenome()
    for i,l in enumerate(open(outCyntenatorAlignment, 'r').readlines()):
        if i == 0:
            assert '#alignment' in l
        elif 'Alignment' in l:
            idSb = l.split(' ')[1]
        else:
            line = l.split()
            if len(line) == 0:
                continue
            assert len(line) == 5
            (_, gn1, s1, gn2, s2) = line
            gn1 = gn1.strip()
            gn2 = gn2.strip()
            s1 = s1[1]
            s2 = s2[1]
            assert s1 == '+' or s1 == '-' or s1 == '.'
            assert s2 == '+' or s2 == '-' or s2 == '.'
            # a '.' means a mismatch in the gene alignment
            if s1 == '.':
                assert gn1 == '-'
            if s2 == '.':
                assert gn2 == '-'

            if s1 == '.' or s2 == '.':
                continue
            else:
                s1 = +1 if s1 == '+' else -1
                s2 = +1 if s2 == '+' else -1
                fam1n = None
                fam2n = None
                fam1 = families.getFamilyByName(gn1, default=None)
                # in our case
                assert fam1.fn == gn1, '%s = %s' % (fam1.fn, gn1)
                if fam1 is not None:
                    fam1n = fam1.fn
                    assert fam1n == gn1
                # in our case
                fam2 = families.getFamilyByName(gn2, default=None)
                assert fam2.fn == gn2, '%s = %s' % (fam2.fn, gn2)
                if fam2 is not None:
                    fam2n = fam2.fn
                    assert fam2n == gn2
                assert fam1 is not None and fam2 is not None
                # no better way than to take genes of genome1 along this mismatch
                if fam1n != fam2n:
                    print >> sys.stderr, 'mismatch in idSb = %s : paired g1fn %s != paired g2fn %s' % (idSb, fam1n, fam2n)
                    # ex :
                    # In an alignment of Cyntenator
                    # +A +D -D -C +E +F +G +H ...
                    # +A +B +C +D +E +F +G +H ...
                    # here the segment (+C +D) was micro-rearranged
                    if includeMicroRearrangedGenesInSbs:
                        # will return  +A +B -D -C +E +F +G +H ...
                        pass
                    else:
                        # will return  +A +B +E +F +G +H ...
                        continue

                ancSbGenome[idSb].append((fam1n, s1, gn1, gn2))

    newAncSbGenome=myLightGenomes.LightGenome()
    cptSB = 0
    for (SBId,ancSb) in ancSbGenome.iteritems():
        for (ancGene,s,tb1Name,tb2Name) in ancSb:
            newAncSbGenome[cptSB].append(myLightGenomes.OGene(ancGene, s))
        cptSB += 1
    (ancSbGenome,newAncSbGenome) = (newAncSbGenome,ancSbGenome)
    del newAncSbGenome
    return ancSbGenome

def launchCyntenator(genome1, genome2, families,
                     # options of cyntenator
                     # threshold is the minimum alignment score
                     threshold=4,
                     # gap is the cost of a gap
                     gap=-2,
                     # mismatch is the cost of a mismatch
                     mismatch=-3,
                     # max nb of same gene in the same alignement
                     coverage=1,
                     # nb of overlapping alignments autorized (unique assignment
                                                               # mean no overlap)
                     filter=0,
                     # options more general
                     tandemGapMax=5,  # minimalLengthForSbs=3,
                     filterType=myDiags.FilterType.InBothGenomes,
                     includeMicroRearrangedGenesInSbs=True,
                     outCyntenatorGenomes="../res/genome.%s.list.cyntenator",
                     pathCyntenatorBin=PATH_CYNTENATOR_BIN,
                     outAlignmentCyntenator='../res/alignment.cyntenator'):

    ((g1_tb, tb2g1, (nCL1, nGL1)), (g2_tb, tb2g2, (nCL2, nGL2))) = myDiags.editGenomes(genome1, genome2, families,
                                                                                       filterType=filterType,
                                                                                       labelWith='FamName',
                                                                                       tandemGapMax=tandemGapMax,
                                                                                       minChromLength=1,
                                                                                       keepOriginal=True)
    # remove genes that are not in the genome in gX_tb (genes not remapped in a anc. gene representative of a tb)
    genes1ToRemove = set([g.n for (chr, chrom) in genome1.iteritems() for (gIdx, g) in enumerate(chrom) if ((chr not in tb2g1) or (gIdx not in tb2g1[chr].old))])
    genes2ToRemove = set([g.n for (chr, chrom) in genome2.iteritems() for (gIdx, g) in enumerate(chrom) if ((chr not in tb2g2) or (gIdx not in tb2g2[chr].old))])
    # print >> sys.stderr, 'nb genes1 removed = %s' % len(genes1ToRemove)
    # print >> sys.stderr, 'nb nGL1 = %s' % nGL1
    # print >> sys.stderr, 'nb genes2 removed = %s' % len(genes2ToRemove)
    # print >> sys.stderr, 'nb nGL2 removed = %s' % nGL2
    assert len(genes1ToRemove) == nGL1
    assert len(genes2ToRemove) == nGL2
    genome1.removeGenes(genes1ToRemove)
    genome2.removeGenes(genes2ToRemove)
    genome1.getGeneNames(checkNoDuplicates=True)
    genome2.getGeneNames(checkNoDuplicates=True)
    assert isinstance(genome1, myLightGenomes.LightGenome)

    for genome in [genome1, genome2]:
        #genome.removeUnofficialChromosomes()
        genome.removeChrsStrictlySmallerThan(1)

    # create Cyntenator input Data
    (outCyntenatorGenome1, outCyntenatorGenome2) = \
        ourGenomesAndFamiliesToCyntenatorGenomes(genome1, genome2, families,
                                                 outCyntenatorGenomes=outCyntenatorGenomes)
    # Launch Cyntenator
    # necessary options
    necOptions = "-t \"(%s %s)\" -h id" % (outCyntenatorGenome1, outCyntenatorGenome2)
    # facultative options
    facOptions = "-thr %s -gap %s -mis %s -coverage %s -filter %s -last " % (threshold, gap, mismatch, coverage, filter)
    bashComand = pathCyntenatorBin + ' ' + necOptions + ' ' + facOptions + ' >%s' % outAlignmentCyntenator
    subProcStdin = None  #  subprocess.PIPE, if you want it
    subProcStderr = None  #  subprocess.PIPE, if you want it
    print >> sys.stderr, os.getcwd()
    print >> sys.stderr, bashComand
    proc = subprocess.Popen(bashComand, shell=True, stdin=subProcStdin, stdout=subprocess.PIPE, stderr=subProcStderr)
    # returns stdout and stderr
    adhoreStdOut, subProcStderr = proc.communicate(subProcStdin)
    # proc.returncode, might be interesting too
    # stdout is a str
    # Use this line to output the output of GRIMM
    print >> sys.stderr, subProcStderr

    # convert synteny blocks of ADHoRe into a genome
    alignmentAsLightGenome = lightGenomeFromPairwiseAlignementOfCyntenator(outAlignmentCyntenator, families,
                                                                           includeMicroRearrangedGenesInSbs=includeMicroRearrangedGenesInSbs)

    return alignmentAsLightGenome

if __name__ == '__main__':
    os.chdir('/home/' + getpass.getuser() + '/Libs/PhylDiag/data')
    genome1 = myLightGenomes.LightGenome('genesST.Homo.sapiens.list.bz2')
    genome2 = myLightGenomes.LightGenome('genesST.Mus.musculus.list.bz2')
    families = myLightGenomes.Families('ancGenes.Euarchontoglires.list.bz2')

    #default
    mismatch = -3
    #mismatch = -1000000
    print launchCyntenator(genome1, genome2, families,
                           threshold=4, gap=-2, mismatch=mismatch,
                           tandemGapMax=5,# minimalLengthForSbs=3,
                           filterType=myDiags.FilterType.InBothGenomes,
                           outCyntenatorGenomes="../res/genome.%s.list.cyntenator",
                           pathCyntenatorBin=PATH_CYNTENATOR_BIN,
                           outAlignmentCyntenator='../res/alignment.cyntenator')

    alignmentAsLightGenome = lightGenomeFromPairwiseAlignementOfCyntenator('../res/alignment.cyntenator', families, includeMicroRearrangedGenesInSbs=True)
    print >> sys.stdout, alignmentAsLightGenome
