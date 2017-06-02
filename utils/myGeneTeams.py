#!/usr/bin/python
# -*- coding: utf-8 -*-

# LibsDyogen version 1.0 (6/11/2015)
# python v2.7 at least is needed
# Copyright © 2015 IBENS/Dyogen : Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : jlucas@ens.fr
# Licences GLP v3 and CeCILL v2

# Wrapper for 'homologyteams'
# http://euler.slu.edu/~goldwasser/homologyteams/
# Identifying Conserved Gene Clusters in the Presence of Homology Families
# Xin He and Michael H. Goldwasser
# Journal of Computational Biology, 12(6-7):638-656, 2005.)

# * Project: homologyteams
# * Authors: Michael H. Goldwasser (goldwamh@slu.edu) and Xin He (xinhe2@uiuc.edu)
# * Version: 1.0 (May 2004)
# *
# * Copyright (C) 2004  Michael H. Goldwasser and Xin He
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# Copyright of the wrapper
# python 2.7
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

# Uses also homologyteams
import os
import sys
import itertools
import collections

import myMapping
import myFile
import myTools
import myGenomes
import myDiags
import myMaths

PATH_HOMOLOGYTEAMS_BIN = "/home/jlucas/Libs/homologyteams-1.1/src/homologyteams"

FilterType = myDiags.FilterType

#
# Extract gene teams in a pairwise comparison of two chromosomes
#################################################################
@myTools.verbose
def extractGtsInPairCompChr(c1, c2, gc1, gc2, gapMax=0, verbose=False):
    # Convert chromosomes into *.cog file format for homologyteams
    families1 = set([f for (f,_) in gc1])
    families2 = set([f for (f,_) in gc2])
    families = dict([(f, ([], [])) for f in families1 & families2])
    if len(families.keys()) == 0:
        return ()

    for i1, (f,_) in enumerate(gc1):
        if f in families and f is not None:
            families[f][0].append(i1)

    for i2, (f,_) in enumerate(gc2):
        if f in families and f is not None:
            families[f][1].append(i2)

    cogFileName = './tmpFile.cog'
    cogFile = open(cogFileName, 'w')
    print >> cogFile, "FAMILY_ID %s %s" % (c1, c2)
    for (f, (i1s, i2s)) in families.iteritems():
        print >> cogFile, myFile.myTSV.printLine([f, ':'.join([str(i1) for i1 in i1s]), ':'.join([str(i2) for i2 in i2s])], delim=' ')
    cogFile.close()
    # The file *.cog is ok

    # Launch homologyteams
    # Requires that the path to the homologyteams binaries has been added in the
    # PATH environment variable.
    tmpGeneTeamsFileName = 'tmpFile.team'

    cmdLine = PATH_HOMOLOGYTEAMS_BIN
    cmdLine += ' -d ' + str(gapMax)
    cmdLine += ' -O witness ' + cogFileName
    cmdLine += ' -W ' + tmpGeneTeamsFileName
    try :
        os.system(cmdLine)
    except:
        raise ValueError("Error, 'homologyteams' is not plugged properly. Read the README file of LibsDyogen for a proper installation.")

    # parse tmpGeneTeamsFileName
    tmpGeneTeamsFile = open(tmpGeneTeamsFileName, 'r')

    def parseWitness(line):
        line = line.strip(' \n')
        line = line.replace(':', ',')
        line = '[' + line + ']'
        # gene team
        # gt = [(f, i), (f, i), ...], with i the index of the gene on the chr
        gt = eval(line)
        gt = [(f,int(i)) for (f,i) in gt]
        return gt

    geneTeams = []
    while True:
        # Read the first line and do nothing
        line = tmpGeneTeamsFile.readline()
        if not line:
            break
        # Read the 2nd line "Witness in the 1st chromosome:\n" and do nothing
        tmpGeneTeamsFile.readline()
        # Read the 3rd line and parse it
        lineGtC1 = tmpGeneTeamsFile.readline()
        gt1 = parseWitness(lineGtC1)
        # read the 4th line "Witness in the 2st chromosome:\n and do nothing
        tmpGeneTeamsFile.readline()
        # Read the 3rd line and parse it
        lineGtC2 = tmpGeneTeamsFile.readline()
        gt2 = parseWitness(lineGtC2)
        # Read the 4th line and do nothing
        tmpGeneTeamsFile.readline()
        # Read the 2nd line "Witness in the 1st chromosome:\n" and do nothing
        geneTeams.append(((c1, gt1), (c2, gt2)))
    tmpGeneTeamsFile.close()

    # remove temporary files
    os.system("rm %s" % tmpGeneTeamsFileName)
    os.system("rm %s" % cogFileName)
    return geneTeams

@myTools.verbose
def extractGtsInPairCompGenomes(g1, g2, families,
                                tandemGapMax=0,
                                gapMax=None,
                                filterType=FilterType.None,
                                minChromLength=0,
                                verbose=True):
    if isinstance(g1, myGenomes.Genome) and isinstance(g2, myGenomes.Genome):
        g1 = g1.intoDict()
        g2 = g2.intoDict()
    elif isinstance(g1, dict) and isinstance(g2, dict):
        pass
    else:
        raise TypeError('g1 and/or g2 must be either myGenomes.Genome or dict')
    #step 1 :filter genomes and rewrite in tandem blocks if needed
    ##############################################################
    # rewrite genomes by family names (ie ancGene names)
    g1_aID = myMapping.labelWithFamID(g1, families)
    g2_aID = myMapping.labelWithFamID(g2, families)
    # genes that are not in ancGene have a aID=None
    print >> sys.stderr, "genome 1 initially contains %s genes" % sum([len(g1[c1]) for c1 in g1])
    print >> sys.stderr, "genome 2 initially contains %s genes" % sum([len(g2[c2]) for c2 in g2])
    # Must be applied on the two genomes, because of the mode inBothGenomes (InCommonAncestor => not only anchor genes are kept but all genes herited from a gene of the LCA)
    #mfilt2origin1 -> mGf2Go1
    ((g1_aID, mGf2Go1, (nCL1, nGL1)), (g2_aID, mGf2Go2, (nCL2, nGL2))) =\
        myDiags.filter2D(g1_aID, g2_aID, filterType, minChromLength)
    print >> sys.stderr, "genome 1 after filterType=%s and minChromLength=%s contains %s genes" %\
        (filterType, minChromLength, sum([len(g1_aID[c1]) for c1 in g1_aID]))
    print >> sys.stderr, "genome 2 after filterType=%s and minChromLength=%s contains %s genes" %\
        (filterType, minChromLength, sum([len(g2_aID[c2]) for c2 in g2_aID]))
    nGD1 = myMapping.nbDup(g1_aID)[0]
    nGD2 = myMapping.nbDup(g2_aID)[0]
    print >> sys.stderr, "genome 1 contains %s gene duplicates (initial gene excluded)" % nGD1
    print >> sys.stderr, "genome 2 contains %s gene duplicates (initial gene excluded)" % nGD2
    (g1_tb, mtb2g1, nGTD1) = myMapping.remapRewriteInTb(g1_aID, tandemGapMax=tandemGapMax, mOld=mGf2Go1)
    (g2_tb, mtb2g2, nGTD2) = myMapping.remapRewriteInTb(g2_aID, tandemGapMax=tandemGapMax, mOld=mGf2Go2)
    print >> sys.stderr, "genome 1 rewritten in tbs, contains %s tbs" % sum([len(g1_tb[c1]) for c1 in g1_tb])
    print >> sys.stderr, "genome 2 rewritten in tbs, contains %s tbs" % sum([len(g2_tb[c2]) for c2 in g2_tb])
    #TODO, optimise next step

    # second level of verbosity
    verbose2 = True if (len(g1) > 500 or len(g2) > 500) else False
    N12s, N12_g = myDiags.numberOfHomologies(g1_tb, g2_tb, verbose=verbose2)
    print >> sys.stderr, "pairwise comparison of genome 1 and genome 2 yields %s hps" % N12_g
    print >> sys.stderr, "genome 1 contains %s tandem duplicated genes (initial gene excluded)" % nGTD1
    print >> sys.stderr, "genome 2 contains %s tandem duplicated genes (initial gene excluded)" % nGTD2
    nDD1 = myMapping.nbDup(g1_tb)[0]
    nDD2 = myMapping.nbDup(g2_tb)[0]
    print >> sys.stderr, "genome 1 contains %s dispersed duplicated tbs (initial tb excluded)" % nDD1
    print >> sys.stderr, "genome 2 contains %s dispersed duplicated tbs (initial tb excluded)" % nDD2
    assert nDD1 + nGTD1 == nGD1
    assert nDD2 + nGTD2 == nGD2

    # step 2 and 3 : build the MHP and extract putative sbs as diagonals
    #################################################################################
    # extract sbs in the tb base
    listOfGts = []
    nbPairwiseComparisons = len(g1_tb.keys())*len(g2_tb.keys())
    listOfPercentage = range(0, 101, 5)[1:]
    print >> sys.stderr, "gene teams extraction",
    for (i, (chr1, chr2)) in enumerate(itertools.product(g1_tb.keys(), g2_tb.keys())):
        tmpListOfGts = extractGtsInPairCompChr(chr1, chr2,
                                               g1_tb[chr1], g2_tb[chr2],
                                               gapMax=gapMax, verbose=True)
        listOfGts.extend(tmpListOfGts)
        progress = int(float(i*100)/nbPairwiseComparisons)
        if progress in listOfPercentage:
            print >> sys.stderr, "%s" % progress + "%",
            listOfPercentage.remove(progress)
    # new line in the print
    print >> sys.stderr, ""


    # setp 4 : statistical validation of gts
    ##################################################
    # TODO

    # format output
    for (i, gt) in enumerate(listOfGts):
        # translate from the tb base to gene base
        ((c1, gt1), (c2, gt2)) = gt
        la = []
        la = collections.defaultdict(lambda: ([], []))
        for (f1, i1) in gt1:
            la[f1][0].append(i1)
        for (f2, i2) in gt2:
            la[f2][1].append(i2)
        # Sort 'la' owing to the first chromosome
        l = []
        for f in la:
            l.append((f, sorted(la[f][0]), sorted(la[f][1])))
        l.sort(key=lambda x: x[1][0])
        la = [f for (f, _, _) in l]
        # all the paralogs of this family in the gene team
        l1 = [i1s for (_, i1s, _) in l]
        l2 = [i2s for (_, _, i2s) in l]
        # DEBUG assertion
        assert len(l1) == len(l2), "len(l1)=%s, len(l2)=%s" % (len(l1), len(l2))

        mtb2gc1 = mtb2g1[c1]
        mtb2gc2 = mtb2g2[c2]
        assert mtb2gc1.__class__.__name__ == mtb2gc2.__class__.__name__ == 'Mapping'
        # gene team format: ([f1, f2, ...], (c1, [f1i1s, f2i1s, ...]), (c2, [f1i2s, f2i2s, ...])
        la = [f for f in la]
        l1 = [[mtb2gc1[i1] for i1 in i1s] for i1s in l1]
        l2 = [[mtb2gc2[i2] for i2 in i2s] for i2s in l2]
        gt = (la, (c1, l1), (c2, l2))
        listOfGts[i] = gt

    # DEBUG assertion error
    for (la, (c1, l1), (c2, l2)) in listOfGts:
       assert len(la) == len(l1) == len(l2), "len(l1)=%s, len(l2)=%s, len(la)=%s\nl1=%s\nl2=%s\nla=%s" % (len(l1), len(l2), len(la), l1, l2, la)

    return listOfGts


def printGtsFile(listOfGts, genome1, genome2, families):
    print >> sys.stderr, "Print gene teams"

    def foo(genomeX, cX, lX, idxAG):
        if isinstance(genomeX, myGenomes.Genome):
            gXs = [genomeX.lstGenes[cX][gIdx].names[0] for gIdxs in lX[idxAG] for gIdx in gIdxs]
        elif isinstance(genomeX, dict):
            gXs = [genomeX[cX][gIdx][0] for gIdxs in lX[idxAG] for gIdx in gIdxs]
        gXs = ' '.join(gXs)
        gXIdxs = [str(gIdx) for gIdxs in lX[idxAG] for gIdx in gIdxs]
        gXIdxs = ' '.join(gXIdxs)
        return (gXs, gXIdxs)

    statsGts = []
    for (idGt, gt) in enumerate(listOfGts):
        (la, (c1, l1), (c2, l2)) = gt
        assert len(la) == len(l1) == len(l2), "len(l1)=%s, len(l2)=%s, len(la)=%s\nl1=%s\nl2=%s\nla=%s" % (len(l1), len(l2), len(la), l1, l2, la)
        nbAG = len(la)
        statsGts.append(nbAG)
        for (idxAG, fId) in enumerate(la):
            (g1s, g1Idxs) = foo(genome1, c1, l1, idxAG)
            (g2s, g2Idxs) = foo(genome2, c2, l2, idxAG)
            print myFile.myTSV.printLine([idGt, families[fId].fn, c1, c2, g1Idxs, g2Idxs, g1s, g2s])

    print >> sys.stderr, "Distribution of the lengths of gene teams\t", myMaths.myStats.syntheticTxtSummary(statsGts)
