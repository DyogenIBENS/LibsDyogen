#!/usr/bin/python
# -*- coding: utf-8 -*-
# PhylDiag v1.02
# needs python 2.7 at least
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

###################################################################
# PhylDiag core algorithm and tools to manage synteny blocks #
###################################################################

# This code is based in major parts on the publication ``PhylDiag : identifying complex cases of conserved synteny that include tandem duplications''
# authors :
# Joseph LUCAS (IBENS, Paris, France, jlucas@ens.fr)
# Matthieu MUFFATO (EBI, Cambridge, United Kingdom, muffato@ebi.ac.uk)
# Lucas TITTMANN (IBENS, Paris, France)
# Hugues ROEST CROLLIUS (IBENS, Paris, France, hrc@ens.fr)

# This code uses personal libraries : myTools and myGenomes

###################################
# Vocabulary (see the publication)
###################################
# g : gene
# sb : synteny block
# tb : tandem block
# s : may refer to 'strand', the transcriptional orientation of a gene or a tb
# hp : homology pack
# sign : may refer to a sign of a homology pack
# MH : matrix of homologies
# MHP: matrix of homology packs
# diagonal : synonym of a sb (a sb is roughly  a daigonal in the MHP)
# diagonal type : either 'slash' ('/', bottom-left to top-right) or 'backslash' ('\', top-left to bottom-right)
# consistent diagonal : either a slash diagonal with hps signs = (+1 or None) or a backslash diagonal with hps signs = (-1 or None)

import os
import sys
import math
import copy
import itertools
import collections
# FIXME not the good way to import

import myTools
import myGenomes
import myProbas
import myFile
import myMaths
import myMapping
import myMultiprocess

import enum

FilterType = enum.Enum('None', 'InCommonAncestor', 'InBothSpecies')


class Diagonal():
    def __init__(self, *args, **kargs):
        if len(args) == 4 and (isinstance(args[0], str) or args[0] is None)\
                and all([isinstance(l, list) for l in args[1:]]):
            # args = [diagType, l1, l2, la]
            self.dt = str(args[0]) if args[0] is not None else None
            # TODO inherit from Diagonal from list and change:
            # self.la into self[0]
            # self.l1 into self[1]
            # self.l2 into self[2]
            # also implement some function to perform easy merges
            # __add__ function that performs the merge operations on all the
            # lists.
            self.l1 = list(args[1])
            self.l2 = list(args[2])
            self.la = list(args[3])
        elif len(args) == 1 and isinstance(args[0], Diagonal):
            # args = [diag]
            diag = args[0]
            # The diagonal type: either '/', '\' or None
            self.dt = str(diag.dt) if diag.dt is not None else None
            # The list of the indexes of genes on the first chromosome
            self.l1 = list(diag.l1)
            # The list of the indexes of genes on the second chromosome
            self.l2 = list(diag.l2)
            # The list of ancestral genes (familyName, ancetsralStrand)
            # ancestralStrand is the strand of the ancestral chromosome if the order
            # of the ancestral synteny block was the same as the order in the first
            # chromosome.
            # la : [..., (aGn, aGs, dist), ...]
            # (ancestral gene name, ancestral gene strand, distance)
            # distance is a measure of the the 2D-spacing between hps in the matrix
            # of homologies that gives an information on the probable spacing between
            # ancestral genes in the ancestor.
            self.la = list(diag.la)
        else:
            raise ValueError("wrong arguments for the constructor")

    def beg(self):
        return (self.l1[0], self.l2[0])

    def end(self):
        return (self.l1[-1], self.l2[-1])

    def minOnG(self, rankGenome):
        assert rankGenome in [1,2]
        return min(self.l1) if rankGenome == 1 else min(self.l2)

    def maxOnG(self, rankGenome):
        assert rankGenome in [1,2]
        return max(self.l1) if rankGenome == 1 else max(self.l2)

    #returns the maximum gap between two hps with the Chebyschev Distance metric (CD)
    def max_g(self):
        #l1s = [(self.l1[i1], self.l1[i1+1]) for i1 in enumerate(self.l1[:-1])]
        l1s = myTools.myIterator.slidingTuple(self.l1)
        max_g1 = max(l1s, key=lambda x: abs(x[1]-x[0]))
        # do not forget the -1
        max_g1 = abs(max_g1[0] - max_g1[1]) - 1
        l2s = myTools.myIterator.slidingTuple(self.l2)
        max_g2 = max(l2s, key=lambda x: abs(x[1]-x[0]))
        # do not forget the -1
        max_g2 = abs(max_g2[0] - max_g2[1]) - 1
        max_g = max(max_g1, max_g2)
        return max_g

    def calculateCharacteristics(self):
        # number of homologies in the sb
        m = len(self.la)
        #(x0,y0) = (l1_tb[0],l2_tb[0])
        #(x1,y1) = (l1_tb[1],l2_tb[1])
        #diagType = '/' if (x1 - x0)*(y1 - y0) else '\\'
        lw1 = [sys.maxint, -sys.maxint]
        lw2 = [sys.maxint, -sys.maxint]
        for idx1 in self.l1:
            lw1[0] = idx1 if idx1 < lw1[0] else lw1[0]
            lw1[1] = idx1 if idx1 > lw1[1] else lw1[1]
        for idx2 in self.l2:
            lw2[0] = idx2 if idx2 < lw2[0] else lw2[0]
            lw2[1] = idx2 if idx2 > lw2[1] else lw2[1]
        # +1 in order to have the nb of genes in the W1
        lw1 = lw1[1] - lw1[0]+1
        # +1 --------------------------------------- W2
        lw2 = lw2[1] - lw2[0]+1
        assert lw1 >= 0, lw1
        assert lw2 >= 0, lw2
        l1_min = min(self.l1)
        l1_max = max(self.l1)
        l2_min = min(self.l2)
        l2_max = max(self.l2)
        if m > 1:
            # Because of the scanning process
            assert self.l1[1] >= self.l1[0]
            max_g = self.max_g()
        else:
            max_g = None
        return (m, max_g, lw1, lw2, l1_min, l1_max, l2_min, l2_max)

    def __repr__(self):
        return "\ndiagType=%s\nl1=%s\nl2=%s\nla=%s\n" % (self.dt, self.l1, self.l2, self.la)


class SyntenyBlock(Diagonal):
    def __init__(self, *args, **kwargs):
        #args -- tuple of anonymous arguments
        #kwargs -- dictionary of named arguments
        if len(args) == 2 and isinstance(args[0], Diagonal) and (isinstance(args[1], float) or isinstance(args[1], int) or args[1] is None):
            # args = [diag, pVal]
            (diag, pVal) = args
            Diagonal.__init__(self, diag)
            self.pVal = pVal
        elif len(args) == 1 and isinstance(args[0], SyntenyBlock):
            # args = [sb]
            sb = args[0]
            Diagonal.__init__(self, sb.dt, sb.l1, sb.l2, sb.la)
            self.pVal = sb.pVal
        else:
            raise ValueError("wrong arguments for the constructor")

    def __repr__(self):
        return Diagonal.__repr__(self) + "pVal=%s\n" % self.pVal

# Diagonal Pseudo Distance
def DPD((x0,y0), (x1,y1)):
    return 2 * max(abs(x1-x0), abs(y1-y0)) - min(abs(x1-x0), abs(y1-y0))

# Chebyshev Distance
def CD((x0, y0), (x1, y1)):
    return max(abs(x1-x0), abs(y1-y0))

# Manhattan Distance
def MD((x0,y0),(x1,y1)):
    return abs(x1-x0) + abs(y1-y0)

# Euclidean Distance
def ED((x0,y0),(x1,y1)):
    return round(math.sqrt((x1-x0)**2 + (y1-y0)**2))

# The frame of all the distances used during the merging process
def framed(f):
    def distance((x0,y0), (x1,y1), diagType):
        if diagType == '/':
            # Only consider the top right part of the distance matrix, (x1>=x0 and y1>=y0)
            if (x1-x0)>=0 and (y1-y0) >= 0: # in i-adhore "...>0", here we allow '=' to take care of close dispered paralogies that can be understood as tandem duplicates
                res = f((x0,y0),(x1,y1))
            else:
                res = sys.maxint
        elif diagType == '\\':
            # Only consider the bottom right part of the distance matrix, (x1>=x0 and y1<=y0)
            if (x1-x0) >= 0 and (y1-y0) <= 0: # in i-adhore "...<0", here we allow '=', idem than before
                res = f((x0,y0),(x1,y1))
            else:
                res = sys.maxint
        else:
            # diagType == None, diag is composed of a single gene
            res = f((x0,y0),(x1,y1))
        return res
    return distance

#
# Generator managing the queue of diagonals for the merging process
####################################################################
class queueWithBackup:
    # gen is a generator
    def __init__(self, gen):
        self.gen = gen
        self.backup = collections.deque()
        self.todofirst = []

    def __iter__(self):
        return self

    # The next returned value comes from either the buffer or from the main generator
    def next(self):
        if len(self.todofirst) > 0:
            return self.todofirst.pop()
        # Returns the next value of the generator
        return self.gen.next()

    # The reinserted element is put on hold, the last reinserted element will be on top of the pile during the next rewind
    def putBack(self, x):
        self.backup.appendleft(x)

    # Recorded elements are put in a prioritary buffer
    def rewind(self):
        self.todofirst.extend(self.backup)
        self.backup = collections.deque()

#
# Merge diagonals if they are separated by a gap less long than 'gapMax' relatively to a distance metric
# inputs:
# listOfDiags : a list containing elements as (l1, l2, la)
#       l1 = [ ..., [i11, i12...], ....] with i1x < i1y for all x<y
#               i1 : gene index of the corresponding tb in the genome 1
#       l2 : [..., [...,ix2...], ...]
#               ix2 : gene index in the genome 2
#       la : [..., (ancGeneIdx,aGs),...]
#               ancGeneIdx is the ancGene number (corresponds to the line index of the gene)
#               aGs is the orientation of the ancetsral gene (usually the same
#               as the corresponding tb in the genome1)
# gapMax : the maximum allowed gap between merged diagonals
# distanceMetric : the metric used for the gap calculation
# outputs:
#       listOfSortedAndMergedDiagonals : same structure than listOfDiags
########################################################################################
# TODO : optimize the search of diags by avoiding considering diags that are on the left of diagA
@myTools.verbose
def mergeSbs(listOfDiags, gapMax, gc2, distanceMetric = 'DPD', verbose = False):
    assert gapMax>=0
    # assert that the diag elements in listOfDiags are either Diagonal objects
    # or SyntenyBlock objects
    assert all(diag.__class__.__name__ == 'Diagonal' for diag in listOfDiags) or\
        all(diag.__class__.__name__ == 'SyntenyBlock' for diag in listOfDiags)

    if distanceMetric == 'DPD':
        print >> sys.stderr, "use Diagonal Pseudo Distance to merge diagonals with a gap up to %s elements" % gapMax
        distance = framed(DPD)
    elif distanceMetric == 'CD':
        print >> sys.stderr, "use Chebyshev Distance to merge diagonals with a gap up to %s elements" % gapMax
        distance = framed(CD)
    elif distanceMetric == 'MD':
        print >> sys.stderr, "use Manhattan Distance to merge diagonals with a gap up to %s elements" % gapMax
        distance = framed(MD)
    elif distanceMetric == 'ED':
        print >> sys.stderr, "use Euclidean Distance to merge diagonals with a gap up to %s elements" % gapMax
        distance = framed(ED)
    else:
        raise ValueError('Must use a distance either DPD (Diagonal PSeudo Distance) or MD (Manhattan Distance), Euclidean Distance (ED) or Chebyshev Distance (CD)')

    print >> sys.stderr, "Number of Diags before DiagMerger = ", len(listOfDiags)
    diagGen = []
    listOfFinishedDiags = []
    nbFusion = 0
    # Add the distance in la : [(ag1,ags1,dist1=0), (ag2,ags2,dist2=1), (ag2,ags3,dist3=3), ...]
    # means that between ag1 and ag2 there is no 2Dgap but between ag2 and ag3
    # there is a gap of 1. The gap depends on the distance metric chosen.
    # This gap gives us an information on the relevance of this adjacency in th
    # synteny block.
    for currGap in range(0, gapMax+1):
        nbFusionCurrGap = 0
        # Sort diagonals by increasing index on the genome A (i.e. by increasing x coordinate of the leftmost hp of the diagonal)
        # This is a key part of the algorithm ! During the whole process we need to keep diagonals sorted by increasing index on the genome A
        if currGap > 0: #After initialisation
            listOfDiags = sorted(list(diagGen) + listOfFinishedDiags, key=lambda diag: diag.beg()[0])
            listOfFinishedDiags=[]
        # assert all([l1_1[0][0] <= l1_2[0][0] for ((_,l1_1,_,_,_,_),(_,l1_2,_,_,_,_)) in myTools.myIterator.slidingTuple(listOfDiags)])
        # assert all([l1[0][0] == min(l1)[0] for (_,l1,_,_,_,_) in listOfDiags])
        # assert all([l1[-1][0] == max(l1)[0] for (_,l1,_,_,_,_) in listOfDiags])
        diagGen = queueWithBackup(d for d in listOfDiags)
        for diagA in diagGen:
            ### DEBUG example
            #if 426 == diagA.beg()[0]+1 and 543 == diagA.beg()[1]+1:
            #       print >> sys.stderr, "Is diag A : [diagType = %s, start =
            #       (%s,%s), end =(%s,%s)] fusionable with" %
            #       (diagA.dt,diagA.beg()[0]+1,diagA.beg()[1]+1,diagA.en()d[0]+1,diagA.end()[1]+1)
            ### DEBUG example
            continueToSearch = True
            fusionableDiags = []
            impossibleToMergeDiags = []
            for diagB in diagGen:
                if not continueToSearch:
                    impossibleToMergeDiags.append(diagB)
                    break
                ### DEBUG example
                #if 427 == diagB.beg()[0]+1 and 541 == diagB.beg()[1]+1:
                #       print >> sys.stderr, "Is diag B : [diagType = %s, start
                #       = (%s,%s), end =(%s,%s)] fusionable with" %
                #       (diagB.dt,diagB.beg()[0]+1,diagB.beg()[1]+1,diagB.end()[0]+1,diagB.end()[1]+1)
                #       print >> sys.stderr, "diagA.end() = (%s,%s)" %
                #       (diagA.end()[0]+1,diagA.end()[1]+1), " ?"
                ### DEBUG example
                # Thanks to the sort at the beginning, we known that the starting hp of diagB is on the right of the starting hp of diagA
                #TODO : change diagGen to put diags that are on the left of diagA in a buffer and avoid considering them
                if diagB.beg()[0] < diagA.end()[0]: # in i-adhore diagB.beg()[0] <= diagA.end()[0]
                    impossibleToMergeDiags.append(diagB)
                    continue
                elif diagB.beg()[0] <= diagA.end()[0] + currGap+1:
                    # diagA.end()[0] < diagB.beg()[0] <= diagA.end()[0] + currGap
                    # Check if diagTypes are compatible
                    if diagA.dt == diagB.dt or diagA.dt == None or diagB.dt == None:
                        # take the known diagType if it is known in at least one of the 2 diags
                        dT = diagA.dt if diagA.dt != None else diagB.dt
                        # Check the distance
                        if distance(diagA.end(),diagB.beg(),dT) == currGap+1 :
                            fusionableDiags.append(diagB)
                        else:
                            impossibleToMergeDiags.append(diagB)
                            continue
                    else:
                        impossibleToMergeDiags.append(diagB)
                        continue
                else:
                    # diagA.end()[0] + currGap < diagB.beg()[0]
                    # stop the loop (remember that the diags are sorted!)
                    # Impossible to merge next diags to diagA
                    continueToSearch = False
                    impossibleToMergeDiags.append(diagB)

            if len(fusionableDiags) > 0 :
                #assert all(x[4][0] >= diagA.end()[0] for x in fusionableDiags)
                #if diagA.dt == '/':len(la)
                #       assert all( (x[0] == '/' or  x[0] == None ) for x in fusionableDiags)
                #       assert all(x[4][1]-diagA.end()[1] >= 0 for x in
                #       fusionableDiags), [(x[4][1],diagA.end()[1]) for x in fusionableDiags]
                #elif diagA.dt == '\\':
                #       assert all( (x[0] == '\\' or x[0] == None ) for x in fusionableDiags)

                #sort fusionableDiags by number of hps, the more hps there are in the diagonal, the more the fusion may yield a long diagonal
                fusionableDiags = sorted(fusionableDiags, key=lambda diag:len(diag.la), reverse=True)
                # if fusionableDiags contains at least 2 diags and if the two
                # first diags of fusionableDiags have the same number of hps, we
                # keep sorting all the diags that have the same number of hps
                if len(fusionableDiags) > 1 and len(fusionableDiags[0].la) == len(fusionableDiags[1].la):
                    maxNbHps = len(fusionableDiags[0].la)
                    tmpBestrFusionableDiags = [diag for diag in fusionableDiags if len(diag.la) == maxNbHps]
                    # score'Non'Diagonality = abs(deltaX - deltaY)
                    scoreNonDiagonality = lambda diag: abs(diag.beg()[0]-diagA.end()[0]) - abs(diag.beg()[1]-diagA.end()[1])
                    # sort from non diagonality score, the best fusionable diag
                    # is the one that is farthest from the diagonal axis
                    tmpBestrFusionableDiags = sorted(tmpBestrFusionableDiags, key=scoreNonDiagonality, reverse=True)
                    fusionableDiags = tmpBestrFusionableDiags + fusionableDiags[len(tmpBestrFusionableDiags):]

                # Sort by priority depending on the choice of the distance :
                # either we want to fuse in priority along the diagonal line or
                # relatively to the proximity on chromosome 1
                #scoreNonDiagonality= lambda x: abs(x[4][0]-diagA[5][0]) - abs(x[4][1]-diagA[5][1])
                #fusionableDiags = sorted(fusionableDiags, key=scoreNonDiagonality, reverse=True)

                #sort by proximity to the diagonal line
                #elif distanceMetric == 'MD': # Doesn't change the result
                #       #sort by proximity on chromosome 1
                #       fusionableDiags = sorted(fusionableDiags, key=lambda x:x[4][0]-diagA[5][0])

                diagToFuse = fusionableDiags.pop(0) # Pop the first element of the list
                # Lists are fused
                if diagToFuse.dt != None and diagA.dt != None:
                    if  diagToFuse.dt == diagA.dt:
                        dt_res = diagToFuse.dt
                    else:
                        # opposed diagTypes
                        raise ValueError
                elif diagToFuse.dt != None and diagA.dt == None:
                    dt_res = diagToFuse.dt
                elif diagA.dt != None and diagToFuse.dt == None:
                    dt_res = diagA.dt
                else:
                    # dt == None and diagA.dt == None:
                    # If the two diagonals are singletons
                    dt_res = '/' if diagA.end()[1] < diagToFuse.beg()[1] else '\\'
                diagA.l1.extend(diagToFuse.l1)
                diagA.l2.extend(diagToFuse.l2)
                # Compute the new dist with currGap
                diagToFuse.la = [(aGN, aGs, dist+diagToFuse.la[-1][2]+currGap) for (aGN, aGs, dist) in diagToFuse.la]
                diagA.la.extend(diagToFuse.la)
                # Previous diagonals without diagTypes have now a diagType and
                # all ancestral genes can now be oriented.
                for (i, (aGN, aGs, dist)) in enumerate(diagA.la):
                    if aGs == None:
                        if dt_res == '/':
                            aGs = gc2[diagA.l2[i]][1]
                        else:
                            assert dt_res == '\\'
                            s_tb2 = gc2[diagA.l2[i]][1]
                            aGs = - s_tb2 if s_tb2 != None else None
                    # update the orientation
                    diagA.la[i] = (aGN, aGs, dist)

                diagA.dt = dt_res
                if diagA.__class__.__name__ == SyntenyBlock:
                    assert diagToFuse.__class__.name__ == SyntenyBlock
                    diagA.pVal = min(diagA.pVal, diagToFuse.pVal)

                nbFusion += 1
                nbFusionCurrGap += 1
                #We still try to merge diags to the new diagA
                diagGen.putBack(diagA) # Needed to update diagA infos in the loop
            else:
                listOfFinishedDiags.append(diagA)
            # Diags are re-sorted owing to the A chromosome index
            for diag in sorted(fusionableDiags + impossibleToMergeDiags, key=lambda diag: diag.l1[0]):
                diagGen.putBack(diag) # diagonals that were not fusionable are recorded to try to merge them after
            diagGen.rewind()
        print >> sys.stderr, "number of merges for currGap=%s :%s" % (currGap, nbFusionCurrGap)
    # Once all merges have been performed for a currGap it is necessary to repeat the merging process with all the diagonals for currGap+1
    print >> sys.stderr, "Total number of merges =", nbFusion
    listOfSortedAndMergedDiagonals = sorted(list(diagGen) + listOfFinishedDiags, key=lambda diag: diag.l1[0])
    print >> sys.stderr, "number of Diags after the merging process =", len(listOfSortedAndMergedDiagonals)

    return listOfSortedAndMergedDiagonals


# TODO merge sbs when they have a small overlap
# issue, the merge may inlclude some other diagonals...
#def mergeOverlappingSbs(sbsInPairComSV, overlapMax=0):
#    return sbsInPairComSV

def buildSetOfVertices(sbsInPairCompWithIds):
    V = set([])
    sbsG1 = collections.defaultdict(list)
    sbsG2 = collections.defaultdict(list)
    for (idsb, (c1, c2), sb) in sbsInPairCompWithIds.iterByOrderedIds():
        # identify all sbs that are projected on chr c1
        sbsG1[c1].append(idsb)
        # identify all sbs that are projected on chr c2
        sbsG2[c2].append(idsb)
        #weight = len(sb.la)
        #V[iDsb] = weight
        V.add(idsb)
    return (V, (sbsG1, sbsG2))


def findOverlapsOnGenome(rankGenome, sbsInPairCompWithIds, sbsGX, overlapMax=0):

    def findOverlapsOnGenomeOrder(rankGenome, order, sbsGXc, N, O, I, id2sb, overlapMax=0):
        assert rankGenome in [1, 2]
        assert order in [+1, -1]
        if order == +1:
            # sort by starting coordinate first and end coordinate second
            sbsGXc.sort(key=lambda iDsb: (id2sb(iDsb).minOnG(rankGenome), id2sb(iDsb).maxOnG(rankGenome)))
            # find overlaps on the right side of sbra
            # sbra: sb range 'a'
            # sbrb: sb range 'b'
            for (ia, iDsba) in enumerate(sbsGXc):
                sba = id2sb(iDsba)
                sba_beg = sba.minOnG(rankGenome) # and not sba.beg()[coord]
                sba_end = sba.maxOnG(rankGenome) # and not sba.end()[coord]
                # Python automatically takes care of limits, if b >= len(l),
                # l[:b] = []
                for (ib, iDsbb) in enumerate(sbsGXc[ia+1:]):
                    sbb = id2sb(iDsbb)
                    sbb_beg = sbb.minOnG(rankGenome) # and not sbb.beg()[coord]
                    sbb_end = sbb.maxOnG(rankGenome) # and not sbb.end()[coord]
                    if sbb_beg <= sba_end:
                        #overlap
                        if sbb_end <= sba_end:
                            #sbb is included in sba on the 'rankGenome'
                            overlap = sbb_end - sbb_beg + 1
                        else:
                            overlap = sba_end - sbb_beg + 1
                            # if overlap == sba_end - sba_beg + 1, sba is
                            # included in sbb on the 'rankGenome'
                        assert 0 < overlap

                        # record inclusions
                        # The inclusion dict is not symetric: sba might be
                        # included in sbb without sbb being included in sba.
                        if overlap == sba_end - sba_beg + 1:
                            # sba is included in sbb
                            I[iDsba].add(iDsbb)
                        elif overlap == sbb_end - sbb_beg + 1:
                            # sbb is included in sba
                            I[iDsbb].add(iDsba)

                        # The overlap dict and the conflict graph N are
                        # symetric.
                        # 'sba has an overlap with sbb' is exactly the same
                        # has 'sbb has an overlap with sbb'
                        # (same with 'conflict' instead of 'overlap')
                        if overlap <= overlapMax:
                            # record allowed overlap
                            O[iDsba].add(iDsbb)
                            O[iDsbb].add(iDsba)
                        else:
                            # record not allowed overlap
                            # sba and sbb are neighbours in the conflict
                            # graph, they are linked by an edge.
                            N[iDsba].add(iDsbb)
                            N[iDsbb].add(iDsba)
                    else:
                        # sba_end < sbb_beg
                        break
        else:
            # order == -1:
            # sort by ending coordinate first and starting coordinate second
            sbsGXc.sort(key=lambda iDsb: (id2sb(iDsb).maxOnG(rankGenome), id2sb(iDsb).minOnG(rankGenome)), reverse=True)
            # find overlaps on the left side of sbra
            # sbra: sb range 'a'
            # sbrb: sb range 'b'
            for (ia, iDsba) in enumerate(sbsGXc):
                sba = id2sb(iDsba)
                sba_beg = sba.minOnG(rankGenome)
                sba_end = sba.maxOnG(rankGenome)
                # Python automatically takes care of limits, if b >= len(l),
                # l[:b] = []
                for (ib, iDsbb) in enumerate(sbsGXc[ia+1:]):
                    sbb = id2sb(iDsbb)
                    sbb_beg = sbb.minOnG(rankGenome)
                    sbb_end = sbb.maxOnG(rankGenome)
                    if sba_beg <= sbb_end:
                        #overlap
                        if sba_beg <= sbb_beg:
                            #sbb is included in sba on the 'rankGenome'
                            overlap = sbb_end - sbb_beg + 1
                        else:
                            overlap = sbb_end - sba_beg + 1
                            # if overlap == sba_end - sba_beg + 1, sba is
                            # included in sbb on the 'rankGenome'
                        assert 0 < overlap

                        # record inclusions
                        # The inclusion dict is not symetric: sba might be
                        # included in sbb without sbb being included in sba.
                        if overlap == sba_end - sba_beg + 1:
                            # sba is included in sbb
                            I[iDsba].add(iDsbb)
                        elif overlap == sbb_end - sbb_beg + 1:
                            # sbb is included in sba
                            I[iDsbb].add(iDsba)

                        # The overlap dict and the conflict graph N are
                        # symetric.
                        # 'sba has an overlap with sbb' is exactly the same
                        # has 'sbb has an overlap with sbb'
                        # (same with 'conflict' instead of 'overlap')
                        if overlap <= overlapMax:
                            # record allowed overlap
                            O[iDsba].add(iDsbb)
                            O[iDsbb].add(iDsba)
                        else:
                            # record not allowed overlap
                            # sba and sbb are neighbours in the conflict
                            # graph, they are linked by an edge.
                            N[iDsba].add(iDsbb)
                            N[iDsbb].add(iDsba)
                    else:
                        # sbb_end < sba_beg
                        break
        if overlapMax == 0:
            assert len(O.keys()) == 0
        return (N, O, I)

    # method "getItemById" of sbsInPairCompWithIds"
    id2sb = sbsInPairCompWithIds.getItemById

    N = collections.defaultdict(set)
    O = collections.defaultdict(set)
    I = collections.defaultdict(set)
    for c in sbsGX.keys():
        (N, O, I) = findOverlapsOnGenomeOrder(rankGenome, +1, sbsGX[c], N, O, I, id2sb, overlapMax=overlapMax)
        if overlapMax > 0:
            # also needed to check for the overlap in the reverse
            # order
            (N, O, I) = findOverlapsOnGenomeOrder(rankGenome, -1, sbsGX[c], N, O, I, id2sb, overlapMax=overlapMax)
    return (N, O, I)


def buildConflictGraph(sbsInPairComp, overlapMax=0):

    def findOverlaps(sbsInPairCompWithIds, sbsG1, sbsG2, overlapMax=0):
        (N1, O1, I1) = findOverlapsOnGenome(1, sbsInPairCompWithIds, sbsG1, overlapMax=overlapMax)
        (N2, O2, I2) = findOverlapsOnGenome(2, sbsInPairCompWithIds, sbsG2, overlapMax=overlapMax)

        def concatenateDictsOfSets(D1, D2):
            newD = collections.defaultdict(set)
            for (key, value) in (list(D1.iteritems()) + list(D2.iteritems())):
                newD[key] = newD[key].union(value)
            return newD

        N = concatenateDictsOfSets(N1, N2)
        O = concatenateDictsOfSets(O1, O2)
        I = concatenateDictsOfSets(I1, I2)
        return (N, O, I)

    sbsInPairCompWithIds = myTools.OrderedDict2dOfLists()
    for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
        new_sb = SyntenyBlock(sb)
        #sbsInPairCompWithIds[c1][c2].append(new_sb)
        sbsInPairCompWithIds.addToLocation((c1, c2), new_sb)
    #sbsInPairCompWithIds.identifyItems()
    (V, (sbsG1, sbsG2)) = buildSetOfVertices(sbsInPairCompWithIds)
    (N, O, I) = findOverlaps(sbsInPairCompWithIds, sbsG1, sbsG2, overlapMax=overlapMax)
    return (sbsInPairCompWithIds, V, N, O, I, (sbsG1, sbsG2))

# Edit sbs in order to have non-overlapping sbs
@myTools.verbose
def filterOverlappingSbs(sbsInPairComp, overlapMax=0, removeSingleHpSbs=True, verbose=False):
    # The general idea for construction non-overlapping sbs come from:
    # "Nonoverlapping local alignments" (Bafna 1996)
    # The problem of finding the best non-overlapping sbs is analog of
    # the "Idependent subset of Rectangles" (IR) problem.
    # 1st build the conflict graph. With axis-parallel rectangles.
    # find overlaps

    # 1st step, build the graph of conflict graph 'G(V,E)'
    # with 'V' a set of weighted Vertices (weights corresponds to the nb of
    # hps for instance)
    # and 'E' a set of edges between these vertices

    # The conflict-graph G = (V,E)
    # List of all vertices (iD of the sbs and its weight), sorted by
    # decreasing weights.
    # Dictionnary of all the overlapping sbs, (it is analof to E, the set of
    # edges in Bafna 1996)

    @myTools.verbose
    def solveConflictGraph(sbsInPairCompWithIds, V, N, verbose=False):
        # as descibed at the easiest algorithm for the 'IR' problem in Bafna
        # 1996:
        # We solve the conflict graph by taking vertices from highest weights to
        # lowest heights. Each time a vetex is taken, remove all its neighbours
        # and update the set of edges.

        # method "getItemById" of sbsInPairCompWithIds"
        id2sb = sbsInPairCompWithIds.getItemById
        id2location = sbsInPairCompWithIds.getItemLocationById
        Vcopy = set(V)
        # sort V by increasing weights
        sortedV = sorted([v for v in Vcopy], key=lambda iDsb: len(id2sb(iDsb).la))

        noOrSmallOverlapIds = set([])
        removedSbsBecauseOfOverlap = set([])
        while len(Vcopy) > 0:
            iD_selectedSb = sortedV.pop()
            # SB Not Removed
            sbnr = id2sb(iD_selectedSb)
            (c1nr, c2nr, _) = id2location(iD_selectedSb)
            Vcopy.remove(iD_selectedSb)
            noOrSmallOverlapIds.add(iD_selectedSb)
            # update the edges
            for iD_removedSb in N[iD_selectedSb]:
                # SB Removed
                sbr = id2sb(iD_removedSb)
                (c1r, c2r, _) = id2location(iD_removedSb)
                print >> sys.stderr, "sb (%s, %s) with %s tbs is removed because of unallowed overlap with sb (%s, %s) with %s tbs" %\
                    (c1r, c2r, len(sbr.la), c1nr, c2nr, len(sbnr.la))
                for iD_needUpdateSb in N[iD_removedSb]:
                    N[iD_needUpdateSb] = N[iD_needUpdateSb] - N[iD_removedSb]
                    # if N[iD_needUpdateSb] is empty, it will be output as a
                    # non-overlapping sb latter
                del N[iD_removedSb]
                removedSbsBecauseOfOverlap.add(iD_removedSb)
                try:
                   Vcopy.remove(iD_removedSb)
                except:
                    pass
            del N[iD_selectedSb]
            sortedV = [idSb for idSb in sortedV if idSb in Vcopy]

        assert len(V) == len(noOrSmallOverlapIds) + len(removedSbsBecauseOfOverlap), "%s = %s + %s" %\
            (len(V), len(noOrSmallOverlapIds), len(removedSbsBecauseOfOverlap))

        # TODO
        # A last step would be to optimise the set using the 't-opt' algorithm
        # on page 9 of "Nonoverlapping local alignments" (Bafna 1996)

        # return the set of ids of sbs that are kept
        return noOrSmallOverlapIds

    @myTools.verbose
    def truncateSbsWithSmallOverlap(sbsInPairCompWithIds, O, verbose=False):
        id2sb = sbsInPairCompWithIds.getItemById
        id2location = sbsInPairCompWithIds.getItemLocationById
        sortedO = list(O.keys())
        # sort by increasing weight
        sortedO.sort(key=lambda idsb: len(id2sb(idsb).la))
        # set of 0X keys
        setO = set(O.keys())
        while len(setO) > 0:
            # iDsbNotTruncated
            idsbnt = sortedO.pop()
            # DEBUG
            sbnt = id2sb(idsbnt)
            (c1nt, c2nt, _) = id2location(idsbnt)
            for idsbt in O[idsbnt]:
                sbt = id2sb(idsbt)
                (c1t, c2t, _) = id2location(idsbt)
                assert c1nt == c1t or c2nt == c2t
                if c1nt == c1t:
                    newl1 = []
                    newl2 = []
                    newla = []
                    for (i, i1) in enumerate(sbt.l1):
                        if ((i1 < sbnt.minOnG(1)) or (sbnt.maxOnG(1) < i1)):
                            newl1.append(i1)
                            newl2.append(sbt.l2[i])
                            newla.append(sbt.la[i])
                    truncation = len(sbt.la) - len(newla)
                    if truncation > 0:
                        print >> sys.stderr, "sb (%s, %s) with %s tbs truncated (-%s tbs) because of sb (%s, %s) with %s tbs" % (c1t, c2t, len(sbt.la), truncation, c1nt, c2nt, len(sbnt.la))
                    sbt.l1 = newl1
                    sbt.l2 = newl2
                    sbt.la = newla
                # not an elif, since c1nt == c1t and c2nt == c2t is possible
                if c2nt == c2t:
                    newl1 = []
                    newl2 = []
                    newla = []
                    for (i, i2) in enumerate(sbt.l2):
                        if ((i2 < sbnt.minOnG(2)) or (sbnt.maxOnG(2) < i2)):
                            newl1.append(sbt.l1[i])
                            newl2.append(i2)
                            newla.append(sbt.la[i])
                    truncation = len(sbt.la) - len(newla)
                    if truncation > 0:
                        print >> sys.stderr, "sb (%s, %s) with %s tbs truncated (-%s tbs) because of sb (%s, %s) with %s tbs" % (c1t, c2t, len(sbt.la), truncation, c1nt, c2nt, len(sbnt.la))
                    sbt.l1 = newl1
                    sbt.l2 = newl2
                    sbt.la = newla
            # FIXME, are you sure for 'len(id2sb(idsb).la) > 0'
            sortedO = [idsb for idsb in sortedO if idsb in setO and len(id2sb(idsb).la) > 0]
            # Need to re-sort since the length of the sb has changed
            sortedO.sort(key=lambda idsb: len(id2sb(idsb).la))
            setO = set(sortedO)
        return sbsInPairCompWithIds
    (sbsInPairCompWithIds, V, N, O, I, (sbsG1, sbsG2)) =\
        buildConflictGraph(sbsInPairComp, overlapMax=overlapMax)
    nbSbsBeforeOverlapFiltering = len(V)
    noOrSmallOverlapIds = solveConflictGraph(sbsInPairCompWithIds, V, N, verbose=verbose)
    # DEBUG assertion
    assert noOrSmallOverlapIds.issubset(V)

    cptRemovedSbBecauseOfNotAllowedOverlap = 0
    for idsb in V - noOrSmallOverlapIds:
            (c1, c2, idx) = sbsInPairCompWithIds.getItemLocationById(idsb)
            sb = sbsInPairCompWithIds[c1][c2][idx]
            print >> sys.stderr, "Removed sb (%s,%s) of %s hps because of overlap" % (c1, c2, len(sb.la))
            cptRemovedSbBecauseOfNotAllowedOverlap += 1
    sbsInPairCompWithIds.keepIds(noOrSmallOverlapIds)

    # Update the small overlap dict
    newO = {}
    for (idsba, idsbbs) in O.iteritems():
        if idsba in noOrSmallOverlapIds:
            newO[idsba] = set([idsbb for idsbb in idsbbs if idsbb in noOrSmallOverlapIds])
    O = newO

    if overlapMax > 0:
        # Need to edit sbs to truncate sbs that have an allowed overlap
        #new_sbsInPairCompWithIds = truncateSbsWithSmallOverlap(new_sbsInPairCompWithIds, O, verbose=verbose)
        sbsInPairCompWithIds = truncateSbsWithSmallOverlap(sbsInPairCompWithIds, O, verbose=verbose)

    removedIdBecauseEmptySb = set([])
    cptRemovedSbAfterTruncation = 0
    for (idsb, (c1, c2), sb) in sbsInPairCompWithIds.iterByOrderedIds():
        if len(sb.la) == 0 or (len(sb.la) == 1 and removeSingleHpSbs):
            #        # Remove all sbs too small, (they are small because of the preceding
    #        # truncation)
    #        assert overlapMax > 0
            removedIdBecauseEmptySb.add(idsb)
            print >> sys.stderr, "Removed sb (%s,%s) because of truncation (less than 1tb in the sb after truncation)" % (c1, c2)
            cptRemovedSbAfterTruncation += 1
    sbsInPairCompWithIds.removeIds(removedIdBecauseEmptySb)

        #assert set(O.keys()) >= noOrSmallOverlapIds
        # FIXME, is it right to do that ?
        #O = dict([(iDsba, iDsbb)  for (iDsba, iDsbb) in O.iteritems() if ((iDsba in noOrSmallOverlapIds) and (iDsbb in noOrSmallOverlapIds))])

        #assert len(noOrSmallOverlapIds) == len([a for a in new_sbsInPairComp])

    new_sbsInPairComp = myTools.Dict2d(list)
    for (idsb, (c1, c2), sb) in sbsInPairCompWithIds.iterByOrderedIds():
        new_sbsInPairComp[c1][c2].append(sb)

    # This print is already in the upstream function
    #print >> sys.stderr, "Nb sbs before overlap-filtering = %s" % nbSbsBeforeOverlapFiltering

    assert cptRemovedSbBecauseOfNotAllowedOverlap == nbSbsBeforeOverlapFiltering - len(noOrSmallOverlapIds)
    nbSbsAfterOverlapFiltering = len(list(new_sbsInPairComp.iteritems2d()))
    # Total nb of hps removed
    nbHpsBeforeOverlapFiltering = sum([len(sb.la) for (_, sb) in sbsInPairComp.iteritems2d()])
    nbHpsAfterOverlapFiltering = sum([len(sb.la) for (_, sb) in new_sbsInPairComp.iteritems2d()])
    print >> sys.stderr, "Nb hps before overlap-filtering = %s" % nbHpsBeforeOverlapFiltering
    print >> sys.stderr, "Nb hps after overlap-filtering = %s" % nbHpsAfterOverlapFiltering
    assert nbSbsAfterOverlapFiltering == nbSbsBeforeOverlapFiltering - (cptRemovedSbAfterTruncation + cptRemovedSbBecauseOfNotAllowedOverlap)
    print >> sys.stderr, "Nb sbs removed because of not allowed overlap = %s" % cptRemovedSbBecauseOfNotAllowedOverlap
    print >> sys.stderr, "Nb sbs removed during truncation = %s" % cptRemovedSbAfterTruncation
    # This print is already in the upstream function
    #print >> sys.stderr, "Nb non-overlapping sbs returned = %s" % nbSbsAfterOverlapFiltering

    return new_sbsInPairComp


def strandProduct(sa, sb):
    if sa != None and sb != None:
        return sa*sb
    else:
        return None

# Return the hm (or the mhp) of two compared chromosomes rewritten with family names instead of gene (or tb) names
# inputs :
#       gcX : a chromosome of the genome X rewritten with the family names : [...,(f,s),...]
# ouputs :
#       M : mh or mhp depending on whether gc1 and gc2 are chromosomes written in genes or in tbs
#               M is a dict of dict : M[i1][i2] = hpSign, this structure is lighter since M is a sparse matrix
#       locG2 : a dict {..., f:[(i2,s2)],...} with i2 the indices of the genes on gc2. For each f locG2 gives all its genes positions in the chromosome 2
# f : family, often the ancGene index
###################################################################################################################################
def homologyMatrix(gc1, gc2):
    #1
    #M={}
    #locG2 = {}
    #for (i2,(f2,s2)) in enumerate(gc2):
    #       if f2 != None:
    #               if f2 not in locG2:
    #                       locG2[f2]=[]
    #                       locG2[f2].append( (i2,s2) ) # LOCalisation on Gc 2
    #               for (i1,(f1,s1)) in enumerate(gc1):
    #                       if f1 == f2:
    #                               if i1 not in  M:
    #                                       M[i1]={}
    #                               M[i1][i2] = strandProduct(s1,s2)

    #2 faster than 1
    locG2 = {}
    for (i2,(f,s2)) in enumerate(gc2):
        if f != None:
            if f not in locG2:
                locG2[f]=[]
            locG2[f].append( (i2,s2) ) # LOCalisation on Gc 2
            #locG2.setdefault(f,[]).append( (i2,s2) ) # LOCalisation on Gc 2, same excecution time but less readable
    M={}
    if not locG2: # if locG2 is empty
        return (M, locG2)

    for (i1,(f,s1)) in enumerate(gc1):
        if f !=None and f in locG2: #TODO : remove by advance the f == None from gc1
            M[i1]={}
            for (i2,s2) in locG2[f]:
                M[i1][i2]= strandProduct(s1,s2)

    #3 (faster than 2), but since we need to return locG2 of the input gc2 for next computations, it is not convenient
    #TODO need to add a boolean to the returned values because locG2 may not correspond to the gc2 chromosome
    # use the dict locG2 on the smallest chromosome
    #switchedGCs = False
    #if len(gc1) < len(gc2):
    #       (gc1,gc2) = (gc2,gc1)
    #       switchedGCs = True
    #       if f != None:
    #               if f not in locG2:
    #                       locG2[f]=[]
    #               locG2[f].append( (i2,s2) ) # LOCalisation on Gc 2
    #M={}
    #for (i1,(f,s1)) in enumerate(gc1):
    #       if f != None and f in locG2: #TODO : remove by advance the f == None from gc1
    #               M[i1]={}
    #               for (i2,s2) in locG2[f]:
    #                       M[i1][i2]= strandProduct(s1,s2)
    #if switchedGCs == True:
    #       (gc1,gc2) = (gc2,gc1)
    #       switchedGCs = False

    return (M, locG2)

#
# At the beginning of a diagonal this function is called to set a diagonal type
# inputs :
#       i1,i2 : the indices of the first hp of the current diagonal in M
#       M : the homology matrix (either mh or mhp)
#       consistentDiagonals : boolean (true if you want consistent diagonals otherwise false)
# output :
#       diagType :
#               either a bottom-left to top-right (slash : '/') diagonal
#               or top-left to bottom right (backslash : '\') diagonal
################################################################################################
def findDiagType(i1, i2, M, consistentSwDType):
    if M[i1][i2] != None and consistentSwDType:
        diagType = '/' if M[i1][i2] == +1 else '\\'
    else:
        if i1+1 in M and i2+1 in M[i1+1] and ((M[i1+1][i2+1] in [1,None]) if consistentSwDType else True):
            diagType = '/'
        elif i1+1 in M and i2-1 in M[i1+1] and ((M[i1+1][i2-1] in [-1,None]) if consistentSwDType else True):
            diagType = '\\'
        else:
            diagType = None
    return diagType


@myTools.verbose
def extractSbsInPairCompChrWrapper(c1, c2, gc1, gc2, gapMax=0, distanceMetric = 'DPD',
                                   consistentSwDType=True, verbose=False):
    print >> sys.stderr, "(PPID = %s, PID = %s) start to extract diagonals on G1[%s]_vs_G2[%s]" %\
        (os.getppid(), os.getpid(), c1, c2)
    listOfDiags = extractSbsInPairCompChr(gc1, gc2, gapMax=0, distanceMetric = 'DPD',
                                          consistentSwDType=True, verbose=False)
    # FIXME diagType=diag[0] could be added as an information on the diagonal here
    return ((c1, c2), listOfDiags)

# Extract sbs in a pairwise comparison of two chromosomes
############################################################
@myTools.verbose
def extractSbsInPairCompChr(gc1, gc2, gapMax=0, distanceMetric='DPD',
                            consistentSwDType=True, verbose=False):
    listOfDiags = []
    (M,locG2) = homologyMatrix(gc1, gc2)
    if not locG2: # if locG2 is empty
        return listOfDiags
    la=[]
    l1=[]
    l2=[]
    diagType = None
    i1_old = None

    #TODO scan M instead of gc1 : impossible since 'dict size cannot change during a loop over its items'
    # scan M from left to right
    for (i1,(f,_)) in enumerate(gc1):
        # f in locG2 means that i1 has at least one homolog on gc2, locG2 never changes, that is why we iterate over locG2
        if f != None and f in locG2:
            i1_old=i1
            # scan M from bottom to top
            for (i2,_) in locG2[f]: # the family name of the gene at index i2 is the same as f
                # When a diagonal is extracted the scanning of M must continue right after the first hp of the extracted diagonal
                i1=i1_old
                # TODO write a recursive function easier understanding
                # While hps can be added to a diagonal
                while i1 in M and i2 in M[i1]:
                    # Here a diagonal is started or hps are added to an already started diagonal
                    f = gc1[i1][0]
                    if len(la) == 0:
                        # first hp of a diagonal, the diagonalType is searched
                        diagType = findDiagType(i1,i2,M,consistentSwDType)

                    # A hp has a unique associated gene family (often a unique ancestral gene)
                    # orientation of tbs on gc1 are used as a reference, if the orientation of a tb on gc1 is unknown, it is possible to infer the ancestral orientation by using the diagType and the orientation of the homologous tb on gc2
                    if gc1[i1][1] != None:
                        ancestralStrand = gc1[i1][1]
                    elif diagType == '/' and gc2[i2][1] != None:
                        ancestralStrand = gc2[i2][1]
                    elif diagType == '\\' and gc2[i2][1] != None:
                        ancestralStrand = -gc2[i2][1]
                    else:
                        ancestralStrand = None

                    # la = [..., (aGn, aGs, dist), ...]
                    la.append((f, ancestralStrand, len(la)+1))
                    l1.append(i1)
                    l2.append(i2)
                    del M[i1][i2]
                    if len(M[i1].keys()) == 0:
                        del M[i1]
                    #diagType == None if sameStrand == False and len(la) == 1, first element of a diagonal which we donnot take care of gene orientation
                    #if (diagType == "/" or diagType == None) and i1+1 in M and i2+1 in M[i1+1] and ((M[i1+1][i2+1] in [+1,None]) if consistentSwDType else True):
                    if diagType == "/" and i1+1 in M and i2+1 in M[i1+1] and ((M[i1+1][i2+1] in [+1,None]) if consistentSwDType else True):
                        i1=i1+1
                        i2=i2+1
                        #assert i2-l2[-1] == +1
                        #assert i1-l1[-1] == +1
                        #assert i2 in M[i1]

                    elif diagType=="\\" and i1+1 in M and i2-1 in M[i1+1] and ((M[i1+1][i2-1] in [-1,None]) if consistentSwDType else True):
                        i1=i1+1
                        i2=i2-1
                        #assert i2-l2[-1] == -1
                        #assert i1-l1[-1] == +1
                        #assert i2 in M[i1]
                    else:
                        # Since no more hps can be added to the current diagonal, the diagonal is recorded
                        #assert len(la) > 0
                        listOfDiags.append(Diagonal(diagType,l1,l2,la))
                        l1=[]
                        l2=[]
                        la=[]
                        diagType=None
                        break # exit the while loop and iter the for loop
    # merging process, fuse diagonals
    if len(listOfDiags) > 0 and gapMax >= 0:
        listOfDiags = mergeSbs(listOfDiags, gapMax, gc2, distanceMetric=distanceMetric, verbose=False)
    return listOfDiags


def crossGeneContent(g1, g2):
    # create a set with all gene names
    geneNames1 = myMapping.setOfGeneNames(g1)
    geneNames2 = myMapping.setOfGeneNames(g2)
    gNsInCommon = geneNames1.intersection(geneNames2)
    gNsInG1notInG2 = geneNames1 - gNsInCommon
    gNsInG2notInG1 = geneNames2 - gNsInCommon
    return (gNsInCommon, gNsInG1notInG2, gNsInG2notInG1)


# Depending on the filterType parameter:
#       - filterType = FilterType.None (genomes are not filtered)
#       - filterType = FilterType.InCommonAncestor (only genes herited from the ancestor are kept)
#       - filterType = FilterType.InBothSpecies (only 'anchor genes', ie genes present in both species, are kept)
# Returns (g1,g2,trans1,trans2)
#      - g1 the genome 'g1' rewritten with ancGenes ID
#       - trans1 = { ..., newi:oldi, ...} with newi, the index of a gene in 'g1' and oldi the index of the same gene in the original genome 'g1'
def filter2D(g1_orig, g2_orig, filterType, minChromLength, keepOriginal=False):
    # Mark genes that are not in the intersection for future removal
    # Marking is done by switching (g,s) into (None,s)
    # warning : modifies g1 and g2
    if keepOriginal:
        g1 = copy.deepcopy(g1_orig)
        g2 = copy.deepcopy(g2_orig)
    else:
        g1 = g1_orig
        g2 = g2_orig

    remapFilterSize = myMapping.remapFilterSize
    remapCoFilterContentAndSize = myMapping.remapCoFilterContentAndSize

    # In all cases filtering on length
    # All genes are conserved
    if filterType == FilterType.None:
        (g1, mG1f2G1o, (nCL1, nGL1)) = remapFilterSize(g1, minChromLength)
        (g2, mG2f2G2o, (nCL2, nGL2)) = remapFilterSize(g2, minChromLength)
    # Only genes herited from an ancestral genes are conserved
    elif filterType == FilterType.InCommonAncestor:
        (g1, mG1f2G1o, (nCL1, nGL1)) = remapCoFilterContentAndSize(g1, set([None]), minChromLength)
        (g2, mG2f2G2o, (nCL2, nGL2)) = remapCoFilterContentAndSize(g2, set([None]), minChromLength)
    # Only conserve genes present in both extant genomes
    elif filterType == FilterType.InBothSpecies:
        mG1f2G1o = None
        mG2f2G2o = None
        while True:
            # after this step genes that have no homolog in the other genome are marked None
            # 'gNsInCommon' is also called the set of 'anchor genes' in bibliography
            (gNsInCommon, gNsInG1notInG2, gNsInG2notInG1) = crossGeneContent(g1, g2)
            removedGNs1 = gNsInG1notInG2 | set([None])
            removedGNs2 = gNsInG2notInG1 | set([None])
            (g1, mG1f2G1o, (nCL1, nGL1)) =\
                remapCoFilterContentAndSize(g1, removedGNs1,
                                            minChromLength,
                                            mOld=mG1f2G1o)
            (g2, mG2f2G2o, (nCL2, nGL2)) =\
                remapCoFilterContentAndSize(g2, removedGNs2,
                                            minChromLength,
                                            mOld=mG2f2G2o)
            hasChanged = (nGL1 > 0) or (nGL2 > 0)
            # If a chromosome has been removed because of the filtering on the length,
            # the filtering is performed once more.
            if not hasChanged:
                break
    else:
        # impossible case
        raise

    # In order to write genomes in files and verify sbs predictions
    # for c in g2:
    #       for (i,(g,s)) in enumerate(g2[c]):
    #               print >> sys.stderr, c, "\t", s, "\t", g2.lstGenes[c][trans2[c][i]].names[0]

    return ((g1, mG1f2G1o, (nCL1, nGL1)), (g2, mG2f2G2o, (nCL2, nGL2)))


# compute the adviced maximal gap length (in genes) for the diagonal merger
# We compute the gap that gives the closest probability p-value(sb(nbHps,gap,N12,N1,N2)) to targetProba in an a MHP(N1,N2,N12).
# Usually this MHP is an average of all the MHPs of the whole genome comparison. Even if considering this average MHP is not very rigorous, we only need reasonable rather than optimal values for recommendedGap.
@myTools.verbose
def recommendedGap(nbHps, targetProba, N12, N1, N2, p_hpSign=None, maxGapThreshold=20, verbose=False):
    firstPrint=True
    tries = []
    for g in range(0,maxGapThreshold+1):
        L = (nbHps-1)*g+nbHps # lengths of the chromosomal windows
        pVal = myProbas.pValue(nbHps,g,L,L,N12,N1,N2,p_hpSign,verbose=verbose)
        if firstPrint == True:
            print >> sys.stderr, "P-values of sbs of %s hps (nbHpsRecommendedGap) hps spaced by a gap 'g' in an average MHP(N1=%s, N2=%s, N12=%s) over all MHPs involved in the whole genome comparison.:" % (nbHps,N1,N2,N12)
            firstPrint=False
        print >> sys.stderr, "nbHps=%s, g=%g, pVal=%s" % (nbHps, g,pVal)
        if pVal != None:
            tries.append(pVal)
        elif pVal == None:
            tries.append(-float('inf'))
        else:
            raise ValueError
    print >> sys.stderr, "The first p-value over targetProba=%s (targetProbaRecommendedGap) sets the recommended gap value" % targetProba
    for (g,pVal) in  enumerate(tries):
        if pVal > targetProba:
            return g
    #min_index, min_value = min(enumerate(tries), key=lambda x:x[1])
    #g = min_index
    #return g


# Number of non-empty value in the mh or the mhp
#@myTools.tictac
@myTools.verbose
def numberOfHomologies(g1, g2, verbose=False):
    nbHomologies = myTools.Dict2d(int)
    listOfPercentage = range(0, 101, 5)[1:]
    nbPairwiseComparisons = len(g1) * len(g2)
    print >> sys.stderr, "pairwise comparison of chromosomes analysis for counting hps",
    for (i, (c1, c2)) in enumerate(itertools.product(g1, g2)):
        gc1 = g1[c1]
        gc2 = g2[c2]
        (Ms, _) = homologyMatrix(gc1, gc2)
        progress = int(float(i*100)/nbPairwiseComparisons)
        if progress in listOfPercentage:
            print >> sys.stderr, "%s" % progress + "%",
            listOfPercentage.remove(progress)
        nbHomologies[c1][c2] = sum([len(Ms[i1]) for i1 in Ms])
    # new line in the print
    print >> sys.stderr, ""
    nbHomologies_g = sum([nbH for nbH in nbHomologies.values2d()])
    #assert nbHomologies_g >= min(sum([len(g1[c]) for c in g1]), sum([len(g2[c]) for c in g2])),"%s,%s" %  (sum([len(g1[c]) for c in g1]), sum([len(g2[c]) for c in g2]))
    #Not needed since the two lineages may have undergone some differential gene losses
    return nbHomologies, nbHomologies_g


# Statistical validation of synteny blocks
# because of collections.Counter
@myTools.minimalPythonVersion((2, 7))
@myTools.verbose
def statisticalValidation(sbsInPairComp, g1_tb, g2_tb, N12s, p_hpSign,
                          pThreshold=0.001, NbOfHomologiesThreshold=50,
                          validateImpossToCalc_mThreshold=3, verbose=False):

    def fNbDiags(sbsInPairCompX):
        return len(list(sbsInPairCompX.iteritems2d()))

    sbsInPairCompRejected = myTools.Dict2d(list)
    sbsInPairCompImpossibleToCalcProba = myTools.Dict2d(list)
    sbsInPairCompStatVal = myTools.Dict2d(list)
    for ((c1, c2), diag) in sbsInPairComp.iteritems2d():
        (m, max_g, lw1, lw2, l1_min, l1_max, l2_min, l2_max) = diag.calculateCharacteristics()
        na = len(g1_tb[c1])  # Nb of Tbs on C1
        nb = len(g2_tb[c2])  # Nb of Tbs on C2
        nab = N12s[c1][c2]  # Nb of homologies in the mhp
        if m <= 1:
            # all diagonals of length 1 are rejected
            sbsInPairCompRejected[c1][c2].append((diag, None))
            continue
        elif m > nab:
            # there are not m hps in the MHP
            sbsInPairCompImpossibleToCalcProba[c1][c2].append((diag, None))
            continue
        elif m > min([lw1, lw2]):
            # there are too many dispersed paralogies in the window
            sbsInPairCompImpossibleToCalcProba[c1][c2].append((diag, None))
            continue

        # This is to avoid excessive time consuming computations
        if m > NbOfHomologiesThreshold:
            p = 0
        else:
            p = myProbas.pValue(m, max_g, lw1, lw2, nab, na, nb, p_hpSign[c1][c2], verbose=verbose)

        if p is None:
            sbsInPairCompImpossibleToCalcProba[c1][c2].append((diag, None))
        elif p < pThreshold:
            sbsInPairCompStatVal[c1][c2].append((diag, p))
        else:
            sbsInPairCompRejected[c1][c2].append((diag, p))
            assert len(diag.l1) == len(diag.l1)
            assert len(diag.la) == len(diag.l2)

    # automatically validate sbs that contain more than validateImpossToCalc_mThreshold hps
    sbsInPairCompImpossibleToCalcProbaStatVal = myTools.Dict2d(list)
    sbsInPairCompImpossibleToCalcProbaRejected = myTools.Dict2d(list)
    for ((c1, c2), (diag, pVal)) in sbsInPairCompImpossibleToCalcProba.iteritems2d():
        m = diag.calculateCharacteristics()[0]
        if m >= validateImpossToCalc_mThreshold:
            # if the diagonal contains more than validateImpossToCalc_mThreshold
            # hps, it is validated in order to avoid to reject long and perfect
            # diagonals because of only one dispersed tandem duplication
            sbsInPairCompImpossibleToCalcProbaStatVal[c1][c2].append((diag, pVal))
        else:
            sbsInPairCompImpossibleToCalcProbaRejected[c1][c2].append((diag, pVal))


    assert fNbDiags(sbsInPairCompImpossibleToCalcProbaStatVal) + fNbDiags(sbsInPairCompImpossibleToCalcProbaRejected) == fNbDiags(sbsInPairCompImpossibleToCalcProba),\
        "%s + %s = %s" % (fNbDiags(sbsInPairCompImpossibleToCalcProbaStatVal), fNbDiags(sbsInPairCompImpossibleToCalcProbaRejected), fNbDiags(sbsInPairCompImpossibleToCalcProba))

    # update the lists
    sbsInPairCompStatVal = sbsInPairCompStatVal + sbsInPairCompImpossibleToCalcProbaStatVal
    sbsInPairCompRejected = sbsInPairCompRejected + sbsInPairCompImpossibleToCalcProbaRejected

    #print >> sys.stderr, fNbDiags(sbsInPairCompStatVal)
    #print >> sys.stderr, fNbDiags(sbsInPairCompRejected)
    #print >> sys.stderr, fNbDiags(sbsInPairComp)

    assert fNbDiags(sbsInPairCompStatVal) + fNbDiags(sbsInPairCompRejected) == fNbDiags(sbsInPairComp),\
        "%s + %s = %s" % (fNbDiags(sbsInPairCompStatVal), fNbDiags(sbsInPairCompRejected), fNbDiags(sbsInPairComp))

    firstPrint = True
    for ((c1, c2), (diag, pVal)) in sbsInPairCompRejected.iteritems2d():
        (m, max_g, lw1, lw2, l1_min, l1_max, l2_min, l2_max) = diag.calculateCharacteristics()
        na = len(g1_tb[c1])  # Nb of Tbs on C1
        nb = len(g2_tb[c2])  # Nb of Tbs on C2
        nab = N12s[c1][c2]  # Nb of homologies in the mhp
        if m>1:
            if firstPrint:
                print >> sys.stderr, "Diagonals of more than 1 hp that have been rejected during the statistical test:"
                firstPrint = False
            # TODO comment se fait-il que certaines diagonales soient en double ici ?
            print >> sys.stderr, "(c1=%s:%s-%s,c2=%s:%s-%s) \t (m=%s, max_g=%s, lw1=%s lw2=%s, nab=%s, na=%s, nb=%s)" %\
                (c1,l1_min,l1_max,c2,l2_min,l2_max,m,max_g,lw1,lw2,nab,na,nb), "\t p=", pVal

    firstPrint = True
    for ((c1, c2), (diag, pVal)) in sbsInPairCompStatVal.iteritems2d():
        (m, max_g, lw1, lw2, l1_min, l1_max, l2_min, l2_max) = diag.calculateCharacteristics()
        if m == 2:
            na = len(g1_tb[c1])  # Nb of Tbs on C1
            nb = len(g2_tb[c2])  # Nb of Tbs on C2
            nab = N12s[c1][c2]  # Nb of homologies in the mhp
            if firstPrint:
                print >> sys.stderr, "Diagonals containing 2 hps that have passed the statistical test:"
                firstPrint = False
            print >> sys.stderr, "(c1=%s:%s-%s,c2=%s:%s-%s) \t (m=%s, max_g=%s, lw1=%s lw2=%s, nab=%s, na=%s, nb=%s)" %\
                (c1,l1_min,l1_max,c2,l2_min,l2_max,m,max_g,lw1,lw2,nab,na,nb), "\t p=", pVal

    def statsDiagsM(sbsInPairCompX, m):
        if len(sbsInPairCompX.keys2d()) == 0:
            return 'RAS'
        (k1Foo, k2Foo) = sbsInPairCompX.keys2d()[0]
        if isinstance(sbsInPairCompX[k1Foo][k2Foo][0], tuple):
            sbsInPairCompXM = [diag for sbs in sbsInPairCompX.values2d() for (diag, pValue) in sbs if len(diag.la) == m]
        else:
            # sbsInPairComp has no pValue
            assert isinstance(sbsInPairCompX[k1Foo][k2Foo][0], Diagonal)
            sbsInPairCompXM = [diag for diags in sbsInPairCompX.values2d() for diag in diags if len(diag.la) == m]
        sbsInPairCompXMGaps = [diag.max_g() for diag in sbsInPairCompXM if len(diag.la) == m]
        diagsXMGaps = collections.Counter(sbsInPairCompXMGaps)
        diagsXMGaps = ["%s:%s" % (length, nb) for (length, nb) in sorted(diagsXMGaps.items())]
        return diagsXMGaps
    print >> sys.stderr, "Over all sbs of 2 hps before stat. val. the distribution of gap maximum is: {%s}" % " ".join(statsDiagsM(sbsInPairComp,2))
    print >> sys.stderr, "Over all rejected sbs of 2 hps the distribution of gap maximum is: {%s}" % " ".join(statsDiagsM(sbsInPairCompRejected,2))
    print >> sys.stderr, "Over all sbs where it is imposs. to calculate the proba. with 2 hps the distribution of gap maximum is : {%s} (due to dispersed paralogies)" % " ".join(statsDiagsM(sbsInPairCompImpossibleToCalcProba,2))
    print >> sys.stderr, "Over all sbs where it is imposs. to calculate the proba. with 2 hps which are finally validated, the distribution of gap maximum is: {%s} (due to dispersed paralogies)" % " ".join(statsDiagsM(sbsInPairCompImpossibleToCalcProbaStatVal,2))
    print >> sys.stderr, "Over all sbs where it is imposs. to calculate the proba. with 2 hps which are finally rejected, the distribution of gap maximum is: {%s} (due to dispersed paralogies)" % " ".join(statsDiagsM(sbsInPairCompImpossibleToCalcProbaRejected,2))

    print >> sys.stderr, "Over all stat. val. sbs of 2 hps the distribution of gap maximum is: {%s}" % " ".join(statsDiagsM(sbsInPairCompStatVal,2))

    def statsDiagLengths(sbsInPairCompX):
        if len(sbsInPairCompX.keys2d()) == 0:
            return 'RAS'
        (k1Foo, k2Foo) = sbsInPairCompX.keys2d()[0]
        if isinstance(sbsInPairCompX[k1Foo][k2Foo][0], tuple):
            sbsInPairCompX_ = [diag for sbs in sbsInPairCompX.values2d() for (diag, pValue) in sbs]
        else:
            # sbsInPairComp has no pValue
            #assert sbsInPairCompX[k1Foo][k2Foo][0].__class__.name == 'Diagonal'
            assert isinstance(sbsInPairCompX[k1Foo][k2Foo][0], Diagonal)
            sbsInPairCompX_ = [diag for diags in sbsInPairCompX.values2d() for diag in diags]
        diagXLengths = collections.Counter([len(diags.la) for diags in sbsInPairCompX_])
        diagXLengths = ["%s:%s" % (length, nb) for (length, nb) in sorted(diagXLengths.items())]
        return diagXLengths
    print >> sys.stderr, "Over all diagonal before the stat. val., distribution of the diag lengths: {%s}" % " ".join(statsDiagLengths(sbsInPairComp))
    print >> sys.stderr, "Over all rejected diagonals, distribution of all diag lengths: {%s}" %  " ".join(statsDiagLengths(sbsInPairCompRejected))
    print >> sys.stderr, "Over all diagonals with p-Value impossible to compute (mainly due to paralogies), distribution of diag lengths: {%s}" %  " ".join(statsDiagLengths(sbsInPairCompImpossibleToCalcProba))
    print >> sys.stderr, "Over all diagonals with p-Value impossible to compute (mainly due to paralogies) which are finally validated, distribution of diag lengths: {%s}" %  " ".join(statsDiagLengths(sbsInPairCompImpossibleToCalcProbaStatVal))
    print >> sys.stderr, "Over all diagonals with p-Value impossible to compute (mainly due to paralogies) which are finally rejected, distribution of diag lengths: {%s}" %  " ".join(statsDiagLengths(sbsInPairCompImpossibleToCalcProbaRejected))
    print >> sys.stderr, "Over all diagonals that passed the stat. val., distribution of diag lengths: {%s}" %  " ".join(statsDiagLengths(sbsInPairCompStatVal))

    for (c1, c2) in sbsInPairCompStatVal.keys2d():
        sbsInPairCompStatVal[c1][c2] = [SyntenyBlock(diag, pVal) for (diag, pVal) in sbsInPairCompStatVal[c1][c2]]
    return sbsInPairCompStatVal


# Split a sb as soon as one of its gap contains another sb
# the other sb may be a micro-inversion or a transposition-like sb
def fIdentifyBreakpointsWithinGaps(sbsInPairComp, removeSingleHpSbs=True):

    # Check if list2 is nested in the gaps of list1
    # list1: contains the indices of the first sb
    # list2: contains the indices of the second sb
    def noOverlapList(list1, list2):
        setList1 = set(list1)
        setList2 = set(list2)
        if len(setList1 & setList2) > 0:
            # TODO: for ancestral genome reconstruction, the fractionation
            # of the segmental duplication should be taken into account to
            # find what was the ancestral gene order.
            # This fractionation should at first sight be similar to the
            # fractionation that follows a WGD.

            #shortestListLength = min(len(setList1), len(setList2))
            #if float(len(setList1 & setList2)) > (0.8 * shortestListLength):
                #print >> sys.stderr, "A probable segmental duplication has been seen (>80% genes in common with the smallest sb of the pair of overlapping sbs)"
            return False
        else:
            return True

    def noOverlapSb(((c1a, c2a), sba), ((c1b, c2b), sbb)):
        sba = id2sb(idsba)
        sbb = id2sb(idsbb)
        if c1a != c1b and c2a != c2b:
            return False
        res = True
        if c1a == c1b:
            res = res and noOverlapList(sba.l1, sbb.l1)
        if c2a == c2b:
            res = res and noOverlapList(sba.l2, sbb.l2)
        return res

    # Split nested lists 'list1' and 'list2' by introducing breakpoints in gaps
    # of the list that contain nested indices of the other list.
    def splitNestedLists(lista, listb):
        alternatingList = [('a', ia) for ia in lista] + [('b', ib) for ib in listb]
        alternatingList.sort(key=lambda x: x[1])
        # DEBUG assertion
        # check that there are no duplicated indices
        #l = [x for (_, x) in alternatingList]
        #assert len(set([x for x in l if l.count(x) > 1])) == 0,\
        #    "Indices of list2 are not nested within the gaps of list1"
        #
        splitRanks_a = []
        splitRanks_b = []
        nbIDa = 0
        nbIDb = 0
        iDold = sys.maxint
        for (iD, _) in alternatingList:
            if iD == 'a':
                if iD != iDold:
                    splitRanks_a.append(nbIDa)
                nbIDa += 1
            else:
                assert iD == 'b'
                if iD != iDold:
                    splitRanks_b.append(nbIDb)
                nbIDb += 1
            iDold = iD
        splitRanks_a.append(nbIDa)
        splitRanks_b.append(nbIDb)
        return (splitRanks_a, splitRanks_b)

        #def adjIdxs(listX):
        #    lAdjIdxs = []
        #    i1Old = sys.maxint
        #    for i1 in listX:
        #        if i1 != i1Old + 1:
        #            # start a new strict diagonal
        #            lAdjIdxs.append([i1])
        #        else:
        #            lAdjIdxs[-1].append(i1)
        #        i1Old = i1
        #    return lAdjIdxs
        #filledList1 = set(range(min(list1), max(list1) +1))
        #gapsList1 = sorted(list(filledList1 - setList1))
        #gapsList1 = adjIdxs(gapsList1)

    def splitSbBySplitRanks(sb, splitRanks):
        newSbs = []
        for (xBeg, xEnd) in myTools.myIterator.slidingTuple(splitRanks):
            newDt = sb.dt
            newl1 = sb.l1[xBeg:xEnd]
            newl2 = sb.l2[xBeg:xEnd]
            newla = sb.la[xBeg:xEnd]
            # Here we assign the former pVal to all sub-sbs. The
            # hypothesis is to consider that since the sub-sbs were
            # neighbours within the former sb, their significances have
            # already been assessed by the former pValue calculation.
            newPval = sb.pVal
            newSbs.append(SyntenyBlock(Diagonal(newDt, newl1, newl2, newla), newPval))
        return newSbs

    def splitNestedSbs(((c1a, c2a), sba), ((c1b, c2b), sbb)):
        # DEBUG assertion
        assert len(sba.l1) == len(sba.l2) == len(sba.la)
        assert len(sbb.l1) == len(sbb.l2) == len(sbb.la)
        #
        if c1a == c1b:
            (splitRanks_1a, splitRanks_1b) = splitNestedLists(sba.l1, sbb.l1)
            (splitRanks_a, splitRanks_b) = (splitRanks_1a, splitRanks_1b)
        if c2a == c2b:
            (splitRanks_2a, splitRanks_2b) = splitNestedLists(sba.l2, sbb.l2)
            # if the diagType is '\', l2 starts from the highest idx on c2 and
            # ends with the lowest idx.
            if sba.dt == '\\':
                splitRanks_2a = [len(sba.l2) - rank for rank in  splitRanks_2a]
            if sbb.dt == '\\':
                splitRanks_2b = [len(sbb.l2) - rank for rank in  splitRanks_2b]
            (splitRanks_a, splitRanks_b) = (splitRanks_2a, splitRanks_2b)
        if c1a == c1b and c2a == c2b:
            # take into account the splits in the two genomes
            # For instance a short transposition
            # C1 = +A+B+C+D+E+F
            # C2 = +A+D+B+C+E+F
            # if gapMax >=2:
            #   sba = [+A+D+E+F]
            #   sbb = [+B+C]
            #   splitRanks_1a = [1]
            #   splitRanks_2a = [2]
            #   splitRanks_1b = []
            #   splitRanks_2b = []
            #   splitRanks_a = [1, 2] # which means a split before element
            #   at idx 1 (D) and another split before element at idx 2 (E)
            #   sba -> [+A] [+D] [+E+F]
            #   Yields 4 sbs: [A], [D], [+E+F], and [+B+C]
            splitRanks_a = sorted(list(set(splitRanks_1a + splitRanks_2a)))
            splitRanks_b = sorted(list(set(splitRanks_1b + splitRanks_2b)))

        #newSbs = splitSbBySplitRanks(sba, splitRanks_a)
        #newSbs.extend(splitSbBySplitRanks(sbb, splitRanks_b))
        return (splitRanks_a, splitRanks_b)


    # Find sbs that are nested within gaps of other sbs
    #   * nested inversions (micro-inversions, nested segments and inversed
    #     diagonalType).
    #   * transpositions (only 'one' nested segment, to avoid considering
    #     segmental duplications followed by fractionation).
    sbsInPairCompWithIds = myTools.OrderedDict2dOfLists()
    for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
        new_sb = SyntenyBlock(sb)
        sbsInPairCompWithIds.addToLocation((c1, c2), new_sb)
    (V, (sbsG1, sbsG2)) = buildSetOfVertices(sbsInPairCompWithIds)
    (N1, O1, I1) = findOverlapsOnGenome(1, sbsInPairCompWithIds, sbsG1, overlapMax=0)
    (N2, O2, I2) = findOverlapsOnGenome(2, sbsInPairCompWithIds, sbsG2, overlapMax=0)
    # method "getItemById" of sbsInPairCompWithIds"
    id2sb = sbsInPairCompWithIds.getItemById
    id2location = sbsInPairCompWithIds.getItemLocationById

    todo = set([])
    id2splitRanks = collections.defaultdict(list)

    # Intra-synteny block rearrangement: either an intra-sb
    # micro-transposition or an intra-sb micro-inversion (inversed
    # diagType). Or a mix of both. To disantangle both phenomenon, it
    # would be necessary to check the rank location in the host sb.
    for idsbb in (set(I1.keys()) & set(I2.keys())):
        for idsba in (I1[idsbb] & I2[idsbb]):
            # sbb is included in sba in both genomes
            sba = id2sb(idsba)
            sbb = id2sb(idsbb)
            (c1a, c2a, _) = id2location(idsba)
            (c1b, c2b, _) = id2location(idsbb)
            assert set([c1a, c2a]) == set([c1b, c2b])
            assert c1a == c1b and c2a == c2b
            #
            if noOverlapSb(((c1a, c2a), sba), ((c1b, c2b), sbb)):
                # FIXME: removed because it is interesting to split the
                # diagonal in the case of an intra-synteny-block transposition
                #if sba.dt == '/' and sbb.dt == '\\' or sba.dt == '\\' and sbb.dt == '/':

                # Return tuples (id du sb splité, splitRanks)
                # This allows to incrementaly update the splitRanks
                #  and then it wil be possible to perform the splits
                #  after the for loops
                (splitRanks_a, splitRanks_b) = splitNestedSbs(((c1a, c2a), sba), ((c1b, c2b), sbb))
                id2splitRanks[idsba] = sorted(list(set(id2splitRanks[idsba] + splitRanks_a)))
                id2splitRanks[idsbb] = sorted(list(set(id2splitRanks[idsbb] + splitRanks_b)))
                todo = todo | set([idsba, idsbb])

    # extra-sb rearrangement, for instance a transposition
    # FIXME: allow only 'one' nested segment, to avoid considering segmental
    # duplications followed by fractionation.
    # Rq: '^' means symmetric difference
    # A ^ B = set with elements in either A or B but not both
    for idsbb in (set(I1.keys()) ^ set(I2.keys())):
        if idsbb in I1:
            idsbas = I1[idsbb]
        else:
            assert idsbb in I2
            idsbas = I2[idsbb]
        for idsba in idsbas:
            # sbb is included in sba only in one genome
            sba = id2sb(idsba)
            sbb = id2sb(idsbb)
            (c1a, c2a, _) = id2location(idsba)
            (c1b, c2b, _) = id2location(idsbb)
            # DEBUG assertion
            assert c1a == c1b or c2a == c2b
            #
            if noOverlapSb(((c1a, c2a), sba), ((c1b, c2b), sbb)):
                # Return tuples (id of the splited sb, splitRanks)
                # This allows to incrementaly update the splitRanks
                #  and then it wil be possible to perform the splits
                #  after the for loops
                (splitRanks_a, splitRanks_b) = splitNestedSbs(((c1a, c2a), sba), ((c1b, c2b), sbb))
                id2splitRanks[idsba] = sorted(list(set(id2splitRanks[idsba] + splitRanks_a)))
                id2splitRanks[idsbb] = sorted(list(set(id2splitRanks[idsbb] + splitRanks_b)))
                todo = todo | set([idsba, idsbb])

    # With 'todo' and (id of the splited sb, splitRanks)s, calculate the newSbs
    # and removed sbs
    # TODO update the set of sbs
    #for (idsb, ) in sbsInPairCompWithIds.iteritemsWithIds:
    #    if idsb in removedSbIds:
    #        pass
    #    else:
    #        # keep
    #for sb in newSbs:
    #    # add sb in sbsInPairCompWithIds
    # END TODO
    newSbs = []
    removedSbIds = set([])
    for idsb in todo:
        (c1, c2, _) = id2location(idsb)
        for newSb in splitSbBySplitRanks(id2sb(idsb), id2splitRanks[idsb]):
            newSbs.append((c1, c2, newSb))
        removedSbIds = removedSbIds | set([idsb])
    sbsInPairCompWithIds.removeIds(removedSbIds)
    intersTmp = set(id for (id, _, _) in sbsInPairCompWithIds.iterByOrderedIds()) & removedSbIds
    # DEBUG assertion
    assert len(intersTmp) == 0
    for (c1, c2, newSb) in newSbs:
        sbsInPairCompWithIds.addToLocation((c1, c2), newSb)

    # also remove sbs of only one hp!
    # TODO print in stderr the reason of the removal
    if removeSingleHpSbs:
        singleHpSbs = set([])
        for (idsb, _, sb) in sbsInPairCompWithIds.iterByOrderedIds():
            assert len(sb.la) > 0
            if len(sb.la) == 1:
                singleHpSbs.add(idsb)
        sbsInPairCompWithIds.removeIds(singleHpSbs)

    # update the small overlap dict ??
    #newO = {}
    #for (idsba, idsbbs) in O.iteritems():
    #    if idsba in noOrSmallOverlapIds:
    #        newO[idsba] = set([idsbb for idsbb in idsbbs if idsbb in noOrSmallOverlapIds])
    #O = newO

    # TODO, translate into a simpler Dict2D!!
    new_sbsInPairComp = myTools.Dict2d(list)
    for (idsb, (c1, c2), sb) in sbsInPairCompWithIds.iterByOrderedIds():
        new_sbsInPairComp[c1][c2].append(sb)
    return new_sbsInPairComp


def extractSbsInPairCompGenomesInTbs(g1_tb, g2_tb, ancGenes,
                                     filterType=FilterType.None,
                                     tandemGapMax=0,
                                     gapMax=None,
                                     distanceMetric='DPD',
                                     pThreshold=0.001,
                                     identifyBreakpointsWithinGaps=False,
                                     nonOverlappingSbs=False,
                                     overlapMax=0,
                                     consistentSwDType=True,
                                     minChromLength=0,
                                     nbHpsRecommendedGap=2,
                                     targetProbaRecommendedGap=0.01,
                                     validateImpossToCalc_mThreshold=3,
                                     multiProcess=False,
                                     verbose=False):

    # second level of verbosity
    #verbose2 = verbose if (len(g1) > 500 or len(g2) > 500) else False
    verbose2 = False
    N12s, N12_g = numberOfHomologies(g1_tb, g2_tb, verbose=verbose2)
    print >> sys.stderr, "pairwise comparison of genome 1 and genome 2 yields %s hps" % N12_g
    # compute the recommended gapMax parameter
    #########################################
    verbose2 = False
    (p_hpSign, p_hpSign_g, (sTBG1, sTBG1_g), (sTBG2, sTBG2_g)) =\
        myProbas.statsHpSign(g1_tb, g2_tb, verbose=verbose2)
    print >> sys.stderr, "genome 1 tb orientation proba = {+1:%.2f%%,-1:%.2f%%,None:%.2f%%} (stats are also calculated for each chromosome)" % (sTBG1_g[+1]*100, sTBG1_g[-1]*100, sTBG1_g[None]*100)
    print >> sys.stderr, "genome 2 tb orientation proba = {+1=%.2f%%,-1:%.2f%%,None:%.2f%%} (stats are also calculated for each chromosome)" % (sTBG2_g[+1]*100, sTBG2_g[-1]*100, sTBG2_g[None]*100)
    print >> sys.stderr, "hp sign proba in the 'global' mhp = {+1:%.2f%%,-1:%.2f%%,None:%.2f%%) (probabilities are also calculated for pairwise mhp)" % (p_hpSign_g[+1]*100, p_hpSign_g[-1]*100, p_hpSign_g[None]*100)
    N1_g = sum([len(g1_tb[c1]) for c1 in g1_tb])
    N2_g = sum([len(g2_tb[c2]) for c2 in g2_tb])
    #Build an average MHP
    #N50 is better than the mean of the chromosome lengths, since it is less sensitive to numerous and very small chromosome lengths
    #N1_N50 = myMaths.myStats.getValueNX(sorted([len(g1_tb[c1]) for c1 in g1_tb]),50) # calculate the N50
    #N2_N50 = myMaths.myStats.getValueNX(sorted([len(g2_tb[c2]) for c2 in g2_tb]),50) # calculate the N50
    #weightedAverage is even better than N50 since it returns a more stable value, not a length of a chromosome of the karyotype, it better reflects the global distribution

    #Waring: if the genome is badly assembled and contain a lot of small contigs, this averaging is not relevant and the recommended gapMax won't be relevant neither
    if len(g1_tb) > 50 or len(g2_tb) > 50:
        print >> sys.stderr, "Warning: one of the two genome seems badly assembled, this may mislead the recommended gapMax calculation"
    N1_weightedAverage = int(myMaths.myStats.getWeightedAverage([len(g1_tb[c1]) for c1 in g1_tb]))
    N2_weightedAverage = int(myMaths.myStats.getWeightedAverage([len(g2_tb[c2]) for c2 in g2_tb]))
    density = float(N12_g)/(N1_g*N2_g)
    # conservation of the density
    N12_weightedAverage = int(density*N1_weightedAverage*N2_weightedAverage)
    gap = recommendedGap(nbHpsRecommendedGap, targetProbaRecommendedGap, N12_weightedAverage, N1_weightedAverage, N2_weightedAverage, p_hpSign=p_hpSign_g, verbose=verbose)
    print >> sys.stderr, "recommended gapMax = %s tbs" % gap
    if gapMax is None:
        gapMax = gap
    print >> sys.stderr, "used gapMax = %s tbs" % gapMax

    # step 2 and 3 : build the MHP and extract putative sbs as diagonals
    #################################################################################
    # extract sbs in the tb base
    sbsInPairComp = myTools.Dict2d(list)
    if multiProcess and (len(g1_tb.keys()) > 1 or len(g2_tb.keys()) > 1):
        # if the multiprocess option is True and if there is more than one pairwise comparison of chromosomes
        tasks = [(c1, c2, g1_tb[c1], g2_tb[c2]) for (c1, c2) in itertools.product([c1 for c1 in g1_tb], [c2 for c2 in g2_tb])]
        resGen = myMultiprocess.multiprocessTasks(extractSbsInPairCompChrWrapper, tasks,
                                                  gapMax=gapMax,
                                                  distanceMetric=distanceMetric,
                                                  consistentSwDType=consistentSwDType,
                                                  verbose=False)
        for ((c1, c2), listOfSbs) in resGen:
            if len(listOfSbs) > 0:
                sbsInPairComp[c1][c2] = listOfSbs
    else:
        nbPairwiseComparisons = len(g1_tb.keys())*len(g2_tb.keys())
        listOfPercentage = range(0, 101, 5)[1:]
        print >> sys.stderr, "synteny block extraction",
        for (i, (c1, c2)) in enumerate(itertools.product(g1_tb.keys(), g2_tb.keys())):
            listOfSbs = extractSbsInPairCompChr(g1_tb[c1], g2_tb[c2],
                                                gapMax=gapMax,
                                                distanceMetric=distanceMetric,
                                                consistentSwDType=consistentSwDType,
                                                verbose=verbose)
            if len(listOfSbs) > 0:
                sbsInPairComp[c1][c2] = listOfSbs
            progress = int(float(i*100)/nbPairwiseComparisons)
            if progress in listOfPercentage:
                print >> sys.stderr, "%s" % progress + "%",
                listOfPercentage.remove(progress)
        # new line in the print
        print >> sys.stderr, ""

    # setp 4 : statistical validation of putative sbs
    ##################################################
    #DEBUG
    #print >> sys.stderr, "sTBG1['Y']=%s" % sTBG1['Y']
    #print >> sys.stderr, "sTBG2['Y']=%s" % sTBG2['Y']
    #print >> sys.stderr, "p_hpSign[('Y','Y')]=%s" % p_hpSign[('Y','Y')]
    # 'SV' stands for statistical validation
    sbsInPairComp = statisticalValidation(sbsInPairComp, g1_tb, g2_tb, N12s, p_hpSign,
                                          pThreshold=pThreshold,
                                          NbOfHomologiesThreshold=50,
                                          validateImpossToCalc_mThreshold=validateImpossToCalc_mThreshold,
                                          verbose=verbose)

    cptLoopIter = 0
    # initialise la condition d'arrêt de la boucle
    atLeastOneNonOverlappingMerge = False
    while True:

        if identifyBreakpointsWithinGaps:
            # Do this task until no more change (stability!)
            nbSbs = len([sb for (_, sb) in sbsInPairComp.iteritems2d()])
            if cptLoopIter == 0:
                print >> sys.stderr, "Nb sbs before identifying breakpoints within gaps = %s" % nbSbs
            while True:
                nbSbsOld = len(list(sbsInPairComp.iteritems2d()))
                sbsInPairComp = fIdentifyBreakpointsWithinGaps(sbsInPairComp)
                nbSbsNew = len(list(sbsInPairComp.iteritems2d()))
                if nbSbsNew == nbSbsOld:
                    break
            nbSbs = len([sb for (_, sb) in sbsInPairComp.iteritems2d()])
            if cptLoopIter == 0:
                print >> sys.stderr, "Nb sbs after identifying breakpoints within gaps = %s" % nbSbs

        if nonOverlappingSbs:
            nbSbs = len([sb for (_, sb) in sbsInPairComp.iteritems2d()])
            if cptLoopIter == 0:
                print >> sys.stderr, "Nb sbs before overlap-truncation-filtering = %s" % nbSbs
            # TODO, compute the variation of the coverage before and after this step
            sbsInPairComp = filterOverlappingSbs(sbsInPairComp, overlapMax=overlapMax, verbose=False)
            # DEBUG, verify that the filtered sbs are not overlapping
            (_, _, N, O, _, _) = buildConflictGraph(sbsInPairComp, overlapMax=0)
            assert len(N) == 0
            assert len(O) == 0
            nbSbs = len([sb for (cc, sb) in sbsInPairComp.iteritems2d()])
            if cptLoopIter == 0:
                print >> sys.stderr, "Nb sbs after overlap-truncation-filtering = %s" % nbSbs

            # mergeNonOverlappingSbs
            for (c1, c2) in sbsInPairComp.keys2d():
                # BM: Before Merge
                nbSbsBM = len([sb for (cc, sb) in sbsInPairComp.iteritems2d()])
                sbsInPairComp[c1][c2] = mergeSbs(sbsInPairComp[c1][c2], gapMax, g2_tb[c2], verbose=False)
                # AM: After Merge
                nbSbsAM = len([sb for (cc, sb) in sbsInPairComp.iteritems2d()])
                assert nbSbsAM <= nbSbsBM
                atLeastOneNonOverlappingMerge = atLeastOneNonOverlappingMerge or (True if nbSbsAM != nbSbsBM else False)
            nbSbs = len([sb for (cc, sb) in sbsInPairComp.iteritems2d()])
            if cptLoopIter == 0:
                print >> sys.stderr, "Nb sbs after merging non-overlapping sbs = %s" % nbSbs

        cptLoopIter += 1

        if atLeastOneNonOverlappingMerge and cptLoopIter <= 30:
            continue
        else:
            if cptLoopIter  > 1:
                # the first merge non-overlapping sbs was usefull, thus other
                # iterations of
                # IdentifyBreakpointsWithinGaps|nonOverlappingSbs|mergeNonOverlappingSbs
                # were performed.
                print >> sys.stderr, "IdentifyBreakpointsWithinGaps|nonOverlappingSbs|mergeNonOverlappingSbs was repeated %s times" % (cptLoopIter - 1)

                # Always finish by identifyBreakpointsWithinGaps followed by
                # nonOverlappingSbs without a merge
                if identifyBreakpointsWithinGaps:
                    # Do this task until no more change (stability!)
                    nbSbs = len([sb for (cc, sb) in sbsInPairComp.iteritems2d()])
                    print >> sys.stderr, "Nb sbs before identifying breakpoints within gaps = %s" % nbSbs
                    while True:
                        nbSbsOld = len(list(sbsInPairComp.iteritems2d()))
                        sbsInPairComp = fIdentifyBreakpointsWithinGaps(sbsInPairComp)
                        nbSbsNew = len(list(sbsInPairComp.iteritems2d()))
                        if nbSbsNew == nbSbsOld:
                            break
                    nbSbs = len([sb for (cc, sb) in sbsInPairComp.iteritems2d()])
                    print >> sys.stderr, "Nb sbs after identifying breakpoints within gaps = %s" % nbSbs
                if nonOverlappingSbs:
                    nbSbs = len([sb for (cc, sb) in sbsInPairComp.iteritems2d()])
                    print >> sys.stderr, "Nb sbs before overlap-truncation-filtering = %s" % nbSbs
                    # TODO, compute the variation of the coverage before and after this step
                    sbsInPairComp = filterOverlappingSbs(sbsInPairComp, overlapMax=overlapMax, verbose=False)
                    # DEBUG, verify that the filtered sbs are not overlapping
                    (_, _, N, O, _, _) = buildConflictGraph(sbsInPairComp, overlapMax=0)
                    assert len(N) == 0
                    assert len(O) == 0
                    nbSbs = len([sb for (cc, sb) in sbsInPairComp.iteritems2d()])
                    print >> sys.stderr, "Nb sbs after overlap-truncation-filtering = %s" % nbSbs
                # No merge non-overlapping sbs
            # leave the while loop
            break

    return sbsInPairComp

# Complete procedure to compute synteny blocks (sbs) between 2 genomes, using homology relationships contained in ancGenes
# Inputs:
#       g1, g2 : myGenomes.Genome (or dict : {...c:[..., (g,s), ...]}) of the two compared species
#       ancGenes : myGenomes.Genome define the homology relationships (usually the ancGenes of the LCA of the two compared species)
#       consistentSwDType :     True  => gene transcription orientations must be consistent with diagonal types ('/' or '\')
#                               False => Do not take care of gene orientation
#       gapMax : distance (in tandemblocks) allowed between two diagonals to be merged together. This allows to go over wrong  annotations. But it can also introduce some Fp
#       minChromLength is the minimal length of the chromosome to be considered
#       distanceMetric : the distance metric (either MD,DPD,ED or CD)
#       pThreshold : the probability threshold under which the sb is considered significant
#       FilterType
#               InCommonAncestor : all genes herited from a gene of the LCA are kept during the filtering of extant genomes
#               InBothGenomes : only 'anchor genes' (genes present in both genomes) are kept
#               None : genomes are not filtered
# Outputs:
#       synteny blocks are yielded as (strand,(c1,l1),(c2,l2),la)
#       c1 : chromosomes of genome 'g1' where there is a synteny block corresponding to the diagonal
#       l1 : the synteny block on genome 'g1'
#       l1 = [..., (i,s), ...] with 'i' the index of the gene and 's' its strand in the genome 1 without species specific genes
#       la = [..., (ancGeneName, ancGeneOrientation, tb width, tb height), ...]]
#           tb -> tandem block associated with the ancGene in this synteny block
#           width : length in tandem duplicates on the 1st genome
#           height: length in tandem duplicates on ths 2nd genome
#########################################################################################################################
@myTools.tictac
@myTools.verbose
def extractSbsInPairCompGenomes(g1, g2, ancGenes,
                                filterType=FilterType.None,
                                tandemGapMax=0,
                                gapMax=None,
                                distanceMetric='DPD',
                                pThreshold=0.001,
                                identifyBreakpointsWithinGaps=False,
                                nonOverlappingSbs=False,
                                overlapMax=0,
                                consistentSwDType=True,
                                minChromLength=0,
                                nbHpsRecommendedGap=2,
                                targetProbaRecommendedGap=0.01,
                                validateImpossToCalc_mThreshold=3,
                                multiProcess=False,
                                verbose=False):
    # TODO, raise a true warning message
    if nonOverlappingSbs is False:
        if overlapMax > 0:
            print >> sys.stderr, "Warning: the maxAllowedGap specified is not used since the nonOverlappingSbs is False"

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
    g1_aID = myMapping.labelWithAncGeneID(g1, ancGenes)
    g2_aID = myMapping.labelWithAncGeneID(g2, ancGenes)
    # genes that are not in ancGene have a aID=None
    print >> sys.stderr, "genome 1 initially contains %s genes" % sum([len(g1[c1]) for c1 in g1])
    print >> sys.stderr, "genome 2 initially contains %s genes" % sum([len(g2[c2]) for c2 in g2])
    # Must be applied on the two genomes, because of the mode inBothGenomes (InCommonAncestor => not only anchor genes are kept but all genes herited from a gene of the LCA)
    #mfilt2origin1 -> mGf2Go1
    ((g1_aID, mGf2Go1, (nCL1, nGL1)), (g2_aID, mGf2Go2, (nCL2, nGL2))) =\
        filter2D(g1_aID, g2_aID, filterType, minChromLength)
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
    print >> sys.stderr, "genome 1 contains %s tandem duplicated genes (initial gene excluded)" % nGTD1
    print >> sys.stderr, "genome 2 contains %s tandem duplicated genes (initial gene excluded)" % nGTD2
    nDD1 = myMapping.nbDup(g1_tb)[0]
    nDD2 = myMapping.nbDup(g2_tb)[0]
    print >> sys.stderr, "genome 1 contains %s dispersed duplicated tbs (initial tb excluded)" % nDD1
    print >> sys.stderr, "genome 2 contains %s dispersed duplicated tbs (initial tb excluded)" % nDD2
    assert nDD1 + nGTD1 == nGD1
    assert nDD2 + nGTD2 == nGD2

    sbsInPairComp = extractSbsInPairCompGenomesInTbs(g1_tb, g2_tb, ancGenes,
                                                     gapMax=gapMax,
                                                     distanceMetric=distanceMetric,
                                                     pThreshold=pThreshold,
                                                     identifyBreakpointsWithinGaps=identifyBreakpointsWithinGaps,
                                                     nonOverlappingSbs=nonOverlappingSbs,
                                                     overlapMax=overlapMax,
                                                     consistentSwDType=consistentSwDType,
                                                     nbHpsRecommendedGap=nbHpsRecommendedGap,
                                                     targetProbaRecommendedGap=targetProbaRecommendedGap,
                                                     validateImpossToCalc_mThreshold=validateImpossToCalc_mThreshold,
                                                     multiProcess=multiProcess,
                                                     verbose=verbose)

    # format output
    for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
        # translate from the tb base to gene base
        mtb2gc1 = mtb2g1[c1]
        mtb2gc2 = mtb2g2[c2]
        assert mtb2gc1.__class__.__name__ == mtb2gc2.__class__.__name__ == 'Mapping'
        #mfilt2origin1 -> mGf2Go1
        #mGf2GoC1 = mGf2Go1[c1]
        #mGf2GoC2 = mGf2Go2[c2]
        la = []
        l1 = []
        l2 = []
        # convert diags from tb to genes before yielding them
        # each hp corresponds to an ancestral gene
        for (indx_HP, aGene) in enumerate(sb.la):
            indx_tb_g1 = sb.l1[indx_HP]
            l1.append([gIndxs for gIndxs in mtb2gc1.new[indx_tb_g1]])
            indx_tb_g2 = sb.l2[indx_HP]
            l2.append([gIdxs for gIdxs in mtb2gc2.new[indx_tb_g2]])
            la.append((ancGenes.lstGenes[None][aGene[0]].names[0], aGene[1], aGene[2]))
        # modify the current object SyntenyBlock, within sbsInPairComp
        sb.l1 = l1
        sb.l2 = l2
        sb.la = la

    return sbsInPairComp


#################
# Print function
#################
# Print diags into [nbHps, c1, c2, la, l1, l2)]
# With:
# c1 and c2 the chromosomes of the pairwise comparison in which the synteny block has been found
# la = [..., (aGName, strand, dist), ...]
# l1 = [..., [g1, ...gN], ...] with [g1, ..., gN] the child tb of the corresponding aG
def printSbsFile(sbsInPairComp, genome1, genome2, sortByDecrLengths=True):
    print >> sys.stderr, "Print synteny blocks"

    def foo(genomeX, cX, lX, idxHp):
        gXs = [genomeX.lstGenes[cX][gIdx].names[0] for gIdx in lX[idxHp]]
        gXs = ' '.join(gXs)
        sXs = [genomeX.lstGenes[cX][gIdx].strand for gIdx in lX[idxHp]]
        for (i, s) in enumerate(sXs):
            if s == +1:
                sXs[i] = '+'
            elif s == -1:
                sXs[i] = '-'
            else:
                raise
        sXs = ''.join(sXs)
        return (gXs, sXs)

    # listOfSbs = [ ..., ((c1, c2), sb), ...]
    listOfSbs = sbsInPairComp.intoList()
    if sortByDecrLengths:
        listOfSbs.sort(key=lambda x: len(x[1].la), reverse=True)

    statsSbs = []
    for (idSb, ((c1, c2), sb)) in enumerate(listOfSbs):
        assert len(sb.la) == len(sb.l1) == len(sb.l2), "len(l1)=%s, len(l2)=%s, len(la)=%s\nl1=%s\nl2=%s\nla=%s" % (len(sb.l1), len(sb.l2), len(sb.la), sb.l1, sb.l2, sb.la)
        nbHps = len(sb.la)
        statsSbs.append(nbHps)
        for (idxHp, (aGname, aGstrand, dist)) in enumerate(sb.la):
            (g1s, s1s) = foo(genome1, c1, sb.l1, idxHp)
            (g2s, s2s) = foo(genome2, c2, sb.l2, idxHp)
            print myFile.myTSV.printLine([idSb, aGname, aGstrand, dist, c1, c2, s1s, s2s, g1s, g2s])

    print >> sys.stderr, "Distribution of the lengths of synteny blocks\t", myMaths.myStats.syntheticTxtSummary(statsSbs)


#################################
# Parser function
#################################
# adding genome1 and genome2 allow to recover l1 and l2, the lists of homologous
# tbs on genome1 and genome2.
def parseSbsFile(fileName, genome1=None, genome2=None):

    def foo(gXs, sXs, genome, chr):
        if genome is None:
            return ([None], [None])
        gXs = gXs.split(' ')
        nbTandemDup = len(sXs)
        assert len(sXs) == len(gXs) == nbTandemDup, "sXs=%s, gXs=%s and nbTandemDup=%s" % (sXs, gXs, nbTandemDup)
        gXs = gXs.split(' ')
        for i in range(nbTandemDup):
            gN = gXs[i]
            gXs[i] = genome.getPosition([gN])
            s = sXs[i]
            if s == '+':
                sXs[i] = +1
            elif s == '-':
                sXs[i] = -1
            else:
                raise
            assert genome.lstGenes[chr][gXs[i]].strand == sXs[i]
        return (gXs, sXs)

    sbsReader =\
        myFile.myTSV.readTabular(fileName,
                                 [int, str, str, int, str, str, str, str, str, str],
                                 delim='\t')
    sbsInPairComp = myTools.Dict2d(list)
    pVal = 0
    idSb_old = 'foo'
    c1_old = 'fooC1'
    c2_old = 'fooC2'
    for (idSb, aGname, aGstrand, dist, c1, c2, s1s, s2s, g1s, g2s) in sbsReader:
        aGstrand = int(dist) if aGstrand != 'None' else None
        dist = int(dist) if dist != 'None' else None
        (g1s, s1s) = foo(g1s, s1s, genome1, c1)
        (g2s, s2s) = foo(g2s, s2s, genome2, c2)
        if idSb != idSb_old:
            if idSb_old == 'foo':
                # first idSb
                l1 = list(zip(g1s, s1s))
                l2 = list(zip(g2s, s2s))
                la = list((aGname, aGstrand, dist))
            else:
                # record the former synteny block that has been parsed
                # TODO find the diagType
                if l1[0] < l1[-1]:
                    if l2[0] < l2[-1]:
                        diagType = '/'
                    elif l2[2] > l2[-1]:
                        diagType = '\\'
                    else:
                        # horizontal synteny block
                        diagType = None
                else:
                    # vertical synteny block
                    diagType = None
                #if c1_old == '10' and c2_old == '14':
                #    print >> sys.stderr, 'Hello !', idSb_old
                #    raw_input()
                sbsInPairComp[c1_old][c2_old].append(SyntenyBlock(Diagonal(diagType, l1, l2, la), pVal))
                # start a new sb
                l1 = list(zip(g1s, s1s))
                l2 = list(zip(g2s, s2s))
                la = list((aGname, aGstrand, dist))
            idSb_old = idSb
            c1_old = c1
            c2_old = c2
        else:
            l1.append(zip(g1s, s1s))
            l2.append(zip(g2s, s2s))
            la.append((aGname, aGstrand, dist))
    # last idSb
    # TODO find the diagType
    if l1[0] < l1[-1]:
        if l2[0] < l2[-1]:
            diagType = '/'
        elif l2[2] > l2[-1]:
            diagType = '\\'
        else:
            # horizontal synteny block
            diagType = None
    else:
        # vertical synteny block
        diagType = None
    sbsInPairComp[c1][c2].append(SyntenyBlock(Diagonal(diagType, l1, l2, la), pVal))
    return sbsInPairComp


###########################################
# Post processing of synteny blocks (sbs)
###########################################

# Build a genome with sbs as contigs and only conserve one ancestral gene by hp
def buildGenomeFromSbs(sbsInPairComp, sbLengthThreshold):
    cptSb=0
    ancSbGenome={}
    for ((c1, c2), (sb, pVal)) in sbsInPairComp.iteritems2d():
        old_ancGene=None
        for ihp,(ancGene, ancStrand, dist) in enumerate(sb.la):
            #assert ancGene != old_ancGene, "%s,%s" % (ancGene, old_ancGene) # Happens when there are paralogous hps in the same diagonal
            if ancGene == old_ancGene:
                continue # This may create diags of size 1 that will need to be removed
            ancSbGenome[cptSb].append((ancGene, ancStrand, dist)) # ancGene synonymous of hp-name
            old_ancGene = ancGene

        if len(ancSbGenome[cptSb]) <= 1: # In every cases sbs of length 1 are removed
            del ancSbGenome[cptSb]
            continue
        elif sbLengthThreshold != None and len(ancSbGenome[cptSb]) <= sbLengthThreshold:
            del ancSbGenome[cptSb]
            continue
        cptSb = cptSb+1
    return ancSbGenome

@myTools.verbose
def computeAncestralCoverageBySbs(g1_tb, g2_tb, ancSbGenome, verbose = False):
    ancGenesInSbs=set([])
    for ancSb in ancSbGenome.values():
        for (ancGene,_,_) in ancSb:
            ancGenesInSbs.add(ancGene)
    print >> sys.stderr, "Nb of ancGenes in synteny blocks (each ancGene can appear at most once)", len(ancGenesInSbs)
    ancGenesInG1Tb = set([])
    ancGenesInG2Tb = set([])
    for tb1 in [tb1 for c1 in g1_tb for tb1 in g1_tb[c1]]:
        (ancGene,s) = tb1
        ancGenesInG1Tb.add(ancGene)
    for tb2 in [tb2 for c2 in g2_tb for tb2 in g2_tb[c2]]:
        (ancGene,s) = tb2
        ancGenesInG2Tb.add(ancGene)
    #{ancGenes that are present in G1 rewritten in tbs and in G2 rewritten in tbs}
    ancGenesInG1TbAndInG2Tb = ancGenesInG1Tb & ancGenesInG2Tb # intersection
    print >> sys.stderr, "coverage LCA_S1-S2 = cardinal({ancGenes in Sb}) / cardinal({ancGenes that are present in G1 rewritten in tbs and in G2 rewritten in tbs})\ncoverage LCA_S1-S2 = ", float(len(ancGenesInSbs)) / len(ancGenesInG1TbAndInG2Tb)
    coverage = float(len(ancGenesInSbs)) / len(ancGenesInG1TbAndInG2Tb)
    return coverage