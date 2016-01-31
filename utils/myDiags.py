#!/usr/bin/python
# -*- coding: utf-8 -*-


#  ____  _           _ ____  _
# |  _ \| |__  _   _| |  _ \(_) __ _  __ _
# | |_) | '_ \| | | | | | | | |/ _` |/ _` |
# |  __/| | | | |_| | | |_| | | (_| | (_| |
# |_|   |_| |_|\__, |_|____/|_|\__,_|\__, |
#              |___/                 |___/
#
# PhylDiag version 2.0 (6/11/2015)
# LibsDyogen version 1.0 (6/11/2015)
# python v2.7 at least is needed
# Copyright © 2015 IBENS/Dyogen : Joseph LUCAS, Lucas TITTMANN and Hugues ROEST CROLLIUS
# mail : jlucas@ens.fr
# Licences GLP v3 and CeCILL v2

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
# cs: conserved segment
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
import myLightGenomes

import time
import multiprocessing.dummy as multiprocessing  # for using threads
import multiprocessing  # for using processes
import enum

if sys.platform != 'win32':
    import extractDiags

from utils.myLightGenomes import OGene

FilterType = enum.Enum('InFamilies', 'InBothGenomes', 'None')

defaultArgsPhylDiag = \
    [("filterType", str, 'InBothGenomes'),
     ("tandemGapMax", int, 10),
     ("gapMax", str, 5),
     ("distinguishMonoGenicDiags", bool, True),
     ('distanceMetric', str, 'CD'),
     ('pThreshold', str, 'None'),
     ('gapMaxMicroInv', str, '1'),
     ('identifyMonoGenicInvs', bool, True),
     ('identifyMicroRearrangements', bool, True),
     ('truncationMax', str, 10),
     ("minChromLength", int, 2),
     ("sameStrand", bool, True),
     ('nbHpsRecommendedGap', int, 2), ('targetProbaRecommendedGap', float, 0.01),
     ('validateImpossToCalc_mThreshold', int, 3),
     # The multiprocess does not seem to work well for most of the data. It work well only for some data
     # with ~ 800 contigs
     ('optimisation', str, 'cython'),
     ('verbose', bool, False)]

defaultKwargsPhylDiagDict = dict((n, (t, v)) for (n, t, v) in defaultArgsPhylDiag)

def defaultKwargsPhylDiag(arguments=None):
    res = dict((n, t(v)) for (n, t, v) in defaultArgsPhylDiag)
    if arguments is None:
        return res
    else:
        res.update(dict((k, v) for (k, v) in arguments.iteritems() if k in defaultKwargsPhylDiagDict.keys()))
        for (argN, tpe) in [('gapMax', int), ('truncationMax', int), ('gapMaxMicroInv', int), ('pThreshold', float), ('optimisation', str)]:
            if arguments[argN] == 'None':
                res[argN] = None
            else:
                try:
                    res[argN] = tpe(arguments[argN])
                except:
                    raise TypeError('%s should be either an int or None, not %s, a %s' % (argN, arguments[argN], type(arguments[argN])))
        res['filterType'] = FilterType[list(FilterType._keys).index(arguments["filterType"])]
    return res

# TODO, write a diagonal class that is made for tbs containing genes
class Diagonal():
    def __init__(self, *args):
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
        assert rankGenome in [1, 2]
        if rankGenome == 1:
            l = self.l1
        else:
            l = self.l2
        if isinstance(l, list) and isinstance(l[0], list):
            res = min([idx for tb in l for idx in tb])
        else:
            res = min(l)
        return res

    def maxOnG(self, rankGenome):
        assert rankGenome in [1, 2]
        if rankGenome == 1:
            l = self.l1
        else:
            l = self.l2
        if isinstance(l, list) and isinstance(l[0], list):
            res = max([idx for tb in l for idx in tb])
        else:
            res = max(l)
        return res

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

    def truncate(self, range1, range2):
        # range1 = (begining, end) on chromosome 1
        # range2 = (begining, end) on chromosome 2
        assert len(range1) == 2 and len(range2) == 2
        new_l1 = []
        new_l2 = []
        new_la = []
        assert len(self.l1) == len(self.l2) == len(self.la)
        # FIXME if the sb is in tbs
        if all([isinstance(item, list) for item in self.l1]) and\
            all([isinstance(item, list) for item in self.l2]):
            for (idxH, aG) in enumerate(self.la):
                i1 = self.l1[idxH][0]
                i2 = self.l2[idxH][0]
                if (range1[0] <= i1 and i1 < range1[1]) and (range2[0] <= i2 and i2 < range2[1]):
                    new_l1.append([i1])
                    new_l2.append([i2])
                    new_la.append(aG)
        else:
            assert all([isinstance(item, int) for item in self.l1]) and\
                all([isinstance(item, int) for item in self.l2]), str(self.l2)
            for (idxH, aG) in enumerate(self.la):
                i1 = self.l1[idxH]
                i2 = self.l2[idxH]
                if (range1[0] <= i1 and i1 < range1[1]) and (range2[0] <= i2 and i2 < range2[1]):
                    new_l1.append(i1)
                    new_l2.append(i2)
                    new_la.append(aG)
                # FIXME, self.lX being lists is not a normal case of the Diagonal
                # object. this condition is for current convenience purpose. It
                # should be replaced y a more general Diagonal object that allow to
                # use self.lX being lists. minOnGX(), beg() and end() should be
                # changed in consequence.
                # if isinstance(tb1, list) and isinstance(tb2, list):
                #    new_tb1 = [g1Idx for g1Idx in tb1 if (range1[0] <= g1Idx and g1Idx <= range1[1])]
                #    new_tb2 = [g2Idx for g2Idx in tb2 if (range2[0] <= g2Idx and g2Idx <= range2[1])]
                #    if len(new_tb1) > 0 and len(new_tb2) > 0:
                #        new_l1.append(new_tb1)
                #        new_l2.append(new_tb2)
                #        new_la.append(aG)
        self.l1 = new_l1
        self.l2 = new_l2
        self.la = new_la

    def __repr__(self):
        return "\ndiagType=%s\nl1=%s\nl2=%s\nla=%s\n" % (self.dt, self.l1, self.l2, self.la)

# Diagonal Pseudo Distance
def DPD((x0, y0), (x1, y1)):
    return 2 * max(abs(x1 - x0), abs(y1 - y0)) - min(abs(x1 - x0), abs(y1 - y0))

# Chebyshev Distance
def CD((x0, y0), (x1, y1)):
    return max(abs(x1 - x0), abs(y1 - y0))

# Manhattan Distance
def MD((x0, y0), (x1, y1)):
    return abs(x1 - x0) + abs(y1 - y0)

# Euclidean Distance
def ED((x0, y0), (x1, y1)):
    return round(math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2))

# not really a distance, but may be used af it is one in several cases
def minimumOfDistancesIn1DGenomes((x1, y1), (x2, y2)):
    return min(abs(x1 - x2), abs(y1 - y2))

# The frame of all the distances used during the merging process
def framed(f):
    def distance((x0, y0), (x1, y1), diagType):
        if diagType == '/':
            # Only consider the top right part of the distance matrix, (x1 >= x0 and y1 >= y0)
            # in i-adhore "...>0", here we allow '=' to take care of close dispered paralogies that can be understood as tandem duplicates
            if (x1-x0) >= 0 and (y1 - y0) >= 0:
                res = f((x0, y0), (x1, y1))
            else:
                res = sys.maxint
        elif diagType == '\\':
            # Only consider the bottom right part of the distance matrix, (x1>=x0 and y1<=y0)
            # in i-adhore "...<0", here we allow '=', idem than before
            if (x1-x0) >= 0 and (y1 - y0) <= 0:
                res = f((x0, y0), (x1, y1))
            else:
                res = sys.maxint
        else:
            # diagType == None, diag is composed of a single gene
            res = f((x0, y0), (x1, y1))
        return res
    return distance

def intervalsOverlap((x1, x2), (y1, y2)):
    """
    check if two 1D intervals overlap
    """
    assert x1 <= x2
    assert y1 <= y2
    # If there is an overlap it means there exists some number C which is in both ranges, i.e.
    # x1 <= C <= x2
    # y1 <= C <= y2
    # Now, if we are allowed to assume that the ranges are well-formed (so that x1 <= x2 and y1 <= y2) then it is sufficient to test
    # x1 <= y2 & y1 <= x2
    # other solution: max(x1, y1) <= min(x2, y2)
    return x1 <= y2 and y1 <= x2

def boundingBoxesOfDiagsOverlap(diagA, diagB):
    """
    check if two bounding boxes of diagonals overlap
    """
    assert isinstance(diagA, Diagonal)
    assert isinstance(diagB, Diagonal)
    res = False
    if intervalsOverlap((diagA.minOnG(1), diagA.maxOnG(1)), (diagB.minOnG(1), diagB.maxOnG(1))):
        res = True
    elif intervalsOverlap((diagA.minOnG(2), diagA.maxOnG(2)), (diagB.minOnG(2), diagB.maxOnG(2))):
        res = True
    return res

def distanceBetweenBoundingBoxesOfDiags(diagA, diagB, distance=CD):
    """
    compute the distance between two bounding boxes of diagonals
    """
    assert isinstance(diagA, Diagonal)
    assert isinstance(diagB, Diagonal)

    if boundingBoxesOfDiagsOverlap(diagA, diagB):
        res = 0
        # pA = None
        # pB = None
    else:
        # invertedOn1 = False
        minA1 = diagA.minOnG(1)
        minB1 = diagB.minOnG(1)
        if minB1 < minA1:
            # invertedOn1 = True
            (diagA, diagB) = (diagB, diagA)
            minA1 = diagA.minOnG(1)
            minB1 = diagB.minOnG(1)
        assert minA1 < minB1
        maxA1 = diagA.maxOnG(1)
        # maxB1 = diagB.maxOnG(1)
        minA2 = diagA.minOnG(2)
        maxA2 = diagA.maxOnG(2)
        minB2 = diagB.minOnG(2)
        maxB2 = diagB.maxOnG(2)
        # since no overlap
        assert maxA1 < minB1
        if maxA2 < minB2:
            pA = (maxA1, maxA2)
            pB = (minB1, minB2)
            res = distance(pA, pB)
        elif maxB2 < minA2:
            pA = (maxA1, minA2)
            pB = (minB1, maxB2)
            res = distance(pA, pB)
        else:
            raise ValueError
        # if invertedOn1:
        #     (pA, pB) = (pB, pA)
    return res

# TODO: use it in mergeDiags()
# def distanceBetweenDiags(diagA, diagB, distance=CD):
#     assert isinstance(diagA, Diagonal)
#     assert isinstance(diagB, Diagonal)
#     inverted = False
#     if diagB.minOnG(1) < diagA.minOnG(1):
#         inverted = True
#         (diagA, diagB) = (diagB, diagA)
#     if diagonalsOverlap(diagA, diagB):
#         res = 0
#     else:
#         d1 = distance(diagA.beg(), diagB.beg())
#         d2 = distance(diagA.beg(), diagB.end())
#     # ensure that diagA.beg() is on the left of diagB.beg()
#     if diagB.beg()[0] <= diagA.end()[0] + gapMax + 1:
#         # Check if diagTypes are compatible
#         if diagA.dt == diagB.dt or diagA.dt == None or diagB.dt == None:
#             # take the known diagType if it is known in at least one of the 2 diags
#             dT = diagA.dt if diagA.dt != None else diagB.dt
#             # Check the distance
#             if distance(diagA.end(), diagB.beg(), dT) == currGap+1:
#                 fusionableDiags.append(diagB)
#             else:
#                 impossibleToMergeDiags.append(diagB)
#                 continue
#         else:
#             impossibleToMergeDiags.append(diagB)
#             continue
#     else:
#         # diagA.end()[0] + currGap < diagB.beg()[0]
#         # stop the loop (remember that the diags are sorted!)
#         # Impossible to merge next diags to diagA
#         continueToSearch = False
#         impossibleToMergeDiags.append(diagB)

class SyntenyBlock(Diagonal):
    def __init__(self, *args):
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

# Merge diagonals if they are separated by a gap less long than 'gapMax' relatively to a distance metric
# inputs:
# listOfDiags : a list containing elements as (l1, l2, la)
#       l1 = [ ..., [i11, i12...], ....] with i1x < i1y for all x<y
#               i1 : gene index of the corresponding tb in the genome1
#       l2 : [..., [...,ix2...], ...]
#               ix2 : gene index in the genome2
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
def mergeDiags(listOfDiags, gapMax, gc2, distanceMetric = 'CD', verbose = False):
    assert gapMax>=0
    # assert that the diag elements in listOfDiags are either Diagonal objects or SyntenyBlock objects
    assert all(diag.__class__.__name__ == 'Diagonal' for diag in listOfDiags) or\
        all(diag.__class__.__name__ == 'SyntenyBlock' for diag in listOfDiags), listOfDiags[0]

    if distanceMetric == 'DPD':
        print >> sys.stderr, "Use Diagonal Pseudo Distance to merge diagonals with a gap up to %s elements" % gapMax
        distance = framed(DPD)
    elif distanceMetric == 'CD':
        print >> sys.stderr, "Use Chebyshev Distance to merge diagonals with a gap up to %s elements" % gapMax
        distance = framed(CD)
    elif distanceMetric == 'MD':
        print >> sys.stderr, "Use Manhattan Distance to merge diagonals with a gap up to %s elements" % gapMax
        distance = framed(MD)
    elif distanceMetric == 'ED':
        print >> sys.stderr, "Use Euclidean Distance to merge diagonals with a gap up to %s elements" % gapMax
        distance = framed(ED)
    else:
        raise ValueError('Must use a distance either DPD (Diagonal PSeudo Distance) or MD (Manhattan Distance), Euclidean Distance (ED) or Chebyshev Distance (CD)')

    print >> sys.stderr, "Nb Diags before DiagMerger = ", len(listOfDiags)
    diagGen = []
    listOfFinishedDiags = []
    nbFusion = 0
    # Add the distance in la : [(ag1,ags1,dist1=0), (ag2,ags2,dist2=1), (ag2,ags3,dist3=3), ...]
    # means that between ag1 and ag2 there is no 2Dgap but between ag2 and ag3
    # there is a gap of 1. The gap depends on the distance metric chosen.
    # This gap gives us an information on the relevance of this adjacency in the
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
                if diagB.beg()[0] < diagA.end()[0]:  # in i-adhore diagB.beg()[0] <= diagA.end()[0]
                    impossibleToMergeDiags.append(diagB)
                    continue
                elif diagB.beg()[0] <= diagA.end()[0] + currGap+1:
                    # diagA.end()[0] < diagB.beg()[0] <= diagA.end()[0] + currGap
                    # Check if diagTypes are compatible
                    if diagA.dt == diagB.dt or diagA.dt == None or diagB.dt == None:
                        # take the known diagType if it is known in at least one of the 2 diags
                        dT = diagA.dt if diagA.dt != None else diagB.dt
                        # Check the distance
                        if distance(diagA.end(), diagB.beg(), dT) == currGap+1:
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
                    if diagToFuse.dt == diagA.dt:
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
                # if len(diagA.la) > 1 and len(diagB.la) > 1:
                #     print >> sys.stderr, "---"
                #     print >> sys.stderr, "diagA : %s-%s %s-%s" % (diagA.minOnG(1), diagA.maxOnG(1), diagA.minOnG(2), diagA.maxOnG(2))
                #     print >> sys.stderr, "diagB : %s-%s %s-%s" % (diagToFuse.minOnG(1), diagToFuse.maxOnG(1), diagToFuse.minOnG(2), diagToFuse.maxOnG(2))
                #     print >> sys.stderr, "ROI : %s-%s %s-%s" % (min(diagA.minOnG(1), diagToFuse.minOnG(1)), max(diagA.maxOnG(1), diagToFuse.maxOnG(1)),                                        min(diagA.minOnG(2), diagToFuse.minOnG(2)), max(diag.maxOnG(2), diagToFuse.maxOnG(2)))
                #     pass
                diagA.l1.extend(diagToFuse.l1)
                diagA.l2.extend(diagToFuse.l2)
                # Compute the new dist with currGap
                diagToFuse.la = [(aGN, aGs, dist+diagA.la[-1][2]+currGap) for (aGN, aGs, dist) in diagToFuse.la]
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
        print >> sys.stderr, "Nb of merges for currGap=%s :%s" % (currGap, nbFusionCurrGap)
    # Once all merges have been performed for a currGap it is necessary to repeat the merging process with all the diagonals for currGap+1
    print >> sys.stderr, "Total nb of merges =", nbFusion
    listOfSortedAndMergedDiagonals = sorted(list(diagGen) + listOfFinishedDiags, key=lambda diag: diag.l1[0])
    print >> sys.stderr, "Nb Diags after the merging process =", len(listOfSortedAndMergedDiagonals)

    return listOfSortedAndMergedDiagonals


# TODO merge sbs when they have a small overlap
# issue, the merge may inlclude some other diagonals...
#def mergeOverlappingSbs(sbsInPairComSV, truncationMax=0):
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

def findOverlapsOnGenome(rankGenome, sbsInPairCompWithIds, sbsGX, truncationMax=0):

    def findOverlapsOnGenomeOrder(rankGenome, order, sbsGXc, N, O, I, id2sb, truncationMax=0):
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
                        if overlap <= truncationMax:
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
                        if overlap <= truncationMax:
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
        if truncationMax == 0:
            assert len(O.keys()) == 0
        return (N, O, I)

    # method "getItemById" of sbsInPairCompWithIds"
    id2sb = sbsInPairCompWithIds.getItemById

    N = collections.defaultdict(set)
    O = collections.defaultdict(set)
    I = collections.defaultdict(set)
    for c in sbsGX.keys():
        (N, O, I) = findOverlapsOnGenomeOrder(rankGenome, +1, sbsGX[c], N, O, I, id2sb, truncationMax=truncationMax)
        if truncationMax > 0:
            # also needed to check for the overlap in the reverse
            # order
            (N, O, I) = findOverlapsOnGenomeOrder(rankGenome, -1, sbsGX[c], N, O, I, id2sb, truncationMax=truncationMax)
    return (N, O, I)


def buildConflictGraph(sbsInPairComp, truncationMax=0):

    def findOverlaps(sbsInPairCompWithIds, sbsG1, sbsG2, truncationMax=0):
        (N1, O1, I1) = findOverlapsOnGenome(1, sbsInPairCompWithIds, sbsG1, truncationMax=truncationMax)
        (N2, O2, I2) = findOverlapsOnGenome(2, sbsInPairCompWithIds, sbsG2, truncationMax=truncationMax)

        def concatenateDictsOfSets(D1, D2):
            newD = collections.defaultdict(set)
            for (key, value) in (list(D1.iteritems()) + list(D2.iteritems())):
                newD[key] = newD[key].union(value)
            return newD

        N = concatenateDictsOfSets(N1, N2)
        O = concatenateDictsOfSets(O1, O2)
        I = concatenateDictsOfSets(I1, I2)
        return (N, O, I)

    if sbsInPairComp.__class__.__name__ == 'Dict2d':
        sbsInPairCompWithIds = myTools.OrderedDict2dOfLists()
        for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
            new_sb = SyntenyBlock(sb)
            #sbsInPairCompWithIds[c1][c2].append(new_sb)
            sbsInPairCompWithIds.addToLocation((c1, c2), new_sb)
        #sbsInPairCompWithIds.identifyItems()
    elif sbsInPairComp.__class__.__name__ == 'OrderedDict2dOfLists':
        # this condition allows to keep former ids
        sbsInPairCompWithIds = sbsInPairComp
    else:
        raise TypeError('sbsInPairComp is either a Dict2D or an OrderedDict2dOfLists')
    (V, (sbsG1, sbsG2)) = buildSetOfVertices(sbsInPairCompWithIds)
    (N, O, I) = findOverlaps(sbsInPairCompWithIds, sbsG1, sbsG2, truncationMax=truncationMax)
    return (sbsInPairCompWithIds, V, N, O, I, (sbsG1, sbsG2))

def isL2NestedInL1(L1, L2):
    assert isinstance(L1, list)
    assert isinstance(L2, list)
    res = True if min(L1) <= min(L2) and max(L2) <= max(L1) else False
    return res

def isThereAnOverlap(L1, L2):
    assert isinstance(L1, list)
    assert isinstance(L2, list)
    res = True if max(L1) >= min(L2) and min(L1) <= max(L2) else False
    return res

def statusOfSbbComparedToSbaProjected((ca, sba), (cb, sbb)):
    res = 'noOverlap'
    if ca == cb:
        if isL2NestedInL1(sba.l1, sbb.l1):
            res = 'nested'
        elif isThereAnOverlap(sba.l1, sbb.l1):
            res = 'overlapping'
    return res

def statusOfSbbComparedToSba(((c1a, c2a), sba), ((c1b, c2b), sbb)):
    assert isinstance(sba, Diagonal)
    assert isinstance(sbb, Diagonal)
    resC1 = statusOfSbbComparedToSbaProjected((c1a, sba), (c1b, sbb))
    resC2 = statusOfSbbComparedToSbaProjected((c2a, sba), (c2b, sbb))
    return (resC1, resC2)

# Edit sbs in order to have non-overlapping sbs
@myTools.verbose
def solveSbsOverlaps(sbsInPairComp, truncationMax=0, removeSingleHpSbs=False, verbose=False):
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
        # sort V by increasing weights
        sortedV = sorted([int(v) for v in V], key=lambda iDsb: len(id2sb(iDsb).la))
        noOrSmallOverlapIds = set([])
        removedSbsBecauseOfOverlap = set([])
        nbRemovedHps = 0
        while len(sortedV) > 0:
            iD_selectedSb = sortedV.pop()
            # SB Not Removed
            sbnr = id2sb(iD_selectedSb)
            (c1nr, c2nr, _) = id2location(iD_selectedSb)
            noOrSmallOverlapIds.add(iD_selectedSb)
            # update the edges
            while len(N[iD_selectedSb]) > 0:
                iD_removedSb = N[iD_selectedSb].pop()
                N[iD_removedSb].remove(iD_selectedSb)
                # SB Removed
                sbr = id2sb(iD_removedSb)
                (c1r, c2r, _) = id2location(iD_removedSb)
                print >> sys.stderr, "sb (%s:%s-%s, %s:%s-%s) with %s tbs is removed because of unallowed overlap with sb (%s:%s-%s, %s:%s-%s) with %s tbs (idxs in tbs)" %\
                    (c1r, sbr.minOnG(1), sbr.maxOnG(1), c2r, sbr.minOnG(2), sbr.maxOnG(2), len(sbr.la),
                     c1nr, sbnr.minOnG(1), sbnr.maxOnG(1), c2nr, sbnr.minOnG(2), sbnr.maxOnG(2), len(sbnr.la))
                nbRemovedHps += len(sbr.la)
                if 'nested' in statusOfSbbComparedToSba(((c1nr, c2nr), sbnr), ((c1r, c2r), sbr)):
                    # sbr is at least partially nested in sbnr
                    print >> sys.stderr, "last overlap is probably a segmental duplication"
                elif 'nested' in statusOfSbbComparedToSba(((c1r, c2r), sbr), ((c1nr, c2nr), sbnr)):
                    # sbnr is at least partially nested in sbr
                    print >> sys.stderr, "last overlap is probably a (strange) segmental duplication"
                for iD_needUpdateSb in N[iD_removedSb]:
                    N[iD_needUpdateSb].remove(iD_removedSb)
                    if len(N[iD_needUpdateSb]) == 0:
                        del N[iD_needUpdateSb]
                del N[iD_removedSb]
                removedSbsBecauseOfOverlap.add(iD_removedSb)
                sortedV.remove(iD_removedSb)
            del N[iD_selectedSb]

        assert len(V) == len(noOrSmallOverlapIds) + len(removedSbsBecauseOfOverlap), "%s = %s + %s" %\
            (len(V), len(noOrSmallOverlapIds), len(removedSbsBecauseOfOverlap))

        # TODO
        # A last step would be to optimise the set using the 't-opt' algorithm
        # on page 9 of "Nonoverlapping local alignments" (Bafna 1996)

        # return the set of ids of sbs that are kept
        return (noOrSmallOverlapIds, nbRemovedHps)

    @myTools.verbose
    def truncateSbsWithSmallOverlap(sbsInPairCompWithIds, O, verbose=False):
        id2sb = sbsInPairCompWithIds.getItemById
        id2location = sbsInPairCompWithIds.getItemLocationById
        sortedO = list(O.keys())
        # sort by increasing weight
        sortedO.sort(key=lambda idsb: len(id2sb(idsb).la))
        # set of 0X keys
        setO = set(O.keys())
        nbRemovedHps = 0
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
                    nbRemovedHps += truncation
                    # if truncation > 1 and len(newla) > 1:
                    #     print >> sys.stderr, "sb (%s:%s-%s, %s:%s-%s) with %s tbs truncated (-%s tbs) because of sb (%s:%s-%s, %s:%s-%s) with %s tbs (idxs in tbs)" %\
                    #                          (c1t, sbt.minOnG(1), sbt.maxOnG(1), c2t, sbt.minOnG(2), sbt.maxOnG(2), len(sbt.la), truncation,
                    #                           c1nt, sbnt.minOnG(1), sbnt.maxOnG(1), c2nt, sbnt.minOnG(2), sbnt.maxOnG(2), len(sbnt.la))
                    #     if c1t == c1nt and c2t == c2nt:
                    #         print >> sys.stderr, "ROI: %s:%s-%s %s:%s-%s" % (c1t, min(sbt.minOnG(1), sbnt.minOnG(1)), max(sbt.maxOnG(1), sbnt.maxOnG(1)),
                    #                                                          c2t, min(sbt.minOnG(2), sbnt.minOnG(2)), max(sbt.maxOnG(2), sbnt.maxOnG(2)))
                    #         pass
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
                    nbRemovedHps += truncation
                    if truncation > 0:
                        print >> sys.stderr, "sb (%s:%s-%s, %s:%s-%s) with %s tbs truncated (-%s tbs) because of sb (%s:%s-%s, %s:%s-%s) with %s tbs (idxs in tbs)" %\
                                             (c1t, sbt.minOnG(1), sbt.maxOnG(1), c2t, sbt.minOnG(2), sbt.maxOnG(2), len(sbt.la), truncation,
                                              c1nt, sbnt.minOnG(1), sbnt.maxOnG(1), c2nt, sbnt.minOnG(2), sbnt.maxOnG(2), len(sbnt.la))
                    sbt.l1 = newl1
                    sbt.l2 = newl2
                    sbt.la = newla
            # FIXME, are you sure for 'len(id2sb(idsb).la) > 0'
            sortedO = [idsb for idsb in sortedO if idsb in setO and len(id2sb(idsb).la) > 0]
            # Need to re-sort since the length of the sb has changed
            sortedO.sort(key=lambda idsb: len(id2sb(idsb).la))
            setO = set(sortedO)
        return (sbsInPairCompWithIds, nbRemovedHps)

    (sbsInPairCompWithIds, V, N, O, I, (sbsG1, sbsG2)) =\
        buildConflictGraph(sbsInPairComp, truncationMax=truncationMax)
    nbSbsBeforeOverlapFiltering = len(V)
    (noOrSmallOverlapIds, nbRemovedHpsBecauseOfNotAllowedOverlap) = solveConflictGraph(sbsInPairCompWithIds, V, N, verbose=verbose)
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

    nbRemovedHpsAfterTruncation = 0
    if truncationMax > 0:
        # Need to edit sbs to truncate sbs that have an allowed overlap
        #new_sbsInPairCompWithIds = truncateSbsWithSmallOverlap(new_sbsInPairCompWithIds, O, verbose=verbose)
        (sbsInPairCompWithIds, nbRemovedHpsAfterTruncation) = truncateSbsWithSmallOverlap(sbsInPairCompWithIds, O, verbose=verbose)

    removedIdBecauseEmptySb = set([])
    cptRemovedSbAfterTruncation = 0
    for (idsb, (c1, c2), sb) in sbsInPairCompWithIds.iterByOrderedIds():
        # Remove all sbs too small, (they are small because of the preceding truncation)
        # assert truncationMax > 0
        if len(sb.la) == 0:
            removedIdBecauseEmptySb.add(idsb)
            print >> sys.stderr, "Removed sb (%s,%s) because of truncation (less than 1tb in the sb after truncation)" % (c1, c2)
            cptRemovedSbAfterTruncation += 1
        elif len(sb.la) == 1 and removeSingleHpSbs:
            removedIdBecauseEmptySb.add(idsb)
            print >> sys.stderr, "Removed sb (%s,%s) because of truncation (1tb in the sb after truncation and removeSingleHpSbs==True)" % (c1, c2)
            nbRemovedHpsAfterTruncation += 1
            cptRemovedSbAfterTruncation += 1

    sbsInPairCompWithIds.removeIds(removedIdBecauseEmptySb)
    # DEBUG
    if __debug__:
        (sbsInPairCompWithIds, V, N, O, I, (sbsG1, sbsG2)) =\
            buildConflictGraph(sbsInPairCompWithIds, truncationMax=0)
        assert len(N) == 0, N
        assert len(O) == 0, O

        #assert set(O.keys()) >= noOrSmallOverlapIds
        # FIXME, is it right to do that ?
        #O = dict([(iDsba, iDsbb)  for (iDsba, iDsbb) in O.iteritems() if ((iDsba in noOrSmallOverlapIds) and (iDsbb in noOrSmallOverlapIds))])
        #assert len(noOrSmallOverlapIds) == len([a for a in new_sbsInPairComp])

    new_sbsInPairComp = myTools.Dict2d(list)
    for (idsb, (c1, c2), sb) in sbsInPairCompWithIds.iterByOrderedIds():
        new_sbsInPairComp[c1][c2].append(sb)

    # This print is already in the upstream function
    #print >> sys.stderr, "Nb sbs before overlap-filtering = %s" % nbSbsBeforeOverlapFiltering
    assert cptRemovedSbBecauseOfNotAllowedOverlap == nbSbsBeforeOverlapFiltering - len(noOrSmallOverlapIds), \
        "%s = %s - %s" % (cptRemovedSbBecauseOfNotAllowedOverlap, nbSbsBeforeOverlapFiltering, len(noOrSmallOverlapIds))
    nbSbsAfterOverlapFiltering = len(list(new_sbsInPairComp.iteritems2d()))
    # Total nb of hps removed
    nbHpsBeforeOverlapFiltering = sum([len(sb.la) for (_, sb) in sbsInPairComp.iteritems2d()])
    nbHpsAfterOverlapFiltering = sum([len(sb.la) for (_, sb) in new_sbsInPairComp.iteritems2d()])
    assert nbHpsAfterOverlapFiltering == nbHpsBeforeOverlapFiltering - (nbRemovedHpsBecauseOfNotAllowedOverlap + nbRemovedHpsAfterTruncation),\
        "%s = %s - (%s + %s)" % (nbHpsAfterOverlapFiltering, nbHpsBeforeOverlapFiltering, nbRemovedHpsBecauseOfNotAllowedOverlap, nbRemovedHpsAfterTruncation)
    print >> sys.stderr, "Nb hps before overlap-filtering = %s" % nbHpsBeforeOverlapFiltering
    assert nbHpsBeforeOverlapFiltering > 0
    print >> sys.stderr, "Nb hps after overlap-filtering = %s (%%%.2f)" % (nbHpsAfterOverlapFiltering, float(nbHpsAfterOverlapFiltering) / float(nbHpsBeforeOverlapFiltering))
    print >> sys.stderr, "Nb hps removed because of not allowed overlap = %s (%%%.2f)" % (nbRemovedHpsBecauseOfNotAllowedOverlap, float(nbRemovedHpsBecauseOfNotAllowedOverlap)/ float(nbHpsBeforeOverlapFiltering))
    print >> sys.stderr, "Nb hps removed because of truncation = %s (%%%.2f)" % (nbRemovedHpsAfterTruncation, float(nbRemovedHpsAfterTruncation)/ float(nbHpsBeforeOverlapFiltering))
    assert nbSbsAfterOverlapFiltering == nbSbsBeforeOverlapFiltering - (cptRemovedSbAfterTruncation + cptRemovedSbBecauseOfNotAllowedOverlap)
    print >> sys.stderr, "Nb sbs removed because of not allowed overlap = %s (%%%.2f)" % (cptRemovedSbBecauseOfNotAllowedOverlap, float(cptRemovedSbBecauseOfNotAllowedOverlap) / float(nbSbsBeforeOverlapFiltering))
    print >> sys.stderr, "Nb sbs removed during truncation = %s (%%%.2f)" % (cptRemovedSbAfterTruncation, float(cptRemovedSbAfterTruncation) / float(nbSbsBeforeOverlapFiltering))
    # This print is already in the upstream function
    #print >> sys.stderr, "Nb non-overlapping sbs returned = %s" % nbSbsAfterOverlapFiltering

    if __debug__:
        # DEBUG, verify that the filtered sbs are not overlapping
        (_, _, N, O, _, _) = buildConflictGraph(new_sbsInPairComp, truncationMax=0)
        assert len(N) == 0, N
        assert len(O) == 0, O

    return (new_sbsInPairComp)


def strandProduct(sa, sb):
    if sa != None and sb != None:
        return sa * sb
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
    # 1 faster
    locG2 = {}
    for (i2, (f, s2)) in enumerate(gc2):
        if f != None:
            if f not in locG2:
                locG2[f] = []
            locG2[f].append((i2, s2))  # LOCalisation on Gc 2
            #locG2.setdefault(f,[]).append( (i2,s2) ) # LOCalisation on Gc 2, same excecution time but less readable
    M = {}
    if not locG2:  # if locG2 is empty
        return (M, locG2)
    for (i1, (f, s1)) in enumerate(gc1):
        if f != None and f in locG2:  #TODO : remove by advance the f == None from gc1
            M[i1]={}
            for (i2, s2) in locG2[f]:
                M[i1][i2] = strandProduct(s1, s2)

    #2 slower than 1
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

    #3 (faster than 1), but since we need to return locG2 of the input gc2 for next computations, it is not convenient
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

    #4 slower than 1
    #familiesLocation1 = collections.defaultdict(list)
    #for (i1, (f, s1)) in enumerate(gc1):
    #    if f is not None:
    #        familiesLocation1[f].append(i1)
    #familiesLocation2 = collections.defaultdict(list)
    #locG2 = collections.defaultdict(list)
    #for (i2, (f, s2)) in enumerate(gc2):
    #    if f is not None:
    #        familiesLocation2[f].append(i2)
    #        locG2[f].append((i2, s2))
    #        M = collections.defaultdict(lambda: collections.defaultdict(list))
    #for f in set(familiesLocation1.keys()) | set(familiesLocation2.keys()):
    #    for (i1, i2) in itertools.product(familiesLocation1[f], familiesLocation2[f]):
    #        M[i1][i2] = strandProduct(gc1[i1][1], gc2[i2][1])

    #5: (slower than 2) method inspired by a discussion with Eric Tannier 10/04/2015
    # M = {}
    # # s: sorted, the strings are sorted following the string '<' operator, the alphabetical order
    # gc1s = sorted(enumerate(gc1), key=lambda x: x[1].n)
    # gc2s = sorted(enumerate(gc2), key=lambda x: x[1].n)
    # # c: compact
    # gc1sc = []
    # oldGn1 = None
    # for (idxG1, (gn1, s1)) in gc1s:
    #     if gn1 != oldGn1:
    #         gc1sc.append((gn1, [(idxG1, s1)]))
    #     else:
    #         gc1sc[-1][1].append((idxG1, s1))
    #     oldGn1 = gn1
    # gc2sc = []
    # oldGn2 = None
    # for (idxG2, (gn2, s2)) in gc2s:
    #     if gn2 != oldGn2:
    #         gc2sc.append((gn2, [(idxG2, s2)]))
    #     else:
    #         gc2sc[-1][1].append((idxG2, s2))
    #     oldGn2 = gn2
    #
    # lenGC1SC = len(gc1sc)
    # lenGC2SC = len(gc2sc)
    # idxGN1 = 0
    # idxGN2 = 0
    # while idxGN1 < lenGC1SC and idxGN2 < lenGC2SC:
    #     gn1 = gc1sc[idxGN1][0]
    #     gn2 = gc2sc[idxGN2][0]
    #     if gn2 == gn1:
    #         lidxGS1 = gc1sc[idxGN1][1]
    #         lidxGS2 = gc2sc[idxGN2][1]
    #         for (idxG1, s1) in lidxGS1:
    #             if idxG1 not in M:
    #                 M[idxG1] = {}
    #             for (idxG2, s2) in lidxGS2:
    #                 M[idxG1][idxG2] = strandProduct(s1, s2)
    #         idxGN1 += 1
    #         idxGN2 += 1
    #     else:
    #         if gn1 < gn2:
    #             idxGN1 += 1
    #         else:
    #             # gn2 < gn1
    #             idxGN2 += 1
    # # needed in some algorithms
    # locG2 = dict(gc2sc)

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
def findDiagType(i1, i2, M, sameStrand):
    if M[i1][i2] != None and sameStrand:
        diagType = '/' if M[i1][i2] == +1 else '\\'
    else:
        if i1+1 in M and i2+1 in M[i1+1] and ((M[i1+1][i2+1] in [1,None]) if sameStrand else True):
            diagType = '/'
        elif i1+1 in M and i2-1 in M[i1+1] and ((M[i1+1][i2-1] in [-1,None]) if sameStrand else True):
            diagType = '\\'
        else:
            diagType = None
    return diagType


def extractDiagsInPairCompChrWrapper((c1, c2, gc1, gc2, sameStrand)):
    # print >> sys.stderr, "(PPID = %s, PID = %s) start to extract diagonals on G1[%s]_vs_G2[%s]" %\
    #     (os.getppid(), os.getpid(), c1, c2)
    (listOfDiags, nbHomologies) = extractDiagsInPairCompChr(gc1, gc2, sameStrand, verbose=False)
    return ((c1, c2), listOfDiags, nbHomologies)

# Extract diags in a pairwise comparison of two chromosomes
############################################################
@myTools.verbose
def extractDiagsInPairCompChr(gc1, gc2, sameStrand=True, verbose=False):
    listOfDiags = []
    (M, locG2) = homologyMatrix(gc1, gc2)
    nbHomologies = sum([len(M[i]) for i in M])
    if not locG2:
        # if locG2 is empty
        return (listOfDiags, nbHomologies)
    la = []
    l1 = []
    l2 = []
    diagType = None

    #TODO scan M instead of gc1 : impossible since 'dict size cannot change during a loop over its items'
    # scan M from left to right
    for (i1, (f, _)) in enumerate(gc1):
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
                        diagType = findDiagType(i1, i2, M, sameStrand)

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
                    #if (diagType == "/" or diagType == None) and i1+1 in M and i2+1 in M[i1+1] and ((M[i1+1][i2+1] in [+1,None]) if sameStrand else True):
                    if diagType == "/" and i1+1 in M and i2+1 in M[i1+1] and ((M[i1+1][i2+1] in [+1,None]) if sameStrand else True):
                        i1 = i1 + 1
                        i2 = i2 + 1
                        #assert i2-l2[-1] == +1
                        #assert i1-l1[-1] == +1
                        #assert i2 in M[i1]

                    elif diagType=="\\" and i1+1 in M and i2-1 in M[i1+1] and ((M[i1+1][i2-1] in [-1,None]) if sameStrand else True):
                        i1 = i1 + 1
                        i2 = i2 - 1
                        #assert i2-l2[-1] == -1
                        #assert i1-l1[-1] == +1
                        #assert i2 in M[i1]
                    else:
                        # Since no more hps can be added to the current diagonal, the diagonal is recorded
                        #assert len(la) > 0
                        diag = Diagonal(diagType, l1, l2, la)
                        listOfDiags.append(diag)
                        l1 = []
                        l2 = []
                        la = []
                        diagType = None
                        break # exit the while loop and iter the for loop
    return (listOfDiags, nbHomologies)


def crossGeneContent(g1, g2):
    # create a set with all gene names
    geneNames1 = g1.getGeneNames(g1, checkNoDuplicates=True)
    geneNames2 = g2.getGeneNames(g2, checkNoDuplicates=True)
    gNsInCommon = geneNames1.intersection(geneNames2)
    gNsInG1notInG2 = geneNames1 - gNsInCommon
    gNsInG2notInG1 = geneNames2 - gNsInCommon
    return (gNsInCommon, gNsInG1notInG2, gNsInG2notInG1)


# Depending on the filterType parameter:
#       - filterType = FilterType.None (genomes are not filtered)
#       - filterType = FilterType.InFamilies (only genes herited from the ancestor are kept)
#       - filterType = FilterType.InBothGenomes (only 'anchor genes', ie genes present in both species, are kept)
# Returns (g1,g2,trans1,trans2)
#      - g1 the genome 'g1' rewritten with ancGenes ID
#       - trans1 = { ..., newi:oldi, ...} with newi, the index of a gene in 'g1' and oldi the index of the same gene in the original genome 'g1'
def filter2D(g1_orig, g2_orig, filterType, minChromLength=1, keepOriginal=False):
    # Mark genes that are not in the intersection for future removal
    # Marking is done by switching (g,s) into (None,s)
    # warning : modifies g1 and g2
    isinstance(g1_orig, myLightGenomes.LightGenome)
    if keepOriginal:
        g1 = myLightGenomes.LightGenome(g1_orig)
        g2 = myLightGenomes.LightGenome(g2_orig)
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
    elif filterType == FilterType.InFamilies:
        (g1, mG1f2G1o, (nCL1, nGL1)) = remapCoFilterContentAndSize(g1, {None}, minChromLength)
        (g2, mG2f2G2o, (nCL2, nGL2)) = remapCoFilterContentAndSize(g2, {None}, minChromLength)
    # Only conserve genes present in both extant genomes
    elif filterType == FilterType.InBothGenomes:
        mG1f2G1o = None
        mG2f2G2o = None
        (nCL1, nGL1) = (0, 0)
        (nCL2, nGL2) = (0, 0)
        while True:
            # after this step genes that have no homolog in the other genome are marked None
            # 'gNsInCommon' is also called the set of 'anchor genes' in bibliography
            geneNames1 = g1.getGeneNames(checkNoDuplicates=False)
            geneNames2 = g2.getGeneNames(checkNoDuplicates=False)
            gNsInCommon = geneNames1 & geneNames2
            gNsToRemove = ((geneNames1 | geneNames2) - gNsInCommon) | {None}
            (g1, mG1f2G1o, (tmp_nCL1, tmp_nGL1)) = remapCoFilterContentAndSize(g1, gNsToRemove, minChromLength,
                                                                       mOld=mG1f2G1o)
            (g2, mG2f2G2o, (tmp_nCL2, tmp_nGL2)) = remapCoFilterContentAndSize(g2, gNsToRemove, minChromLength,
                                                                       mOld=mG2f2G2o)
            nCL1 += tmp_nCL1
            nGL1 += tmp_nGL1
            nCL2 += tmp_nCL2
            nGL2 += tmp_nGL2

            hasChanged = (tmp_nGL1 > 0) or (tmp_nGL2 > 0)
            # If a chromosome has been removed because of the filtering on the length,
            # the filtering is performed once more.
            if not hasChanged:
                break
        # there may be duplicates which family is in both genomes
        assert g1.getGeneNames(checkNoDuplicates=False) == g2.getGeneNames(checkNoDuplicates=False)
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
    for g in range(0, maxGapThreshold+1):
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
@myTools.tictac
@myTools.verbose
def numberOfHomologies(g1, g2, verbose=False):
   nbHomologies = myTools.Dict2d(int)
   totalNbComps = len(g1) * len(g2)
   print >> sys.stderr, "pairwise comparison of chromosomes analysis for counting hps"
   progressBar = myTools.ProgressBar(totalNbComps)
   for (cptComp, (c1, c2)) in enumerate(itertools.product(g1, g2)):
       gc1 = g1[c1]
       gc2 = g2[c2]
       (Ms, _) = homologyMatrix(gc1, gc2)
       progressBar.printProgressIn(sys.stderr, cptComp)
       nbHomologies[c1][c2] = sum([len(Ms[i1]) for i1 in Ms])
   # new line in the print
   nbHomologies_g = sum([nbH for nbH in nbHomologies.values2d()])
   # assert nbHomologies_g >= min(sum([len(g1[c]) for c in g1]), sum([len(g2[c]) for c in g2])),"%s,%s" %  (sum([len(g1[c]) for c in g1]), sum([len(g2[c]) for c in g2]))
   # Not needed since the two lineages may have undergone some differential gene losses
   return nbHomologies, nbHomologies_g


# Statistical validation of synteny blocks
# because of collections.Counter
@myTools.minimalPythonVersion((2, 7))
@myTools.verbose
def statisticalValidation(diagsInPairComp, g1_tb, g2_tb, N12s, p_hpSign,
                          pThreshold=1.0,
                          NbOfHomologiesThreshold=50,
                          validateImpossToCalc_mThreshold=3,
                          considerMonogenicSb=False,
                          validateSbsWithGapMaxZero=True,
                          verbose=False):

    def fNbDiags(diagsInPairCompX):
        return len(list(diagsInPairCompX.iteritems2d()))

    diagsInPairCompRejected = myTools.Dict2d(list)
    diagsInPairCompImpossibleToCalcProba = myTools.Dict2d(list)
    sbsInPairCompStatVal = myTools.Dict2d(list)
    for ((c1, c2), diag) in diagsInPairComp.iteritems2d():
        (m, max_g, lw1, lw2, l1_min, l1_max, l2_min, l2_max) = diag.calculateCharacteristics()
        assert m > 0
        na = len(g1_tb[c1])  # Nb of Tbs on C1
        nb = len(g2_tb[c2])  # Nb of Tbs on C2
        nab = N12s[c1][c2]  # Nb of homologies in the mhp
        if m == 1 and not considerMonogenicSb:
            # all diagonals of length 1 are rejected
            diagsInPairCompRejected[c1][c2].append((diag, None))
            continue
        elif m > nab:
            # there are not m hps in the MHP
            diagsInPairCompImpossibleToCalcProba[c1][c2].append((diag, None))
            continue
        elif m > min([lw1, lw2]):
            # there are too many dispersed paralogies in the window
            diagsInPairCompImpossibleToCalcProba[c1][c2].append((diag, None))
            continue

        # This is to avoid excessive time consuming computations
        if m == 1 and considerMonogenicSb:
            p = 1
        elif m >= 2 and max_g == 0 and validateSbsWithGapMaxZero:
            p = 0
        elif m > NbOfHomologiesThreshold:
            p = 0
        else:
            p = myProbas.pValue(m, max_g, lw1, lw2, nab, na, nb, p_hpSign[c1][c2], verbose=verbose)

        if p is None:
            diagsInPairCompImpossibleToCalcProba[c1][c2].append((diag, None))
        elif p <= pThreshold:
            sbsInPairCompStatVal[c1][c2].append(SyntenyBlock(diag, p))
        else:
            diagsInPairCompRejected[c1][c2].append((diag, p))
            assert len(diag.l1) == len(diag.l1)
            assert len(diag.la) == len(diag.l2)

    # automatically validate sbs that contain more than validateImpossToCalc_mThreshold hps
    sbsInPairCompImpossibleToCalcProbaStatVal = myTools.Dict2d(list)
    diagsInPairCompImpossibleToCalcProbaRejected = myTools.Dict2d(list)
    for ((c1, c2), (diag, pVal)) in diagsInPairCompImpossibleToCalcProba.iteritems2d():
        m = diag.calculateCharacteristics()[0]
        if m >= validateImpossToCalc_mThreshold:
            # if the diagonal contains more than validateImpossToCalc_mThreshold
            # hps, it is validated in order to avoid to reject long and perfect
            # diagonals because of only one dispersed tandem duplication
            sbsInPairCompImpossibleToCalcProbaStatVal[c1][c2].append(SyntenyBlock(diag, pVal))
        else:
            diagsInPairCompImpossibleToCalcProbaRejected[c1][c2].append((diag, pVal))


    assert fNbDiags(sbsInPairCompImpossibleToCalcProbaStatVal) + fNbDiags(diagsInPairCompImpossibleToCalcProbaRejected) == fNbDiags(diagsInPairCompImpossibleToCalcProba),\
        "%s + %s = %s" % (fNbDiags(sbsInPairCompImpossibleToCalcProbaStatVal), fNbDiags(diagsInPairCompImpossibleToCalcProbaRejected), fNbDiags(diagsInPairCompImpossibleToCalcProba))

    # update the lists
    sbsInPairCompStatVal = sbsInPairCompStatVal + sbsInPairCompImpossibleToCalcProbaStatVal
    diagsInPairCompRejected = diagsInPairCompRejected + diagsInPairCompImpossibleToCalcProbaRejected

    assert fNbDiags(sbsInPairCompStatVal) + fNbDiags(diagsInPairCompRejected) == fNbDiags(diagsInPairComp),\
        "%s + %s = %s" % (fNbDiags(sbsInPairCompStatVal), fNbDiags(diagsInPairCompRejected), fNbDiags(diagsInPairComp))

    firstPrint = True
    for ((c1, c2), (diag, pVal)) in diagsInPairCompRejected.iteritems2d():
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
                (c1, l1_min, l1_max, c2, l2_min, l2_max, m, max_g, lw1, lw2, nab, na, nb), "\t p=", pVal

    firstPrint = True
    for ((c1, c2), sb) in sbsInPairCompStatVal.iteritems2d():
        (m, max_g, lw1, lw2, l1_min, l1_max, l2_min, l2_max) = sb.calculateCharacteristics()
        if m == 2:
            na = len(g1_tb[c1])  # Nb of Tbs on C1
            nb = len(g2_tb[c2])  # Nb of Tbs on C2
            nab = N12s[c1][c2]  # Nb of homologies in the mhp
            if firstPrint:
                print >> sys.stderr, "Diagonals containing 2 hps that have passed the statistical test:"
                firstPrint = False
            print >> sys.stderr, "(c1=%s:%s-%s,c2=%s:%s-%s) \t (m=%s, max_g=%s, lw1=%s lw2=%s, nab=%s, na=%s, nb=%s)" %\
                (c1, l1_min, l1_max, c2, l2_min, l2_max, m, max_g, lw1, lw2, nab, na, nb), "\t p=", sb.pVal

    def statsDiagsM(diagsInPairCompX, m):
        if len(diagsInPairCompX.keys2d()) == 0:
            return 'RAS'
        assert all([len(diagsInPairCompX[k1Foo][k2Foo]) > 0 for (k1Foo, k2Foo) in diagsInPairCompX.keys2d()])
        (k1Foo, k2Foo) = diagsInPairCompX.keys2d()[0]
        if isinstance(diagsInPairCompX[k1Foo][k2Foo][0], tuple):
            diagsInPairCompXM = [diag for sbs in diagsInPairCompX.values2d() for (diag, pValue) in sbs if len(diag.la) == m]
        else:
            # diagsInPairComp has no pValue
            assert isinstance(diagsInPairCompX[k1Foo][k2Foo][0], Diagonal)
            diagsInPairCompXM = [diag for diags in diagsInPairCompX.values2d() for diag in diags if len(diag.la) == m]
        diagsInPairCompXMGaps = [diag.max_g() for diag in diagsInPairCompXM if len(diag.la) == m]
        diagsXMGaps = collections.Counter(diagsInPairCompXMGaps)
        diagsXMGaps = ["%s:%s" % (length, nb) for (length, nb) in sorted(diagsXMGaps.items())]
        return diagsXMGaps
    print >> sys.stderr, "Over all sbs of 2 hps before stat. val. the distribution of gap maximum is: {%s}" % " ".join(statsDiagsM(diagsInPairComp,2))
    print >> sys.stderr, "Over all rejected sbs of 2 hps the distribution of gap maximum is: {%s}" % " ".join(statsDiagsM(diagsInPairCompRejected,2))
    print >> sys.stderr, "Over all sbs where it is imposs. to calculate the proba. with 2 hps the distribution of gap maximum is : {%s} (due to dispersed paralogies)" % " ".join(statsDiagsM(diagsInPairCompImpossibleToCalcProba,2))
    print >> sys.stderr, "Over all sbs where it is imposs. to calculate the proba. with 2 hps which are finally validated, the distribution of gap maximum is: {%s} (due to dispersed paralogies)" % " ".join(statsDiagsM(sbsInPairCompImpossibleToCalcProbaStatVal,2))
    print >> sys.stderr, "Over all sbs where it is imposs. to calculate the proba. with 2 hps which are finally rejected, the distribution of gap maximum is: {%s} (due to dispersed paralogies)" % " ".join(statsDiagsM(diagsInPairCompImpossibleToCalcProbaRejected,2))

    print >> sys.stderr, "Over all stat. val. sbs of 2 hps the distribution of gap maximum is: {%s}" % " ".join(statsDiagsM(sbsInPairCompStatVal,2))

    def statsDiagLengths(diagsInPairCompX):
        if len(diagsInPairCompX.keys2d()) == 0:
            return 'RAS'
        (k1Foo, k2Foo) = diagsInPairCompX.keys2d()[0]
        if isinstance(diagsInPairCompX[k1Foo][k2Foo][0], tuple):
            diagsInPairCompX_ = [diag for sbs in diagsInPairCompX.values2d() for (diag, pValue) in sbs]
        else:
            # diagsInPairComp has no pValue
            # assert diagsInPairCompX[k1Foo][k2Foo][0].__class__.name == 'Diagonal'
            assert isinstance(diagsInPairCompX[k1Foo][k2Foo][0], Diagonal)
            diagsInPairCompX_ = [diag for diags in diagsInPairCompX.values2d() for diag in diags]
        diagXLengths = collections.Counter([len(diags.la) for diags in diagsInPairCompX_])
        diagXLengths = ["%s:%s" % (length, nb) for (length, nb) in sorted(diagXLengths.items())]
        return diagXLengths
    print >> sys.stderr, "Over all diagonal before the stat. val., distribution of the diag lengths: {%s}" % " ".join(statsDiagLengths(diagsInPairComp))
    print >> sys.stderr, "Over all rejected diagonals, distribution of all diag lengths: {%s}" % " ".join(statsDiagLengths(diagsInPairCompRejected))
    print >> sys.stderr, "Over all diagonals with p-Value impossible to compute (mainly due to paralogies), distribution of diag lengths: {%s}" % " ".join(statsDiagLengths(diagsInPairCompImpossibleToCalcProba))
    print >> sys.stderr, "Over all diagonals with p-Value impossible to compute (mainly due to paralogies) which are finally validated, distribution of diag lengths: {%s}" % " ".join(statsDiagLengths(sbsInPairCompImpossibleToCalcProbaStatVal))
    print >> sys.stderr, "Over all diagonals with p-Value impossible to compute (mainly due to paralogies) which are finally rejected, distribution of diag lengths: {%s}" % " ".join(statsDiagLengths(diagsInPairCompImpossibleToCalcProbaRejected))
    print >> sys.stderr, "Over all diagonals that passed the stat. val., distribution of diag lengths: {%s}" % " ".join(statsDiagLengths(sbsInPairCompStatVal))


    tmp = myTools.Dict2d(list)
    for ((c1, c2), (diag, sb)) in diagsInPairCompRejected.iteritems2d():
        tmp[c1][c2].append(diag)
    diagsInPairCompRejected = tmp

    return (sbsInPairCompStatVal, diagsInPairCompRejected)


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

def splitSbBySplitRanks(sb, splitRanks):
    newSbs = []
    for (xBeg, xEnd) in myTools.myIterator.slidingTuple(splitRanks):
        newDt = sb.dt
        newl1 = sb.l1[xBeg:xEnd]
        newl2 = sb.l2[xBeg:xEnd]
        newla = sb.la[xBeg:xEnd]
        if newla[0][2] > 1:
            newla = [(aG, aGs, (dist - newla[0][2] + 1)) for (aG, aGs, dist) in newla]
        assert newla[0][2] == 1
        # Here we assign the former pVal to all sub-sbs. The
        # hypothesis is to consider that since the sub-sbs were
        # neighbours within the former sb, their significances have
        # already been assessed by the former pValue calculation.
        if isinstance(sb, SyntenyBlock):
            newPval = sb.pVal
        else:
            assert isinstance(sb, Diagonal)
            newPval = 0
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

        assert splitRanks_a[0] == 0 and splitRanks_a[-1] == len(sba.la)
        assert splitRanks_b[0] == 0 and splitRanks_b[-1] == len(sbb.la)

    #newSbs = splitSbBySplitRanks(sba, splitRanks_a)
    #newSbs.extend(splitSbBySplitRanks(sbb, splitRanks_b))
    return (splitRanks_a, splitRanks_b)

# Split a sb as soon as one of its gap contains another sb
# the other sb may be a micro-inversion or a transposition-like sb
# FIXME : it seems that breakpoints in gaps due to mono-genic conserved segments are not identified
def fIdentifyBreakpointsWithinGaps(sbsInPairComp, removeSingleHpSbs=False):
    # Find sbs that are nested within gaps of other sbs
    #   * nested inversions (micro-inversions, nested segments and inverted
    #     diagonalType).
    #   * transpositions (only 'one' nested segment, to avoid considering
    #     segmental duplications followed by fractionation).
    sbsInPairCompWithIds = myTools.OrderedDict2dOfLists()
    for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
        new_sb = SyntenyBlock(sb)
        sbsInPairCompWithIds.addToLocation((c1, c2), new_sb)
    id2sb = sbsInPairCompWithIds.getItemById
    id2location = sbsInPairCompWithIds.getItemLocationById

    (V, (sbsG1, sbsG2)) = buildSetOfVertices(sbsInPairCompWithIds)
    (N1, O1, I1) = findOverlapsOnGenome(1, sbsInPairCompWithIds, sbsG1, truncationMax=0)
    (N2, O2, I2) = findOverlapsOnGenome(2, sbsInPairCompWithIds, sbsG2, truncationMax=0)

    # Intra-synteny block rearrangement: either an intra-sb
    # micro-transposition or an intra-sb micro-inversion (inverted
    # diagType). Or a mix of both. To disentangle both phenomenon, it
    # would be necessary to check the rank location in the host sb.
    todo = set([])
    id2splitRanks = collections.defaultdict(list)
    for idsbb in (set(I1.keys()) & set(I2.keys())):
        for idsba in (I1[idsbb] & I2[idsbb]):
            # sbb is included in sba in both genomes
            sba = id2sb(idsba)
            sbb = id2sb(idsbb)
            (c1a, c2a, _) = id2location(idsba)
            (c1b, c2b, _) = id2location(idsbb)
            assert {c1a, c2a} == {c1b, c2b}
            assert c1a == c1b and c2a == c2b
            #
            if noOverlapSb(((c1a, c2a), sba), ((c1b, c2b), sbb)):
                # FIXME: removed because it is interesting to split the
                # diagonal in the case of an intra-synteny-block transposition
                #if sba.dt == '/' and sbb.dt == '\\' or sba.dt == '\\' and sbb.dt == '/':

                (splitRanks_al1, splitRanks_bl1) = splitNestedLists(sba.l1, sbb.l1)
                (splitRanks_al2, splitRanks_bl2) = splitNestedLists(sba.l2, sbb.l2)
                if len(splitRanks_al1) == len(splitRanks_al2) and len(splitRanks_bl1)  == len(splitRanks_bl2):
                    # this last condition was for avoiding strange cases
                    # +A+B+C+D+E
                    # +A-D+B-C+E
                    # where -C-D is seen has an included diagonal

                    # Return tuples (id of the splited sb, splitRanks)
                    # This allows to incrementally update the splitRanks
                    #  and then it wil be possible to perform the splits
                    #  after the for loops
                    (splitRanks_a, splitRanks_b) = splitNestedSbs(((c1a, c2a), sba), ((c1b, c2b), sbb))

                    id2splitRanks[idsba] = sorted(list(set(id2splitRanks[idsba] + splitRanks_a)))
                    id2splitRanks[idsbb] = sorted(list(set(id2splitRanks[idsbb] + splitRanks_b)))
                    todo = todo | {idsba, idsbb}

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
                # Return tuples (id of the spited sb, splitRanks)
                # This allows to incrementally update the splitRanks
                #  and then it wil be possible to perform the splits
                #  after the for loops
                (splitRanks_a, splitRanks_b) = splitNestedSbs(((c1a, c2a), sba), ((c1b, c2b), sbb))
                id2splitRanks[idsba] = sorted(list(set(id2splitRanks[idsba] + splitRanks_a)))
                id2splitRanks[idsbb] = sorted(list(set(id2splitRanks[idsbb] + splitRanks_b)))
                todo = todo | {idsba, idsbb}

    # With 'todo' and (id of the splited sb, splitRanks)s, calculate the newSbs and removed sbs
    removedSbIds = set([])
    for idsb in todo:
        (c1, c2, _) = id2location(idsb)
        for newSb in splitSbBySplitRanks(id2sb(idsb), id2splitRanks[idsb]):
            sbsInPairCompWithIds.addToLocation((c1, c2), newSb)
        removedSbIds = removedSbIds | {idsb}
    sbsInPairCompWithIds.removeIds(removedSbIds)

    assert len(set(sbsInPairCompWithIds.orderedIds) & removedSbIds) == 0

    # TODO print in stderr the reason of the removal
    if removeSingleHpSbs:
        # remove sbs of only one hp
        singleHpSbs = set([])
        for (idsb, _, sb) in sbsInPairCompWithIds.iterByOrderedIds():
            assert len(sb.la) > 0
            if len(sb.la) == 1:
                singleHpSbs.add(idsb)
        sbsInPairCompWithIds.removeIds(singleHpSbs)

    new_sbsInPairComp = myTools.Dict2d(list)
    for (idsb, (c1, c2), sb) in sbsInPairCompWithIds.iterByOrderedIds():
        new_sbsInPairComp[c1][c2].append(sb)
    return new_sbsInPairComp

# TODO: split this function in two modules:
# 1) IdentifyMicroInversions nested in sbs gaps
def fIdentifyMonoGcSsAndIdentifyMonoGinvsNestedInSbsGaps(sbsInPairComp, putativeMicroInversionsInPairComp,
                                            gapMaxMicroInv=0, identifyMonogInvs=True):
    # 1) Give unique id to each sb and diag
    idsSbs = set()
    idsPutativeMicroInversions = set()
    sbsAndDiagsInPairCompWithIds = myTools.OrderedDict2dOfLists()
    id = None
    for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
        new_sb = SyntenyBlock(sb)
        id = sbsAndDiagsInPairCompWithIds.addToLocation((c1, c2), new_sb)
        idsSbs.add(id)
    assert id in {None, sbsAndDiagsInPairCompWithIds.maxId}
    if id is None:
        id = -1
    for ((c1, c2), diag) in putativeMicroInversionsInPairComp.iteritems2d():
        id += 1
        new_diag = Diagonal(diag)
        sbsAndDiagsInPairCompWithIds.addToLocationWithId((c1, c2), new_diag, id)
        idsPutativeMicroInversions.add(id)
    assert len(idsSbs & idsPutativeMicroInversions) == 0
    id2sb = sbsAndDiagsInPairCompWithIds.getItemById
    id2location = sbsAndDiagsInPairCompWithIds.getItemLocationById

    # 2) build the overlap graph
    (V, (sbsG1, sbsG2)) = buildSetOfVertices(sbsAndDiagsInPairCompWithIds)
    (N1, O1, I1) = findOverlapsOnGenome(1, sbsAndDiagsInPairCompWithIds, sbsG1, truncationMax=0)
    (N2, O2, I2) = findOverlapsOnGenome(2, sbsAndDiagsInPairCompWithIds, sbsG2, truncationMax=0)

    # 3) identify micro-inversions nested in sbs gaps
    todo = set()
    id2splitRanks = collections.defaultdict(list)
    idsIdentifiedMicroInv = set()
    # Intra-synteny block microInv
    for idsbb in (set(I1.keys()) & set(I2.keys())):
        sbb = id2sb(idsbb)
        for idsba in (I1[idsbb] & I2[idsbb]):
            # sbb is included in sba in both genomes
            sba = id2sb(idsba)
            (c1a, c2a, _) = id2location(idsba)
            (c1b, c2b, _) = id2location(idsbb)
            assert c1a == c1b and c2a == c2b
            if idsbb in idsPutativeMicroInversions:
                if noOverlapSb(((c1a, c2a), sba), ((c1b, c2b), sbb)):
                    if sba.dt == '/' and sbb.dt == '\\' or sba.dt == '\\' and sbb.dt == '/':
                        if len(sbb.la) == 1:
                            assert sbb.la[0][1] is not None, sbb.la[0][1]
                            if not identifyMonogInvs:
                                continue
                        (splitRanks_a, splitRanks_b) = splitNestedSbs(((c1a, c2a), sba), ((c1b, c2b), sbb))
                        if len(splitRanks_a[1:-1]) == 1 and len(splitRanks_b[1:-1]) == 0:
                            # if sbb is nested in a gap of sba, with a gapMaxMicroInv
                            sbb_l = min(sbb.l1)
                            sbb_r = max(sbb.l1)
                            sbb_b = min(sbb.l2)
                            sbb_t = max(sbb.l2)
                            sba_ll = min([(sbb_l - idx, idx) for idx in sba.l1 if idx < sbb_l])[1]
                            sba_rr = min([(idx - sbb_r, idx) for idx in sba.l1 if sbb_r < idx])[1]
                            sba_bb = min([(sbb_b - idx, idx) for idx in sba.l2 if idx < sbb_b])[1]
                            sba_tt = min([(idx - sbb_t, idx) for idx in sba.l2 if sbb_t < idx])[1]
                            assert all(sba_xx > 0 for sba_xx in {sba_rr, sba_tt}), "%s" % [sba_rr, sba_tt]
                            assert all(sba_xx >= 0 for sba_xx in {sba_ll, sba_bb})
                            distance1 = CD((sbb_l, sbb_b), (sba_ll, sba_bb))
                            distance2 = CD((sbb_r, sbb_t), (sba_rr, sba_tt))
                            distance = max(distance1, distance2)
                            gap = distance - 1
                            assert gap >= 0, gap
                            # FIXME, gapMaxMicroInv=0 yields the same result as gapMaxMicroInv=None or -1 ...
                            if gap <= gapMaxMicroInv:
                                # its a micro-inversion perfectly nested
                                # a is splited in one point
                                idsIdentifiedMicroInv.add(idsbb)
                                splitRanks_b = splitRanks_b[1:-1]
                                assert len(splitRanks_b) == 0, "%s" % str((splitRanks_b, len(sbb.la), len(sbb.l1), len(sbb.l2)))
                                id2splitRanks[idsba] = sorted(list(set(id2splitRanks[idsba] + splitRanks_a)))
                                # id2splitRanks[idsbb] = sorted(list(set(id2splitRanks[idsbb] + splitRanks_b)))
                                todo = todo | {idsba}
    print >> sys.stderr, "%s diags identified as micro-inversions" % len(idsIdentifiedMicroInv)

    # With 'todo' and (id of the splited sb, splitRanks)s, calculate the newSbs and removed sbs
    # 4) splits sbs when a gap contains a micro-inversion
    removedSbIds = set()
    newSbIds = set()
    for idsb in todo:
        (c1, c2, _) = id2location(idsb)
        for newSb in splitSbBySplitRanks(id2sb(idsb), id2splitRanks[idsb]):
            newSbId = sbsAndDiagsInPairCompWithIds.addToLocation((c1, c2), newSb)
            newSbIds.add(newSbId)
        removedSbIds = removedSbIds | {idsb}
    sbsAndDiagsInPairCompWithIds.removeIds(removedSbIds)
    assert len((idsSbs - removedSbIds) & newSbIds) == 0
    assert len(idsIdentifiedMicroInv & idsSbs) == 0
    assert len(set(sbsAndDiagsInPairCompWithIds.orderedIds) & removedSbIds) == 0
    idsSbs = (idsSbs - removedSbIds) | newSbIds

    id2sb = sbsAndDiagsInPairCompWithIds.getItemById
    id2location = sbsAndDiagsInPairCompWithIds.getItemLocationById

    # transform diags into sbs for identified microInv
    for idsb in idsIdentifiedMicroInv:
        diag = id2sb(idsb)
        (c1, c2, idx) = id2location(idsb)
        sbsAndDiagsInPairCompWithIds[c1][c2][idx] = SyntenyBlock(diag, 0)
    idsSbs |= idsIdentifiedMicroInv
    idsPutativeMicroInversions -= idsIdentifiedMicroInv
    idsDiagsThatAreNotMicroInv = idsPutativeMicroInversions

    # TODO, translate into a simpler Dict2D
    new_sbsInPairComp = myTools.Dict2d(list)
    diagsThatAreNotSbsInPairComp = myTools.Dict2d(list)
    for (idsb, (c1, c2), sb) in sbsAndDiagsInPairCompWithIds.iterByOrderedIds():
        assert all(isinstance(idx, int) for idx in sb.l1)
        assert all(isinstance(idx, int) for idx in sb.l2)
        if idsb in idsSbs:
            new_sbsInPairComp[c1][c2].append(sb)
        else:
            assert idsb in idsDiagsThatAreNotMicroInv
            diagsThatAreNotSbsInPairComp[c1][c2].append(sb)

    return (new_sbsInPairComp, diagsThatAreNotSbsInPairComp)

def fIdentifyInversionsAtSbsExtremities(sbsInPairComp, putativeMicroInversionsInPairComp,
                                        gapMaxMicroInv=0, identifyMonogInvs=True, hyperSensitive=False):
    # 1) Give unique id to each sb and diag
    idsSbs = myTools.Dict2d(set)
    idsOfPutativeMicroInversions = set()
    sbsAndDiagsInPairCompWithIds = myTools.OrderedDict2dOfLists()
    id = None
    for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
        new_sb = SyntenyBlock(sb)
        id = sbsAndDiagsInPairCompWithIds.addToLocation((c1, c2), new_sb)
        idsSbs[c1][c2].add(id)
    assert id in {None, sbsAndDiagsInPairCompWithIds.maxId}
    if id is None:
        id = -1
    for ((c1, c2), diag) in putativeMicroInversionsInPairComp.iteritems2d():
        id += 1
        new_diag = Diagonal(diag)
        sbsAndDiagsInPairCompWithIds.addToLocationWithId((c1, c2), new_diag, id)
        idsOfPutativeMicroInversions.add(id)
    assert len(set([id for (_, id) in idsSbs.iteritems2d()]) & idsOfPutativeMicroInversions) == 0
    id2sb = sbsAndDiagsInPairCompWithIds.getItemById
    id2location = sbsAndDiagsInPairCompWithIds.getItemLocationById

    # 3) identify micro-inversions at extremities of sbs
    idsIdentifiedMicroInv = set()
    for idsbb in idsOfPutativeMicroInversions:
        sbb = id2sb(idsbb)
        (c1b, c2b, _) = id2location(idsbb)
        for idsba in idsSbs[c1b][c2b]:
            # sbb is included in sba in both genomes
            sba = id2sb(idsba)
            (c1a, c2a, _) = id2location(idsba)
            if (c1a, c2a) == (c1b, c2b):
                if hyperSensitive:
                    if sba.dt == '/' and sbb.dt == '\\' or sba.dt == '\\' and sbb.dt == '/' or sbb.dt == None:
                        if 1 <= distanceBetweenBoundingBoxesOfDiags(sbb, sba, distance=minimumOfDistancesIn1DGenomes) <= gapMaxMicroInv + 1:
                            # if the distance from one infinite line around the bounding box is less than gapMaxMicroInv
                            if len(sbb.la) == 1:
                                if not identifyMonogInvs:
                                    continue
                            idsIdentifiedMicroInv.add(idsbb)
                else:
                # The next commented functioning is less sensible but more specific
                    if sba.dt == '/' and sbb.dt == '\\' or sba.dt == '\\' and sbb.dt == '/':
                        if 1 <= distanceBetweenBoundingBoxesOfDiags(sbb, sba, distance=CD) <= gapMaxMicroInv + 1:
                            # if the distance from one corner is less than gapMaxMicroInv
                            if len(sbb.la) == 1:
                                if not identifyMonogInvs:
                                    continue
                            idsIdentifiedMicroInv.add(idsbb)

    # transform diags into sbs for identified microInv
    for idsb in idsIdentifiedMicroInv:
        diag = id2sb(idsb)
        (c1, c2, idx) = id2location(idsb)
        sbsAndDiagsInPairCompWithIds[c1][c2][idx] = SyntenyBlock(diag, 0)
    idsSbs = set([id for (_, id) in idsSbs.iteritems2d()])
    idsSbs |= idsIdentifiedMicroInv
    idsOfPutativeMicroInversions -= idsIdentifiedMicroInv
    idsDiagsThatAreNotSbs = idsOfPutativeMicroInversions

    # TODO, translate into a simpler Dict2D
    new_sbsInPairComp = myTools.Dict2d(list)
    diagsThatAreNotSbsInPairComp = myTools.Dict2d(list)
    for (idsb, (c1, c2), sb) in sbsAndDiagsInPairCompWithIds.iterByOrderedIds():
        assert all(isinstance(idx, int) for idx in sb.l1)
        assert all(isinstance(idx, int) for idx in sb.l2)
        if idsb in idsSbs:
            new_sbsInPairComp[c1][c2].append(sb)
        else:
            assert idsb in idsDiagsThatAreNotSbs
            diagsThatAreNotSbsInPairComp[c1][c2].append(sb)

    return (new_sbsInPairComp, diagsThatAreNotSbsInPairComp)

def doDistinguishMonoGenicDiags(diagsInPairComp):
    """
    Filter diags of length 1 that might be due to annotation errors or dispersed duplications maybe followed by gaped
    tandem duplications
    :param diagsInPairComp:
    :return:
    """
    monogenicDiagsInPairComp = myTools.Dict2d(list)
    for c1 in diagsInPairComp.keys():
        for c2 in diagsInPairComp[c1].keys():
            lDiags = []
            for diag in diagsInPairComp[c1][c2]:
                (m, max_g, lw1, lw2, l1_min, l1_max, l2_min, l2_max) = diag.calculateCharacteristics()
                # m==1 and (lw1 == 1 or lw2 == 1) corresponds to a dispersed duplication followed by gaped tandem duplications
                if m > 1 and lw1 > 1 and lw2 > 1:
                    lDiags.append(diag)
                else:
                    monogenicDiagsInPairComp[c1][c2].append(diag)
            if len(lDiags) > 0:
                diagsInPairComp[c1][c2] = lDiags
            else:
                del diagsInPairComp[c1][c2]
        if len(diagsInPairComp[c1].keys()) == 0:
            del diagsInPairComp[c1]
    return (diagsInPairComp, monogenicDiagsInPairComp)

@myTools.verbose
def doIdentifyMicroInversions(sbsInPairComp, putativeMicroInvsInPairComp, gapMaxMicroInv, identifyMonogInvs, verbose=False):
    # WARNING: the identification of "monogenic micro-inversion" may in fact hide a tandem duplication with a revert
    # orientation followed by the deletion of the initial gene.
    # Example:
    # 1) initial genome: +A+B+C
    # 2) tandem duplication: +A+B-B+C
    # 3) gene deletion: A-B+C
    # This looks like an inversion but it is not an inversion.
    if identifyMonogInvs:
        print >> sys.stderr, 'WARNING: returned mono-genic micro-inversions may be a hidden tandem duplication, with a ' \
                             'revert orientation, followed by the deletion of the initial gene.'

    print >> sys.stderr, "Nb sbs before identifying micro inversions within sbs gaps = %s" % len(sbsInPairComp.items2d())
    (sbsInPairComp, diagsNotSbsInPairComp) = fIdentifyMonoGcSsAndIdentifyMonoGinvsNestedInSbsGaps(sbsInPairComp,
                                                                                     putativeMicroInvsInPairComp,
                                                                                     gapMaxMicroInv=gapMaxMicroInv,
                                                                                     identifyMonogInvs=identifyMonogInvs)
    print >> sys.stderr, "Nb sbs after identifying micro inversions within sbs gaps = %s" % len(sbsInPairComp.items2d())
    print >> sys.stderr, "Nb sbs before identifying micro inversions at sbs extremities = %s" % len(sbsInPairComp.items2d())

    # TODO identify mono-genic sb when it is close to two projections of sbs-extremities, or two chromosome extremities,
    # or one projection of sb and one chromosome extremity

    # example:
    # G1: ABCDE-K-ZYX----123
    # G2: K-321------ABCDE-ZYK
    (sbsInPairComp, diagsNotSbsInPairComp) = fIdentifyInversionsAtSbsExtremities(sbsInPairComp,
                                                                                 diagsNotSbsInPairComp,
                                                                                 gapMaxMicroInv=gapMaxMicroInv,
                                                                                 identifyMonogInvs=identifyMonogInvs)
    print >> sys.stderr, "Nb sbs after identifying micro inversions at sbs extremities = %s" % len(sbsInPairComp.items2d())
    return (sbsInPairComp, diagsNotSbsInPairComp)

def loopIdentifyBreakpointsWithinGaps(sbsInPairComp):
    # Do this task until no more change (stability!)
    while True:
        nbSbsOld = len(sbsInPairComp.items2d())
        sbsInPairComp = fIdentifyBreakpointsWithinGaps(sbsInPairComp)
        nbSbsNew = len(sbsInPairComp.items2d())
        assert nbSbsNew >= nbSbsOld
        if nbSbsNew == nbSbsOld:
            break
    return sbsInPairComp

def doMergeAllDiags(diagsInPairComp, gapMax, g2_tb, distanceMetric):
    for (c1, c2) in diagsInPairComp.keys2d():
        listOfDiags = diagsInPairComp[c1][c2]
        if len(listOfDiags) > 0:
            diagsInPairComp[c1][c2] = mergeDiags(listOfDiags, gapMax, g2_tb[c2], distanceMetric, verbose=False)
    return diagsInPairComp

def loop_Imi_Ibwg_Om_Mnosbs(sbsInPairComp, diagsInPairCompRejected,
                            gapMaxMicroInv, identifyMonoGenicInvs, distinguishMonoGenicDiags,
                            identifyMicroRearrangements, truncationMax, gapMax,
                            g2_tb, distanceMetric,
                            monogenicDiagsInPairComp=None,
                            loopIterMax=4, mergeAfterSolveOverlaps=True):
    """
    excecute loopIterMax times the pipeline Ibwg|Om|Mnosbs
    with:
    * Ibwg: Identify Breakpoints Within Gaps
    * Om: Overlap Max, the operation that solves overlaps by truncating and removing overlaping sbs
    * Mnosbs: Merge non-overlapping sbs
    The whole process always finishes by Ibwg|Om without Mnosbs
    """

    def isfirstOrLastIter(cptIter):
        if cptIter == 0 or cptIter == loopIterMax:
            return True
        else:
            return False

    cptLoopIter = 0
    while cptLoopIter <= loopIterMax:
        if cptLoopIter  == loopIterMax:
            print >> sys.stderr, "identifyMonoGSbAndMicroInvs|identifyMicroRearrangements|nonOverlappingSbs|mergeNonOverlappingSbs was repeated %s times" % (cptLoopIter - 1)

        # identify micro-inversions
        # TODO: Put that into the loop
        if gapMaxMicroInv is not None:
            if cptLoopIter == 0:
                if identifyMonoGenicInvs and distinguishMonoGenicDiags:
                    assert monogenicDiagsInPairComp is not None
                    putativeMicroInversionsInPairComp = diagsInPairCompRejected + monogenicDiagsInPairComp
                    typeOfMicroInv = ''
                else:
                    typeOfMicroInv = ' non-monogenic'
                    putativeMicroInversionsInPairComp = diagsInPairCompRejected
                    assert not distinguishMonoGenicDiags or all([len(sb.la) > 1 for (_, sb) in putativeMicroInversionsInPairComp.iteritems2d()])
            else:
                putativeMicroInversionsInPairComp = diagsNotSbsInPairComp
            if isfirstOrLastIter(cptLoopIter): print >> sys.stderr, "Nb sbs before identifying%s micro-inversions = %s" % (typeOfMicroInv, len(sbsInPairComp.items2d()))
            (sbsInPairComp, diagsNotSbsInPairComp) = doIdentifyMicroInversions(sbsInPairComp, putativeMicroInversionsInPairComp, gapMaxMicroInv, identifyMonoGenicInvs, verbose=False)
            # remark, doIdentifyMicroInversions may yield mono-genic sbs even if identifyMonoGenicInvs because of the identifyMicroRearrangements that may split a sb in two sub-sbs with one of them of one gene.
            if isfirstOrLastIter(cptLoopIter): print >> sys.stderr, "Nb sbs after identifying%s micro-inversions = %s" % (typeOfMicroInv, len(sbsInPairComp.items2d()))

        # identify breakpoints within gaps
        if identifyMicroRearrangements:
            if isfirstOrLastIter(cptLoopIter): print >> sys.stderr, "Nb sbs before identifying breakpoints within gaps = %s" % len(sbsInPairComp.items2d())
            sbsInPairComp = loopIdentifyBreakpointsWithinGaps(sbsInPairComp)
            if isfirstOrLastIter(cptLoopIter): print >> sys.stderr, "Nb sbs after identifying breakpoints within gaps = %s" % len(sbsInPairComp.items2d())

        # solve overlaps (truncation and removal) + merge non-overlapping sbs
        if truncationMax is not None:
            # solve overlaps (truncation and removal)
            if isfirstOrLastIter(cptLoopIter): print >> sys.stderr, "Nb sbs before solve overlaps = %s" % len(sbsInPairComp.items2d())
            # TODO, compute the variation of the coverage before and after this step
            sbsInPairComp= solveSbsOverlaps(sbsInPairComp, truncationMax=truncationMax, verbose=(True if cptLoopIter == 0 else False))
            if isfirstOrLastIter(cptLoopIter): print >> sys.stderr, "Nb sbs after solve overlaps = %s" % len(sbsInPairComp.items2d())

            # merge non-overlapping sbs
            # Always finish by post processes without a merge !
            if mergeAfterSolveOverlaps and cptLoopIter < loopIterMax:
                if isfirstOrLastIter(cptLoopIter): print >> sys.stderr, "Nb sbs before merging non-overlapping sbs = %s" % len(sbsInPairComp.items2d())
                sbsInPairComp = doMergeAllDiags(sbsInPairComp, gapMax, g2_tb, distanceMetric)
                if isfirstOrLastIter(cptLoopIter): print >> sys.stderr, "Nb sbs after merging non-overlapping sbs = %s" % len(sbsInPairComp.items2d())
        cptLoopIter += 1

    if cptLoopIter > 1 and cptLoopIter < loopIterMax:
        # the first merge non-overlapping sbs was useful, thus other iterations of
        # identifyMicroRearrangements|nonOverlappingSbs|mergeNonOverlappingSbs were performed.
        print >> sys.stderr, "identifyMonoGSbAndMicroInvs|identifyMicroRearrangements|nonOverlappingSbs|mergeNonOverlappingSbs was repeated %s times" % (cptLoopIter - 1)

    return sbsInPairComp

def postProcessDiags(diagsInPairComp, distinguishMonoGenicDiags,
                     pThreshold,  g1_tb, g2_tb, N12s, p_hpSign, validateImpossToCalc_mThreshold,
                     identifyMonoGenicInvs, gapMaxMicroInv, identifyMicroRearrangements, truncationMax,
                     gapMax, distanceMetric):

    # distinguish mono-genic diagonals
    if distinguishMonoGenicDiags:
        print >> sys.stderr, "Total nb diags = %s" % len(diagsInPairComp.items2d())
        (diagsInPairComp, monogenicDiagsInPairComp) = doDistinguishMonoGenicDiags(diagsInPairComp)
        print >> sys.stderr, "Nb non-monogenic diags = %s" % len(diagsInPairComp.items2d())
        print >> sys.stderr, "Nb monogenic diags = %s" % len(monogenicDiagsInPairComp.items2d())
        assert all([len(sb.la) > 1 for (_, sb) in diagsInPairComp.iteritems2d()])

    # statistical validation of diags (putative sbs) into sbs
    typeOfDiags = ' non-monogenic' if distinguishMonoGenicDiags else ''
    if pThreshold is None:
        for (c1, c2) in diagsInPairComp.keys2d():
            diagsInPairComp[c1][c2] = [SyntenyBlock(diag, None) for diag in diagsInPairComp[c1][c2]]
        sbsInPairCompStatVal = diagsInPairComp
        diagsInPairCompRejected = myTools.Dict2d(list)
    else:
        print >> sys.stderr, "Nb%s diags before statistical validation = %s" % (typeOfDiags, len(diagsInPairComp.items2d()))
        (sbsInPairCompStatVal, diagsInPairCompRejected) = statisticalValidation(diagsInPairComp, g1_tb, g2_tb, N12s, p_hpSign,
                                                                                pThreshold=pThreshold,
                                                                                NbOfHomologiesThreshold=50,
                                                                                validateImpossToCalc_mThreshold=validateImpossToCalc_mThreshold,
                                                                                # diagsInPairComp contains no monogenic diag
                                                                                considerMonogenicSb=False,
                                                                                validateSbsWithGapMaxZero=True,
                                                                                verbose=False)
        print >> sys.stderr, "Nb%s diags not validated as sbs = %s" % (typeOfDiags, len(diagsInPairCompRejected.items2d()))
        print >> sys.stderr, "Nb%s diags%s validated as sbs = %s" % (typeOfDiags, ' stat' if pThreshold else '', len(sbsInPairCompStatVal.items2d()))
    sbsInPairComp = sbsInPairCompStatVal

    assert not (identifyMonoGenicInvs and not distinguishMonoGenicDiags), "if 'identifyMonoGenicInvs' is set to True, 'distinguishMonoGenicDiags' should also be set to True"

    totalNbOfHps = sum([len(sb.la) for (_, sb) in sbsInPairComp.iteritems2d()])

    # TODO: identify micro-inversions in the loop (especially identify mono-genic sbs)
    sbsInPairComp = loop_Imi_Ibwg_Om_Mnosbs(sbsInPairComp, diagsInPairCompRejected,
                                            gapMaxMicroInv, identifyMonoGenicInvs, distinguishMonoGenicDiags,
                                            identifyMicroRearrangements, truncationMax, gapMax,
                                            g2_tb, distanceMetric,
                                            monogenicDiagsInPairComp=monogenicDiagsInPairComp,
                                            loopIterMax=4, mergeAfterSolveOverlaps=True)

    # TODO: add an option to identify all one:one orthology (Best reciprocal hit) in mono-genic diagonal, as a mono-genic sb

    return (sbsInPairComp, diagsInPairComp)

def extractSbsInPairCompGenomesInTbs(g1_tb, g2_tb,
                                     gapMax=None,
                                     distanceMetric='CD',
                                     distinguishMonoGenicDiags=True,
                                     pThreshold=None,
                                     gapMaxMicroInv=0,
                                     identifyMonoGenicInvs=False,
                                     identifyMicroRearrangements=True,
                                     truncationMax=None,
                                     sameStrand=True,
                                     nbHpsRecommendedGap=2,
                                     targetProbaRecommendedGap=0.01,
                                     validateImpossToCalc_mThreshold=3,
                                     optimisation=None,
                                     verbose=False):

    # step 2 and 3 : build the MHP and extract strict and consistent diagonals
    #################################################################################
    diagsInPairComp = myTools.Dict2d(list)
    N12s = myTools.Dict2d(int)
    N12_g = 0
    print >> sys.stderr, "synteny block extraction"
    tic = time.time()
    if optimisation == 'cython' and (len(g1_tb.keys()) > 1 or len(g2_tb.keys()) > 1):
        (p_hpSign, p_hpSign_g, N12s, N12_g, diagsInPairComp) = \
            extractDiags.extractDiagsInPairCompChr(g1_tb, g2_tb, sameStrand, distanceMetric)
    else:
        if optimisation == None or (len(g1_tb.keys()) <= 1 and len(g2_tb.keys()) <= 1):
            totalNbComps = len(g1_tb) * len(g2_tb)
            progressBar = myTools.ProgressBar(totalNbComps)
            currCompNb = 0
            for c1 in g1_tb.keys():
                for c2 in g2_tb.keys():
                    (listOfDiags, N12) = extractDiagsInPairCompChr(g1_tb[c1], g2_tb[c2], sameStrand, verbose=verbose)
                    if len(listOfDiags) > 0:
                        diagsInPairComp[c1][c2] = listOfDiags
                    N12s[c1][c2] = N12
                    N12_g += N12
                    currCompNb += 1
                    progressBar.printProgressIn(sys.stderr, currCompNb)
        elif optimisation == 'multiProcess':
            # if the multiprocess option is True and if there is more than one pairwise comparison of chromosomes
            pool = multiprocessing.Pool()
            tasks = [(c1, c2, g1_tb[c1], g2_tb[c2], sameStrand) for (c1, c2) in itertools.product([c1 for c1 in g1_tb], [c2 for c2 in g2_tb])]
            for ((c1, c2), listOfDiags, N12) in pool.map(extractDiagsInPairCompChrWrapper, tasks, chunksize=(len(tasks)/4)):
                if len(listOfDiags) > 0:
                    diagsInPairComp[c1][c2] = listOfDiags
                N12s[c1][c2] = N12
                N12_g += N12
            tac = time.time()
            print >> sys.stderr, "Multiprocessing(extractDiagsInPairCompChr) was executed in %ss" % (tac - tic)
        else:
            raise ValueError('optimisation should be in [\'cython\', \'multiprocess\', None]')

        # second level of verbosity
        verbose2 = False
        (p_hpSign, p_hpSign_g, (sTBG1, sTBG1_g), (sTBG2, sTBG2_g)) =\
                myProbas.statsHpSign(g1_tb, g2_tb, verbose=verbose2)
        print >> sys.stderr, "genome1 tb orientation proba = {+1:%.2f%%,-1:%.2f%%,None:%.2f%%} (stats are also calculated for each chromosome)" % (sTBG1_g[+1]*100, sTBG1_g[-1]*100, sTBG1_g[None]*100)
        print >> sys.stderr, "genome2 tb orientation proba = {+1=%.2f%%,-1:%.2f%%,None:%.2f%%} (stats are also calculated for each chromosome)" % (sTBG2_g[+1]*100, sTBG2_g[-1]*100, sTBG2_g[None]*100)
        print >> sys.stderr, "hp sign proba in the 'global' mhp = {+1:%.2f%%,-1:%.2f%%,None:%.2f%%) (probabilities are also calculated for pairwise mhp)" % (p_hpSign_g[+1]*100, p_hpSign_g[-1]*100, p_hpSign_g[None]*100)

    tac = time.time()
    print >> sys.stderr, "%s(extractDiagsInPairCompChr) was executed in %ss" % (optimisation, (tac - tic))

    print >> sys.stderr, "Nb strict and consistent diags = %s" % len(diagsInPairComp.items2d())
    assert all([len(diagsInPairComp[k1Foo][k2Foo]) > 0 for (k1Foo, k2Foo) in diagsInPairComp.keys2d()])

    # Compute a recommended gapMax parameter value
    ###############################################
    if gapMax is None:
        #Build an average MHP
        #weightedAverage is even better than N50 since it returns a more stable value, not a length of a chromosome of the karyotype, it better reflects the global distribution
        #Waring: if the genome is badly assembled and contain a lot of small contigs, this averaging is not relevant and the recommended gapMax won't be relevant neither
        if len(g1_tb) > 50 or len(g2_tb) > 50:
            print >> sys.stderr, "Warning: one of the two genome seems badly assembled, this may mislead the recommended gapMax calculation"
        N1_weightedAverage = int(myMaths.myStats.getWeightedAverage([len(g1_tb[c1]) for c1 in g1_tb]))
        N2_weightedAverage = int(myMaths.myStats.getWeightedAverage([len(g2_tb[c2]) for c2 in g2_tb]))
        N1_g = sum([len(g1_tb[c1]) for c1 in g1_tb])
        N2_g = sum([len(g2_tb[c2]) for c2 in g2_tb])
        density = float(N12_g)/(N1_g*N2_g)
        # conservation of the density
        N12_weightedAverage = int(density*N1_weightedAverage*N2_weightedAverage)
        gap = recommendedGap(nbHpsRecommendedGap, targetProbaRecommendedGap, N12_weightedAverage, N1_weightedAverage, N2_weightedAverage, p_hpSign=p_hpSign_g, verbose=verbose)
        print >> sys.stderr, "recommended gapMax = %s tbs" % gap
        gapMax = gap
    print >> sys.stderr, "used gapMax = %s tbs" % gapMax

    # merge strict and consistent diagonals into gaped consistent diagonals
    ########################################################################
    print >> sys.stderr, "Nb diags before merging diags = %s" % len(diagsInPairComp.items2d())
    diagsInPairComp = doMergeAllDiags(diagsInPairComp, gapMax, g2_tb, distanceMetric)
    print >> sys.stderr, "Nb diags after merging diags = %s" % len(diagsInPairComp.items2d())

    (sbsInPairComp, diagsInPairComp) = postProcessDiags(diagsInPairComp, distinguishMonoGenicDiags,
                                                        pThreshold,  g1_tb, g2_tb, N12s, p_hpSign, validateImpossToCalc_mThreshold,
                                                        identifyMonoGenicInvs, gapMaxMicroInv, identifyMicroRearrangements, truncationMax,
                                                        gapMax, distanceMetric)
    return sbsInPairComp


def editGenomes(g1, g2, families,
                filterType=FilterType.InBothGenomes,
                labelWith='FamID',
                tandemGapMax=0,
                minChromLength=2,
                keepOriginal=False):
    assert labelWith in {'FamID', 'FamName'}
    nCini1 = len(g1.keys())
    nCini2 = len(g2.keys())
    nGini1 = sum([len(chrom1) for chrom1 in g1.values()])
    nGini2 = sum([len(chrom2) for chrom2 in g2.values()])
    #step 1 :filter genomes and rewrite in tandem blocks if needed
    ##############################################################
    # rewrite genomes by family names (ie ancGene names)
    if labelWith == 'FamID':
        g1_fID = myMapping.labelWithFamID(g1, families)
        g2_fID = myMapping.labelWithFamID(g2, families)
    else:
        assert labelWith == 'FamName'
        g1_fID = myMapping.labelWithFamNames(g1, families)
        g2_fID = myMapping.labelWithFamNames(g2, families)
    # genes that are not in ancGene have a aID=None
    nGiniInFam1 = len([fID for chrom1 in g1_fID.values() for (fID, _) in chrom1 if fID is not None])
    nGiniInFam2 = len([fID for chrom2 in g2_fID.values() for (fID, _) in chrom2 if fID is not None])
    print >> sys.stderr, "genome1 initially contains %s chromosomes" % nCini1
    print >> sys.stderr, "genome2 initially contains %s chromosomes" % nCini2
    print >> sys.stderr, "genome1 initially contains %s genes (%s genes are in families, %.2f%%)" % (nGini1, nGiniInFam1, (100 * float(nGiniInFam1) / float(nGini1)))
    print >> sys.stderr, "genome2 initially contains %s genes (%s genes are in families, %.2f%%)" % (nGini2, nGiniInFam2, (100 * float(nGiniInFam2) / float(nGini2)))
    # Must be applied on the two genomes, because of the mode inBothGenomes (InFamilies => not only anchor genes are kept but all genes herited from a gene of the LCA)
    #mfilt2origin1 -> mGf2Go1
    ((g1_fID, mGf2Go1, (nCL1, nGL1)), (g2_fID, mGf2Go2, (nCL2, nGL2))) =\
        filter2D(g1_fID, g2_fID, filterType, minChromLength, keepOriginal=keepOriginal)
    print >> sys.stderr, "genome1 after filterType=%s and minChromLength=%s contains %s genes" %\
        (filterType, minChromLength, sum([len(g1_fID[c1]) for c1 in g1_fID]))
    print >> sys.stderr, "genome2 after filterType=%s and minChromLength=%s contains %s genes" %\
        (filterType, minChromLength, sum([len(g2_fID[c2]) for c2 in g2_fID]))
    nGD1 = myMapping.nbDup(g1_fID)[0]
    nGD2 = myMapping.nbDup(g2_fID)[0]
    (g1_tb, mtb2g1, nGTD1) = myMapping.remapRewriteInTb(g1_fID, tandemGapMax=tandemGapMax, mOld=mGf2Go1)
    (g2_tb, mtb2g2, nGTD2) = myMapping.remapRewriteInTb(g2_fID, tandemGapMax=tandemGapMax, mOld=mGf2Go2)
    print >> sys.stderr, "genome1 rewritten in tbs, contains %s tbs" % sum([len(g1_tb[c1]) for c1 in g1_tb])
    print >> sys.stderr, "genome2 rewritten in tbs, contains %s tbs" % sum([len(g2_tb[c2]) for c2 in g2_tb])
    nDD1 = myMapping.nbDup(g1_tb)[0]
    nDD2 = myMapping.nbDup(g2_tb)[0]
    print >> sys.stderr, "genome1 contains %s gene duplicates (initial gene excluded)" % nGD1
    print >> sys.stderr, "genome1 contains %s tandem duplicated genes (initial gene excluded)" % nGTD1
    print >> sys.stderr, "genome1 contains %s dispersed duplicated tbs (initial tb excluded)" % nDD1
    print >> sys.stderr, "genome2 contains %s gene duplicates (initial gene excluded)" % nGD2
    print >> sys.stderr, "genome2 contains %s tandem duplicated genes (initial gene excluded)" % nGTD2
    print >> sys.stderr, "genome2 contains %s dispersed duplicated tbs (initial tb excluded)" % nDD2
    assert nDD1 + nGTD1 == nGD1
    assert nDD2 + nGTD2 == nGD2

    # conservation law genes
    def nbOfGenesInAGenomeInTbs(g_tb, mtb2g):
        nbGenes = 0
        for (chr, chrom) in g_tb.iteritems():
            nbGenes += sum(len(mtb2g[chr][itb]) for (itb, _) in enumerate(chrom))
        return nbGenes
    assert nGini1 == nGL1 + nbOfGenesInAGenomeInTbs(g1_tb, mtb2g1)
    assert nGini2 == nGL2 + nbOfGenesInAGenomeInTbs(g2_tb, mtb2g2)
    # conservation law chromosomes
    assert nCini1 == nCL1 + len(g1_tb.keys())
    assert nCini2 == nCL2 + len(g2_tb.keys())

    return ((g1_tb, mtb2g1, (nCL1, nGL1)), (g2_tb, mtb2g2, (nCL2, nGL2)))



# Complete procedure to compute synteny blocks (sbs) between 2 genomes, using homology relationships contained in ancGenes
# Inputs:
#       g1, g2 : myGenomes.Genome (or dict : {...c:[..., (g,s), ...]}) of the two compared species
#       ancGenes : myGenomes.Genome define the homology relationships (usually the ancGenes of the LCA of the two compared species)
#       sameStrand :     True  => gene transcription orientations must be consistent with diagonal types ('/' or '\')
#                               False => Do not take care of gene orientation
#       gapMax : distance (in tandemblocks) allowed between two diagonals to be merged together. This allows to go over wrong  annotations. But it can also introduce some Fp
#       minChromLength is the minimal length of the chromosome to be considered
#       distanceMetric : the distance metric (either MD,DPD,ED or CD)
#       pThreshold : the probability threshold under which the sb is considered significant
#       FilterType
#               InFamilies : all genes herited from a gene of the LCA are kept during the filtering of extant genomes
#               InBothGenomes : only 'anchor genes' (genes present in both genomes) are kept
#               None : genomes are not filtered
# Outputs:
#       synteny blocks are yielded as (strand,(c1,l1),(c2,l2),la)
#       c1 : chromosomes of genome 'g1' where there is a synteny block corresponding to the diagonal
#       l1 : the synteny block on genome 'g1'
#       l1 = [..., (i,s), ...] with 'i' the index of the gene and 's' its strand in the genome1 without species specific genes
#       la = [..., (ancGeneName, ancGeneOrientation, tb width, tb height), ...]]
#           tb -> tandem block associated with the ancGene in this synteny block
#           width : length in tandem duplicates on the 1st genome
#           height: length in tandem duplicates on ths 2nd genome
#########################################################################################################################
@myTools.tictac
@myTools.verbose
def extractSbsInPairCompGenomes(g1, g2, families,
                                filterType=FilterType.InBothGenomes,
                                tandemGapMax=0,
                                gapMax=None,
                                distanceMetric='CD',
                                distinguishMonoGenicDiags=True,
                                pThreshold=None,
                                gapMaxMicroInv=0,
                                identifyMonoGenicInvs=True,
                                identifyMicroRearrangements=True,
                                truncationMax=None,
                                sameStrand=True,
                                minChromLength=2,
                                nbHpsRecommendedGap=2,
                                targetProbaRecommendedGap=0.01,
                                validateImpossToCalc_mThreshold=3,
                                optimisation=None,
                                verbose=False):
    isinstance(g1, myLightGenomes.LightGenome)
    isinstance(g2, myLightGenomes.LightGenome)
    isinstance(families, myLightGenomes.Families)

    if isinstance(g1, myGenomes.Genome) and isinstance(g2, myGenomes.Genome):
        g1 = g1.intoDict()
        g2 = g2.intoDict()
    elif isinstance(g1, dict) and isinstance(g2, dict):
        pass
    else:
        raise TypeError('g1 and/or g2 must be either myGenomes.Genome or dict')

    ((g1_tb, mtb2g1, (nCL1, nGL1)), (g2_tb, mtb2g2, (nCL2, nGL2))) = editGenomes(g1, g2, families,
                                                                                 filterType=filterType,
                                                                                 labelWith='FamID',
                                                                                 tandemGapMax=tandemGapMax,
                                                                                 minChromLength=minChromLength,
                                                                                 keepOriginal=False)

    sbsInPairComp = extractSbsInPairCompGenomesInTbs(g1_tb, g2_tb,
                                                     gapMax=gapMax,
                                                     distanceMetric=distanceMetric,
                                                     distinguishMonoGenicDiags=distinguishMonoGenicDiags,
                                                     pThreshold=pThreshold,
                                                     identifyMonoGenicInvs=identifyMonoGenicInvs,
                                                     gapMaxMicroInv=gapMaxMicroInv,
                                                     identifyMicroRearrangements=identifyMicroRearrangements,
                                                     truncationMax=truncationMax,
                                                     sameStrand=sameStrand,
                                                     nbHpsRecommendedGap=nbHpsRecommendedGap,
                                                     targetProbaRecommendedGap=targetProbaRecommendedGap,
                                                     validateImpossToCalc_mThreshold=validateImpossToCalc_mThreshold,
                                                     optimisation=optimisation,
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
            assert isinstance(indx_tb_g1, int), indx_tb_g1
            l1.append([gIndxs for gIndxs in mtb2gc1.new[indx_tb_g1]])
            indx_tb_g2 = sb.l2[indx_HP]
            l2.append([gIdxs for gIdxs in mtb2gc2.new[indx_tb_g2]])

            # FIXME la.append((ancGenes.lstGenes[None][aGene[0]].names[0], aGene[1], aGene[2]))
            la.append((families.getFamNameByID(aGene[0]), aGene[1], aGene[2]))

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
def printSbsFile(sbsInPairComp, genome1, genome2, sortByDecrLengths=True, stream=sys.stdout):
    print >> sys.stderr, "Print synteny blocks"

    def foo(genomeX, cX, lX, idxHp, reverseOrder=False):
        assert reverseOrder in {True, False, None}, reverseOrder
        tb = lX[idxHp] if (reverseOrder in {+1, None}) else list(reversed(lX[idxHp]))
        gXs = [genomeX[cX][gIdx].n for gIdx in tb]
        gXs = ' '.join(gXs)
        sXs = [genomeX[cX][gIdx].s for gIdx in tb]
        for (i, s) in enumerate(sXs):
            if s == +1:
                sXs[i] = '+'
            elif s == -1:
                sXs[i] = '-'
            else:
                raise
        sXs = ''.join(sXs)
        return (gXs, sXs)

    if isinstance(sbsInPairComp, myTools.Dict2d):
        # listOfSbs = [ ..., ((c1, c2), sb), ...]
        listOfSbs = sbsInPairComp.intoList()
    else:
        assert isinstance(sbsInPairComp, myTools.OrderedDict2dOfLists)
        # listOfSbs = [ ..., ((c1, c2), sb, id), ...]
        listOfSbs = sbsInPairComp.intoList()

    if sortByDecrLengths:
        listOfSbs.sort(key=lambda x: len(x[1].la), reverse=True)

    statsSbs = []
    if isinstance(sbsInPairComp, myTools.OrderedDict2dOfLists):
        generatorSbs = ((id, (c1, c2), sb) for ((c1, c2), sb, id) in listOfSbs)
    else:
        assert isinstance(sbsInPairComp, myTools.Dict2d)
        generatorSbs = ((idSb, (c1, c2), sb) for (idSb, ((c1, c2), sb)) in enumerate(listOfSbs))

    for (idSb, (c1, c2), sb) in generatorSbs:
        assert len(sb.la) == len(sb.l1) == len(sb.l2), "len(l1)=%s, len(l2)=%s, len(la)=%s\nl1=%s\nl2=%s\nla=%s" % (len(sb.l1), len(sb.l2), len(sb.la), sb.l1, sb.l2, sb.la)
        nbHps = len(sb.la)
        statsSbs.append(nbHps)
        for (idxHp, (aGname, aGstrand, dist)) in enumerate(sb.la):
            (g1s, s1s) = foo(genome1, c1, sb.l1, idxHp)
            reverseOrderOnG2 = True if sb.dt in {'/' or None} else False
            (g2s, s2s) = foo(genome2, c2, sb.l2, idxHp, reverseOrder=reverseOrderOnG2)
            print >> stream, myFile.myTSV.printLine([idSb, aGname, aGstrand, dist, c1, c2, s1s, s2s, g1s, g2s])

    print >> sys.stderr, "Distribution of the lengths of synteny blocks:", myMaths.myStats.syntheticTxtSummary(statsSbs)

def findDiagTypeFromFileInfos(l1, l2, s1s=None, s2s=None):
    assert len(l1) == len(l2)
    # single gene diagonal
    if len(l1) == 1:
        assert len(l2) == 1 and s1s is not None and s2s is not None
        assert len(l1[0]) == len(s1s), '%s %s' % (len(l1[0]), len(s1s))
        # find orientations of genes
        if all([s == +1 for s in s1s]):
            s1 = +1
        elif all([s == -1 for s in s1s]):
            s1 = -1
        else:
            s1 = None
        if all([s == +1 for s in s2s]):
            s2 = +1
        elif all([s == -1 for s in s2s]):
            s2 = -1
        else:
            s2 = None
        if s1 == +1:
            if s2 == +1 or s2 == None:
                diagType = '/'
            elif s2 == -1:
                diagType = '\\'
            else:
                diagType = None
        elif s1 == -1:
            if s2 == -1 or s2 == None:
                diagType = '/'
            elif s2 == +1:
                diagType = '\\'
            else:
                diagType = None
        else:
            diagType = None
    else:
        # FIXME another idea instead of using 'min'
        # if myMaths.mean(l1[0]) < myMaths.mean(l1[-1]):
        # but this should be coherent with the rewriting in tbs

        # since the rewriting in tbs use the localisation of the first gene
        if min(l1[0]) < min(l1[-1]):
            if min(l2[0]) < min(l2[-1]):
                diagType = '/'
            elif min(l2[0]) > min(l2[-1]):
                diagType = '\\'
            else:
                # horizontal synteny block
                diagType = None
        else:
            # vertical synteny block
            diagType = None
    return diagType


#################################
# Parser function
#################################
# adding genome1 and genome2 allow to recover l1 and l2, the lists of homologous
# tbs on genome1 and genome2.
def parseSbsFile(fileName, genome1=None, genome2=None, withIds=False):
    assert isinstance(genome1, myLightGenomes.LightGenome)
    assert isinstance(genome2, myLightGenomes.LightGenome)

    def foo(gXs, sXs, genome, chr):
        if not genome.withDict:
            genome.computeDictG2P()
        if genome is None:
            return ([None], [None])
        gXs = gXs.split(' ')
        nbTandemDup = len(sXs)
        assert len(sXs) == len(gXs) == nbTandemDup, "sXs=%s, gXs=%s and nbTandemDup=%s" % (sXs, gXs, nbTandemDup)
        new_sXs = []
        new_gXs = []
        for i in range(nbTandemDup):
            gN = gXs[i]
            genePos = genome.getPosition(gN, default=None)
            if genePos is None:
                raise ValueError("gene %s is not in %s genome" % (gN, genome.name))
            chromosome = genePos.c
            assert chr == chromosome, "chr = %s(%s) and chromosome = %s(%s)" % (chr, type(chr), chromosome, type(chromosome))
            index = genePos.idx
            #gXs[i] = genome.lstGenes[chromosome][index].names[0]
            new_gXs.append(index)
            s = sXs[i]
            if s == '+':
                # insert at index 'i'
                assert len(new_sXs) == i
                new_sXs.append(+1)
            elif s == '-':
                # insert at index 'i'
                assert len(new_sXs) == i
                new_sXs.append(-1)
            else:
                raise
            #FIXME
            #assert genome.lstGenes[chromosome][index].strand == new_sXs[i], "len(genome.lstGenes[chr])=%s, index=%s" % (len(genome.lstGenes[chr]), gXs[i].index)
        assert len(new_gXs) == len(new_sXs)
        assert len(new_gXs) == len(gXs) and len(new_sXs) == len(sXs)
        return (new_gXs, new_sXs)

    sbsReader =\
        myFile.myTSV.readTabular(fileName,
                                 [int, str, str, int, str, str, str, str, str, str],
                                 delim='\t')
    if withIds:
        sbsInPairCompWithIds = myTools.OrderedDict2dOfLists()
    else:
        sbsInPairComp = myTools.Dict2d(list)
    pVal = 0
    idSb_old = 'foo'
    c1_old = 'fooC1'
    c2_old = 'fooC2'
    for (idSb, aGname, aGstrand, dist, c1, c2, s1s, s2s, g1s, g2s) in sbsReader:
        # iteration over each line of the input file that corresponds to tandem blocks
        aGstrand = int(aGstrand) if aGstrand != 'None' else None
        dist = int(dist) if dist != 'None' else None
        assert genome1
        (g1s, s1s) = foo(g1s, s1s, genome1, c1)
        (g2s, s2s) = foo(g2s, s2s, genome2, c2)
        if idSb != idSb_old:
            if idSb_old == 'foo':
                # first idSb
                l1 = [g1s]
                l2 = [g2s]
                la = [(aGname, aGstrand, dist)]
                s1s_old = s1s
                s2s_old = s2s
            else:
                # record the former synteny block that has been parsed

                #assert l1[0] <= l1[-1], "%s <= %s" % (l1[0], l1[-1])
                if len(la) == 1:
                    diagType = findDiagTypeFromFileInfos(l1, l2, s1s_old, s2s_old)
                else:
                    diagType = findDiagTypeFromFileInfos(l1, l2)
                #if c1_old == '10' and c2_old == '14':
                newSb = SyntenyBlock(Diagonal(diagType, l1, l2, la), pVal)
                if withIds:
                    sbsInPairCompWithIds.addToLocationWithId((c1_old, c2_old), newSb, idSb_old)
                else:
                    sbsInPairComp[c1_old][c2_old].append(newSb)
                # start a new sb
                l1 = [g1s]
                l2 = [g2s]
                la = [(aGname, aGstrand, dist)]
                s1s_old = s1s
                s2s_old = s2s
            idSb_old = idSb
            c1_old = c1
            c2_old = c2
        else:
            l1.append(g1s)
            l2.append(g2s)
            la.append((aGname, aGstrand, dist))

    # last idSb
    # find the diagType
    if len(la) == 1:
        diagType = findDiagTypeFromFileInfos(l1, l2, s1s_old, s2s_old)
    else:
        diagType = findDiagTypeFromFileInfos(l1, l2)
    newSb = SyntenyBlock(Diagonal(diagType, l1, l2, la), pVal)
    if withIds:
        sbsInPairCompWithIds.addToLocationWithId((c1, c2), newSb, idSb)
        return sbsInPairCompWithIds
    else:
        sbsInPairComp[c1][c2].append(newSb)
        return sbsInPairComp


###########################################
# Post processing of synteny blocks (sbs)
###########################################

# Build a genome with sbs as contigs and only conserve one ancestral gene by hp
def buildGenomeFromSbs(sbsInPairComp, sbLengthThreshold=None):
    cptSb = 0
    ancSbGenome = myLightGenomes.LightGenome()
    for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
        old_ancGene = None
        for (ihp, (ancGeneName, ancStrand, dist)) in enumerate(sb.la):
            # Happens when there are paralogous hps in the same diagonal
            # assert ancGene != old_ancGene, "%s,%s" % (ancGene, old_ancGene)
            if ancGeneName == old_ancGene:
                # Be waware that this may create diags of size 1
                continue
            # ancGene synonymous of hp-name
            # ancSbGenome[cptSb].append((ancGeneName, ancStrand, dist))
            ancSbGenome[cptSb].append(OGene(ancGeneName, ancStrand))
            old_ancGene = ancGeneName
        # # In every cases sbs of length 1 are removed
        # if len(ancSbGenome[cptSb]) <= 1:
        #     del ancSbGenome[cptSb]
        #     continue
        if sbLengthThreshold != None and len(ancSbGenome[cptSb]) <= sbLengthThreshold:
            del ancSbGenome[cptSb]
            continue
        else:
            cptSb = cptSb + 1
    return ancSbGenome

def lengthProjectionOnGenome(rankGenome, sb, c, genome, correctLens=True, meanInterGeneLen=None):
    if meanInterGeneLen is None:
        meanInterGeneLen = 0.0
    assert isinstance(sb, SyntenyBlock)
    assert isinstance(genome, myGenomes.Genome)
    lGeneCoord = []
    sblX = sb.l1 if rankGenome == 1 else sb.l2
    for tbX in sblX:
        for gIdx in tbX:
            #g1Pos = genome1.getPosition([g1n]).pop()
            #chromosome = g1Pos.chromosome
            #assert c1 == chromosome
            #index = g1Pos.index
            c = myGenomes.commonChrName(c)
            g = genome.lstGenes[c][gIdx]
            assert isinstance(g, myGenomes.Gene)
            # beg = 5' extremity end=3' extremity
            # assert (g.strand in {+1, None} and g.beginning < g.end) or (g.strand == -1 and g.end < g.beginning)
            lGeneCoord.extend([g.beginning, g.end])
    minOnG = min(lGeneCoord)
    maxOnG = max(lGeneCoord)
    lengthG = maxOnG - minOnG
    if correctLens:
        # see Nadeau & Taylor 1984 815-816
        n = len(sblX)
        r = lengthG
        if n >= 2:
            # this consider the local density of markers
            m = r * float(n+1)/float(n-1)
        elif n == 1:
            # add for each extremity (twice), the mean intergene length / 2
            m = r + 2 * float(meanInterGeneLen) / 2
        else:
            raise ValueError('A sb should at least contain one marker')
        lengthG = m
    return lengthG

# TODO
# def computeMeanInterGeneLenInChrs(sbsInPairCompWithIds, )
#     # preprocess data
#     # record all mapped homologs
#     setOfMappedHomologs1 = set()
#     setOfMappedHomologs2 = set()
#     for (id, (c1, c2), sb) in _sbsInPairCompWithIds.iterByOrderedIds():
#         sblX = sb.l1
#         for tbX in sblX:
#             for gIdx in tbX:
#                 g = genome1L[c1][gIdx]
#                 setOfMappedHomologs1.add(g.n)
#         sblX = sb.l2
#         for tbX in sblX:
#             for gIdx in tbX:
#                 g = genome2L[c2][gIdx]
#                 setOfMappedHomologs2.add(g.n)
#
#     (meanInterHomologLenInChr1, densityInG1) = genome1.computeMeanInterGeneLen(setOfMappedHomologs1)
#     (meanInterHomologLenInChr2, densityInG2) = genome2.computeMeanInterGeneLen(setOfMappedHomologs2)
#     return meanInterGeneLenInChrs

# genome1 and genome2 must be myGenomes.Genomes
def getSbsMeanLengths(genome1, genome2, sbsInPairComp, lengthUnit='Mb', meanInterHomologsLenInChrs=None):
    assert lengthUnit in ['Mb', 'gene']
    assert isinstance(genome1, myGenomes.Genome)
    assert isinstance(genome2, myGenomes.Genome)
    if meanInterHomologsLenInChrs:
        assert len(meanInterHomologsLenInChrs) == 2
        (meanInterGeneLenInChr1, meanInterGeneLenInChr2) = meanInterHomologsLenInChrs

    if isinstance(sbsInPairComp, myTools.OrderedDict2dOfLists):
        _sbsInPairCompWithIds = sbsInPairComp
    else:
        assert isinstance(sbsInPairComp, myTools.Dict2d)
        _sbsInPairCompWithIds = myTools.OrderedDict2dOfLists()
        for ((c1, c2), sb) in sbsInPairComp.iteritems2d():
            _sbsInPairCompWithIds.addToLocation((c1, c2), sb)

    # compute sb lengths
    sbId2Length = {}
    if lengthUnit == 'gene':
        for (id, (c1, c2), sb) in _sbsInPairCompWithIds.iterByOrderedIds():
            nbAncGenes = len(sb.la)
            sbId2Length[id] = nbAncGenes
    elif lengthUnit == 'Mb':
        for (id, (c1, c2), sb) in _sbsInPairCompWithIds.iterByOrderedIds():
            if meanInterHomologsLenInChrs:
                lengthG1 = lengthProjectionOnGenome(1, sb, c1, genome1, correctLens=True, meanInterGeneLen=meanInterGeneLenInChr1[c1])
                lengthG2 = lengthProjectionOnGenome(2, sb, c2, genome2, correctLens=True, meanInterGeneLen=meanInterGeneLenInChr2[c2])
            else:
                lengthG1 = lengthProjectionOnGenome(1, sb, c1, genome1, correctLens=True)
                lengthG2 = lengthProjectionOnGenome(2, sb, c2, genome2, correctLens=True)
            averageSbLength = float(lengthG1 + lengthG2) / 2.0
            # in megabases
            averageSbLength = averageSbLength / 1000000
            sbId2Length[id] = averageSbLength
    return sbId2Length

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
