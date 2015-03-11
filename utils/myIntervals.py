# -*- coding: utf-8 -*-
# LibsDyogen
# python 2.7
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

import collections
from utils.myLightGenomes import OGene, LightGenome

# Region
# n = gene name
# ht = 'h' if the region is close to the 'h'ead of the gene else ht = 't' the region is close to the 't'ail of the gene
#   ht = None if the gene has no orientation
Region = collections.namedtuple("Region", ['n', 'ht'])

# returns the region either on the left (-1) or the right (+1), depending on leftOrRight, of the oriented gene 'og'
def regionFromGene(og, leftOrRight):
    # leftOrRight = +1 means that we want the region on the right of og
    assert leftOrRight in {+1, -1}
    assert isinstance(og, OGene)
    if og.s == +1:
        return Region(og.n, 'h' if leftOrRight == +1 else 't')
    elif og.s == -1:
        return Region(og.n, 't' if leftOrRight == +1 else 'h')
    else:
        assert og.s == None
        raise ValueError('Since og has no orientation it is impossible to find the region')

class Adjacency(tuple):
    def __new__(cls, g1n, g2n):
        assert isinstance(g1n, str) and isinstance(g2n, str)
        assert g1n != g2n
        # this relationship of order works for if both g1n and g2n are strings as well as if they are both integers
        if g1n < g2n:
            return tuple.__new__(cls, (g1n, g2n))
        else:
            return tuple.__new__(cls, (g2n, g1n))
    def reversed(self):
        (g1n, g2n) = self
        # we do not want this function to return an OAdjacency Object or to change the current object since all
        # adjacencies are recorded in only one fashion using the relation og1.n < og2.n.
        # It could be miscunfusing to return Adjacencies that do not ensure that og1.n < og2.n.
        return (g2n, g1n)
    def __repr__(self):
        #return self.__class__.__name__ + tuple.__repr__(self)
        return 'Adj' + tuple.__repr__(self)

# class used to have only one repentant of one adjacency given the two gene names
# example:
# g1n = 'G1'
# g2n = 'G2'
# #(works also if g1n and g2n are both integers)
# g1s = +1
# g2s = -1
# adj1 = OAdjacency((g1n, g1s), (g2n, g2s))
# adj2 = OAdjacency((g2n, -g2s), (g1n, -g1s))
# adj1 == adj2
# > return True
class OAdjacency(tuple):
    # cannot use __init__ since tuple is an immutable class
    # (cf http://stackoverflow.com/questions/1565374/subclassing-python-tuple-with-multiple-init-arguments)
    def __new__(cls, og1, og2):
        if not isinstance(og1, OGene):
            og1 = OGene(og1)
        if not isinstance(og2, OGene):
            og2 = OGene(og2)
        assert og1.n != og2.n
        # this relationship of order works for if both g1n and g2n are strings as well as if they are both integers
        if og1.n < og2.n:
            return tuple.__new__(cls, (og1, og2))
        else:
            return tuple.__new__(cls, (og2.reversed(), og1.reversed()))
    def reversed(self):
        (og1, og2) = self
        # we do not want this function to return an OAdjacency Object or to change the current object since all
        # adjacencies are recorded in only one fashion using the relation og1.n < og2.n.
        # It could be miscunfusing to return Adjacencies that do not ensure that og1.n < og2.n.
        return ((og2.reversed(), og1.reversed()))
    def __repr__(self):
        #return self.__class__.__name__ + tuple.__repr__(self)
        return 'Adj' + tuple.__repr__(self)


def findConsideredFlankingGenesIdxs(chrom, x, ignoredGeneNames=set()):
    # index of the gene at the left of the intergene
    idxGl = x - 1
    # index of the gene at the right of the intergene
    idxGr = x
    while 0 <= idxGl and chrom[idxGl].n in ignoredGeneNames:
        idxGl = idxGl - 1
    while idxGr < len(chrom) and chrom[idxGr].n in ignoredGeneNames:
        idxGr = idxGr + 1
    if idxGl < 0:
        # if the breakpoint has no considered gene on its left (it is at the chrom extremity of considered genes)
        idxGl = None
    if len(chrom) <= idxGr:
        # if the breakpoint has no considered gene on its right (it is at the chrom extremity of considered genes)
        idxGr = None
    return (idxGl, idxGr)


def regionsFromAdjacency(adj):
    assert isinstance(adj, OAdjacency)
    (og1, og2) = adj
    assert og1.s in {+1, -1} and og2.s in {+1, -1}
    # 'Head or Tail of the 1st gene', if it is equal to 'h' it means that the breakpoint occurred near
    # the head of the 1st gene
    g1ht = 'h' if og1.s == +1 else 't'
    # 'Head or Tail of the 2nd gene', (same principle as above)
    g2ht = 't' if og2.s == +1 else 'h'
    return (Region(og1.n, g1ht), Region(og2.n, g2ht))

def intergeneFromRegion(region, genomeWithDict):
    assert isinstance(genomeWithDict, LightGenome)
    assert isinstance(region, Region), region
    assert genomeWithDict.withDict
    (gn, ht) = region
    assert isinstance(gn, int) or isinstance(gn, str)
    assert ht in ['h', 't']
    # localisation of gene gn: 'c' is the chromosome and 'idxG' is the index of the gene on 'c'
    gpos = genomeWithDict.getPosition(gn)
    if gpos is not None:
        (c, idxG) = gpos
    else:
        raise ValueError("Gene %s is not in the genome" % gn)
    chrom = genomeWithDict[c]
    g = chrom[idxG]
    assert gn == g.n
    if g.s == +1:
        # the breakpoint location on 'chrom'
        x = idxG if ht == 't' else idxG + 1
    else:
        assert g.s == -1
        # the breakpoint location on 'chrom'
        x = idxG if ht == 'h' else idxG + 1
    return (c, x)

def intergeneFromAdjacency(adjacency, genomeWithDict, default=None):
    assert isinstance(adjacency, OAdjacency)
    assert isinstance(genomeWithDict, LightGenome)
    assert genomeWithDict.withDict
    (og1, og2) = adjacency
    g1pos = genomeWithDict.getPosition(og1.n, default=None)
    grpos = genomeWithDict.getPosition(og2.n, default=None)
    if g1pos is None or grpos is None:
        raise ValueError('At least one gene of the adjacency is not in genomeWithDict')
    if g1pos.c != grpos.c:
        return default
    assert g1pos.c in genomeWithDict and grpos.c in genomeWithDict
    # The two genes are located on the same chromosome
    if abs(g1pos.idx - grpos.idx) != 1:
        return default
    # The two genes are direct neighbours
    if g1pos.idx < grpos.idx:
        (ogl, ogr) = (og1, og2)
        # localisation of the intergene
        (c, x) = (g1pos.c, g1pos.idx + 1)
        # FIXME
        # check orientations!
        # assert (ogl, ogr) == (genomeWithDict[g1pos.c][g1pos.idx], genomeWithDict[grpos.c][grpos.idx]),\
        #     "%s == %s" % ((ogl, ogr), (genomeWithDict[g1pos.c][g1pos.idx], genomeWithDict[grpos.c][grpos.idx]))
    else:
        (ogl, ogr) = (og2.reversed(), og1.reversed())
        # localisation of the intergene
        (c, x) = (grpos.c, grpos.idx + 1)
        # FIXME
        # assert (ogl, ogr) == (genomeWithDict[grpos.c][grpos.idx], genomeWithDict[g1pos.c][g1pos.idx]),\
        #     "%s == %s" % ((ogl, ogr), (genomeWithDict[grpos.c][grpos.idx], genomeWithDict[g1pos.c][g1pos.idx]))
    return (c, x)


# TODO
# def adjacencyFromIntergene(intergene, genomeWithDict, default=None):
#     assert isinstance(genomeWithDict, LightGenome)
