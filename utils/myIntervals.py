# -*- coding: utf-8 -*-
# LibsDyogen
# python 2.7
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

import collections
from myLightGenomes import OGene, LightGenome
import myTools
import sys

# GeneExtremity
# n = gene name
# ht = 'h' if the gene extremity is close to the 'h'ead of the gene else ht = 't' the gene extremity is close to the 't'ail of the gene
#   ht = None if the gene has no orientation
GeneExtremity = collections.namedtuple("GeneExtremity", ['n', 'ht'])


# returns the gene extremity either on the left (-1) or the right (+1), depending on leftOrRight, of the oriented gene 'og'
def geneExtremityFromGene(og, leftOrRight):
    # leftOrRight = +1 means that we want the gene extremity on the right of og
    assert leftOrRight in {+1, -1, None}
    assert isinstance(og, OGene)
    if leftOrRight is None:
        geneExtr = GeneExtremity(og.n, None)
    else:
        if og.s == +1:
            geneExtr = GeneExtremity(og.n, 'h' if leftOrRight == +1 else 't')
        elif og.s == -1:
            geneExtr = GeneExtremity(og.n, 't' if leftOrRight == +1 else 'h')
        else:
            assert og.s == None
            geneExtr = GeneExtremity(og.n, None)
    return geneExtr

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


def geneExtremitiesFromAdjacency(adj):
    assert isinstance(adj, OAdjacency)
    (og1, og2) = adj
    assert og1.s in {+1, -1} and og2.s in {+1, -1}
    # 'Head or Tail of the 1st gene', if it is equal to 'h' it means that the breakpoint occurred near
    # the head of the 1st gene
    g1ht = 'h' if og1.s == +1 else 't'
    # 'Head or Tail of the 2nd gene', (same principle as above)
    g2ht = 't' if og2.s == +1 else 'h'
    return (GeneExtremity(og1.n, g1ht), GeneExtremity(og2.n, g2ht))

def intergeneFromGeneExtremity(geneExtr, genomeWithDict):
    assert isinstance(genomeWithDict, LightGenome)
    assert isinstance(geneExtr, GeneExtremity), geneExtr
    assert genomeWithDict.withDict
    (gn, ht) = geneExtr
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
    else:
        (ogl, ogr) = (og2.reversed(), og1.reversed())
        # localisation of the intergene
        (c, x) = (grpos.c, grpos.idx + 1)
    return (c, x)


def analyseGenomeIntoChromExtremities(genome, oriented=True):
    isinstance(genome, LightGenome)
    res = set()
    for chrom in genome.values():
        if oriented:
            chromGeneExtrLeft = geneExtremityFromGene(chrom[0], -1)
            chromGeneExtrRight = geneExtremityFromGene(chrom[-1], +1)
        else:
            assert oriented == False
            chromGeneExtrLeft = chrom[0].n
            chromGeneExtrRight = chrom[-1].n
        res.update({chromGeneExtrLeft, chromGeneExtrRight})
    return res

def analyseGenomeIntoAdjacencies(genome, oriented=True):
    isinstance(genome, LightGenome)
    setAdjs = set()
    for chrom in genome.values():
        assert len(chrom) > 0
        if len(chrom) > 1:
            for (og1, og2) in myTools.myIterator.slidingTuple(chrom):
                og1 = OGene(str(og1.n), og1.s)
                og2 = OGene(str(og2.n), og2.s)
                assert og1.n != og2.n
                if oriented:
                    # oriented adjacencies
                    setAdjs.add(OAdjacency(og1, og2))
                else:
                    # unoriented adjacencies
                    setAdjs.add(Adjacency(og1.n, og2.n))
    return setAdjs

# tp = nb True Positives
# tn = nb True Negatives
# fp = nb False Positives
# fn = nb False Negatives
# sn = sensitivity
# sp = specificity
Efficiency = collections.namedtuple('efficiency', ('tp', 'tn', 'fp', 'fn', 'sn', 'sp'))

@myTools.verbose
def compareGenomes(genome, refGenome, mode='adjacency', oriented=True, verbose=False):
    """
    :param genome: studied genome
    :param refGenome: the reference genome
    :param mode: either 'adjacency' or 'chromExtremity'
    :param oriented: use oriented genes if True, else use unoriented genes
    :param verbose: print infos in sys.stderr
    :return: Efficiency(Tp, Tn, Fp, Fn, sensitivity, specificity)
    """
    assert mode in {'adjacency', 'chromExtremity'}
    isinstance(genome, LightGenome)
    isinstance(refGenome, LightGenome)

    print >> sys.stderr, "nb of chroms = %s" % len(genome.keys())
    print >> sys.stderr, "nb of monogenic-chroms = %s" % sum([1 for chrom in genome.values() if len(chrom) == 1])
    print >> sys.stderr, "nb of reference chroms = %s" % len(refGenome.keys())
    print >> sys.stderr, "nb of reference monogenic-chroms = %s" % sum([1 for chrom in refGenome.values() if len(chrom) == 1])

    setGeneNames = genome.getGeneNames(asA=set, checkNoDuplicates=False)
    # R: reference
    setGeneNamesR = refGenome.getGeneNames(asA=set, checkNoDuplicates=True)
    print >> sys.stderr, "nb gene names = %s" % len(setGeneNames)
    print >> sys.stderr, "nb gene names in reference = %s" % len(setGeneNamesR)
    print >> sys.stderr, "nb of gene names in both = %s" % len(setGeneNamesR & setGeneNames)
    # this one is empti for the uniquely Amniota_mm_gg simu, whereas it contains mainy values in the whole simualtion
    print >> sys.stderr, "gene names that are in setGeneNames - setGeneNamesR =", setGeneNames - setGeneNamesR
    print >> sys.stderr, "gene names that are in setGeneNamesR - setGeneNames =", setGeneNamesR - setGeneNames

    setRes = set()
    refSetRes = set()
    if mode == 'adjacency':
        itemType = 'OAdjacency' if oriented else 'Adjacency'
        setRes = analyseGenomeIntoAdjacencies(genome, oriented=oriented)
        refSetRes = analyseGenomeIntoAdjacencies(refGenome, oriented=oriented)
    else:
        assert mode == 'chromExtremity'
        itemType = 'OChromosomeExtremity' if oriented else "chromosomeExtremity"
        setRes = analyseGenomeIntoChromExtremities(genome, oriented=oriented)
        refSetRes = analyseGenomeIntoChromExtremities(refGenome, oriented=oriented)

    print >> sys.stderr, "nb of %s = %s" % (itemType, len(setRes))
    print >> sys.stderr, "nb of %s in ref = %s" % (itemType, len(refSetRes))

    # True positive adjs
    sTp = setRes & refSetRes
    Tp = len(sTp)
    Tn = None
    sFp = setRes - refSetRes
    Fp = len(sFp)
    sFn = refSetRes - setRes
    Fn = len(sFn)
    assert len(refSetRes) == Tp + Fn
    sensitivity = float(Tp) / float(Tp + Fn)
    specificity = float(Tp) / float(Tp + Fp)

    print >> sys.stderr, "Tp=%s" % Tp
    print >> sys.stderr, "Fp=%s" % Fp
    print >> sys.stderr, "Fn=%s" % Fn
    print >> sys.stderr, "sensitivity=%s" % sensitivity
    print >> sys.stderr, "specificity=%s" % specificity

    return (Efficiency(Tp, Tn, Fp, Fn, sensitivity, specificity), (sTp, sFp, sFn))

# TODO
# def adjacencyFromIntergene(intergene, genomeWithDict, default=None):
#     assert isinstance(genomeWithDict, LightGenome)

# def buildGenomeFromListAdjs():
#     genome = LightGenome()
#     chr = '0'
#     if not ggInGenome and not gdInGenome: # Nouveau contig etant donne que les deux genes n'ont pas ete vus
#             genome.append(collections.deque([gg,gd]))
#         #print >> sys.stderr, "Nouveau contig de 2 genes : deltaTicTac = %s" % ticTac()
#
#         elif ggInGenome and not gdInGenome:
#             if gg_gIdx == len(genome[gg_cIdx])-1: #gg est a droite d'un contig
#                 genome[gg_cIdx].append(gd) # on ajoute gd a droite du contig
#             elif gg_gIdx == 0: #gg est a gauche d'un contig
#                 genome[gg_cIdx].appendleft(Gene(gd.name,-gd.strand)) # on ajoute gd a gauche du contig en le renversant
#             else:
#                 print >> sys.stderr, "pb de non linearite : %s n'est pas sur un bord de contig dans l'adjacence (%s,%s) etudiee " % (gg.name,gg.name,gd.name)
#             #print >> sys.stderr, "ajout de gd au contig de gg : deltaTicTac = %s" % ticTac()
#
#         elif gdInGenome and not ggInGenome:
#             if gd_gIdx == 0: #gd est a gauche d'un contig
#                 genome[gd_cIdx].appendleft(gg) # on ajoute gg a gauche du contig
#             elif gd_gIdx == len(genome[gd_cIdx])-1: #gd est a droite d'un contig
#                 genome[gd_cIdx].append(Gene(gg.name,-gg.strand)) # on ajoute gg a droite du contig en le retournant
#             else:
#                 print >> sys.stderr, "pb de non linearite : %s n'est pas sur un bord de contig dans l'adjacence (%s,%s) etudiee " % (gd.name,gg.name,gd.name)
#             #print >> sys.stderr, "ajout de gg au contig de gd : deltaTicTac = %s" % ticTac()
#
#         elif ggInGenome and gdInGenome:
#             if gg_cIdx == gd_cIdx:
#                 print >> sys.stderr, "pb de linearite : %s et %s appartiennent au meme contig " % (gg.name,gd.name)
#             else:
#                 assert gg_cIdx != gd_cIdx
#                 if gg_gIdx == len(genome[gg_cIdx])-1: # gg est a droite d'un contig
#                     if gd_gIdx == 0: # gd est a gauche d'un contig
#                         genome[gg_cIdx].extend(genome[gd_cIdx]) # ajoute le contig de gd a celui de gg
#                         del genome[gd_cIdx]
#                         gd_cIdx = gg_cIdx
#                     elif gd_gIdx == len(genome[gd_cIdx])-1:
#                         genome[gd_cIdx].reverse()
#                         tmp = [Gene(g.name,-g.strand) for g in genome[gd_cIdx]]
#                         genome[gg_cIdx].extend(tmp)
#                         del genome[gd_cIdx]
#                         gd_cIdx = gg_cIdx
#                     else:
#                         print >> sys.stderr, "pb de linearite : %s n'est pas sur le bord d'un contig" % gd
#                         continue
#                     #assert gd_gIdx != 0 and gd_gIdx != len(genome[gd_cIdx])-1 and gg_gIdx != 0 and gg_gIdx != len(genome[gg_cIdx])-1 # Onverifie que gg et gd ne sont plus a une extremite d'un contig
#                     #print >> sys.stderr, "Fusion, gg a droite d'un contig : deltaTicTac = %s" % ticTac()
#
#                 elif gg_gIdx == 0: # gg est a gauche d'un contig
#                     if gd_gIdx == 0: # gd est a gauche d'un contig
#                         genome[gg_cIdx].reverse()
#                         tmp = collections.deque(Gene(g.name,-g.strand) for g in genome[gg_cIdx])
#                         genome[gg_cIdx] = tmp
#                         genome[gg_cIdx].extend(genome[gd_cIdx])
#                         del genome[gd_cIdx]
#                         gd_cIdx = gg_cIdx
#                     elif gd_gIdx == len(genome[gd_cIdx])-1:
#                         genome[gd_cIdx].extend(genome[gg_cIdx])
#                         del genome[gg_cIdx]
#                         gg_cIdx = gd_cIdx
#                     else:
#                         print >> sys.stderr, "pb de linearite : %s n'est pas sur le bord d'un contig" % gd
#                         continue
#                     #assert gd_gIdx != 0 and gd_gIdx != len(genome[gd_cIdx])-1 and gg_gIdx != 0 and gg_gIdx != len(genome[gg_cIdx])-1 # Onverifie que gg et gd ne sont plus a une extremite d'un contig
#                     #print >> sys.stderr, "Fusion, gg a gauche d'un contig : deltaTicTac = %s" % ticTac()
#         else:
#             print >> sys.stderr, "Cas impossible"
#             raise # on leve une erreur