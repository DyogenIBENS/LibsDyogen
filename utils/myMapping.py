# -*- coding: utf-8 -*-
# LibsDyogen
# python 2.7
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

# This file contains classes and methods for 1D management of chromosomes

import myTools
import myGenomes
import myLightGenomes
from myLightGenomes import OGene as OGene
import collections


def setOfGeneNames(genome):
    assert isinstance(genome, myLightGenomes.LightGenome)
    geneNames = set()
    for chrom in genome.values():
        geneNames.update(geneName for (geneName, _) in chrom)
    return geneNames


# because of collections.Counter
@myTools.minimalPythonVersion((2, 7))
def nbDup(g_fID_filt):
    assert isinstance(g_fID_filt, myLightGenomes.LightGenome)
    # create a long chromosome with all the chromosome gene names
    newF = []
    for chr in g_fID_filt.values():
        chr = [g for (g, _) in chr]
        newF = newF + chr
    dupCounter = collections.Counter(newF)
    nbGeneDup = sum([duplicates-1 for duplicates in dupCounter.values()])
    return (nbGeneDup, dupCounter)


# Rewrite genomes as a list of family Ids
def labelWithFamID(genome, families):
    assert isinstance(families, myLightGenomes.Families)
    assert isinstance(genome, myLightGenomes.LightGenome)
    assert type(genome.values()[0]) == list
    assert isinstance(genome.values()[0][0], tuple)
    assert all(len(chrom) > 0 for chrom in genome.values())
    newGenome = myLightGenomes.LightGenome()
    for c in genome.keys():
        #assert len(genome[c]) >=1
        for (g, s) in genome[c]:
            fid = families.getFamID(g, default=None)
            newGenome[c].append(OGene(fid, s))
    return newGenome


# Rewrite genomes as a list of family names
def labelWithFamNames(genome, families):
    assert isinstance(families, myLightGenomes.Families)
    assert isinstance(genome, myLightGenomes.LightGenome)
    assert type(genome.values()[0]) == list
    assert isinstance(genome.values()[0][0], tuple)
    assert all(len(chrom) > 0 for chrom in genome.values())
    newGenome = myLightGenomes.LightGenome()
    for c in genome.keys():
        #assert len(genome[c]) >=1
        for g in genome[c]:
            fn = families.getFamNameByName(g.n, default=None)
            newGenome[c].append(OGene(fn, g.s))
    return newGenome


# Allow to change the 1D mapping
class Mapping(object):
    # from new to old:
    # newMapping format: [[0,1,2], [3,6], [5]]
    # means that [0,1,2,3,4,5,6] -> [a,b,c] with a= the collapsed [0,1,2],
    # b = [3,6] and c = [5].
    # from old to new:
    # newMapping format: [0, 0, 0, 1, 2, 1], means that [0,1,2,3,4,5,6] -> [[0,1,2], [3,6], [5]]
    # if newMapping = [[0],[1],[3]], it means that gene the 3rd gene (idx=2) was
    # deleted during the mapping.
    def __init__(self, newMapping, oldToNew=False):
        # Forward version of the mapping: newIdx -> oldIdx
        self.new = []
        # Reverse version of the mapping: oldIdx -> newIdx
        # Need to use a dict since some oldIdx may have been removed.
        self.old = {}
        newToOld = not oldToNew
        if newToOld:
            self.new = newMapping

            # DEBUG check consistency
            allOldIdxs = reduce(lambda x, y: x + y, newMapping)
            if len(allOldIdxs) != len(set(allOldIdxs)):
                raise ValueError('Their is a duplicated oldIdx in the newMapping argument')

            for newIdx, oldIdxs in enumerate(newMapping):
                for oldIdx in oldIdxs:
                    self.old[oldIdx] = newIdx
        else:
            allNewIdxs = set(newMapping)
            self.old = dict(enumerate(newMapping))
            self.new = [[] for _ in range(len(set(allNewIdxs)))]
            for oldIdx, newIdx in enumerate(newMapping):
                self.new[newIdx].append(oldIdx)

    def __iter__(self):
        for (newIdx, oldIdxs) in enumerate(self.new):
            yield (newIdx, oldIdxs)

    def __get__(self, key):
        return self.new[key]

    def __str__(self):
        return str(self.new)

    def __repr__(self):
        return str(self.new)

    # list-like index access
    def __getitem__(self, key):
        return self.new[key]

    # addition like method MappingC<-A = MappingC<-B + MappingB<-A
    def __add__(self, other):
            assert other.__class__.__name__ == 'Mapping'
            res = Mapping(self.new)
            res.drawFrom(other)
            return res

    # Draw from something: "puiser à partir de qch"
    def drawFrom(self, oldMapping):
        assert oldMapping.__class__.__name__ == 'Mapping'
        newMapping = []
        for oldIdxs in self.new:
            newMapping.append([])
            for oldIdx in oldIdxs:
                for veryOldIdx in oldMapping.new[oldIdx]:
                    newMapping[-1].append(veryOldIdx)
        self = self.__init__(newMapping)


# remap genes following instructions in genomeMapping
# while remapping and collapsing genes into blocks, find the best consensus
# orientation for the block.
def remap(genome, genomeMapping, assertCollapsedGenesHaveSameName=True):
    if isinstance(genome, myLightGenomes.LightGenome):
        assert type(genomeMapping) == dict
        # if genome is a dict or a defaultdict
        newGenome = myLightGenomes.LightGenome()
        for c in genomeMapping:
            newGenome[c] = []
            chrMapping = genomeMapping[c]
            assert chrMapping.__class__.__name__ == 'Mapping'
            for (newIdx, oldIdxs) in chrMapping:
                if assertCollapsedGenesHaveSameName:
                    # Assert that all the names of the collapsed old genes are the
                    # same.
                    assert len(set([genome[c][oldIdx].n for oldIdx in oldIdxs])) == 1,\
                        "At least one collapsed gene is not in the same family as one of its co-collapsed gene."
                name = genome[c][oldIdxs[0]].n
                if len(oldIdxs) == 1:
                    s = genome[c][oldIdxs[0]].s
                elif len(oldIdxs) >= 2:
                    firstOldS = genome[c][oldIdxs[0]].s
                    # Set a None orientation if no consensus
                    # idxsOfChainedTbs[0] = min(idxsOfChainedTbs)
                    if (firstOldS == +1 or firstOldS == -1)\
                            and all([genome[c][oldIdx].s == firstOldS
                                     for oldIdx in oldIdxs]):
                        s = firstOldS
                    else:
                        s = None
                newGenome[c].append(OGene(name, s))
    # Deprecated
    # elif genome.__class__.__name__ == myGenomes.Genome:
    #     assert set(genome.lstGenes.keys()) == set(genomeMapping.keys())
    #     # FIXME no empty constructor
    #     newGenome = myGenomes.Genome()
    #     newGenome.lstGenes = collections.defaultdict(list)
    #     for c in genomeMapping:
    #         chrom = genome.lstGenes[c]
    #         # chrom == genome.lstGenes[c]
    #         chrMapping = genomeMapping[c]
    #         assert chrMapping.__class__.__name__ == 'Mapping'
    #         for (newIdx, oldIdxs) in chrMapping:
    #             if assertCollapsedGenesHaveSameName:
    #                 # Assert that all the names of the collapsed old genes are the
    #                 # same.
    #                 assert len(set([chrom[oldIdx].names[0] for oldIdx in oldIdxs])) == 1,\
    #                     "At least one collapsed gene is not in the same family as one of its co-collapsed gene."
    #             name = chrom[oldIdxs[0]].names[0]
    #             if len(oldIdxs) == 1:
    #                 s = chrom[oldIdxs[0]].strand
    #             elif len(oldIdxs) >= 2:
    #                 firstOldS = chrom[oldIdxs[0]].strand
    #                 # Set a None orientation if no consensus
    #                 # idxsOfChainedTbs[0] = min(idxsOfChainedTbs)
    #                 if firstOldS == +1 or firstOldS == -1\
    #                         and all([chrom[oldIdx].strand == firstOldS
    #                                  for oldIdx in oldIdxs]):
    #                     s = firstOldS
    #                 else:
    #                     s = None
    #             #FIXME not sure for beg and end...
    #             beg = chrom[oldIdxs[0]].beginning
    #             end = chrom[oldIdxs[-1]].end
    #             #FIXME fix the name list, it should maybe be update with all the names of the co-collapsed genes.
    #             names = chrom[oldIdxs[0]].names
    #             newGenome.lstGenes[c].append(myGenomes.Gene(c, beg, end, s, names))
    #     newGenome.init()
    #     #TODO
    #     #elif type(genome) == list:
    else:
        raise TypeError('Not a known type of genome')
    return newGenome

def remapCoFilterContentAndSize(genome, removedNames, minChromLength, mOld=None):
    assert type(removedNames) == set
    assert isinstance(genome, myLightGenomes.LightGenome)
    (newGenome, mfGC2old, (nbChrLoss, nbGeneLoss)) = remapFilterGeneContent(genome, removedNames)
    assert len(newGenome.keys()) == len(mfGC2old.keys())
    (newGenome, mfS2fGC, (nbChrLossb, nbGeneLossb)) = remapFilterSize(newGenome, minChromLength)
    nbChrLoss += nbChrLossb
    nbGeneLoss += nbGeneLossb
    # Update the mapping
    mf2old = {}
    for c in mfS2fGC:
        mf2old[c] = mfS2fGC[c] + mfGC2old[c]
        if mOld is not None:
            mf2old[c] = mf2old[c] + mOld[c]
    #for c in mf2old:
    #    mf2old[c].drawFrom(mfGC2old[c])
    #    mfilt2old[c].drawFrom(mOld[c])
    #asert len(newGenome.keys()) == len(mfilt2old.keys())
    ##DEBUG Loop assertion
    #for (c, chrom) in genome.iteritems():
    #    maxOldIdx = max([max(oldIdxs) for oldIdxs in mfilt2old[c].new])
    #    assert maxOldIdx < len(chrom)
    return (newGenome, mf2old, (nbChrLoss, nbGeneLoss))


# TODO FIXME
def remapFilterGeneLength(genome, minGeneLengthInBp, mOld=None):
    assert type(genome, myLightGenomes.LightGenome)
    (mf2old, (nbChrLoss, nbGeneLoss)) = mapFilterGeneLength(genome, minGeneLengthInBp)
    newGenome = remap(genome, mf2old, assertCollapsedGenesHaveSameName=False)
    return (newGenome, mf2old, (nbChrLoss, nbGeneLoss))

# TODO FIXME
def mapFilterGeneLength(genome, minGeneLengthInBp, mOld=None):
    nbChrLoss = 0
    nbGeneLoss = 0
    mfilt2old = {}
    for (c, chrom) in genome:
        newIdx = []
        for (i, gene) in enumerate(chrom):
            # FIXME: access gene length
            if gene.length < minGeneLengthInBp:
                nbGeneLoss += 1
                continue
            else:
                newIdx.append([i])
        if len(newIdx) > 0:
            mfilt2old[c] = Mapping(newIdx)
        else:
            nbChrLoss += 1
    return (mfilt2old, (nbChrLoss, nbGeneLoss))


def remapFilterSize(genome, minChromLength, mOld=None):
    assert isinstance(genome, myLightGenomes.LightGenome)
    (mf2old, (nbChrLoss, nbGeneLoss)) = mapFilterSize(genome, minChromLength)
    newGenome = remap(genome, mf2old, assertCollapsedGenesHaveSameName=False)
    if mOld is not None:
        for c in mf2old:
            mf2old[c] = mf2old[c] + mOld[c]
    # Conservation law of genes
    #nbGenesInOldG = sum(len(chrom) for chrom in genome.values())
    #nbGenesInNewG = sum(len(chrom) for chrom in newGenome.values())
    #assert nbGenesInOldG == nbGenesInNewG + nbGeneLoss
    return (newGenome, mf2old, (nbChrLoss, nbGeneLoss))


def remapFilterGeneContent(genome, removedNames, mOld=None):
    assert type(removedNames) == set
    assert isinstance(genome, myLightGenomes.LightGenome)
    (mf2old, (nbChrLoss, nbGeneLoss)) = mapFilterGeneContent(genome, removedNames, mOld=mOld)
    #DEBUG Loop assertion
    #for (c, chrom) in genome.iteritems():
    #    nbGenesInOldC = len(chrom)
    #    maxIndexOldGenesInCMapping = max([max(oldIdxs) for oldIdxs in mfilt2old[c].new])
    #    assert maxIndexOldGenesInCMapping < nbGenesInOldC
    newGenome = remap(genome, mf2old, assertCollapsedGenesHaveSameName=False)
    if mOld is not None:
        for c in mf2old:
            mf2old[c] = mf2old[c] + mOld[c]
    # Conservation law of genes
    #nbGenesInOldG = sum(len(chrom) for chrom in genome.values())
    #nbGenesInNewG = sum(len(chrom) for chrom in newGenome.values())
    #assert nbGenesInOldG == nbGenesInNewG + nbGeneLoss
    return (newGenome, mf2old, (nbChrLoss, nbGeneLoss))


def remapRewriteInTb(genome_fID, tandemGapMax=0, mOld=None):
    assert isinstance(genome_fID, myLightGenomes.LightGenome)
    # Take care to not give mOld=mOld, since mOld is not a mapping
    # from genome_fID.
    (mf2old, (nbTandemDup)) = mapRewriteInTb(genome_fID, tandemGapMax=tandemGapMax)
    newGenome = remap(genome_fID, mf2old, assertCollapsedGenesHaveSameName=True)
    if mOld is not None:
        for c in mf2old:
            mf2old[c] = mf2old[c] + mOld[c]
    # Conservation law of genes
    #nbGenesInOldG = sum(len(chrom) for chrom in genome_fID.values())
    #nbGenesInNewG = sum(len(chrom) for chrom in newGenome.values())
    #assert nbGenesInOldG == nbGenesInNewG + nbTandemDup
    return (newGenome, mf2old, (nbTandemDup))


# Returns a mapping corresponding to the removal of too small chromosomes.
# Warning, it could be interesting to update an old mapping associated with
# 'genome'.
# 'mOld' is an Old mapping. At the end of the construction of the new mapping this
# old mapping can be used to tranfer the information of the old mapping to the
# new mapping.
def mapFilterSize(genome, minChromLength, mOld=None):
    assert isinstance(genome, myLightGenomes.LightGenome)
    nbChrLoss = 0
    nbGeneLoss = 0
    mfilt2old = {}
    for c in genome.keys():
        nbGenes = len(genome[c])
        if nbGenes < minChromLength:
            nbChrLoss += 1
            nbGeneLoss += nbGenes
        else:
            mfilt2old[c] = Mapping([[i] for (i, _) in enumerate(genome[c])])
    return (mfilt2old, (nbChrLoss, nbGeneLoss))


# Returns a mapping corresponding to the removal a set of gene names from a
# genome.
# Warning, it could be interesting to update an old mapping associated with
# 'genome'.
# 'mOld' is an Old mapping. At the end of the construction of the new mapping this
# old mapping can be used to tranfer the information of the old mapping to the
# new mapping.
def mapFilterGeneContent(genome, removedNames, mOld=None):
    assert isinstance(genome, myLightGenomes.LightGenome)
    nbGeneLoss = 0
    nbChrLoss = 0
    mfilt2old = {}
    for c in genome:
        for (oldIdx, (g, s)) in enumerate(genome[c]):
            if g in removedNames:
                nbGeneLoss += 1
            else:
                if c not in mfilt2old:
                    mfilt2old[c] = []
                mfilt2old[c].append([oldIdx])
        if c in mfilt2old:
            mfilt2old[c] = Mapping(mfilt2old[c])
        else:
            nbChrLoss += 1
    #TODO
    #elif myGenomes
    #elif list
    return (mfilt2old, (nbChrLoss, nbGeneLoss))


# Inputs :
#       genome_fID = [..., (fID,s), ...] with 'fID' the line number of the gene in the ancGene ans 's' the strand of the gene in the genome
# Outputs :
#       mtb2g : a mapping corresponding to the rewritting process
#####################################################################################################################################
# 'mOld' is an Old mapping. At the end of the construction of the new mapping this
# old mapping can be used to tranfer the information of the old mapping to the
# new mapping.
def mapRewriteInTb(genome_fID, tandemGapMax=0):
    assert isinstance(genome_fID, myLightGenomes.LightGenome)
    # the rewritten genome
    # Need to keep dictionaries because there is often a gap in chromosome notations
    tb2g = {}
    # TODO be sure that this step is consistent with the 2D distance metric chosen.
    # For instance, if the DPD is chosen, on vertical and horizontal lines, the
    # distances are not consistent with the 1D distance metric.
    nbTandemDup = 0
    tandemDistMax = tandemGapMax + 1
    #TODO next loops could be optimised
    combinator = myTools.myCombinator()
    for c, chrom_tb in genome_fID.iteritems():
        tb2g[c] = []

        # print >> sys.stderr, "Length in tbs before tandem gap = %s" % len(chrom_tb)
        for (i, fID) in enumerate(chrom_tb):
            isAlone = True
            for dist in range(1, min(tandemDistMax + 1, len(chrom_tb) - i)):
                if chrom_tb[i + dist] == chrom_tb[i]:
                    pair = (i + dist, i)
                    # Add a link between the elements of a pair
                    combinator.addLink(pair)
                    isAlone = False
            if isAlone:
                combinator.addLink((i, i))
        combinator.reduce()

        tbChains = list(combinator)
        # print >> sys.stderr, "Nb of chains of at least 1 tb = %s" % len(tbChains)
        # print >> sys.stderr, "Nb of chains of at least 2 tbs = %s" % len([a for a in tbChains if len(a) >=2])
        # print >> sys.stderr, "Chains of at least 2 tbs = %s" % [a for a in tbChains if len(a) >= 2]
        # Sort neighbourhood by the increasing smallest index of the neighbourhood.
        # TODO List could be yielded sorted (improve the myCombinator class)
        for tbChain in tbChains:
            tbChain.sort()
        #print >> sys.stderr, "len(tbChains)=%s" % len(tbChains)
        # Since what precedes the next script line is equivalent to:
        # chainsOfTbs.sort(lambda x: min(x))
        tbChains.sort(key=lambda x: x[0])

        firstChainTbIdxs = []
        otherChainTbIdxs = []
        for tbChain in tbChains:
            if len(tbChain) == 1:
                firstChainTbIdx = tbChain[0]
                otherChainTbIdxs.append([])
            elif len(tbChain) >= 2:
                firstChainTbIdx = tbChain[0]
                otherChainTbIdxs.append([])
                for tbIdx in tbChain[1:]:
                    otherChainTbIdxs[-1].append(tbIdx)
            else:
                raise
            firstChainTbIdxs.append(firstChainTbIdx)
        assert len(firstChainTbIdxs) == len(otherChainTbIdxs)

        for (firstTb, otherTbs) in zip(firstChainTbIdxs, otherChainTbIdxs):
            tb2g[c].append([])
            for oldTbIdx in [firstTb] + otherTbs:
                # lis1 + list2 returns the concatenation of the two lists
                tb2g[c][-1].append(oldTbIdx)
        assert len(tb2g[c]) == len(firstChainTbIdxs), len(tb2g[c])

        nbTandemDup +=\
            sum([len(gTandemDuplicates)-1 for gTandemDuplicates in tb2g[c]])
        # DEBUG assertion
        listIsSorted = lambda l: all([i] <= l[i+1] for i in range(len(l)-1))
        assert listIsSorted(tb2g[c])
        combinator.reset()

    mtb2g = {}
    for (c, newMapC) in tb2g.iteritems():
        mtb2g[c] = Mapping(newMapC)

    # #DEBUG assertion
    # nbOffIDGenes = sum([len(chrom) for chrom in genome_fID.values()])
    # nbOfTbs = sum([len(chrom) for chrom in tb2g.values()])
    # assert nbTandemDup == nbOffIDGenes - nbOfTbs, "%s == %s - %s" % (nbTandemDup, nbOffIDGenes, nbOfTbs)

    return (mtb2g, (nbTandemDup))

@myTools.deprecated
def mapRewriteInTbOld(genome_fID, tandemGapMax=0):
    nbTandemDup = 0
    # the rewritten genome
    # Need to keep dictionnaries because there is often a gap in chromosome notations
    tmp_genome_tb = {}
    tb2g = {}
    # Number of Genic Tandem Duplication global
    for c in genome_fID:
        tb2g[c] = []
        # temp values
        tmp_genome_tb[c] = []
        old_fID = 'foo'  # Not None, since None is used to mark genes to be removed
        for i, (fID, s) in enumerate(genome_fID[c]):
            # if the current gene is not a paralog of the previous gene
            if fID != old_fID:
                # start a new TB
                tmp_genome_tb[c].append(fID)
                # add the current gene index to the new tb component of tb2g
                tb2g[c].append([i])
            else:
                # fID == old_fID:
                # The index of the current gene is added to the last tb
                tb2g[c][-1].append(i)
                nbTandemDup += 1
            old_fID = fID

    # TODO be sure that this step is consistent with the 2D distance metric chosen.
    # For instance, if the DPD is chosen, on vertical and horizontal lines, the
    # distances are not consistent with the 1D distance metric.
    if tandemGapMax > 0:
        nbTandemDup = 0
        tandemDistMax = tandemGapMax + 1
        #TODO next loops could be optimised
        tmp_tb2g = {}
        combinator = myTools.myCombinator()
        for c, chrom_tb in tmp_genome_tb.iteritems():
            # print >> sys.stderr, "Length in tbs before tandem gap = %s" % len(chrom_tb)
            for (i, fID) in enumerate(chrom_tb):
                # Range starts at 2 because since genomes are already rewritten
                # in tbs there is no adjacent tbs.
                isAlone = True
                for dist in range(2, min(tandemDistMax + 1, len(chrom_tb) - i)):
                    if chrom_tb[i + dist] == chrom_tb[i]:
                        pair = (i + dist, i)
                        # Add a link between the elements of a pair
                        combinator.addLink(pair)
                        isAlone = False
                if isAlone:
                    combinator.addLink((i, i))
            combinator.reduce()

            tbChains = list(combinator)
            # print >> sys.stderr, "Nb of chains of at least 1 tb = %s" % len(tbChains)
            # print >> sys.stderr, "Nb of chains of at least 2 tbs = %s" % len([a for a in tbChains if len(a) >=2])
            # print >> sys.stderr, "Chains of at least 2 tbs = %s" % [a for a in tbChains if len(a) >= 2]
            # Sort neighbourhood by the increasing smallest index of the neighbourhood.
            # TODO List could be yielded sorted (improve the myCombinator class)
            for tbChain in tbChains:
                tbChain.sort()
            #print >> sys.stderr, "len(tbChains)=%s" % len(tbChains)
            # Since what precedes the next script line is equivalent to:
            # chainsOfTbs.sort(lambda x: min(x))
            tbChains.sort(key=lambda x: x[0])

            firstChainTbIdxs = []
            otherChainTbIdxs = []
            for tbChain in tbChains:
                if len(tbChain) == 1:
                    firstChainTbIdx = tbChain[0]
                    otherChainTbIdxs.append([])
                elif len(tbChain) >= 2:
                    firstChainTbIdx = tbChain[0]
                    otherChainTbIdxs.append([])
                    for tbIdx in tbChain[1:]:
                        otherChainTbIdxs[-1].append(tbIdx)
                else:
                    raise
                firstChainTbIdxs.append(firstChainTbIdx)
            assert len(firstChainTbIdxs) == len(otherChainTbIdxs)

            tmp_tb2g[c] = []
            for (firstTb, otherTbs) in zip(firstChainTbIdxs, otherChainTbIdxs):
                tmp_tb2g[c].append([])
                for oldTbIdx in [firstTb] + otherTbs:
                    # lis1 + list2 returns the concatenation of the two lists
                    tmp_tb2g[c][-1] = tmp_tb2g[c][-1] + tb2g[c][oldTbIdx]
            assert len(tmp_tb2g[c]) == len(firstChainTbIdxs), len(tmp_tb2g[c])

            nbTandemDup +=\
                sum([len(gTandemDuplicates)-1 for gTandemDuplicates in tmp_tb2g[c]])
            # DEBUG assertion
            listIsSorted = lambda l: all([i] <= l[i+1] for i in range(len(l)-1))
            assert listIsSorted(tmp_tb2g[c])
            combinator.reset()
        tb2g = tmp_tb2g

    mtb2g = {}
    for (c, newMapC) in tb2g.iteritems():
        mtb2g[c] = Mapping(newMapC)

    # DEBUG assertion
    # nbOffIDGenes = sum([len(chrom) for chrom in genome_fID.values()])
    # nbOfTbs = sum([len(chrom) for chrom in genome_tb.values()])
    # assert N_GTD_g == nbOffIDGenes - nbOfTbs, "%s == %s - %s" % (N_GTD_g, nbOffIDGenes, nbOfTbs)

    return (mtb2g, (nbTandemDup))


def calcTandemDup(genome, family, allowedGap=0, calcCDF=True):
    allowedDist = allowedGap
    genome_fID = labelWithFamID(genome, family)

    (genomeFilt, gf2gfID, _) = remapFilterGeneContent(genome_fID, set())

    familySizes = collections.defaultdict(int)
    for chrom in genomeFilt:
        for gene in genomeFilt[chrom]:
            familySizes[gene.n] += 1

    nGeneDupl = 0
    for familySize in familySizes.values():
        nGeneDupl += familySize - 1

    dupsInDist = 0
    for genes in genomeFilt.itervalues():
        for iGene in xrange(len(genes) - 1 - allowedDist):
            correctedDist = [genes[iGene].n != genes[iGene + dist].n for dist in xrange(1, allowedDist)]
            if all(correctedDist) and genes[iGene].n == genes[iGene + allowedDist + 1].n:
                dupsInDist += 1

    if calcCDF:
        # recDist = recent distance
        for recDist in xrange(allowedDist):
            for iGene in xrange(len(genes) - 1 - recDist):
                correctedDist = [genes[iGene].n != genes[iGene + dist].n for dist in xrange(1, recDist + 1)]
                if all(correctedDist) and genes[iGene].n == genes[iGene + recDist + 1].n:
                    dupsInDist += 1

    return (dupsInDist, nGeneDupl)
