# -*- coding: utf-8 -*-

# LibsDyogen version 1.0 (6/11/2015)
# python v2.7 at least is needed
# Copyright © 2015 IBENS/Dyogen : Joseph LUCAS, Lucas TITTMANN and Hugues ROEST CROLLIUS
# mail : jlucas@ens.fr
# Licences GLP v3 and CeCILL v2

# This file contains classes and methods for 1D management of chromosomes
import sys

import myTools
import myLightGenomes
from myLightGenomes import OGene as OGene
import collections

# because of collections.Counter
@myTools.minimalPythonVersion((2, 7))
def nbDup(g_fID):
    assert isinstance(g_fID, myLightGenomes.LightGenome)
    # create a long chromosome with all the chromosome gene names
    newF = []
    for chrom in g_fID.values():
        chrom = [gene.n for gene in chrom]
        newF = newF + chrom
    dupCounter = collections.Counter(newF)
    nbGeneDup = sum([duplicates-1 for duplicates in dupCounter.values()])
    return (nbGeneDup, dupCounter)

# Rewrite genomes as a list of family Ids
def labelWithFamID(genome, families=None):
    """
    This function keeps the original copy
    """
    assert isinstance(genome, myLightGenomes.LightGenome)
    # there is at least one chromosome
    assert len(genome.values()) > 0
    assert isinstance(genome.values()[0], list)
    if families is None:
        return genome
    assert isinstance(families, myLightGenomes.Families)
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
def labelWithFamNames(genome, families, keepGnOfGenesNotInFamilies=False):
    assert isinstance(families, myLightGenomes.Families)
    assert isinstance(genome, myLightGenomes.LightGenome)
    assert isinstance(genome.values()[0], list)
    assert all(len(chrom) > 0 for chrom in genome.values()), " ".join(['%s:%s' % (c,len(chrom)) for c, chrom in genome.iteritems()])
    newGenome = myLightGenomes.LightGenome()
    for c in genome.keys():
        for g in genome[c]:
            fn = families.getFamNameByName(g.n, default=None)
            if fn is None and keepGnOfGenesNotInFamilies is True:
                newGenome[c].append(OGene(g.n, g.s))
            else:
                # either
                # 1) fn is None and keepGnOfGenesNotInFamilies is False -> the new gene name is None
                # 2) fn is not None -> gene is in family
                newGenome[c].append(OGene(fn, g.s))
    return newGenome

def labelWithOrthologNames(genome, families):
    """
    Rewrite genomes as a list of orthologs
    works with the gene name nomenclature of MagSimus1
    :param genome: myLightGenomes.LightGenome
    :param families: myLightGenomes.Families
    :return newGenome:
    """
    assert isinstance(families, myLightGenomes.Families)
    assert isinstance(genome, myLightGenomes.LightGenome)
    assert type(genome.values()[0]) == list
    assert all(len(chrom) > 0 for chrom in genome.values())

    newGenome = myLightGenomes.LightGenome()
    for c in genome.keys():
        for g in genome[c]:
            family = families.getFamilyByName(g.n, default=None)
            if family is None:
                newGenome[c].append(OGene(None, g.s))
            elif g.n == family.fn:
                # positional ortholog
                newGenome[c].append(OGene(g.n, g.s))
            elif family is not None:
                # paralog
                newGenome[c].append(OGene(None, g.s))
            else:
                raise ValueError
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
    """
    This function conserves the genome
    """
    if isinstance(genome, myLightGenomes.LightGenome):
        assert type(genomeMapping) == dict
        # if genome is a dict or a defaultdict
        newGenome = myLightGenomes.LightGenome()
        for c in genomeMapping:
            newGenome[c] = myLightGenomes.Chromosome()
            chrMapping = genomeMapping[c]
            assert chrMapping.__class__.__name__ == 'Mapping'
            for (newIdx, oldIdxs) in chrMapping:
                if assertCollapsedGenesHaveSameName:
                    # Assert that all the names of the collapsed old genes are the same.
                    assert len(set([genome[c][oldIdx].n for oldIdx in oldIdxs])) == 1,\
                        "At least one collapsed gene is not in the same family as one of its co-collapsed gene."
                name = genome[c][oldIdxs[0]].n
                if len(oldIdxs) == 1:
                    s = genome[c][oldIdxs[0]].s
                elif len(oldIdxs) >= 2:
                    firstOldS = genome[c][oldIdxs[0]].s
                    if (firstOldS == +1 or firstOldS == -1) and \
                            all([genome[c][oldIdx].s == firstOldS for oldIdx in oldIdxs]):
                        s = firstOldS
                    else:
                        # Set a None orientation if no consensus
                        s = None
                newGenome[c].append(OGene(name, s))
    else:
        raise TypeError('Not a known type of genome')
    return newGenome

def remapCoFilterContentAndSize(genome, removedNames, minChromLength, mOld=None):
    assert type(removedNames) == set
    assert isinstance(genome, myLightGenomes.LightGenome)
    (newGenome, mfGC2old, (nbChrLoss, nbGeneLoss)) = remapFilterGeneContent(genome, removedNames)
    assert len(newGenome.keys()) == len(mfGC2old.keys())
    (newGenome, mfS2fGC, (nbChrLossSize, nbGeneLossSize)) = remapFilterSize(newGenome, minChromLength)
    nbChrLoss += nbChrLossSize
    nbGeneLoss += nbGeneLossSize
    # Update the mapping
    mf2old = {}
    for c in mfS2fGC:
        mf2old[c] = mfS2fGC[c] + mfGC2old[c]
        if mOld is not None:
            mf2old[c] = mf2old[c] + mOld[c]
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
    return (mfilt2old, (nbChrLoss, nbGeneLoss))

# 'mOld' is an Old mapping. At the end of the construction of the new mapping this
# old mapping can be used to tranfer the information of the old mapping to the
# new mapping.
def mapRewriteInTb(genome_fID, tandemGapMax=0):
    """
    :param genome_fID: [..., (fID,s), ...] with 'fID' the line number of the gene in the ancGene ans 's' the strand of the gene in the genome
    :param tandemGapMax:
    :return: mtb2g : a mapping corresponding to the rewritting process
    """
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
                if chrom_tb[i + dist].n == chrom_tb[i].n:
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
                    if chrom_tb[i + dist].n == chrom_tb[i].n:
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

def calcGeneDist(pos1, pos2, ifNotOnSameChr = -1):
    # Returns the distance between two Genes
    # pos1 = (chromosome, positionOnChromosome)
    # Returns ifNotOnSameChr if Genes are on different chromosome
    # Returns None if either of the positions is None
    if pos1 is None or pos2 is None:
        return(None)
    if pos1[0] != pos2[0]:
        return(ifNotOnSameChr)
    else:
        return(abs(pos1[1] - pos2[1]))

def calcDupDist(genome, family):
    # Returns tuple: (dict of distances of duplications, genome remapped with family names)
    # the former has a list of distances of a gene to its next duplication, with
    # chromosomes in order of the remapped chromosome (latter output)
    # if there are no duplications, the list in dict of distances is empty
    # if the next duplication is on another chromosome, it is coded as -1
    genome_fID = labelWithFamID(genome, family)
    (genomeFilt, gf2gfID, _) = remapFilterGeneContent(genome_fID, {None})

    posDict = {}
    distDict = {}
    for (chr, genes) in genomeFilt.iteritems():
        for (pos, gene) in enumerate(genes):
            if gene.n in posDict:
                posDict[gene.n].append((chr, pos))
                newDist = calcGeneDist(posDict[gene.n][-2], posDict[gene.n][-1])
                distDict[gene.n].append(newDist)
            else:
                posDict[gene.n] = [(chr, pos)]
                distDict[gene.n] = []

    return (distDict, genomeFilt)

def calcTotalDupFromDist(distDict):
    # Given the dict of distances from calcDupDist, return number of duplication in genome
    return(sum([len(distances) for distances in distDict.itervalues()]))

def calcTandemDupFromDist(distDict, gapMax, cumulated = True):
    # Given the dict of distances from calcDupDist and a gapMax, return number of tandemDuplications
    distMax = gapMax+1
    if cumulated:
        return(len([distance for distances in distDict.itervalues() for distance in distances if distance != -1 and distance <= distMax]))
    else:
        return (len([distance for distances in distDict.itervalues() for distance in distances if distance != -1 and distance == distMax]))

def getAncFamNames(genomeOrFamDesc, famAnc):
    if isinstance(genomeOrFamDesc, myLightGenomes.LightGenome):
        geneNameList = [gene.n for genes in genomeOrFamDesc.itervalues() for gene in genes]
        out = [famAnc.getFamNameByName(gene, default=None) for gene in geneNameList]
    elif isinstance(genomeOrFamDesc, myLightGenomes.Families):
        out = [famAnc.getFamNameByName(next(iter(recFam.dns)), default=None) for recFam in genomeOrFamDesc]
    else:
        raise ValueError('Function not implemented for type ' + str(type(genomeOrFamDesc)))
    return(out)

def calcNumberOfGeneDeletions(genomeOrFamDesc, family):
    geneNameList = getAncFamNames(genomeOrFamDesc, family)
    nGeneDeleted = len(set([recFamily.fn for recFamily in family]) - set(geneNameList))
    return (nGeneDeleted)

def calcNumberOfGeneBirths(genomeOrFamDesc, family):
    geneNameList = getAncFamNames(genomeOrFamDesc, family)
    nGeneBirths = geneNameList.count(None)
    return (nGeneBirths)

def calcNumberOfGeneDuplications(genomeOrFamDesc, family):
    geneNameList = getAncFamNames(genomeOrFamDesc, family)
    descOfAncGenes = [gene for gene in geneNameList if gene is not None]
    # nGeneDups = sum([nbParalogs - 1 for nbParalogs in collections.Counter(descOfAncGenes).values()]) would give the same result
    nGeneDups = len(descOfAncGenes) - len(set(descOfAncGenes))
    return (nGeneDups)

@myTools.deprecated
def calcTandemDup(genome, family, allowedGap=0, cumulated=True):

    def calcTandemDupWithStrictGap(genomeFilt, allowedGap):
        # Helper function to calcTandemDup
        # Returns on filtered genome the exact number of tandemDups with given gap
        dupsInDist = 0
        for genes in genomeFilt.itervalues():
            for iGene in xrange(len(genes) - 1 - allowedGap):
                correctedDist = [genes[iGene].n != genes[iGene + gap].n for gap in xrange(1, allowedGap+1)]
                if all(correctedDist) and genes[iGene].n == genes[iGene + allowedGap + 1].n:
                    dupsInDist += 1
        return(dupsInDist)
    # DEPRECATED: Used as a test function as result is verified by hand
    # Use calcTandemDupFromDist in combination with calcDupDist for faster results
    #
    # Function to calculate the number of tandemDups at given allowedGap
    # returns ( nTandemDup, nAllDup)
    # If cumulated is True, returns number of tandemDuplications with
    # gap <= allowedGap

    genome_fID = labelWithFamID(genome, family)
    (genomeFilt, gf2gfID, _) = remapFilterGeneContent(genome_fID, {None})

    familySizes = collections.defaultdict(int)
    for chrom in genomeFilt:
        for gene in genomeFilt[chrom]:
            familySizes[gene.n] += 1

    nGeneDupl = 0
    for familySize in familySizes.values():
        nGeneDupl += familySize - 1

    dupsInDist = calcTandemDupWithStrictGap(genomeFilt, allowedGap)

    if cumulated:
        for recentGap in xrange(allowedGap):
            dupsInDist += calcTandemDupWithStrictGap(genomeFilt, recentGap)

    return (dupsInDist, nGeneDupl)
