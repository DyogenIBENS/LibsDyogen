# -*- coding: utf-8 -*-

# LibsDyogen version 1.0 (6/11/2015)
# python v2.7 at least is needed
# Copyright © 2015 IBENS/Dyogen : Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : jlucas@ens.fr
# Licences GLP v3 and CeCILL v2

# To use this file install homology teams, cf the README of LibsDyogen
import copy
import sys
import os
import collections
import itertools

import myTools
import myDiags
import myMapping
import mySvgDrawer
import myLightGenomes
from mySvgDrawer import Point

# see the css joint file "HomologyGroup" ranges from 0 to 44 (included)
HomologyGroupRange = range(3, 44+1)
# see the css joint file "NoHomologyInWindow" ranges from 0 to 14 (included)
NoHomologyInWindowRange = range(3, 14+1)

FilterType = myDiags.FilterType

# parse the user input (text) for the chromosome range and asses if this query is consistent with the genome data
# g2gtb is a dictionnary to convert indices from gene coordinates into tb coordinates
def parseChrRange(text, genome, g2gtb=None):
    if len(text.split(":")) == 2:
        chr = text.split(":")[0]
        if chr not in genome.keys():
            raise ValueError("chr %s not in genome" % chr)
    else:
        raise ValueError('range not formated as expected : \"chr:deb-fin\"')

    range = text.split(":")[-1].split("-")
    if len(range) == 2 and range[0].isdigit and (range[1].isdigit() or range[1] == '~'):
        if g2gtb is not None:
            assert isinstance(g2gtb, dict)
            targetIdxG1 = int(range[0]) - 1
            if targetIdxG1 in g2gtb[chr]:
                range[0] = g2gtb[chr][targetIdxG1]
            else:
                print >> sys.stderr, "Warning, gene %s at %s:%s is not part of a tandem block, thus we take the nearest tandem block" % \
                                     (genome[chr][targetIdxG1].n, chr, targetIdxG1 + 1)
                idxGeneThatIsInTargetTb1 = min(g2gtb[chr].keys(), key=lambda x: abs(x-targetIdxG1))
                print >> sys.stderr, "Warning: abs(idxGeneThatIsInTargetTb1 - targetIdxG1) = %s" % abs(idxGeneThatIsInTargetTb1 - targetIdxG1)
                range[0] = g2gtb[chr][idxGeneThatIsInTargetTb1]
            if range[1] == '~':
                range[1] = max(idxTb for idxTb in g2gtb[chr])
            else:
                targetIdxG2 = int(range[1]) - 1
                if targetIdxG2 in g2gtb[chr]:
                    range[1] = g2gtb[chr][targetIdxG2] + 1
                else:
                    print >> sys.stderr, "Warning, gene %s at %s:%s is not part of a tandem block, thus we take the nearest tandem block" % \
                                         (genome[chr][targetIdxG2].n, chr, targetIdxG2 + 1)
                    idxGeneThatIsInTargetTb2 = min(g2gtb[chr].keys(), key=lambda x: abs(x-targetIdxG2))
                    print >> sys.stderr, "Warning: abs(idxGeneThatIsInTargetTb2 - targetIdxG2) = %s" % abs(idxGeneThatIsInTargetTb2 - targetIdxG2)
                    range[1] = g2gtb[chr][idxGeneThatIsInTargetTb2] + 1
        else:
            range = (int(range[0])-1, int(range[1]) if range[1] != '~' else len(genome[chr]))
        if range[0] < 0 or len(genome[chr]) < range[1] or range[1] <= 0 or range[1] <= range[0]:
            raise ValueError("range %s is incoherent for chr %s. FYI chr %s contains %s elements. Be sure that beginning < end of range"\
                % ([range[0]+1, range[1]], chr, chr, len(genome[chr])))
    else:
        raise ValueError
    return (chr, range)


def TbComputeHomologyInformations(chrom1_tb, chrom2_tb):
    # Build MHP = { i1_tb : {i2_tb : hpSign, ...}...} with hpSign
    # the hpSign (s1*s2) of the homology at the (i1,i2) coordinate
    (MHP, locG2) = myDiags.homologyMatrix(chrom1_tb, chrom2_tb)

    tbsOnChr2ThatHaveHomologies = set([])
    for tb1 in MHP:
        tbsOnChr2ThatHaveHomologies |= set(MHP[tb1].keys())

    # Build TBNoHomologiesInWindowX = [..., i_tb, ... ] with i_tb the index of
    # the tb of chromosomeX_tb with no homology in the window
    TbNoHomologiesInWindowC1 = []
    for (i1_tb, _) in enumerate(chrom1_tb):
        if i1_tb not in MHP:
            TbNoHomologiesInWindowC1.append(i1_tb)
    TbNoHomologiesInWindowC2 = []
    for (i2_tb, _) in enumerate(chrom2_tb):
        if i2_tb not in tbsOnChr2ThatHaveHomologies:
            TbNoHomologiesInWindowC2.append(i2_tb)
    # Build homologyGroupsInWindow : [..., ([tbC1_4,tbC1_46,tbC1_80],[tbC2_2,tbC2_7]), ...]
    #   a list of 2-uples of homologous tbs indices to be coloured in the same color
    #      with tbCX_X : [..., i, ...] list of indices of genes in the tb
    locC1_tbIndx = collections.defaultdict(list)
    locC2_tbIndx = collections.defaultdict(list)
    TbHomologyGroupsInWindow = []
    listOfHPsCoordinates = []
    for i1_tb in MHP:
        for i2_tb in MHP[i1_tb]:
            listOfHPsCoordinates.append((i1_tb, i2_tb))
    for (i1_tb, i2_tb) in listOfHPsCoordinates:
        locC2_tbIndx[i1_tb].append(i2_tb)
        locC1_tbIndx[i2_tb].append(i1_tb)
    HomologyGroup1 = set()
    HomologyGroup2 = set()
    for (i1_tb, i2_tb) in listOfHPsCoordinates:
        if i1_tb not in HomologyGroup1 or i2_tb not in HomologyGroup2:
            TbHomologyGroupsInWindow.append(([], []))
            HomologyGroup1.add(i1_tb)
            for i2_tbb in locC2_tbIndx[i1_tb]:
                TbHomologyGroupsInWindow[-1][1].append(i2_tbb)
                HomologyGroup2.add(i2_tbb)
                for i1_tbb in locC1_tbIndx[i2_tbb]:
                    TbHomologyGroupsInWindow[-1][0].append(i1_tbb)
                    HomologyGroup1.add(i1_tbb)
                del locC1_tbIndx[i2_tbb]
            del locC2_tbIndx[i1_tb]
    del HomologyGroup1
    del HomologyGroup2

    return (MHP, (TbNoHomologiesInWindowC1, TbNoHomologiesInWindowC2), TbHomologyGroupsInWindow)

def editGenomes(genome1, genome2, families, filterType, minChromLength, tandemGapMax):
    g1_fID = myMapping.labelWithFamID(genome1, families)
    g2_fID = myMapping.labelWithFamID(genome2, families)

    # Must be applied on the two genomes
    if len((g1_fID.getGeneNames(checkNoDuplicates=False) & g2_fID.getGeneNames(checkNoDuplicates=False)) - {None}) == 0:
        print >> sys.stderr, "Warning, no homologs"
    ((g1_fID_filt, Gf2GfID1, (nCL1, nGL1)),
     (g2_fID_filt, Gf2GfID2, (nCL2, nGL2))) = myDiags.filter2D(g1_fID, g2_fID,
                                                               filterType,
                                                               minChromLength,
                                                               keepOriginal=True)
    (g1_tb, Gtb2GfID1, nGTD1) = myMapping.remapRewriteInTb(g1_fID_filt,
                                                             tandemGapMax=tandemGapMax,
                                                             mOld=Gf2GfID1)
    (g2_tb, Gtb2GfID2, nGTD2) = myMapping.remapRewriteInTb(g2_fID_filt,
                                                             tandemGapMax=tandemGapMax,
                                                             mOld=Gf2GfID2)
    return ((g1_tb, g1_fID, Gtb2GfID1), (g2_tb, g2_fID, Gtb2GfID2))

def computeHomologyInformations(chr1, chr2, (g1, g1_tb, g1_fID, Gtb2GfID1), (g2, g2_tb, g2_fID, Gtb2GfID2)):
    #Focus on the chromosome of the window
    chrom1_ = g1[chr1]
    chrom2_ = g2[chr2]
    c1_fID = g1_fID[chr1]
    c2_fID = g2_fID[chr2]
    chrom1_tb = g1_tb[chr1] if chr1 in g1_tb else []
    chrom2_tb = g2_tb[chr2] if chr2 in g2_tb else []

    Ctb2CfID1 = Gtb2GfID1[chr1]
    Ctb2CfID2 = Gtb2GfID2[chr2]

    # Build genesRemovedDuringFilteringCX = [..., i, ...] the list of genes that
    # have been removed during the filtering process
    genesRemovedDuringFilteringC1 = [i1 for (i1, (anc, _)) in enumerate(c1_fID) if anc is None]
    genesRemovedDuringFilteringC2 = [i2 for (i2, (anc, _)) in enumerate(c2_fID) if anc is None]

    (TbHpSign, (TbNoHomologiesInWindowC1, TbNoHomologiesInWindowC2), TbHomologyGroupsInWindow) = \
        TbComputeHomologyInformations(chrom1_tb, chrom2_tb)

    genesHomologiesHpSign = collections.defaultdict(lambda: collections.defaultdict(int))
    for i1_tb in TbHpSign:
        for i2_tb in TbHpSign[i1_tb]:
            for (i1, i2) in itertools.product([ii1 for ii1 in Ctb2CfID1[i1_tb]],
                                              [ii2 for ii2 in Ctb2CfID2[i2_tb]]):
                s1 = chrom1_[i1][1]
                s2 = chrom2_[i2][1]
                genesHomologiesHpSign[i1][i2] = s1*s2

    # Build genesNoHomologiesInWindow1 = [..., [i5,i6,i7], ...] list of tbs with [i5,i6,i7] a tb of three genes whose
    # indices are i5,i6 and i7
    genesNoHomologiesInWindowC1 = []
    genesNoHomologiesInWindowC2 = []
    for i1_tb in TbNoHomologiesInWindowC1:
        genesNoHomologiesInWindowC1.append([i1 for i1 in Ctb2CfID1[i1_tb]])
    for i2_tb in TbNoHomologiesInWindowC2:
        genesNoHomologiesInWindowC2.append([i2 for i2 in Ctb2CfID2[i2_tb]])

    # Build genesHomologyGroupsInWindow : [..., ([tbC1_4,tbC1_46,tbC1_80],[tbC2_2,tbC2_7]), ...]
    #   a list of 2-uples of homologous tbs indices to be coloured in the same color
    #      with tb_CX_X : [..., i, ...] list of indices of genes in the tb
    genesHomologyGroupsInWindow = []
    for (tbs1, tbs2) in TbHomologyGroupsInWindow:
        genesHomologyGroupsInWindow.append(([], []))
        for i1_tb_group in tbs1:
            genesHomologyGroupsInWindow[-1][0].append([i1 for i1 in Ctb2CfID1[i1_tb_group]])
        for i2_tb_group in tbs2:
            genesHomologyGroupsInWindow[-1][1].append([i2 for i2 in Ctb2CfID2[i2_tb_group]])

    return ((genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2),
            genesHomologiesHpSign,
            (genesNoHomologiesInWindowC1, genesNoHomologiesInWindowC2),
            genesHomologyGroupsInWindow)

def editGenomesAndComputeHomologyInformations(chr1, chr2, genome1, genome2, families,
                                     filterType, minChromLength, tandemGapMax):

    ((g1_tb, g1_fID, Gtb2GfID1), (g2_tb, g2_fID, Gtb2GfID2)) = editGenomes(genome1, genome2, families,
                                                                           filterType, minChromLength, tandemGapMax)
    ((genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2),
     genesHomologiesHpSign,
     (genesNoHomologiesInWindowC1, genesNoHomologiesInWindowC2),
     genesHomologyGroupsInWindow) = computeHomologyInformations(chr1, chr2,
                                                                (genome1, g1_tb, g1_fID, Gtb2GfID1),
                                                                (genome2, g2_tb, g2_fID, Gtb2GfID2))

    return ((genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2),
            genesHomologiesHpSign,
            (genesNoHomologiesInWindowC1, genesNoHomologiesInWindowC2),
            genesHomologyGroupsInWindow)

# Generator of levels for colors or gray indices within a palette:
# farIdxs may be an int. The more this int is high, the more neighbour color will be different
class LevelIdxGenerator():
    def __init__(self, farIdxs=None, inGrays=False):
        if inGrays:
            self.availableLevels = NoHomologyInWindowRange
        else:
            self.availableLevels = HomologyGroupRange

        if farIdxs is not None:
            new_availableLevels = []
            for unit in range(0, farIdxs):
                for i in range(0, len(self.availableLevels) / farIdxs + 1):
                    if i*farIdxs + unit < len(self.availableLevels):
                        new_availableLevels.append(self.availableLevels[(i*farIdxs) + unit])
            assert len(new_availableLevels) == len(self.availableLevels)
            self.availableLevels = new_availableLevels
        self.currIdx = 0

    def getLevel(self, differentFrom=set()):
        if len(differentFrom.intersection(self.availableLevels)) == len(self.availableLevels):
            print >> sys.stderr, "Warning: too many colors too avoid, thus a non-optimal choice of the color is made"
            self.currIdx = 0
        else:
            while self.currIdx in differentFrom:
                if 0 <= self.currIdx < len(self.availableLevels) - 1:
                    self.currIdx += 1
                else:
                    self.currIdx = 0
        level = self.availableLevels[self.currIdx]
        if self.currIdx < len(self.availableLevels) - 2:
            self.currIdx += 1
        else:
            self.currIdx = 0
        return level


def neighboursLevels(chromosome, i):
    # level of the 'L'eft neighbour
    levelL = chromosome[i-1].SVGclass if 0 <= i-1 <= len(chromosome)-1 else None
    # level of the 'R'ight neighbour
    levelR = chromosome[i+1].SVGclass if 0 <= i+1 <= len(chromosome)-1 else None

    def convertSVGclassIntoInt(strSVGclass):
        if strSVGclass is not None:
            try:
                return int(strSVGclass[-2:])
            except:
                if strSVGclass == 'SpeciesSpecificGenes':
                    return None
                else:
                    return int(strSVGclass[-1])

    levelL = convertSVGclassIntoInt(levelL)
    levelR = convertSVGclassIntoInt(levelR)
    neighboursLevels = {levelL, levelR}
    neighboursLevels = set([l for l in neighboursLevels if l is not None])
    return neighboursLevels


def drawChromosomeFromRawInformations(genesStrands,
                                      # coloredName2Idxs = {colorName: {..., gIdx, ...}, ...}
                                      coloredName2Idxs=None,
                                      name2color=None,
                                      symbolsInGenes=None,
                                      lengthGenes=0.7,
                                      halfIntergeneLengths=0.15,
                                      desiredLength=None,
                                      # In a MH, greyIdToGeneIdxs are the genes of a chromosome with no homolog in the compared chromosome
                                      # {greyId: {idxG, ...}, ...}, each list is colored with the same level of grey
                                      greyIdToGeneIdxs=None,
                                      # In MH, white genes are species specific genes that in non-ancestral gene families
                                      whiteGenesIdxs=None,
                                      greyLevelsGenerator=None,
                                      geneWidth=None):
    if greyIdToGeneIdxs:
        greyGidxs = set([gIdx for gIdxs in greyIdToGeneIdxs.values() for gIdx in gIdxs])
    if whiteGenesIdxs and greyIdToGeneIdxs:
        assert len(set(whiteGenesIdxs) &  greyGidxs) == 0
    if coloredName2Idxs:
        coloredGidx = set([gIdx for gIdxs in coloredName2Idxs.values() for gIdx in gIdxs])
    if whiteGenesIdxs and coloredName2Idxs:
        assert len(set(whiteGenesIdxs) & coloredGidx) == 0
    if greyIdToGeneIdxs and coloredName2Idxs:
        assert len(greyGidxs & coloredGidx) == 0
    if whiteGenesIdxs and greyIdToGeneIdxs and coloredName2Idxs:
        assert whiteGenesIdxs | greyGidxs | coloredGidx ==  set(range(len(genesStrands)))

    differentGeneLengths = False
    differentHalfIntergeneLengths = False
    if isinstance(lengthGenes, list):
        assert len(lengthGenes) == len(genesStrands)
        differentGeneLengths = True
    if isinstance(halfIntergeneLengths, list):
        assert len(halfIntergeneLengths) == 2 * len(genesStrands)
        differentHalfIntergeneLengths = True

    def giveGreyLevelsTo(chromosome, greyIdToGeneIdxs):
        for (greyId, geneIdxs) in greyIdToGeneIdxs.iteritems():
            nLevels = reduce(lambda x, y: x | y, [neighboursLevels(chromosome, gIdx) for gIdx in geneIdxs])
            for gIdx in geneIdxs:
                # Choose a level different from the direct neighbours
                grey = greyLevelsGenerator.getLevel(differentFrom=nLevels)
                chromosome[gIdx].SVGclass = "NoHomologyInWindow%s" % grey
        return chromosome

    chromosomeItems = []

    # homologousTbs = [homolog1[tb1=[gene1Idx, ...], tb2=[geneAIdx, ...]],
    #                  homolog2[tb1'=[gene1'Idx, ...], tb2'=[geneA'Idx, ...]],
    #                 ... ]
    if not greyLevelsGenerator:
        greyLevelsGenerator = LevelIdxGenerator(farIdxs=None, inGrays=True)
    if not name2color:
        homologsColorsGenerator = LevelIdxGenerator(farIdxs=None)

    # create chromosomes
    if differentGeneLengths:
        assert isinstance(lengthGenes, list)
        cumGeneLengths = sum(lengthGenes)
    else:
        cumGeneLengths = len(genesStrands) * lengthGenes
    if differentHalfIntergeneLengths:
        cumIntergeneLengths = sum(halfIntergeneLengths)
    else:
        cumIntergeneLengths = (2 * halfIntergeneLengths) * len(genesStrands)
    # width = lengthGenes * 0.7
    # # stroke_width  = 0.05 * lengthGenes
    # stroke_width  = 0.05 * width
    chromLength = cumGeneLengths + cumIntergeneLengths
    if desiredLength is not None:
        assert isinstance(desiredLength, float) or isinstance(desiredLength, int)
        factorLengh = desiredLength / float(chromLength)
        chromLength *= factorLengh
        if differentGeneLengths:
            lengthGenes = [lengthGene * factorLengh for lengthGene in lengthGenes]
        else:
            lengthGenes *= factorLengh
        if differentHalfIntergeneLengths:
            halfIntergeneLengths = [halfIntergeneLength * factorLengh for halfIntergeneLength in halfIntergeneLengths]
        else:
            halfIntergeneLengths *= factorLengh

    if not geneWidth:
        geneWidth = chromLength / 250.0
    geneStrokeWidth = geneWidth / 10.0
    cx = 0
    for (i, s) in enumerate(genesStrands):
        if differentGeneLengths:
            lengthGene = lengthGenes[i]
        else:
            lengthGene = lengthGenes
        if differentHalfIntergeneLengths:
            (intergeneLeft, intergeneRight) = (halfIntergeneLengths[2*i], halfIntergeneLengths[2*i+1])
        else:
            (intergeneLeft, intergeneRight) = (halfIntergeneLengths, halfIntergeneLengths)
        symbol = symbolsInGenes[i] if symbolsInGenes is not None else None
        chromosomeItems.append(mySvgDrawer.Gene(Point(cx + intergeneLeft, 0),
                                                Point(cx + intergeneLeft + lengthGene, 0),
                                                strand=s, width=geneWidth, stroke_width=geneStrokeWidth, SVGclass=None,
                                                text=symbol))
        cx += intergeneLeft + lengthGene + intergeneRight

    if coloredName2Idxs:
        for (coloredName, idxGs) in coloredName2Idxs.iteritems():
            if name2color:
                color = name2color[coloredName]
            else:
                nLevels = reduce(lambda x, y: x | y, [neighboursLevels(chromosomeItems, idxG) for idxG in idxGs])
                color = homologsColorsGenerator.getLevel(differentFrom=nLevels)
            for idxG in idxGs:
                # for idxG in idxTb:
                chromosomeItems[idxG].SVGclass = "HomologGroup%s" % color

    # give grey levels to genes that have no homology in the window
    # this step has to be done after colouring
    if greyIdToGeneIdxs:
        chromosomeItems = giveGreyLevelsTo(chromosomeItems, greyIdToGeneIdxs)

    if whiteGenesIdxs:
        for i in whiteGenesIdxs:
            chromosomeItems[i].SVGclass = "SpeciesSpecificGenes"

    # add a line that goes through all genes of the chromosome
    chromosomeItems = [mySvgDrawer.Line(Point(0, 0), Point(chromLength, 0))] + chromosomeItems
    if desiredLength is not None:
        return (chromosomeItems, factorLengh)
    else:
        return chromosomeItems

# TODO: implement all the basic classes of myLigthGenomes.LightGenome into myLightGenome.Chromosome for more modularity
# then in myLigthGenomes.LightGenome, call the classes of myLigthGenomes.Chromosome whenever it is possible
def svgItemsChromosome(chromosome,
                       families=None,
                       name2color=None,
                       lengthGenes=1.0,
                       genesCoordinates=None,
                       halfIntergeneLengths=0.1,
                       desiredLength=None,
                       geneNamesWithSameGreyLevels=None,
                       greyLevelsGenerator=None,
                       whiteGeneNames=None,
                       symbolsInGenes=None,
                       colorsGenerator=None,
                       geneWidth=None):
    genomeWithOnlyChrom = myLightGenomes.LightGenome()
    chr = '0'
    genomeWithOnlyChrom[chr] = chromosome
    if lengthGenes is not None:
        assert genesCoordinates is None
        if isinstance(lengthGenes, list):
            lengthGenes = {chr: lengthGenes}
        else:
            assert isinstance(lengthGenes, float)
    if genesCoordinates is not None:
        assert lengthGenes is None
        genesCoordinates = {chr: genesCoordinates}
    if isinstance(halfIntergeneLengths, list):
        halfIntergeneLengths = {chr: halfIntergeneLengths}
    if isinstance(symbolsInGenes, list):
        symbolsInGenes = {chr: symbolsInGenes}
    res = svgItemsLightGenome(genomeWithOnlyChrom,
                               families=families,
                               name2color=name2color,
                               lengthGenes=lengthGenes,
                               genesCoordinates=genesCoordinates,
                               halfIntergeneLengths=halfIntergeneLengths,
                               desiredLength=desiredLength,
                               geneNamesWithSameGreyLevels=geneNamesWithSameGreyLevels,
                               greyLevelsGenerator=greyLevelsGenerator,
                               whiteGeneNames=whiteGeneNames,
                               symbolsInGenes=symbolsInGenes,
                               colorsGenerator=colorsGenerator,
                               geneWidth=geneWidth)
    if desiredLength is not None:
        assert  isinstance(res, tuple), '%s' % res
        assert len(res) == 2
        (genomeItems, factorLength) = res
        # FIXME, factorLengthShould be a dict over chrs keys
        # factorLength = desiredLength / float(chromLength)
        return (genomeItems[chr], factorLength)
    else:
        assert len(res) == 1
        genomeItems = res
        return genomeItems[chr]

def genesCoordinates2GeneAndHalfIntergeneLengths(listOfGenesCoordinates):
    chrGeneLengths = []
    chrHalfIntergeneLengths = []
    beg_old = 0
    end_old = 0
    for (idxG, (beg, end)) in enumerate(listOfGenesCoordinates):
        # DEBUG
        print >> sys.stderr, str((beg_old, end_old, beg, end))
        # 1) take the coordinates (beg, end) of the genes and edit the coordinates if they overlap
        if idxG > 0: assert beg_old < end_old
        assert beg < end
        #(beg, end) = (end, beg) if end < beg else (beg, end)
        if end_old <= beg:
            if idxG > 0: assert beg_old < end_old <= beg < end
        else:
            assert beg < end_old
            # overlap!
            print >> sys.stderr, 'Warning: Overlapping genes extremities beg_right=%s < end_left=%s' % (beg, end_old)
            print >> sys.stderr, 'This should not happen if genes were indexed on there 5\' extremities and if they do not overlap'
            # The only solution to get out of this mess is to truncated the former gene
            if end_old < end:
                if (beg_old < beg < end_old < end):
                    # FIXME arbitrary choice : could be (beg_old, end_old - 0.0000000001, end_old, end)
                    (beg_old, end_old, beg, end) = (beg_old, beg - 0.0000000001, beg, end)
                else:
                    assert (beg < beg_old < end_old < end), '%s' % str((beg_old, beg, end_old, end))
                    (beg_old, end_old, beg, end) = (beg_old, end_old - 0.0000000001, end_old, end)
            else:
                print >> sys.stderr, 'Warning: genes are nested!'
                assert beg_old < beg < end <= end_old
                # the current gene is nested in the previous beg_old
                (beg_old, end_old, beg, end) = (beg_old, beg - 0.0000000001, end, end_old)
            # edit the former gene length
            if idxG > 0: assert beg_old < end_old <= beg < end
            print >> sys.stderr, 'Warning: because of overlapping genes, the gene coordinates have been edited'
            listOfGenesCoordinates[idxG - 1] = (beg_old, end_old)
            listOfGenesCoordinates[idxG] = (beg, end)
            chrGeneLengths[idxG - 1] = abs(beg_old - end_old)

        # 2) compute gene lengths
        chrGeneLengths.append(end - beg)

        # 3) compute halfIntergenes
        if idxG == 0:
            leftTelomereLength = beg
            chrHalfIntergeneLengths.append(leftTelomereLength)
        else:
            assert 1 <= idxG < len(listOfGenesCoordinates)
            assert end_old <= beg
            leftIntergeneLength = abs(end_old - beg)
            halfLeftIntergeneLength = leftIntergeneLength / 2.0
            chrHalfIntergeneLengths.extend([halfLeftIntergeneLength, halfLeftIntergeneLength])
            if idxG == len(listOfGenesCoordinates) - 1:
                rightTelomereLength = leftTelomereLength
                chrHalfIntergeneLengths.append(rightTelomereLength)
        (beg_old, end_old) = (beg, end)

    # DEBUG ASSERTION
    _cumChrLength = 0
    assert len(chrGeneLengths) == len(listOfGenesCoordinates)
    for (i, geneLength) in enumerate(chrGeneLengths):
        (beg, end) = listOfGenesCoordinates[i]
        assert beg < end
        if i == 0:
            leftTelomereLength = chrHalfIntergeneLengths[0]
            _cumChrLength += leftTelomereLength
            assert abs(beg - _cumChrLength) < 0.0001
            _cumChrLength += geneLength
            assert abs(end - _cumChrLength) < 0.0001
        elif 1 <= i < len(chrGeneLengths):
            # add length of ith intergene
            leftInterGeneLength = sum(chrHalfIntergeneLengths[i*2-1:i*2+1])
            _cumChrLength += leftInterGeneLength
            assert abs(beg - _cumChrLength) < 0.0000001, '%s ~= %s' % (beg, _cumChrLength)
            _cumChrLength += geneLength
            assert abs(end - _cumChrLength) < 0.0000001, '%s ~= %s' % (end, _cumChrLength)
    return (chrGeneLengths, chrHalfIntergeneLengths)


def svgItemsLightGenome(genome,
                        families=None,
                        # name2color = {famName: color, ...}
                        # or {geneName: color, ...} if no families
                        name2color=None,
                        # lengthGenes may also be a dict: {chr1: [lengthGene1, ....], ...}
                        lengthGenes=1.0,
                        # genesCoordinates: dict {chr1: [(begG1, endG1), ...], ...}, begG1 might be a length in bp
                        # !!! WARNING: this object might be edited if genes overlap !!!
                        genesCoordinates=None,
                        # halfIntergeneLengths may also be a dict {chr1: [telomere1, 1sthalfIntergene1, 2ndhalfIntergene1, 1sthalfIntergene2, ...], ...}
                        halfIntergeneLengths=0.1,
                        # the desired length of the genome
                        desiredLength=None,
                        # In a MH, geneNamesWithSameGreyLevels are the genes of a chromosome with no homolog in the compared chromosome
                        # geneNamesWithSameGreyLevels = [..., {geneName, ...}, ...]
                        geneNamesWithSameGreyLevels=None,
                        # In MH, white genes are species specific genes that in non-ancestral gene families
                        greyLevelsGenerator=None,
                        # whiteGeneNames = {..., geneName, ...}
                        whiteGeneNames=None,
                        # symbols in genes before and after gene names: symbolsInGenes={chr: [..., symbolStr, ...], ...}
                        symbolsInGenes=None,
                        # if name2color is not provided, this generator defines the colors
                        colorsGenerator=None,
                        geneWidth=None):
    assert isinstance(genome, myLightGenomes.LightGenome)
    assert families is None or isinstance(families, myLightGenomes.Families)
    assert (lengthGenes is not None) or (genesCoordinates is not None) and not (lengthGenes and genesCoordinates)
    if genesCoordinates:
        assert isinstance(genesCoordinates, dict)
        assert lengthGenes is None
    if geneNamesWithSameGreyLevels:
        if len(geneNamesWithSameGreyLevels) > 0:
            assert isinstance(geneNamesWithSameGreyLevels[0], set)

    if whiteGeneNames or geneNamesWithSameGreyLevels:
        genome.computeDictG2Ps()
    whiteGenesIdxs = collections.defaultdict(set)
    if whiteGeneNames:
        for whiteGeneName in whiteGeneNames:
            whiteGenesPositions = genome.getPositions(whiteGeneName, default=None)
            if whiteGenesPositions is not None:
                for whiteGenePosition in whiteGenesPositions:
                    whiteGenesIdxs[whiteGenePosition.c].add(whiteGenePosition.idx)
    idxsOfGreyGenes = collections.defaultdict(set)
    greyIdToGeneIdxs = collections.defaultdict(lambda: collections.defaultdict(set))
    if geneNamesWithSameGreyLevels:
        greyIdToGeneIdxs = collections.defaultdict(lambda: collections.defaultdict(set))
        for (iGrey, geneNamesWithSameGreyLevel) in enumerate(geneNamesWithSameGreyLevels):
            for greyGeneName in geneNamesWithSameGreyLevel:
                greyGenesPositions = genome.getPositions(greyGeneName, default=None)
                if greyGenesPositions is not None:
                    for greyGenePosition in greyGenesPositions:
                        greyIdToGeneIdxs[greyGenePosition.c][iGrey].add(greyGenePosition.idx)
                        idxsOfGreyGenes[greyGenePosition.c].add(greyGenePosition.idx)

    # rewrite with families
    if families:
        genome_fam = myMapping.labelWithFamNames(genome, families, keepGnOfGenesNotInFamilies=True)
        #geneNamesOfColoredOrGreyGenes = genome_fam.getOwnedFamilyNames(families, asA=set)
        geneNamesOfColoredOrGreyGenes = genome_fam.getGeneNames(checkNoDuplicates=False, asA=set)
    else:
        genome_fam = genome
        geneNamesOfColoredOrGreyGenes = genome_fam.getGeneNames(checkNoDuplicates=True, asA=set)

    genome_fam.computeDictG2Ps()
    coloredName2Idxs = collections.defaultdict(lambda: collections.defaultdict(set))
    for geneName in geneNamesOfColoredOrGreyGenes:
        genesPositions = genome_fam.getPositions(geneName, default=None)
        for genePosition in genesPositions:
            if genePosition.idx not in idxsOfGreyGenes[genePosition.c] | whiteGenesIdxs[genePosition.c]:
                coloredName2Idxs[genePosition.c][geneName].add(genePosition.idx)

    if not name2color:
        if not colorsGenerator:
            colorsGenerator = LevelIdxGenerator(farIdxs=5)
        name2color = {}
        for geneName in geneNamesOfColoredOrGreyGenes:
            name2color[geneName] = colorsGenerator.getLevel()

    genomeItems = collections.OrderedDict()
    for (chr, chrom) in genome.iteritems():
        genesStrands = [s for (_, s) in chrom]
        if genesCoordinates is None:
            chrGeneLengths = lengthGenes[chr] if isinstance(lengthGenes, dict) else lengthGenes
            chrHalfIntergeneLengths = halfIntergeneLengths[chr] if isinstance(halfIntergeneLengths, dict) else halfIntergeneLengths
        else:
            # This step may edit the genesCoordinates[chr] if genes overlap
            (chrGeneLengths, chrHalfIntergeneLengths) = genesCoordinates2GeneAndHalfIntergeneLengths(genesCoordinates[chr])

        res = drawChromosomeFromRawInformations(genesStrands,
                                                coloredName2Idxs=coloredName2Idxs[chr],
                                                name2color=name2color,
                                                symbolsInGenes=symbolsInGenes[chr] if symbolsInGenes else None,
                                                lengthGenes=chrGeneLengths,
                                                halfIntergeneLengths=chrHalfIntergeneLengths,
                                                desiredLength=desiredLength,
                                                greyIdToGeneIdxs=greyIdToGeneIdxs[chr],
                                                whiteGenesIdxs=whiteGenesIdxs[chr],
                                                greyLevelsGenerator=greyLevelsGenerator,
                                                geneWidth=geneWidth)
        if desiredLength is not None:
            assert isinstance(res, tuple)
            assert len(res) == 2
            (genomeItems[chr], factorLength) = res
        else:
            genomeItems[chr] = res

    if desiredLength is not None:
        return (genomeItems, factorLength)
    else:
        return genomeItems

def svgItemsTicks(range1, sizeCell, widthTicks=None, sizeText=None, rankGenome=1):
    # height = 2 * sizeCell
    nx = range1[1] - range1[0]
    width = nx * sizeCell
    if not widthTicks:
        widthTicks = float(width)/1000
    if not sizeText:
        sizeText = widthTicks*10

    listOfTickItems = []
    for (i, ni) in enumerate(range(range1[0], range1[1])):
        rankI = ni + 1
        cx = i * sizeCell + sizeCell / 2.0
        if rankGenome == 1:
            if rankI % 10 == 0:
                # ticks
                listOfTickItems.append(mySvgDrawer.Line(Point(cx, 0.0), Point(cx, sizeCell/2.0), width=widthTicks))
            if rankI % 50 == 0:
                # text + large ticks
                listOfTickItems.append(mySvgDrawer.Line(Point(cx, 0.0), Point(cx, sizeCell/2.0), width=widthTicks * 3))
                listOfTickItems.append(mySvgDrawer.Text(Point(cx, sizeCell), str(rankI), text_anchor="middle", size=sizeText))

        else:
            assert rankGenome == 2
            if rankI % 10 == 0:
                listOfTickItems.append(mySvgDrawer.Line(Point(cx, sizeCell * 3.0/2.0), Point(cx, 2.0 * sizeCell), width=widthTicks))
            if rankI % 50 == 0:
                listOfTickItems.append(mySvgDrawer.Text(Point(cx, sizeCell), str(rankI), text_anchor="middle", size=sizeText))
                listOfTickItems.append(mySvgDrawer.Line(Point(cx, sizeCell * 3.0/2.0), Point(cx, 2.0 * sizeCell), width=widthTicks * 3))
    return listOfTickItems

def svgItemsMatrixFromRawData(range1, range2,
                              coordsHomologiesAndSymbols,
                              sizeCell,
                              coordsOfHomologiesPerDiags=None,
                              diagNames=None,
                              diagColorGenerator=None,
                              drawSmallLines=True,
                              drawLargeLines=True,
                              drawHomologySymbols=True,
                              scaleFactorRectangles=1.0):
    assert isinstance(coordsHomologiesAndSymbols, dict)
    assert isinstance(scaleFactorRectangles, float) or isinstance(scaleFactorRectangles, int)
    sizeText = float(sizeCell*0.9)
    listOfMatrixItems = []
    if not diagColorGenerator:
        diagColorGenerator = LevelIdxGenerator(farIdxs=7)

    nx = range1[1] - range1[0]
    # Nb of vertical lines (x varies) in the matrix
    nbLinesX = nx + 1

    ny = range2[1] - range2[0]
    # Nb of horizontal lines (y varies) in the matrix
    nbLinesY = ny + 1

    width = nx * sizeCell
    height = ny * sizeCell

    smallestVisibleWidthOnFirefox = 1.0 / 100.0
    # conversion either with inkscape or rsvg-convert
    smallestVisibleWidthForGoodConversionIntoPDF = 1.0 / 10.0
    smallestVisibleWidth=smallestVisibleWidthForGoodConversionIntoPDF
    strokeWidthHomologies = 0 if sizeCell * scaleFactorRectangles < 10 * smallestVisibleWidth else smallestVisibleWidth
    # draw Diagonals first because they are on the background
    print >> sys.stderr, "Nb of diagonals showed = ", len(coordsOfHomologiesPerDiags)

    # sort diagonals to give colors according to localisations of diagonals

    # FIXME Do not sort!!  otherwise the diagNames do not correspond any more!!!
    # coordsOfHomologiesPerDiags.sort(key=lambda x: x[0])

    invYaxis = lambda y: height - (y + sizeCell)

    for (idiag, diag) in enumerate(coordsOfHomologiesPerDiags):
        # choose a color different from the neighbours
        color = diagColorGenerator.getLevel()
        for (i, j) in diag:
            cx_s = i * sizeCell
            cy_s = j * sizeCell
            if nx >= 300 or ny >= 300:
                listOfMatrixItems.append(mySvgDrawer.Rectangle(Point(cx_s, invYaxis(cy_s)),
                                                               sizeCell, sizeCell, fill_opacity=0.90,
                                                               svgClass="HomologGroup%s" % color,
                                                               bigger=scaleFactorRectangles,
                                                               strokeWidth=strokeWidthHomologies))
            else:
                listOfMatrixItems.append(mySvgDrawer.Rectangle(Point(cx_s, invYaxis(cy_s)),
                                                               sizeCell, sizeCell, fill_opacity=0.90,
                                                               svgClass="HomologGroup%s" % color))
        # draw rectangles around diagonals
        min_i = min(diag, key=lambda x: x[0])[0]
        max_i = max(diag, key=lambda x: x[0])[0]
        min_j = min(diag, key=lambda x: x[1])[1]
        max_j = max(diag, key=lambda x: x[1])[1]
        cx_s = min_i * sizeCell
        cy_s = max_j * sizeCell

        homologyWidth = (sizeCell * scaleFactorRectangles)
        strokeWidthBoundingBoxDiags = smallestVisibleWidth
        listOfMatrixItems.append(mySvgDrawer.Rectangle(Point(cx_s, invYaxis(cy_s)),
                                                       (max_j-min_j)*sizeCell + sizeCell, (max_i-min_i)*sizeCell + sizeCell,
                                                       stroke='black', fill='none', strokeWidth=strokeWidthBoundingBoxDiags))
        # need to do it after to be on top layer
        if diagNames is not None:
            widthText = min(float(width)/100, float(height)/100)
            # write the id of the diagonal if any
            listOfMatrixItems.append(mySvgDrawer.Text(Point(cx_s, invYaxis(cy_s)), str(diagNames[idiag]), size=widthText))

    if drawSmallLines:
        for i in range(nbLinesX):
            cxLine = i * sizeCell
            listOfMatrixItems.append(mySvgDrawer.Line(Point(cxLine, height),
                                                      Point(cxLine, 0), width=sizeCell*0.01))
        for j in range(nbLinesY):
            cyLine = j * sizeCell
            listOfMatrixItems.append(mySvgDrawer.Line(Point(0, invYaxis(cyLine - sizeCell)),
                                                      Point(width, invYaxis(cyLine - sizeCell)), width=sizeCell*0.01))
    if drawLargeLines:
        for i in range(nx):
            cxLine = i * sizeCell
            if (i + range1[0] + 1) % 50 == 0:
                listOfMatrixItems.append(mySvgDrawer.Line(Point(cxLine + sizeCell/2.0, height),
                                                          Point(cxLine + sizeCell/2.0, 0), width=sizeCell*0.05))
        for j in range(ny):
            cyLine = j * sizeCell
            if (j + range2[0] + 1) % 50 == 0:
                listOfMatrixItems.append(mySvgDrawer.Line(Point(0.0, invYaxis(cyLine) + sizeCell/2.0),
                                                          Point(width, invYaxis(cyLine) + sizeCell/2.0), width=sizeCell*0.05))

    # fill homologies
    nonZeroValues = []
    allCoordsOfHomologiesInDiags = set([dot for diag in coordsOfHomologiesPerDiags for dot in diag])
    for i1 in coordsHomologiesAndSymbols:
        for i2 in coordsHomologiesAndSymbols[i1]:
            cx = i1 * sizeCell
            cy = i2 * sizeCell
            if drawHomologySymbols:
                nonZeroValues.append((i1, i2))
                s = coordsHomologiesAndSymbols[i1][i2]
                assert s == +1 or s == -1 or s is None, "s=%s" % s
                assocValue = (("+" if s == +1 else "-") if s is not None else '?')
                listOfMatrixItems.append(mySvgDrawer.Text(Point(cx + sizeCell/2.0, invYaxis(cy + sizeText*0.16 - sizeCell / 2.0)),
                                                          assocValue, text_anchor="middle", size=sizeText))
            else:
                # represent homologies not in diags with a black rectangle
                if (i1, i2) not in allCoordsOfHomologiesInDiags:
                    listOfMatrixItems.append(mySvgDrawer.Rectangle(Point(cx, invYaxis(cy)),
                                                                   sizeCell, sizeCell, fill=(0, 0, 0), fill_opacity=0.90,
                                                                   bigger=scaleFactorRectangles,
                                                                   strokeWidth=strokeWidthHomologies))
    if nx < 20 and ny < 20:
        for (i1, i2) in itertools.product(range(nx), range(ny)):
            if (i1, i2) not in nonZeroValues:
                cx = i1 * sizeCell
                cy = i2 * sizeCell
                listOfMatrixItems.append(mySvgDrawer.Text(Point(cx + sizeCell/2.0, invYaxis(cy + sizeText*0.16 - sizeCell / 2.0 )),
                                                          "0", text_anchor="middle", size=sizeText,
                                                          fill=(200, 200, 200), stroke=None))
    return (listOfMatrixItems, (width, height))

def drawWholeGenomeHomologyMatrices(genome1, genome2, families,
                                    inSbsInPairComp=None,
                                    filterType=FilterType.InBothGenomes,
                                    minChromLength=0,
                                    tandemGapMax=0,
                                    scaleFactorRectangles=4.0,
                                    maxWidth=100, maxHeight=100,
                                    outputFileName=None,
                                    fillCompsWithSbs=False):

    WHM = prepareWholeGenomeHomologyMatrices(genome1, genome2, families,
                                             inSbsInPairComp=inSbsInPairComp,
                                             filterType=filterType,
                                             minChromLength=minChromLength,
                                             tandemGapMax=tandemGapMax,
                                             scaleFactorRectangles=scaleFactorRectangles,
                                             maxWidth=maxWidth, maxHeight=maxHeight,
                                             fillCompsWithSbs=fillCompsWithSbs)

    scene = mySvgDrawer.Scene(name='homology_matrix', width=maxWidth, height=maxHeight)
    for item in WHM:
            scene.add(item)
    if outputFileName:
        scene.write_svg(filename=str(outputFileName))
    return scene.strarray()

def writeSVGFileForPairwiseCompOfGenomes(genomeName1,
                                         genomeName2,
                                         HMs,
                                         chromosomesRewrittenInTbs,
                                         filterType,
                                         tandemGapMax,
                                         gapMax,
                                         distanceMetric,
                                         gapMaxMicroInv,
                                         identifyMicroRearrangements,
                                         truncationMax,
                                         nbSbs,
                                         outImageFileName,
                                         switchOnDirectView):
    # KEEP
    # copy the css style sheet
    # dirNameImage = os.path.dirname(outImageFileName)
    # dirNameImage = dirNameImage if dirNameImage != "" else "."
    # print >> sys.stderr, "cp %s/styleForHomologyMatrixWithSBs.css %s/" % (os.path.dirname(os.path.realpath(sys.argv[0])), dirNameImage)
    # os.system("cp %s/styleForHomologyMatrixWithSBs.css %s/" % (os.path.dirname(os.path.realpath(sys.argv[0])), dirNameImage))

    # Header of the SVG file
    height = 100
    width = 100
    var = ['<?xml version="1.0" encoding="utf-8" standalone="no"?>\n',
           # KEEP
           # Warning : requests the css file
           #'<?xml-stylesheet type="text/css" href="styleForHomologyMatrixWithSBs.css" ?>\n',
           '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n',
           "<svg height=\"100%%\" version=\"1.1\" viewBox=\"0 0 %s %s\" width=\"100%%\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n" % (width, height),
           '<defs>\n',
              '<style type="text/css">\n',
                 '*{stroke-linecap:square;stroke-linejoin:round;}\n',
              '</style>\n',
           '</defs>\n']

    # Title
    title = \
        "%s, f=%s, tgm=%s tbs, gm=%s%s, gmmi=%s, ibwg=%s, tm=%s, %s sbs" % \
        ('MHPs' if chromosomesRewrittenInTbs else 'MHs',
         filterType,
         tandemGapMax,
         gapMax,
         distanceMetric,
         gapMaxMicroInv,
         identifyMicroRearrangements,
         truncationMax,
         nbSbs)
    var += ['<svg x="5" y="0" viewBox="5 0 95 5" width="95" height="5">\n',
            ''.join(mySvgDrawer.Text(mySvgDrawer.Point(float(5 + 95)/2.0, 5.0/2.0), title,
                                     text_anchor='middle', fontWeight=300, size=2).strarray()),
            '</svg>\n']

    # Add legends (genomes names)
    var += [
        # X-axis genome1
        '<svg x="5" y="95" viewBox="0 0 95 5" width="95" height="5">\n',
        ''.join(mySvgDrawer.Text(mySvgDrawer.Point(float(5 + 95)/2, 5.0/2.0),
                                 str(genomeName1),
                                 text_anchor='middle', fontWeight=300, size=2).strarray()),
        '</svg>\n',
        # Y-axis genome2
        '<svg x="0" y="0" viewBox="0 0 5 95" width="5" height="95">\n',
        ''.join(mySvgDrawer.Text(mySvgDrawer.Point(5.0/2.0, float(5 + 95)/2.0),
                                 str(genomeName2),
                                 text_anchor='middle', size=2, fontWeight=300,
                                 transform='rotate(-90, %s, %s)' % (5.0/2.0, float(5 + 95)/2.0)).strarray()),
        '</svg>\n'
    ]

    # Add the matrix of homologies (packs)
    # little transformation : viewBox = "the part of the forthcoming images that we want to see", width and height = the
    # width and height of the image that will be printed on the screen. This instructions takes a viewBox of the
    # forthcoming images and display it in a image of the specified width and height
    var += ['<svg preserveAspectRatio="xMidYMid meet" x="5" y="5" viewBox="0 0 100 100" width="90" height="90" >\n']
    for line in HMs:
        # or line.find("<svg")>=0 or line.find("</svg")>=0:
        if line.find("<?xml") >= 0 or line.find("<!DOCTYPE") >= 0:
            continue
        else:
            var += line
    var += ["</svg>\n"]
    var += ["</svg>\n"]

    file = open(outImageFileName, 'w')
    file.writelines(var)
    file.close()
    if switchOnDirectView:
        os.system("%s %s" % ('firefox', outImageFileName))


def writeSVGFileForPairwiseCompOfChrs(genomeName1, chr1, range1,
                                      genomeName2, chr2, range2,
                                      pictureHM, widthHM, heightHM,
                                      chromosomesRewrittenInTbs,
                                      filterType,
                                      tandemGapMax,
                                      gapMax,
                                      distanceMetric,
                                      gapMaxMicroInv,
                                      identifyMicroRearrangements,
                                      truncationMax,
                                      nbSbs,
                                      outImageFileName,
                                      switchOnDirectView):
    # 0.02 and 0.1 are in order to not truncate strokes at the borders of the scene, and to see text that may be larger than expected
    widthHMPicture = widthHM + widthHM * 0.04
    heightHMPicture = heightHM + heightHM * 0.04
    scene = mySvgDrawer.Scene(width=widthHMPicture, height=heightHMPicture)
    mySvgDrawer.translateItems(pictureHM, Point(widthHM * 0.02, heightHM * 0.02))
    for item in pictureHM:
        scene.add(item)
    pictureHM = scene.strarray()

    # KEEP
    # copy the css style sheet
    # dirNameImage = os.path.dirname(outImageFileName)
    # dirNameImage = dirNameImage if dirNameImage != "" else "."
    # print >> sys.stderr, "cp %s/styleForHomologyMatrixWithSBs.css %s/" % (os.path.dirname(os.path.realpath(sys.argv[0])), dirNameImage)
    # os.system("cp %s/styleForHomologyMatrixWithSBs.css %s/" % (os.path.dirname(os.path.realpath(sys.argv[0])), dirNameImage))

    # Header of the SVG file
    height = 100
    width = 100
    # this defines the width and height of the bounding box
    widthViewBox = widthHMPicture
    heightViewBox = heightHMPicture
    # these coordinates are in the scale defined by widthViewBox and heightViewBox
    xViewBox = 0
    yViewBox = 0


    var = ['<?xml version="1.0" encoding="utf-8" standalone="no"?>\n',
           # KEEP
           # Warning : requests the css file
           #'<?xml-stylesheet type="text/css" href="styleForHomologyMatrixWithSBs.css" ?>\n',
           '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n',
           "<svg height=\"100%%\" version=\"1.1\" viewBox=\"0 0 %s %s\" width=\"100%%\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n" % (width, height),
           '<defs>\n',
              '<style type="text/css">\n',
                 '*{stroke-linecap:square;stroke-linejoin:round;}\n',
              '</style>\n',
           '</defs>\n'
           '<g style="fill-opacity:1.0; stroke:black; stroke-width:1;">\n']

    # Title
    title = \
        "%s, f=%s, tgm=%s tbs, gm=%s%s, gmmi=%s, ibwg=%s, tm=%s, %s sbs" % \
        ('MHP' if chromosomesRewrittenInTbs else 'MH',
         filterType,
         tandemGapMax,
         gapMax,
         distanceMetric,
         gapMaxMicroInv,
         identifyMicroRearrangements,
         truncationMax,
         nbSbs)
    var += ['<svg x="5" y="0" viewBox="5 0 95 5" width="95" height="5">\n',
            ''.join(mySvgDrawer.Text(mySvgDrawer.Point(float(5 + 95)/2.0, 5.0/2.0), title,
                                     text_anchor='middle', fontWeight=300, size=2).strarray()),
            '</svg>\n']

    # Add legends (genomes names and ranges)
    var += [
        # X-axis genome1
        '<svg x="5" y="95" viewBox="0 0 95 5" width="95" height="5">\n',
        ''.join(mySvgDrawer.Text(mySvgDrawer.Point(float(5 + 95)/2, 5.0/2.0),
                                 str(genomeName1) + " chr" + str(chr1) + ":" + str(range1[0]+1) + "-" + str(range1[1]),
                                 text_anchor='middle', fontWeight=300, size=2).strarray()),
        '</svg>\n',
        # Y-axis genome2
        '<svg x="0" y="0" viewBox="0 0 5 95" width="5" height="95">\n',
        ''.join(mySvgDrawer.Text(mySvgDrawer.Point(5.0/2.0, float(5 + 95)/2.0),
                                 str(genomeName2) + " chr" + str(chr2) + ":" + str(range2[0]+1) + "-" + str(range2[1]),
                                 text_anchor='middle', size=2, fontWeight=300,
                                 transform='rotate(-90, %s, %s)' % (5.0/2.0, float(5 + 95)/2.0)).strarray()),
        '</svg>\n'
    ]

    # Add the matrix of homologies (packs)
    # little transformation : viewBox = "the part of the forthcoming images that we want to see", width and height = the
    # width and height of the image that will be printed on the screen. This instructions takes a viewBox of the
    # forthcoming images and display it in a image of the specified width and height

    #var += ['<svg preserveAspectRatio="xMidYMid meet" x="5" y="5" viewBox="0 0 100 100" width="90" height="90" >\n']
    # var += ['<svg preserveAspectRatio="xMidYMid meet" x="5" y="5" viewBox="%s %s %s %s" width="90" height="90" >\n' % (xViewBox, yViewBox, widthViewBox, heightViewBox)]
    var += ['<svg x="5" y="5" viewBox="%s %s %s %s" width="90" height="90" >\n' % (xViewBox, yViewBox, widthViewBox, heightViewBox)]
    for line in pictureHM:
        # or line.find("<svg")>=0 or line.find("</svg")>=0:
        if line.find("<?xml") >= 0 or line.find("<!DOCTYPE") >= 0:
            continue
        else:
            var += line
    var += ["</svg>\n"]
    var += ["</g>\n", "</svg>\n"]

    file = open(outImageFileName, 'w')
    file.writelines(var)
    file.close()
    if switchOnDirectView:
        os.system("%s %s" % ('firefox', outImageFileName))


def prepareWholeGenomeHomologyMatrices(genome1, genome2, families,
                                       inSbsInPairComp=None,
                                       filterType=FilterType.InBothGenomes,
                                       minChromLength=0,
                                       tandemGapMax=0,
                                       outputFileName=None,
                                       scaleFactorRectangles=4.0,
                                       maxWidth=100, maxHeight=100,
                                       # for illustrating Mazowita2006
                                       fillCompsWithSbs=False,
                                       withIds=False):

    assert isinstance(genome1, myLightGenomes.LightGenome)
    assert isinstance(genome2, myLightGenomes.LightGenome)
    assert isinstance(families, myLightGenomes.Families)
    assert inSbsInPairComp is None or isinstance(inSbsInPairComp, myTools.Dict2d)

    if inSbsInPairComp is not None:
        if isinstance(inSbsInPairComp, myTools.OrderedDict2dOfLists):
            withIds = True
            sbsInPairComp = inSbsInPairComp
        else:
            withIds = False
            assert isinstance(inSbsInPairComp, myTools.Dict2d)
            sbsInPairComp = myTools.OrderedDict2dOfLists()
            for ((c1, c2), sb) in inSbsInPairComp.iteritems2d():
                sbsInPairComp.addToLocation((c1, c2), sb)
    else:
        sbsInPairComp = myTools.OrderedDict2dOfLists()

    cx = len(genome1.keys())
    # whole number of genes on the x-axis genome (genome1)
    wnx = sum([len(chrom) for chrom in genome1.values()])

    cy = len(genome2.keys())
    # whole number of genes on the y-axis genome (genome2)
    wny = sum([len(chrom) for chrom in genome2.values()])

    drawChromosomes = True if max(wny, wnx) <= 300 else False

    # the size of the components of the matrices is chosen using the smallest and more restricting dimension
    # the smallest dimension is the one that contains more genes and chromosome borders compared to its size

    # room for chromosome names
    lchrNames = min(float(maxWidth)/50, float(maxHeight)/50)
    # wnx units for homologies
    # cx+1 for boundaries between chromosomes
    # same for wny and and cy+1
    sizeCase = min(float(maxWidth - lchrNames) / (wnx + (cx+1)), float(maxHeight - lchrNames) / (wny + (cy+1)))
    width  = (wnx + (cx + 1)) * sizeCase + lchrNames
    height = (wny + (cy + 1)) * sizeCase + lchrNames

    minVisibleStrokeWidthOnFirefox = 1.0 / 100.0
    # conversion either with inkscape or rsvg-convert
    minVisibleStrokeForGoodConversionIntoPDF = 1.0 / 10.0
    minVisibleStroke = minVisibleStrokeForGoodConversionIntoPDF

    ((g1_tb, g1_fID, Gtb2GfID1), (g2_tb, g2_fID, Gtb2GfID2)) = editGenomes(genome1, genome2, families,
                                                                           filterType, minChromLength, tandemGapMax)
    # sort chromosome by decreasing lengths
    sortedChrs1 = [c for (c, nbGenes) in sorted(g1_tb.iteritems(), key=lambda x: len(x[1]), reverse=True)]
    sortedChrs2 = [c for (c, nbGenes) in sorted(g2_tb.iteritems(), key=lambda x: len(x[1]), reverse=True)]

    # for each pairwise comparison prepare the homology matrix
    cumulatedX = lchrNames
    cumulatedY = height - lchrNames
    listOfItems = []
    progressBar = myTools.ProgressBar(len(sortedChrs1) * len(sortedChrs2))

    listOfChrNames = []
    lineWidthBetweenChrs = minVisibleStroke
    nbCompWithSb = 0
    for (i1, c1) in enumerate(sortedChrs1):
        # vertical line separating chromosome comparisons on the x-axis
        listOfItems.append(mySvgDrawer.Line(Point(cumulatedX,                  0),
                                            Point(cumulatedX, height - lchrNames), width=lineWidthBetweenChrs))
        nx = len(genome1[c1])
        listOfChrNames.append(mySvgDrawer.Text(Point(cumulatedX + float(nx * sizeCase)/2, height - float(lchrNames)/2),
                                               c1, size=float(lchrNames)/3, text_anchor='middle'))
        range1 = (0, nx)
        for (i2, c2) in enumerate(sortedChrs2):
            ny = len(genome2[c2])
            if i1 == 0:
                # horizontal line separating chromosome comparisons on the y-axis
                listOfItems.append(mySvgDrawer.Line(Point(lchrNames, cumulatedY),
                                                    Point(width,     cumulatedY), width=lineWidthBetweenChrs))
                listOfChrNames.append(mySvgDrawer.Text(Point(float(lchrNames)/2, cumulatedY - float(ny * sizeCase)/2),
                                                       c2, size=float(lchrNames)/3, text_anchor='middle'))
            range2 = (0, ny)

            listOfMatrixItems = []
            # for illustrating Mazowita2006
            if fillCompsWithSbs:
                if len(sbsInPairComp[c1][c2]) > 0:
                    listOfMatrixItems += [mySvgDrawer.Rectangle(Point(0, 0), (ny * sizeCase), (nx * sizeCase), fill=(0,0,0), fill_opacity=1.0)]
                    nbCompWithSb += 1
            else:
                listOfDiagsWithIds = [(sb, id) for (sb, id) in sbsInPairComp.getItemsAndIdsByLocation((c1, c2))]
                listOfMatrixItems += svgItemsHMChrom1Chrom2(genome1[c1], genome2[c2],
                                                            families,
                                                            listOfDiagsWithIds,
                                                            symbolsInChrom1=None, symbolsInChrom2=None,
                                                            range1=range1, range2=range2,
                                                            drawChromosomes=drawChromosomes,
                                                            drawLargeLines=False,
                                                            drawSmallLines=False,
                                                            drawTicks=drawChromosomes,
                                                            drawHomologySymbols=False,
                                                            scaleFactorRectangles=scaleFactorRectangles,
                                                            withSbIds=withIds,
                                                            maxWidth=(nx * sizeCase), maxHeight=(ny * sizeCase))
            listOfMatrixItems = mySvgDrawer.translateItems(listOfMatrixItems, Point(cumulatedX, cumulatedY - (ny * sizeCase)))
            cumulatedY -= ny * sizeCase + sizeCase
            listOfItems += listOfMatrixItems
            progressBar.printProgressIn(sys.stderr, i1 + i2)
        cumulatedY = height - lchrNames
        cumulatedX += nx * sizeCase + sizeCase
    print >> sys.stderr, '## Nb pair. comp. with sbs = %s' % nbCompWithSb
    nbComps = len(sortedChrs1) * len(sortedChrs2)
    print >> sys.stderr, '## Nb chrs1 = %s, Nb chrs2 = %s' % (len(sortedChrs1), len(sortedChrs2))
    print >> sys.stderr, '## Nb pair. comp. = %s' % nbComps
    print >> sys.stderr, '## prop. comp. with sbs = %.2f' % (float(nbCompWithSb) / float(nbComps))
    # last horizontal line
    listOfItems.append(mySvgDrawer.Line(Point(lchrNames, 0),
                                        Point(width,     0), width=float(width)/2000))
    # last vertical line
    listOfItems.append(mySvgDrawer.Line(Point(width, 0),
                                        Point(width, height - lchrNames), width=float(width)/2000))
    # add chromosome names
    #listOfChrNames = mySvgDrawer.translateItems(listOfChrNames, Point(float(lchrNames)/2, float(lchrNames)/2))
    listOfItems += listOfChrNames
    # whole genome homology matrix
    WHM = listOfItems
    return WHM

def svgItemsHMChrom1Chrom2(chrom1, chrom2,
                           families,
                           listOfDiagsWithIds,
                           symbolsInChrom1=None, symbolsInChrom2=None,
                           range1=None, range2=None,
                           diagNames=None,
                           drawChromosomes=None,
                           drawSmallLines=None,
                           drawLargeLines=True,
                           drawTicks=True,
                           drawHomologySymbols=None,
                           scaleFactorRectangles=1.0,
                           withSbIds=False,
                           maxWidth=100,
                           maxHeight=100):

    assert isinstance(chrom1, myLightGenomes.Chromosome)
    assert isinstance(chrom2, myLightGenomes.Chromosome)
    assert isinstance(listOfDiagsWithIds, list)
    if len(listOfDiagsWithIds) > 0:
        assert isinstance(listOfDiagsWithIds[0][0], myDiags.Diagonal)
    if not range1:
        range1 = (0, len(chrom1))
    if not range2:
        range2 = (0, len(chrom2))
    if not diagNames:
        if withSbIds:
            diagNames = [id for (sb, id) in listOfDiagsWithIds]
        else:
            diagNames = None

    coordsHomologiesPerDiags = []
    symbolsInHomologiesInDiags = collections.defaultdict(lambda: collections.defaultdict(int))
    for (sb, id) in listOfDiagsWithIds:
        coordsHomologiesPerDiags.append([])
        for idxHp, aG in enumerate(sb.la):
            assert isinstance(aG, tuple)
            (aGn, aGs, dist) = aG
            for (gi1, gi2) in itertools.product(sb.l1[idxHp], sb.l2[idxHp]):
                assert range1[0] <= gi1 <= range1[1], '%s <= %s <= %s' % (range1[0], gi1, range1[1])
                assert range2[0] <= gi2 <= range2[1], '%s <= %s <= %s' % (range2[0], gi2, range2[1])
                idxRoiG1 = gi1 - range1[0]
                idxRoiG2 = gi2 - range2[0]
                coordsHomologiesPerDiags[-1].append((idxRoiG1, idxRoiG2))
                assert aGs in {+1, -1, None}, "%s" % aGs
                symbolsInHomologiesInDiags[idxRoiG1][idxRoiG2] = aGs

    (MHP, locG2) = myDiags.homologyMatrix(chrom1, chrom2)

    # add two very ad hoc families for visualisation purposes
    whiteGeneNames = chrom1.getGeneNamesNotInFamilies(families) | chrom2.getGeneNamesNotInFamilies(families)
    # a ^ b : symmetric difference: either in a or in b but not in both
    geneNamesInEither1or2ButNotBoth = set()
    familiesConservedInOnlyOneGenome = chrom1.getOwnedFamilyNames(families) ^ chrom2.getOwnedFamilyNames(families)
    for g in chrom1 + chrom2:
        famName = families.getFamNameByName(g.n, default=None)
        if famName in familiesConservedInOnlyOneGenome:
            geneNamesInEither1or2ButNotBoth.add(g.n)
    assert len(geneNamesInEither1or2ButNotBoth) < len(chrom1) + len(chrom2)

    # cluster the genes of the same family
    tmp = collections.defaultdict(set)
    for gn in geneNamesInEither1or2ButNotBoth:
        famName = families.getFamNameByName(gn, default=None)
        if famName is not None:
            tmp[famName].add(gn)
    geneNamesWithSameGreys = tmp.values()

    nx = len(chrom1)
    if range1:
        assert nx == range1[1] - range1[0], "%s = %s - %s" % (nx, range1[0], range1[1])
    ny = len(chrom2)
    if range2:
        assert ny == range2[1] - range2[0], "%s = %s - %s" % (ny, range2[0], range2[1])
    if drawChromosomes is None:
        drawChromosomes = max(nx, ny) < 300
    if drawSmallLines is None:
        drawSmallLines = max(nx, ny) < 300
    if drawHomologySymbols is None:
        drawHomologySymbols = max(nx, ny) < 150

    # the size of the components of the matrix is chosen using the smallest and more restricting dimension
    # (contains more genes comparing to its size)
    tmp_marginChrom = 1 if drawChromosomes else 0
    tmp_marginTicks = 3 if drawTicks else 0
    ncx = nx + tmp_marginChrom + tmp_marginTicks
    ncy = ny + tmp_marginChrom + tmp_marginTicks
    sizeCase = float(min(float(maxWidth) / ncx, float(maxHeight) / ncy))
    marginChrom = tmp_marginChrom * sizeCase
    marginTicks = tmp_marginTicks * sizeCase
    finalWidth = sizeCase * ncx
    finalHeight = sizeCase * ncy
    assert maxHeight - finalHeight > -0.0001
    assert maxWidth - finalWidth > -0.0001
    (listOfMatrixItems, (width, height)) = svgItemsMatrixFromRawData(range1, range2,
                                                                     MHP,
                                                                     sizeCase,
                                                                     coordsOfHomologiesPerDiags=coordsHomologiesPerDiags,
                                                                     diagNames=diagNames,
                                                                     diagColorGenerator=LevelIdxGenerator(farIdxs=7),
                                                                     drawSmallLines=drawSmallLines,
                                                                     drawLargeLines=drawLargeLines,
                                                                     drawHomologySymbols=drawHomologySymbols,
                                                                     scaleFactorRectangles=scaleFactorRectangles)
    # width = nx * sizeCase
    # height = ny * sizeCase
    assert abs(width - nx * sizeCase) < 0.001 and abs(height - ny * sizeCase) < 0.001
    if drawChromosomes:
        formerWidth = width
        formerHeight = height
        (listOfMatrixItems, (width, height)) = addChromosomesAlongMH(listOfMatrixItems, chrom1, chrom2, families,
                                                                     sizeCase,
                                                                     formerWidth, formerHeight,
                                                                     geneNamesWithSameGreys=geneNamesWithSameGreys,
                                                                     whiteGeneNames=whiteGeneNames,
                                                                     symbolsInChrom1=symbolsInChrom1, symbolsInChrom2=symbolsInChrom2)
        assert abs(width - (formerWidth + marginChrom)) < 0.001 and abs(height - (formerHeight + marginChrom)) < 0.001
    if drawTicks:
        formerWidth = width
        formerHeight = height
        (listOfMatrixItems, (width, height)) = addTicksAlongMH(listOfMatrixItems, chrom1, chrom2, sizeCase,
                                                                     range1, range2,
                                                                     marginChrom,
                                                                     marginTicks,
                                                                     formerWidth, formerHeight)
        assert abs(width - ((nx * sizeCase) + marginChrom + marginTicks)) < 0.001 and abs(height - ((ny * sizeCase) + marginChrom + marginTicks)) < 0.001
    assert abs(finalHeight - ((ny * sizeCase) + marginChrom + marginTicks)) < 0.001
    assert abs(finalWidth - ((nx * sizeCase) + marginChrom + marginTicks)) < 0.001
    assert abs(height - finalHeight) < 0.001
    assert abs(width - finalWidth) < 0.001
    if not drawChromosomes and not drawTicks:
        assert finalWidth == maxWidth or finalHeight == maxHeight
    return listOfMatrixItems

def addChromosomesAlongMH(listOfMatrixItems, chrom1, chrom2, families,
                          sizeCase,
                          formerWidth, formerHeight,
                          geneNamesWithSameGreys=None,
                          whiteGeneNames=None,
                          symbolsInChrom1=None, symbolsInChrom2=None):
    marginChroms = 1 * sizeCase
    listOfMatrixItems = mySvgDrawer.translateItems(listOfMatrixItems, Point(marginChroms, 0))

    geneNamesNotInColor = set()
    if geneNamesWithSameGreys:
        geneNamesNotInColor |= set([gn for geneNamesWithSameGrey in geneNamesWithSameGreys for gn in geneNamesWithSameGrey])
    if whiteGeneNames:
        geneNamesNotInColor |= whiteGeneNames

    familyName2color = {}
    colorGenerator = LevelIdxGenerator()
    coloredGeneNames1 = set()
    for (gn, _) in chrom1:
        if gn not in geneNamesNotInColor:
            famName = families.getFamNameByName(gn)
            if famName:
                coloredGeneNames1.add(gn)
                familyName2color[famName] = colorGenerator.getLevel()
    coloredGeneNames2 = set()
    for (gn, _) in chrom2:
        if gn not in geneNamesNotInColor:
            famName = families.getFamNameByName(gn)
            if famName:
                coloredGeneNames2.add(gn)
                if famName not in familyName2color:
                    familyName2color[famName] = colorGenerator.getLevel()
    geneWidth = sizeCase * 0.75

    def assertWhiteGreyAndColourGeneNamesAreComplementary(whiteGeneNames, geneNamesWithSameGreys, coloredGeneNames, chrom):
        greyGeneNames = set([gn for gns in geneNamesWithSameGreys for gn in gns])
        geneNamesChrom = set([g.n for g in chrom])
        whiteGeneNames = whiteGeneNames & geneNamesChrom
        greyGeneNames = greyGeneNames & geneNamesChrom
        assert len(coloredGeneNames &  whiteGeneNames) == 0
        assert len(coloredGeneNames &  greyGeneNames) == 0
        assert whiteGeneNames | greyGeneNames | coloredGeneNames == geneNamesChrom
    assertWhiteGreyAndColourGeneNamesAreComplementary(whiteGeneNames, geneNamesWithSameGreys, coloredGeneNames1, chrom1)
    assertWhiteGreyAndColourGeneNamesAreComplementary(whiteGeneNames, geneNamesWithSameGreys, coloredGeneNames2, chrom2)

    chromosome1Items = svgItemsChromosome(chrom1,
                                          families=families,
                                          name2color=familyName2color,
                                          lengthGenes=sizeCase * 0.90,
                                          halfIntergeneLengths=sizeCase * 0.05,
                                          geneNamesWithSameGreyLevels=geneNamesWithSameGreys,
                                          whiteGeneNames=whiteGeneNames,
                                          greyLevelsGenerator=LevelIdxGenerator(inGrays=True),
                                          symbolsInGenes=symbolsInChrom1,
                                          geneWidth=geneWidth)
    chromosome2Items = svgItemsChromosome(chrom2,
                                          families=families,
                                          name2color=familyName2color,
                                          lengthGenes=sizeCase * 0.90,
                                          halfIntergeneLengths=sizeCase * 0.05,
                                          geneNamesWithSameGreyLevels=geneNamesWithSameGreys,
                                          whiteGeneNames=whiteGeneNames,
                                          greyLevelsGenerator=LevelIdxGenerator(inGrays=True),
                                          symbolsInGenes=symbolsInChrom2,
                                          geneWidth=geneWidth)
    chromosome1Items = mySvgDrawer.translateItems(chromosome1Items, Point(marginChroms, formerHeight + marginChroms - sizeCase/2.0))
    chromosome2Items = mySvgDrawer.rotateItems(chromosome2Items, Point(0.0, 0.0), -90)
    chromosome2Items = mySvgDrawer.translateItems(chromosome2Items, Point(sizeCase/2.0, formerHeight))
    listOfMatrixItems += chromosome1Items + chromosome2Items
    newWidth = formerWidth + marginChroms
    newHeight = formerHeight + marginChroms
    return (listOfMatrixItems, (newWidth, newHeight))


def addTicksAlongMH(listOfMatrixItems, chrom1, chrom2, sizeCase,
                    range1, range2,
                    marginChrom,
                    marginTicks,
                    formerWidth, formerHeight):
    newHeight = formerHeight + marginTicks
    newWidth = formerWidth + marginTicks
    invYaxis = lambda y: newHeight - (sizeCase + y)
    widthTicks = max(formerWidth, formerHeight) / 500.0
    sizeText = widthTicks * 7.5
    listOfMatrixItems = mySvgDrawer.translateItems(listOfMatrixItems, Point(marginTicks, 0.0))
    x_itemTicks = svgItemsTicks(range1, sizeCase, widthTicks=widthTicks, sizeText=sizeText, rankGenome=1)
    x_itemTicks = mySvgDrawer.translateItems(x_itemTicks, Point(marginTicks + marginChrom, invYaxis(marginTicks - sizeCase)))
    y_itemTicks = svgItemsTicks(range2, sizeCase, widthTicks=widthTicks, sizeText=sizeText, rankGenome=2)
    y_itemTicks = mySvgDrawer.rotateItems(y_itemTicks, Point(0.0,0.0), -90)
    y_itemTicks = mySvgDrawer.translateItems(y_itemTicks, Point(sizeCase, invYaxis(marginTicks + marginChrom - sizeCase)))
    listOfMatrixItems += x_itemTicks + y_itemTicks
    return (listOfMatrixItems, (newWidth, newHeight))

def editSbsInROIWithGenomesInGenes(inSbsInPairComp,
                                   genome1WithOnlyOneChr1, genome2WithOnlyOneChr2,
                                   chr1, chr2,
                                   range1, range2,
                                   considerAllPairComps):
    assert inSbsInPairComp is not None
    assert isinstance(genome1WithOnlyOneChr1, myLightGenomes.LightGenome)
    assert isinstance(genome2WithOnlyOneChr2, myLightGenomes.LightGenome)

    if isinstance(inSbsInPairComp, myTools.OrderedDict2dOfLists):
        withSbIds = True
        sbsInPairComp = inSbsInPairComp
    else:
        assert isinstance(inSbsInPairComp, myTools.Dict2d)
        withSbIds = False
        sbsInPairComp = myTools.OrderedDict2dOfLists()
        for ((c1, c2), sb) in  inSbsInPairComp.items2d():
            sbsInPairComp.addToLocation((c1, c2), sb)

    assert isinstance(sbsInPairComp, myTools.OrderedDict2dOfLists)
    if not considerAllPairComps:
        # need to change the indexes of diags!
        for (_, sb) in sbsInPairComp.iteritems2d():
            sb.l1 = [[range1[0] + i1 for i1 in tb1] for tb1 in sb.l1]
            sb.l2 = [[range2[0] + i2 for i2 in tb2] for tb2 in sb.l2]

    # truncate sbs to fit the ROI
    listOfDiags = []
    for (sb, id) in sbsInPairComp.getItemsAndIdsByLocation((chr1, chr2)):
        newl1 = []
        newl2 = []
        newla = []
        for (idxHp, aG) in enumerate(sb.la):
            tb1 = []
            tb2 = []
            for i1g in sb.l1[idxHp]:
                if (range1[0] <= i1g < range1[1]):
                    tb1.append(i1g)
            for i2g in sb.l2[idxHp]:
                if (range2[0] <= i2g < range2[1]):
                    tb2.append(i2g)
            if len(tb1) > 0 and len(tb2) > 0:
                newl1.append(tb1)
                newl2.append(tb2)
                newla.append(aG)
        if len(newla) > 0:
            newSb = myDiags.SyntenyBlock(myDiags.Diagonal(sb.dt, newl1, newl2, newla), sb.pVal)
            listOfDiags.append((newSb, id))
    return (listOfDiags, withSbIds)

# load precomputed sbs if any, otherwise compute sbs
def editSbsInROIWithGenomesInTbs(sbsInPairComp,
                                 g1_tb, g2_tb,
                                 chr1, chr2,
                                 range1, range2,
                                 considerAllPairComps,
                                 g2tb1, g2tb2):
    assert isinstance(g1_tb, myLightGenomes.LightGenome)
    assert isinstance(g2_tb, myLightGenomes.LightGenome)

    assert len(sbsInPairComp.intoList()) > 0
    if isinstance(sbsInPairComp, myTools.OrderedDict2dOfLists):
        withSbIds = True
        newSbsInPairComp = sbsInPairComp
    else:
        assert isinstance(sbsInPairComp, myTools.Dict2d)
        withSbIds = False
        newSbsInPairComp = myTools.OrderedDict2dOfLists()
        for ((c1, c2), sb) in  sbsInPairComp.items2d():
            newSbsInPairComp.addToLocation((c1, c2), sb)

    listOfDiags = []
    for (sb, id) in newSbsInPairComp.getItemsAndIdsByLocation((chr1, chr2)):
        assert len(sb.l1) > 0
        if isinstance(sb.l1[0], list):
            # change the sb.lX structure from list of lists to list of ints
            sb.l1 = [g2tb1[tb[0]] for tb in sb.l1]
            sb.l2 = [g2tb2[tb[0]] for tb in sb.l2]
            if len(sb.l1) > 0:
                assert isinstance(sb.l1[0], int)

        # need to change the indexes of diags!
        if not considerAllPairComps:
            xoffset = range1[0]
            yoffset = range2[0]
        else:
            xoffset = 0
            yoffset = 0
        if len(sb.l1) > 0:
                assert isinstance(sb.l1[0], int)
        sb.l1 = [[xoffset + i1] for i1 in sb.l1]
        sb.l2 = [[yoffset + i2] for i2 in sb.l2]

       # truncate sbs to fit the ROI
        if (range1[0] <= sb.minOnG(1) and sb.maxOnG(1) <= range1[1]) \
                and (range2[0] <= sb.minOnG(2) and sb.maxOnG(2) <= range2[1]):
            # sb is perfectly included in the ROI
            listOfDiags.append((sb, id))
        elif (sb.maxOnG(1) < range1[0] or range1[1] < sb.minOnG(1)) \
                or (sb.maxOnG(2) < range2[0] or range2[1] < sb.minOnG(2)):
            # sb is not in the ROI
            continue
        else:
            # sb is partially included in the ROI
            sb.truncate(range1, range2)
            if len(sb.la) > 0:
                listOfDiags.append((sb, id))
    return (listOfDiags, withSbIds)

def reduceFamiliesToFamiliesInFilteredGs(genome1, genome2, families, filterType=myDiags.FilterType.InBothGenomes):
    fnsInG1 = genome1.getOwnedFamilyNames(families)
    fnsInG2 = genome2.getOwnedFamilyNames(families)
    familiesInFiltered2DGs = myLightGenomes.Families()
    if filterType == myDiags.FilterType.InBothGenomes:
        for fnInGs in  fnsInG1 & fnsInG2:
            fam = families.getFamilyByName(fnInGs)
            assert fam is not None
            familiesInFiltered2DGs.addFamily(fam)
    elif filterType == myDiags.FilterType.InFamilies:
        for fnInGs in  fnsInG1 | fnsInG2:
            fam = families.getFamilyByName(fnInGs)
            assert fam is not None
            familiesInFiltered2DGs.addFamily(fam)
    return familiesInFiltered2DGs

def reduceFamiliesToFamiliesInFilteredChroms(chrom1, chrom2, families, filterType=myDiags.FilterType.InBothGenomes):
    # FIXME: should use a Chrom class
    genome1 = myLightGenomes.LightGenome()
    genome1['0'] = chrom1
    genome2 = myLightGenomes.LightGenome()
    genome2['0'] = chrom2
    genome1.computeDictG2Ps()
    genome2.computeDictG2Ps()

    fnsInG1 = genome1.getOwnedFamilyNames(families)
    fnsInG2 = genome2.getOwnedFamilyNames(families)
    familiesInFiltered2DGs = myLightGenomes.Families()
    if filterType == myDiags.FilterType.InBothGenomes:
        for fnInGs in  fnsInG1 & fnsInG2:
            fam = families.getFamilyByName(fnInGs)
            assert fam is not None
            familiesInFiltered2DGs.addFamily(fam)
    elif filterType == myDiags.FilterType.InFamilies:
        for fnInGs in  fnsInG1 | fnsInG2:
            fam = families.getFamilyByName(fnInGs)
            assert fam is not None
            familiesInFiltered2DGs.addFamily(fam)
    return familiesInFiltered2DGs

def homologyMatrixViewer(genome1, genome2, families, CDF1, CDF2,
                         convertGenicToTbCoordinates=False,
                         filterType=FilterType.InBothGenomes,
                         minChromLength = 2,
                         tandemGapMax=0,
                         distanceMetric='CD',
                         gapMax=None,
                         distinguishMonoGenicDiags=True,
                         pThreshold=None,
                         gapMaxMicroInv=0,
                         identifyMonoGenicInvs=False,
                         identifyMicroRearrangements=True,
                         truncationMax=None,
                         sameStrand=True,
                         validateImpossToCalc_mThreshold=3,
                         nbHpsRecommendedGap=2,
                         targetProbaRecommendedGap=0.01,
                         chromosomesRewrittenInTbs=False,
                         drawAllInformations=False,
                         doDrawChromosomes=False,
                         scaleFactorRectangles=2.0,
                         considerAllPairComps=True,
                         switchOnDirectView=False,
                         optimisation=None,
                         inSbsInPairComp=None,
                         outSyntenyBlocksFileName="./syntenyBlocksDrawer.txt",
                         outImageFileName="./homologyMatrix.svg",
                         # the number of gene name caracters to write in gene symbols
                         nbCaractersForGeneNamesInSymlbols=1,
                         showGeneNames=False,
                         verbose=True):

    assert isinstance(genome1, myLightGenomes.LightGenome)
    assert isinstance(genome2, myLightGenomes.LightGenome)
    assert isinstance(families, myLightGenomes.Families)

    kwargs = {'gapMax': gapMax,
              'distinguishMonoGenicDiags': distinguishMonoGenicDiags,
              'gapMaxMicroInv': gapMaxMicroInv,
              'distanceMetric': distanceMetric,
              'identifyMonoGenicInvs': identifyMonoGenicInvs,
              'identifyMicroRearrangements': identifyMicroRearrangements,
              'truncationMax': truncationMax,
              'sameStrand': sameStrand,
              'pThreshold': pThreshold,
              'nbHpsRecommendedGap': nbHpsRecommendedGap,
              'targetProbaRecommendedGap': targetProbaRecommendedGap,
              'validateImpossToCalc_mThreshold': validateImpossToCalc_mThreshold,
              'optimisation': optimisation,
              'verbose': True}

    assert distanceMetric == 'DPD' or distanceMetric == 'MD' or distanceMetric == 'CD' or distanceMetric == 'ED'
    assert (convertGenicToTbCoordinates and chromosomesRewrittenInTbs) or not convertGenicToTbCoordinates

    maxWidthHm = 100
    maxHeightHm = 100

    if not chromosomesRewrittenInTbs:
        # define the ROI (Region Of Interest)
        # rangeX = (idx of the first element (start at 0), index of the last element + 1)
        (chr1, range1) = parseChrRange(CDF1, genome1)
        (chr2, range2) = parseChrRange(CDF2, genome2)
        genome1OnlyChr1 = myLightGenomes.LightGenome()
        genome2OnlyChr2 = myLightGenomes.LightGenome()
        for chrom in genome1.values():
            assert isinstance(chrom, myLightGenomes.Chromosome), type(chrom)
        genome1OnlyChr1[chr1] = genome1[chr1][range1[0]:range1[1]]
        genome2OnlyChr2[chr2] = genome2[chr2][range2[0]:range2[1]]
        for chrom in genome1OnlyChr1.values():
            assert isinstance(chrom, myLightGenomes.Chromosome), type(chrom)
        symbolsInChrom1 = [g.n[0:nbCaractersForGeneNamesInSymlbols] for g in genome1OnlyChr1[chr1]]
        symbolsInChrom2 = [g.n[0:nbCaractersForGeneNamesInSymlbols] for g in genome2OnlyChr2[chr2]]
        nbSpeciesSpecificGenes1 = len([gn for gn in genome1OnlyChr1.getGeneNames(asA=list, checkNoDuplicates=False)
                                       if families.getFamilyByName(gn, default=None) is None])
        print >> sys.stderr, "the ROI1 contains %s genes (%s species specific genes)" % (len(genome1OnlyChr1[chr1]), nbSpeciesSpecificGenes1)
        nbSpeciesSpecificGenes2 = len([gn for gn in genome2OnlyChr2.getGeneNames(asA=list, checkNoDuplicates=False) if families.getFamilyByName(gn, default=None) is None])
        print >> sys.stderr, "the ROI2 contains %s genes (%s species specific genes)" % (len(genome2OnlyChr2[chr2]), nbSpeciesSpecificGenes2)

        if considerAllPairComps:
            comparedGenome1 = genome1
            comparedGenome2 = genome2
            familiesForComp = families
        else:
            # extract diagonals in the ROI without considering other pairwise comparisons
            comparedGenome1 = genome1OnlyChr1
            comparedGenome2 = genome2OnlyChr2
            familiesForComp = reduceFamiliesToFamiliesInFilteredChroms(genome1[chr1], genome2[chr2], families, filterType=filterType)

        if inSbsInPairComp is None:
            inSbsInPairComp = myDiags.extractSbsInPairCompGenomes(comparedGenome1,
                                                                  comparedGenome2,
                                                                  familiesForComp,
                                                                  tandemGapMax=tandemGapMax,
                                                                  minChromLength=minChromLength,
                                                                  filterType=filterType,
                                                                  **kwargs)

        (listOfDiagsWithIds, withSbIds)= editSbsInROIWithGenomesInGenes(inSbsInPairComp,
                                                           comparedGenome1, comparedGenome2,
                                                           chr1, chr2,
                                                           range1, range2,
                                                           considerAllPairComps)
        # FIXME : not sure for the next two lines but seems to work to show homologies
        genome1OnlyChr1 = myMapping.labelWithFamNames(genome1OnlyChr1, familiesForComp)
        genome2OnlyChr2 = myMapping.labelWithFamNames(genome2OnlyChr2, familiesForComp)
        svgItemsHM = svgItemsHMChrom1Chrom2(genome1OnlyChr1[chr1], genome2OnlyChr2[chr2],
                                            familiesForComp,
                                            listOfDiagsWithIds,
                                            symbolsInChrom1=symbolsInChrom1, symbolsInChrom2=symbolsInChrom2,
                                            range1=range1, range2=range2,
                                            diagNames=None,
                                            drawChromosomes=None,
                                            drawSmallLines=None,
                                            drawLargeLines=True,
                                            drawTicks=True,
                                            scaleFactorRectangles=scaleFactorRectangles,
                                            maxWidth=maxWidthHm,
                                            maxHeight=maxHeightHm,
                                            withSbIds=withSbIds)
    else:
        assert chromosomesRewrittenInTbs
        ((g1_tb, mtb2g1, (nCL1, nGL1)), (g2_tb, mtb2g2, (nCL2, nGL2))) = \
            myDiags.editGenomes(genome1, genome2, families,
                                filterType=filterType, labelWith='FamName', tandemGapMax=tandemGapMax,
                                minChromLength=minChromLength, keepOriginal=True)

        if not convertGenicToTbCoordinates:
            (chr1, range1) = parseChrRange(CDF1, g1_tb)
            (chr2, range2) = parseChrRange(CDF2, g2_tb)
        else:
            mg2tb1 = dict((c, m.old) for (c, m) in mtb2g1.iteritems())
            mg2tb2 = dict((c, m.old) for (c, m) in mtb2g2.iteritems())
            (chr1, range1) = parseChrRange(CDF1, genome1, g2gtb=mg2tb1)
            (chr2, range2) = parseChrRange(CDF2, genome2, g2gtb=mg2tb2)

        # TChr: Truncated Chromosome
        genome1tbOnlyTChr1 = myLightGenomes.LightGenome()
        genome2tbOnlyTChr2 = myLightGenomes.LightGenome()
        genome1tbOnlyTChr1[chr1] = g1_tb[chr1][range1[0]:range1[1]]
        genome2tbOnlyTChr2[chr2] = g2_tb[chr2][range2[0]:range2[1]]
        symbolsInChromTb1 =  [str(len(tb)) for tb in mtb2g1[chr1].new]
        symbolsInChromTb2 =  [str(len(tb)) for tb in mtb2g2[chr2].new]

        print >> sys.stderr, "the ROI1 contains %s genes" % (sum(len(mtb2g1[chr1][itb]) for (itb, _) in enumerate(genome1tbOnlyTChr1[chr1])))
        print >> sys.stderr, "the ROI2 contains %s genes" % (sum(len(mtb2g2[chr2][itb]) for (itb, _) in enumerate(genome2tbOnlyTChr2[chr2])))
        #Focus on the chromosome of the window, just give simple name to the chromosome of interest
        tb2g1 = mtb2g1[chr1]
        g2tb1 = mtb2g1[chr1].old
        tb2g2 = mtb2g2[chr2]
        g2tb2 = mtb2g2[chr2].old

        if considerAllPairComps:
            comparedGenome1 = g1_tb
            comparedGenome2 = g2_tb
            familiesForComp = families
        else:
            # extract diagonals in the ROI without considering other pairwise comparisons
            comparedGenome1 = genome1tbOnlyTChr1
            comparedGenome2 = genome2tbOnlyTChr2
            familiesForComp = reduceFamiliesToFamiliesInFilteredChroms(genome1[chr1], genome2[chr2], families,
                                                                       filterType=filterType)

        if inSbsInPairComp is None:
            inSbsInPairComp = myDiags.extractSbsInPairCompGenomesInTbs(comparedGenome1,
                                                                       comparedGenome2,
                                                                       **kwargs)
        (listOfDiagsWithIds, withSbIds) = editSbsInROIWithGenomesInTbs(inSbsInPairComp,
                                                           comparedGenome1, comparedGenome2,
                                                           chr1, chr2,
                                                           range1, range2,
                                                           considerAllPairComps,
                                                        g2tb1, g2tb2)

        svgItemsHM = svgItemsHMChrom1Chrom2(genome1tbOnlyTChr1[chr1], genome2tbOnlyTChr2[chr2],
                                            familiesForComp,
                                            listOfDiagsWithIds,
                                            symbolsInChrom1=symbolsInChromTb1, symbolsInChrom2=symbolsInChromTb2,
                                            range1=range1, range2=range2,
                                            diagNames=None,
                                            drawChromosomes=None,
                                            drawSmallLines=None,
                                            drawLargeLines=True,
                                            drawTicks=True,
                                            scaleFactorRectangles=scaleFactorRectangles,
                                            maxWidth=maxWidthHm,
                                            maxHeight=maxHeightHm,
                                            withSbIds=withSbIds)

    # write a simple file with all diagonals into output file
    with open(outSyntenyBlocksFileName, 'w') as f:
        print >> f, "Mode : %s" % 'Genic scale' if chromosomesRewrittenInTbs is False else 'Tandem Blocks scale'
        print >> f, "chromosome %s de %s\t%s\t%s\tchromosome %s de %s\t%s\t%s\t%s" % (chr1, genome1.name, 'beginC1', 'endC1', chr2, genome2.name, 'beginC2', 'endC2', 'length in families')
        print >> f, "c1\tbeg1\tend1\tc2\tbeg2\tend2\thps\tpVal"

        for (sb, id) in listOfDiagsWithIds:
            if isinstance(sb.l1[0], list) and isinstance(sb.l2[0], list):
                minl1 = min(sb.l1[0])
                maxl1 = max(sb.l1[-1])
                minl2 = min(sb.l2[0])
                maxl2 = max(sb.l2[-1])
            elif isinstance(sb.l2[0], int) and isinstance(sb.l2[0], int):
                minl1 = min(sb.l1)
                maxl1 = max(sb.l1)
                minl2 = min(sb.l2)
                maxl2 = max(sb.l2)
            # indices of genes start at 1 to be coherent with the output image
            print >> f, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chr1, minl1 + 1, maxl1 + 1, chr2, minl2 + 1, maxl2 + 1, len(sb.la), sb.pVal)

    nbSbs = len(listOfDiagsWithIds)
    writeSVGFileForPairwiseCompOfChrs(genome1.name, chr1, range1,
                                      genome2.name, chr2, range2,
                                      svgItemsHM, maxWidthHm, maxHeightHm,
                                      chromosomesRewrittenInTbs,
                                      filterType,
                                      tandemGapMax,
                                      gapMax,
                                      distanceMetric,
                                      gapMaxMicroInv,
                                      identifyMicroRearrangements,
                                      truncationMax,
                                      nbSbs,
                                      outImageFileName,
                                      switchOnDirectView)

if __name__ == '__main__':
    # genome1 = myLightGenomes.LightGenome('/home/jlucas/Libs/PhylDiag/data/genesST.Homo.sapiens.light2.list.bz2')
    # genome2 = myLightGenomes.LightGenome('/home/jlucas/Libs/PhylDiag/data/genesST.Mus.musculus.light.list.bz2')
    families = myLightGenomes.Families('/home/jlucas/Libs/PhylDiag/data/ancGenes.Amniota.list.bz2')
    #genome1 = myLightGenomes.LightGenome('/home/jlucas/Libs/PhylDiag/data/genesST.Homo.sapiens.list.bz2')
    genome1 = myLightGenomes.LightGenome('/home/jlucas/Libs/PhylDiag/data/genesST.Homo.sapiens.list.bz2')
    genome2 = myLightGenomes.LightGenome('/home/jlucas/Libs/PhylDiag/data/genesST.Gallus.gallus.list.bz2')
    #genome2 = myLightGenomes.LightGenome('/home/jlucas/Libs/PhylDiag/data/genesST.Mus.musculus.list.bz2')
    genome1.removeUnofficialChromosomes()
    genome2.removeUnofficialChromosomes()
    #families = myLightGenomes.Families('/home/jlucas/Libs/PhylDiag/data/ancGenes.Euarchontoglires.list.bz2')

    prepareWholeGenomeHomologyMatrices(genome1, genome2, families, inSbsInPairComp=None, maxWidth=100, maxHeight=100,
                                filterType=FilterType.None,
                                minChromLength=0,
                                tandemGapMax=0,
                                outputFileName='toto.svg',
                                scaleFactorRectangles=10)

#Find diags with the more paralogs
@myTools.deprecated
def searchInterestingDiags(listOfDiags, range1, range2):
    ecart=-sys.maxint-1
    diag_=None
    for diag in listOfDiags:
        l1=diag[0][1]
        l2=diag[1][1]
        la=diag[2]
        if len(l1) != len(l2):
            print >> sys.stderr, "assymetric (on genome1 and genome2) SB of ", len([anc for (anc, _, _, _) in la]), " HPs"
            print >> sys.stderr, "diag on G1 from %s to %s" % (min(l1[0][0]+range1[0], l1[-1][0]+range1[0]),
                                                               max(l1[0][0]+range1[0], l1[-1][0]+range1[0]))
            print >> sys.stderr, "diag on G2 from %s to %s" % (min(l2[0][0]+range2[0], l2[-1][0]+range2[0]),
                                                               max(l2[0][0]+range2[0], l2[-1][0]+range2[0]))
            if ecart < abs(len(l1) - len(la)):
                ecart = abs(len(l1) - len(la))
                diag_=diag
            if ecart < abs(len(l2) - len(la)):
                ecart = abs(len(l2) - len(la))
                diag_=diag
    diag = diag_ if diag_ != None else listOfDiags[0]
    l1 = diag[0][1]
    l2 = diag[1][1]
    la = diag[2]
    print >> sys.stderr, "The most assymetric diag"
    print >> sys.stderr, " SB of ", len([anc for (anc, _, _, _) in la]), " HPs"
    print >> sys.stderr, "diag on G1 = ", [gene[0] + range1[0] for gene in l1]
    print >> sys.stderr, "diag on G2 = ", [gene[0] + range2[0] for gene in l2]

    minIndiceG1 = sys.maxint
    maxIndiceG1 = -sys.maxint-1
    minIndiceG2 = sys.maxint
    maxIndiceG2 = -sys.maxint-1
    for diag in listOfDiags :
        ((_, l1), (_, l2), _) = diag
        for (i1, _) in l1:
            minIndiceG1 = i1 if i1 < minIndiceG1 else minIndiceG1
            maxIndiceG1 = i1 if i1 > maxIndiceG1 else maxIndiceG1
        for (i2, _) in l2:
            minIndiceG2 = i2 if i2 < minIndiceG2 else minIndiceG2
            maxIndiceG2 = i2 if i2 > maxIndiceG2 else maxIndiceG2
    print >> sys.stderr, 'Indices entre lesquels se trouvents les diagonales sur le G1 = [%s,%s]' % (minIndiceG1+range1[0],
                                                                                                     maxIndiceG1+range1[0])
    print >> sys.stderr, 'Indices entre lesquels se trouvents les diagonales sur le G2 = [%s,%s]' % (minIndiceG2+range2[0],
                                                                                                     maxIndiceG2+range2[0])

    return

# ask the user for the desired chromosome ranges
@myTools.deprecated
def chooseChrsAndRanges(genome1, genome2, families, distanceMetric = 'DPD'):
    while True:
        try:
            (chr1, range1) = parseChrRange(raw_input("chr1:deb1-fin1 = "), genome1)
            break
        except ValueError:
            print >> sys.stderr, "You need to write something as chr1:deb1-fin1 with chr1 a chr of G1 and deb1 and fin1 indices of the first and last gene (indices start at 1)"

    while True:
        try:
            (chr2, range2) = parseChrRange(raw_input("chr2:deb2-fin2 = "), genome2)
            break
        except ValueError:
            print >> sys.stderr, "You need to write something as chr2:deb2-fin2 with chr2 a chr of G2 and deb2 and fin2 indices of the first and last gene (indices start at 1)"

    chrom1 ={}
    chrom2 ={}
    chrom1[chr1] = genome1[chr1][range1[0]:range1[1]]
    chrom2[chr2] = genome2[chr2][range2[0]:range2[1]]
    listOfDiags = myDiags.extractSbsInPairCompGenomes(chrom1, chrom2, families,
                                                      gapMax=gapMax, sameStrand=sameStrand,
                                                      filterType=filterType, minChromLength=minChromLength,
                                                      distanceMetric=distanceMetric, verbose=verbose)
    listOfDiags = list(listOfDiags)
    print >> sys.stderr, "pairwise comparison of the two chromosomes yields", len(listOfDiags), "diagonals."
    if len(listOfDiags) == 0:
        print >> sys.stderr, "There is no diag in the considered region of interest"
    else:
        #recherche des diagonales interessantes a afficher
        print >> sys.stderr, range1
        print >> sys.stderr, range2
        searchInterestingDiags(listOfDiags, range1, range2)

    print >> sys.stderr, "Do you want to chose a new region of interest ?"
    if raw_input('y/n ? ') == 'y':
        return chooseChrsAndRanges(genome1, genome2, families)
    else:
        return (chrom1, chrom2, range1, range2)