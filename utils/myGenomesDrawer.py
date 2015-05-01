# -*- coding: utf-8 -*-
# LibsDyogen
# python 2.7
# Copyright © 2015 IBENS/Dyogen Joseph LUCAS, and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

import sys
import os
import collections
import random
import itertools
import myTools
import myDiags
import myMapping
import mySvgDrawer
import myLightGenomes
from mySvgDrawer import Point

# parse the user input (text) for the chromosome range and asses if this query is consistent with the genome data
# g2gtb is a dictionnary to convert indices from gene coordinates into tb coordinates
def parseChrRange(text, genome, g2gtb=None):
    if len(text.split(":")) == 2:
        chr = text.split(":")[0]
        if chr not in genome.keys():
            print >> sys.stderr, "chr %s not in genome" % chr
            raise ValueError
    else:
        print >> sys.stderr, "range not formated as expected : \"chr:deb-fin\""
        raise ValueError

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
            print >> sys.stderr, \
                "range %s is incoherent for chr %s. FYI chr %s contains %s elements. Be sure that beginning < end of range" % ([range[0]+1, range[1]], chr, chr, len(genome[chr]))
            raise ValueError
    else:
        raise ValueError
    return (chr, range)


def TbComputeHomologyInformations(chrom1_tb, chrom2_tb):
    ###
    # Build MHP = { i1_tb : {i2_tb : hpSign, ...}...} with hpSign
    # the hpSign (s1*s2) of the homology at the (i1,i2) coordinate
    ###
    (MHP, locG2) = myDiags.homologyMatrix(chrom1_tb, chrom2_tb)

    tbsOnChr2ThatHaveHomologies = set([])
    for tb1 in MHP:
        tbsOnChr2ThatHaveHomologies |= set(MHP[tb1].keys())

    ###
    # Build TBNoHomologiesInWindowX = [..., i_tb, ... ] with i_tb the index of
    # the tb of chromosomeX_tb with no homology in the window
    ###
    TbNoHomologiesInWindowC1 = []
    for (i1_tb, _) in enumerate(chrom1_tb):
        if i1_tb not in MHP:
            TbNoHomologiesInWindowC1.append(i1_tb)
    TbNoHomologiesInWindowC2 = []
    for (i2_tb, _) in enumerate(chrom2_tb):
        if i2_tb not in tbsOnChr2ThatHaveHomologies:
            TbNoHomologiesInWindowC2.append(i2_tb)
    ###
    # Build homologyGroupsInWindow : [..., ([tbC1_4,tbC1_46,tbC1_80],[tbC2_2,tbC2_7]), ...]
    #   a list of 2-uples of homologous tbs indices to be coloured in the same color
    #      with tbCX_X : [..., i, ...] list of indices of genes in the tb
    ###
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


def genesComputeHomologyInformations(chr1, chr2, chrom1, chrom2, families,
                                     filterType,
                                     minChromLength, tandemGapMax):
    c1_fID = myMapping.labelWithFamID(chrom1, families)
    c2_fID = myMapping.labelWithFamID(chrom2, families)
    # Must be applied on the two genomes
    if len((c1_fID.getGeneNames(checkNoDuplicates=False) & c2_fID.getGeneNames(checkNoDuplicates=False)) - {None}) == 0:
        print >> sys.stderr, "Warning, no homologs"
    ((c1_fID_filt, Cf2CfID1, (nCL1, nGL1)),
     (c2_fID_filt, Cf2CfID2, (nCL2, nGL2))) = myDiags.filter2D(c1_fID, c2_fID,
                                                               filterType,
                                                               minChromLength,
                                                               keepOriginal=True)
    (chrom1_tb, Ctb2CfID1, nGTD1) = myMapping.remapRewriteInTb(c1_fID_filt,
                                                             tandemGapMax=tandemGapMax,
                                                             mOld=Cf2CfID1)
    (chrom2_tb, Ctb2CfID2, nGTD2) = myMapping.remapRewriteInTb(c2_fID_filt,
                                                             tandemGapMax=tandemGapMax,
                                                             mOld=Cf2CfID2)
    # Ctb2CfID1 = {}
    # for c in Ctb2Cf1:
    #     # see Mapping class addition
    #     Ctb2CfID1[c] = Ctb2Cf1[c] + Cf2CfID1[c]
    # Ctb2CfID2 = {}
    # for c in Ctb2Cf2:
    #     # see Mapping class addition
    #     Ctb2CfID2[c] = Ctb2Cf2[c] + Cf2CfID2[c]



    #Focus on the chromosome of the window
    chrom1_ = chrom1[chr1]
    chrom2_ = chrom2[chr2]
    c1_fID = c1_fID[chr1]
    c2_fID = c2_fID[chr2]
    # c1_fID_filt = c1_fID_filt[chr1]
    # c2_fID_filt = c2_fID_filt[chr2]
    # Cf2CfID1 = Cf2CfID1[chr1]
    # Cf2CfID2 = Cf2CfID2[chr2]
    chrom1_tb = chrom1_tb[chr1] if chr1 in chrom1_tb else []
    chrom2_tb = chrom2_tb[chr2] if chr2 in chrom2_tb else []
    # Ctb2Cf1 = Ctb2Cf1[chr1]
    # Ctb2Cf2 = Ctb2Cf2[chr2]

    print >> sys.stderr, Ctb2CfID1.keys()
    Ctb2CfID1 = Ctb2CfID1[chr1]
    Ctb2CfID2 = Ctb2CfID2[chr2]

    ###
    # Build genesRemovedDuringFilteringCX = [..., i, ...] the list of genes that
    # have been removed during the filtering process
    ###
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

    ###
    # Build genesNoHomologiesInWindow1 = [..., [i5,i6,i7], ...] list of tbs with [i5,i6,i7] a tb of three genes whose indices are i5,i6 and i7
    ###
    genesNoHomologiesInWindowC1 = []
    genesNoHomologiesInWindowC2 = []
    for i1_tb in TbNoHomologiesInWindowC1:
        genesNoHomologiesInWindowC1.append([i1 for i1 in Ctb2CfID1[i1_tb]])
    for i2_tb in TbNoHomologiesInWindowC2:
        genesNoHomologiesInWindowC2.append([i2 for i2 in Ctb2CfID2[i2_tb]])

    ###
    # Build genesHomologyGroupsInWindow : [..., ([tbC1_4,tbC1_46,tbC1_80],[tbC2_2,tbC2_7]), ...]
    #   a list of 2-uples of homologous tbs indices to be coloured in the same color
    #      with tb_CX_X : [..., i, ...] list of indices of genes in the tb
    ###
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


# Generator of levels for colors or gray indices within a palette:
# farIdxs may be an int. The more this int is high, the more neighbour color will be different
class levelIdxGenerator():
    def __init__(self, farIdxs=None, grays=False):
        if grays is not None:
            # see the css joint file "HomologyGroup" ranges from 0 to 44 (included)
            self.firstLevelIdx = 3
            self.lastLevelIdx = 44
        else:
            # see the css joint file "NoHomologyInWindow" ranges from 0 to 14 (included)
            self.firstLevelIdx = 3
            self.lastLevelIdx = 14

        if farIdxs is None:
            self.availableLevels = range(self.firstLevelIdx, self.lastLevelIdx+1)
        else:
            assert farIdxs < self.lastLevelIdx
            self.availableLevels = []
            tmp = range(0, farIdxs)
            random.shuffle(tmp)
            for i in tmp:
                self.availableLevels += range(self.firstLevelIdx + i, self.lastLevelIdx + 1, farIdxs)
            print >> sys.stderr, self.availableLevels
        self.currIdx = 0
        assert len(self.availableLevels) == (self.lastLevelIdx+1) - self.firstLevelIdx

    def getLevel(self, differentFrom=set()):
        if len(differentFrom.intersection(self.availableLevels)) == len(self.availableLevels):
            print >> sys.stderr, "Warning: too many colors too avoid, thus a non-optimal choice of the color is made"
        else:
            while self.currIdx in differentFrom:
                if self.currIdx < self.lastLevelIdx - self.firstLevelIdx:
                    self.currIdx += 1
                else:
                    self.currIdx = 0
        level = self.availableLevels[self.currIdx]
        self.currIdx += 1
        if self.currIdx > self.lastLevelIdx - self.firstLevelIdx:
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


def prepareChromosome(genesStrands,
                      familyName2Idxs=None,
                      familyName2color=None,
                      tbsWithNoHomolog=None,
                      genesRemovedDuringFiltering=None,
                      noHomologGraysGenerator=None,
                      symbolsInGenes=None,
                      lengthGene=1):

    def giveGreyLevelsTo(chromosome, tbsWithNoHomologyInWindow=None):
        if tbsWithNoHomologyInWindow:
            for tb in tbsWithNoHomologyInWindow:
                # Choose a level different from the direct neighbours
                nLevels = reduce(lambda x, y: x | y, [neighboursLevels(chromosome, i) for i in tb])
                grey = noHomologGraysGenerator.getLevel(differentFrom=nLevels)
                for i in tb:
                    chromosome[i].SVGclass = "NoHomologyInWindow%s" % grey
        return chromosome

    chromosomeItems = []

    # homologousTbs = [homolog1[tb1=[gene1Idx, ...], tb2=[geneAIdx, ...]],
    #                  homolog2[tb1'=[gene1'Idx, ...], tb2'=[geneA'Idx, ...]],
    #                 ... ]
    if not noHomologGraysGenerator:
        noHomologGraysGenerator = levelIdxGenerator(farIdxs=None, grays=True)
    if not familyName2color:
        homologsColorsGenerator = levelIdxGenerator(farIdxs=None)

    # create chromosomes
    for (i, s) in enumerate(genesStrands):
        cx = i * lengthGene
        symbol = symbolsInGenes[i] if symbolsInGenes is not None else None
        chromosomeItems.append(mySvgDrawer.Gene(Point(cx, 0),
                                                Point(cx + lengthGene, 0),
                                                strand=s, width=lengthGene*0.7, stroke_width=0.05*lengthGene, SVGclass=None,
                                                text=symbol))

    # give grey levels to genes that have no homology in the window
    chromosomeItems = giveGreyLevelsTo(chromosomeItems, tbsWithNoHomolog)

    if genesRemovedDuringFiltering:
        for i in genesRemovedDuringFiltering:
            chromosomeItems[i].SVGclass = "SpeciesSpecificGenes"

    if familyName2Idxs:
        for (famName, idxGs) in familyName2Idxs.iteritems():
            if familyName2color:
                color = familyName2color[famName]
            else:
                nLevels = reduce(lambda x, y: x | y, [neighboursLevels(chromosomeItems, idxG) for idxG in idxGs])
                color = homologsColorsGenerator.getLevel(differentFrom=nLevels)
            for idxG in idxGs:
                # for idxG in idxTb:
                chromosomeItems[idxG].SVGclass = "HomologGroup%s" % color
    return chromosomeItems

def drawLightGenome(genome,
                    families=None,
                    familyName2color=None,
                    filterType='InFamilies',
                    tandemGapMax=None,
                    lengthGene=1,
                    homologsColorsGenerator=None):
    assert filterType in {None, 'InFamilies'}
    assert isinstance(genome, myLightGenomes.LightGenome)
    assert families is None or isinstance(families, myLightGenomes.Families)

    # rewrite with families
    if families:
        keepGnOfGenesNotInFamilies = True if filterType is None else False
        genome_fam = myMapping.labelWithFamNames(genome, families,
                                                 keepGnOfGenesNotInFamilies=keepGnOfGenesNotInFamilies)
    else:
        genome_fam = genome

    if filterType == 'InFamilies':
        # remove all genes that are not in families
        (genome_famf, Cfamf2C, (nbChrLoss, nbGeneLoss)) = \
            myMapping.remapFilterGeneContent(genome_fam, removedNames={None}, mOld=None)
    else:
        Cfamf2C = {}
        for c in genome.keys():
            Cfamf2C[c] = myMapping.Mapping([[i] for i, _ in enumerate(genome[c])])
        genome_famf = genome
    print >> sys.stderr, genome_famf

    # rewrite in tandem blocks
    if tandemGapMax:
        (genome_tb, Ctb2Cfamf, nGTD) = myMapping.remapRewriteInTb(genome_famf, tandemGapMax=tandemGapMax, mOld=None)
    else:
        Ctb2Cfamf = {}
        for c in genome.keys():
            Ctb2Cfamf[c] = myMapping.Mapping([[i] for i, _ in enumerate(genome_famf[c])])
        genome_tb = genome_famf
    Ctb2Cfam = {}
    for c in Ctb2Cfamf:
        # see Mapping class addition
        Ctb2Cfam[c] = Ctb2Cfamf[c] + Cfamf2C[c]

    if not homologsColorsGenerator:
        homologsColorsGenerator = levelIdxGenerator(farIdxs=5)
    tmp_familyName2color = {}
    familyName2Idxs = collections.defaultdict(lambda: collections.defaultdict(list))
    genome_tb.computeDictG2Ps()
    for famName in sorted(genome_tb.getOwnedFamilyNames(families, asA=set)):
        tmp_familyName2color[famName] = homologsColorsGenerator.getLevel()
        famPositions = genome_tb.getPositions(famName, default=None)
        for tbPos in famPositions:
            for i in Ctb2Cfam[tbPos.c][tbPos.idx]:
                familyName2Idxs[tbPos.c][famName].append(i)
    if familyName2color:
        assert set(tmp_familyName2color.keys()) <= set(familyName2color.keys())
    else:
        familyName2color = tmp_familyName2color

    genomeItems = collections.OrderedDict()
    genome_tb.computeDictG2Ps()
    for (chr, chrom) in genome.iteritems():
        genesStrands = [s for (_, s) in chrom]

        # genesRemovedDuringFilteringCX = [..., i, ...]
        # the list of genes that have been removed during the filtering process
        genesRemovedDuringFiltering = [i1 for (i1, (anc, _)) in enumerate(genome_fam[chr]) if anc is None]

        # symbolsInGenes : [ 4,5,1,1,6,2, ...]
        # number of genes in each TB of C
        symbolsInGenes = []
        if tandemGapMax is None:
            for (i_tb, g) in enumerate(genome_fam[chr]):
                if g.n is not None:
                    symbolsInGenes.append(g.n)
                else:
                    symbolsInGenes.append('')
        else:
            for (i_tb, g) in enumerate(genome_fam[chr]):
                tbSize = len(Ctb2Cfam[chr][i_tb])
                symbolsInGenes.append("%s%s" % (tbSize, g.n))

        # add a line that goes through all genes of the chromosome
        genomeItems[chr] = []
        genomeItems[chr].append(mySvgDrawer.Line(Point(0, 0), Point(lengthGene * len(genome[chr]), 0)))
        genomeItems[chr].extend(prepareChromosome(genesStrands,
                                                  familyName2Idxs=familyName2Idxs[chr],
                                                  familyName2color=familyName2color,
                                                  tbsWithNoHomolog=None,
                                                  genesRemovedDuringFiltering=genesRemovedDuringFiltering,
                                                  noHomologGraysGenerator=None,
                                                  symbolsInGenes=symbolsInGenes,
                                                  lengthGene=lengthGene))
    return genomeItems

def drawChromosomes(genesStrandsC1, tbWithNoHomologyInWindowC1, genesRemovedDuringFilteringC1,
                    genesStrandsC2, tbWithNoHomologyInWindowC2, genesRemovedDuringFilteringC2,
                    homologyGroupsInWindow, closeColorsGenerator,
                    symbolsInGenes, sizeCase, height):

    chromosome1 = prepareChromosome(genesStrandsC1, None, None, tbWithNoHomologyInWindowC1, genesRemovedDuringFilteringC1,
                                    symbolsInGenes=symbolsInGenes[0] if symbolsInGenes else None,
                                    lengthGene=sizeCase)

    chromosome2 = prepareChromosome(genesStrandsC2, None, None, tbWithNoHomologyInWindowC2, genesRemovedDuringFilteringC2,
                                    symbolsInGenes=symbolsInGenes[1] if symbolsInGenes else None,
                                    lengthGene=sizeCase)

    # give a color to each gene using homology relationships
    for (tbs1, tbs2) in homologyGroupsInWindow:
        # Choose a level different from the direct neighbours
        neighboursLevelsOnBothsChrs = reduce(lambda x, y: x | y, [neighboursLevels(chromosome1, i1) for tb1 in tbs1 for i1 in tb1]) | \
                                      reduce(lambda x, y: x | y, [neighboursLevels(chromosome2, i2) for tb2 in tbs2 for i2 in tb2])
        color = closeColorsGenerator.getLevel(differentFrom=neighboursLevelsOnBothsChrs)
        for tb1 in tbs1:
            for i1 in tb1:
                chromosome1[i1].SVGclass = "HomologGroup%s" % color
        for tb2 in tbs2:
            for i2 in tb2:
                chromosome2[i2].SVGclass = "HomologGroup%s" % color

    # insert at the beginning to draw first
    chromosome1.insert(0, mySvgDrawer.Line(Point(0, 0), Point(len(genesStrandsC1) * sizeCase, 0)))
    chromosome2.insert(0, mySvgDrawer.Line(Point(0, 0), Point(len(genesStrandsC2) * sizeCase, 0)))
    return (chromosome1, chromosome2)

def drawMatrix(nx, ny, (begC1, endC1), (begC2, endC2), hpSigns, diagsIndices, sizeCell, width, height,
               diagColorGenerator=None, scaleFactorRectangles=1.0):
    print >> sys.stderr, scaleFactorRectangles
    assert isinstance(scaleFactorRectangles, float)
    sizeText = float(sizeCell*0.9)
    listOfMatrixItems = []
    if not diagColorGenerator:
        diagColorGenerator = levelIdxGenerator(farIdxs=True)

    nbLinesX = nx+1 #Nb of vertical lines (x varies) in the matrix
    nbLinesY = ny+1 #Nb of horizontal lines (y varies) in the matrix

    # draw Diagonals first because they are on the background
    print >> sys.stderr, "Nb of diagonals showed = ", len(diagsIndices)

    # sort diagonals to give colors according to localisations of diagonals
    diagsIndices.sort(key=lambda x: x[0])
    for diag in diagsIndices:
        # choose a color different from the neighbours
        color = diagColorGenerator.getLevel()
        for (i, j) in diag:
            cx_s = i*sizeCell
            cy_s = j*sizeCell
            if nx >= 300 or ny >= 300:
                listOfMatrixItems.append(mySvgDrawer.Rectangle(Point(cx_s, height-(cy_s+sizeCell)),
                                                               sizeCell, sizeCell, fill_opacity=0.90, svgClass="HomologGroup%s" % color,
                                                               bigger=scaleFactorRectangles))
            else:
                listOfMatrixItems.append(mySvgDrawer.Rectangle(Point(cx_s, height-(cy_s+sizeCell)),
                                                               sizeCell, sizeCell, fill_opacity=0.90, svgClass="HomologGroup%s" % color))
                #sizeCell(scgDrw.Rectangle((cx_s,height-(cy_s+sizeCell)), sizeCell, sizeCell, fill_opacity=0.90, svgClass = "chromgrp%s" % color ))
                #listOfMatrixItems.append(mySvgDrawer.Rectangle((cx_s,height-(cy_s+sizeCell)), sizeCell, sizeCell, fill_opacity=fill_opacity=0.90, color=(color,color,color))
        # draw rectangles around diagonals
        min_i = min(diag, key=lambda x: x[0])[0]
        max_i = max(diag, key=lambda x: x[0])[0]
        min_j = min(diag, key=lambda x: x[1])[1]
        max_j = max(diag, key=lambda x: x[1])[1]
        cx_s = min_i*sizeCell
        cy_s = max_j*sizeCell

        scaleFactorBoundingBoxDiags = 1.5 * scaleFactorRectangles if nx >= 300 or ny >= 300 else 1.0
        listOfMatrixItems.append(mySvgDrawer.Rectangle(Point(cx_s, height-(cy_s+sizeCell)),
                                                       (max_j-min_j)*sizeCell + sizeCell, (max_i-min_i)*sizeCell + sizeCell,
                                                       stroke='black', fill='none', strokeWidth=0.2*sizeCell*scaleFactorBoundingBoxDiags))

    # tick lines
    widthTicks = min(float(width)/1000, float(height)/1000)
    sizeTextTicks = widthTicks*10
    for (i, ni) in enumerate(range(begC1, endC1)):
        cx = i*sizeCell
        if ni % 10 == 0:
            # ticks
            listOfMatrixItems.append(mySvgDrawer.Line(Point(sizeCell/2+cx, height + sizeCell/2),
                                                      Point(sizeCell/2+cx, height), width=widthTicks))
        if ni % 50 == 0:
            cyText = height - max(sizeCell/2, sizeTextTicks/2)
            cxx = cx + sizeCell/2
            if nx > 750 or ny > 750:
                if ni % 100 == 0:
                    # TODO FIXME, too high for long chromosomes
                    listOfMatrixItems.append(mySvgDrawer.Text(Point(cxx, cyText + sizeCell), str(ni), text_anchor="middle", size=sizeTextTicks))
                    #listOfMatrixItems.append(mySvgDrawer.Text(Point(cxx, cyText), str(ni), text_anchor="middle", size=sizeTextTicks))
                    listOfMatrixItems.append(mySvgDrawer.Line(Point(cxx, ny * sizeCell),
                                                              Point(cxx, sizeCell), width=sizeCell*0.1))
            else:
                listOfMatrixItems.append(mySvgDrawer.Text(Point(cxx, cyText + sizeCell), str(ni), text_anchor="middle", size=sizeTextTicks))
                listOfMatrixItems.append(mySvgDrawer.Line(Point(cxx, height),
                                                          Point(cxx, height - ny * sizeCell), width=sizeCell*0.1))

    for (j, nj) in enumerate(range(begC2, endC2)):
        cy = j*sizeCell
        if nj % 10 == 0:
            listOfMatrixItems.append(mySvgDrawer.Line(Point(sizeCell/2 - sizeCell, height - (sizeCell/2+cy)),
                                                      Point(0, (height - (sizeCell/2+cy))), width=widthTicks))
        if nj % 50 == 0:
            cxText = - max(sizeCell/2, sizeTextTicks/2)
            cyy = height - (sizeCell/2 + cy)
            if nx > 750 or ny > 750:
                if nj % 100 == 0:
                    # TODO FIXME, not visible for big chromosomes
                    #listOfMatrixItems.append(mySvgDrawer.Text((cxx, cyy), str(nj), text_anchor="middle", size=sizeTextTicks, transform="translate(%s) rotate(90,%s,%s)" % (sizeCell,cxx+sizeCell,cyy)))
                    listOfMatrixItems.append(mySvgDrawer.Text(Point(0, 0),  str(nj),
                                                              text_anchor="middle", size=sizeTextTicks, transform="translate(%s,%s) rotate(-90)" % (cxText, cyy - 6*sizeCell)))
                    # listOfMatrixItems.append(mySvgDrawer.Text(Point(cxText, cyText), str(nj),
                    #                       text_anchor="middle", size=sizeTextTicks, transform="rotate(90,%s,%s)" % (cxText, cyText)))
                    listOfMatrixItems.append(mySvgDrawer.Line(Point(0, cyy),
                                                              Point(width-sizeCell, cyy), width=sizeCell*0.1))
            else:
                #listOfMatrixItems.append(mySvgDrawer.Text((cxx, cyy), str(nj), text_anchor="middle", size=sizeTextTicks, transform="translate(%s) rotate(90,%s,%s)" % (sizeCell,cxx+sizeCell,cyy)))
                # FIXME: why - 6*sizeCell ?
                listOfMatrixItems.append(mySvgDrawer.Text(Point(0, 0),  str(nj),
                                                          text_anchor="middle", size=sizeTextTicks, transform="translate(%s,%s) rotate(90)" % (cxText, cyy - 6*sizeCell)))
                listOfMatrixItems.append(mySvgDrawer.Line(Point(0, cyy),
                                                          Point(nx * sizeCell, cyy), width=sizeCell*0.1))

    if nx < 300 and ny < 300:
        for i in range(nbLinesX):
            cxLine = i*sizeCell
            listOfMatrixItems.append(mySvgDrawer.Line(Point(cxLine, height),
                                                      Point(cxLine, height-((nbLinesY-1)*sizeCell)), width=sizeCell*0.01))
        for j in range(nbLinesY):
            cyLine = j*sizeCell
            listOfMatrixItems.append(mySvgDrawer.Line(Point(0, height-(cyLine)),
                                                      Point((nbLinesX-1)*sizeCell, height-(cyLine)), width=sizeCell*0.01))

        # fill homologies with +1, -1 or ? or 0
        nonZeroValues = []
        for i1 in hpSigns:
            for i2 in hpSigns[i1]:
                nonZeroValues.append((i1, i2))
                #s = hpSigns[i1][i2][0][1]
                s = hpSigns[i1][i2]
                cx = i1*sizeCell + float(sizeCell)/2
                cy = i2*sizeCell + float(sizeCell)/2
                assert s == +1 or s == -1 or s is None, "s=%s" % s
                assocValue = (("+" if s == +1 else "-") if s is not None else '?')
                listOfMatrixItems.append(mySvgDrawer.Text(Point(cx, height-(cy+sizeText*0.16)),
                                                          assocValue, text_anchor="middle", size=sizeText))
        if nx < 20 and ny < 20:
            for (i1, i2) in itertools.product(range(nx), range(ny)):
                if (i1, i2) not in nonZeroValues:
                    cx = i1*sizeCell + float(sizeCell)/2
                    cy = i2*sizeCell + float(sizeCell)/2
                    listOfMatrixItems.append(mySvgDrawer.Text(Point(cx, height-(cy+sizeText*0.16)),
                                                              "0", text_anchor="middle", size=sizeText, fill=(200, 200, 200), stroke=None))
    else:
        # represent homologies with a black rectangle
        for i1 in hpSigns:
            for i2 in hpSigns[i1]:
                if (i1, i2) not in [dot for diag in diagsIndices for dot in diag]:
                    cx_s = i1*sizeCell
                    cy_s = i2*sizeCell
                    listOfMatrixItems.append(mySvgDrawer.Rectangle(Point(cx_s, height-(cy_s+sizeCell)),
                                                                   sizeCell, sizeCell, fill=(0, 0, 0), fill_opacity=0.90, bigger=scaleFactorRectangles))
        print >> sys.stderr, "Warning : some supplementary informations are not displayed because one of the two dimension of the window is > 300"

    return listOfMatrixItems


# draw either the mh or the mhp, if draw mode is 'writeinTB'
# inputs :
#       genesStrandsCX = [+1, -1, ...] of length = to nX
#       genesRemovedDuringFilteringC1 = [..., i, ...] with i the index of the gene removed during the filtering process (CX : Chromosome X)
#      tbWithNoHomologyInWindowC1 = [..., [i6,i7,i8], ...] list of tbs with no homologies in the window, inside the index of genes. If draw mode is 'writeinTB' : tbWithNoHomologyInWindowC1 = [..., i6, ...] : just the index of the TB
#       hpSigns = { i1 : {i2 : s, ...}...} with hpSign the hp sign (s1*s2) of the homology at the (i1,i2) coordinate (i1-th gene on C1 and i2-th gene on C2)
#      homologyGroupsInWindow = [..., ([tbC1_4,tbC1_46,tbC1_80],[tbC2_2,tbC2_7]), ...] a list of 2-uples of homologous tbs indices to be coloured in the same color
#              with tb_CX_X : [..., i, ...] list of indices of genes in the tb
#       diagIndices = [..., [...,(i16,j16),...], ...] list of diagonals with diagonals = list of all the points of the diagonal
# output :
#       string with the svg drawing of the mhp (or the mh)
def drawHomologyMatrix(((begC1, endC1), (begC2, endC2)), (genesStrandsC1, genesStrandsC2),
                       (genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2),
                       (tbWithNoHomologyInWindowC1, tbWithNoHomologyInWindowC2),
                       hpSigns, homologyGroupsInWindow, diagsIndices,
                       outputFileName=None, maxWidth=100, maxHeight=100, symbolsInGenes=None, scaleFactorRectangles=1):
    # example
    # genesStrandsC1 = [1,-1, None, ...]
    # tbWithNoHomologyInWindowC1=[0,2,3,6,10]
    # tbWithNoHomologyInWindowC2=[0,7]
    # homologyGroupsInWindow=[([1],[5]),([4],[4]),([5],[3]),([9],[1]),([7],[6]),([8],[2])]
    # diagsIndices=[[(4,4),(5,3),(8,2),(9,1)]]
    # genesRemovedDuringFilteringC1=[]
    # genesRemovedDuringFilteringC2=[]
    # symbolsInGenes=[[1,1,1,1,1,1,1,2,2,1,1],[1,1,2,1,1,1,2,1]]
    # hpSigns a collections.defaultdict() with hpSigns[i1][i2]=aV
    assert isinstance(genesStrandsC1, list) and isinstance(genesStrandsC2, list)
    assert isinstance(genesRemovedDuringFilteringC1, list) and isinstance(genesRemovedDuringFilteringC2, list)
    assert isinstance(tbWithNoHomologyInWindowC1, list) and isinstance(tbWithNoHomologyInWindowC2, list)
    assert isinstance(homologyGroupsInWindow, list)
    assert isinstance(hpSigns, dict)
    assert isinstance(diagsIndices, list)
    # For the print on screen
    begC1 = begC1 + 1
    begC2 = begC2 + 1
    endC1 = endC1 + 1
    endC2 = endC2 + 1
    # nx = number of genes on the first genome
    # ny = number of genes on the second genome
    nx = len(genesStrandsC1)
    ny = len(genesStrandsC2)

    # the size of the components of the matrix is chosen using the smallest and more restricting dimension
    # (contains more genes comparing to its size)
    sizeCase = float(min(float(maxWidth) / (nx + 3), float(maxHeight) / (ny + 3)))  # +3 for margins

    width = float(nx+3) * sizeCase
    height = float(ny+3) * sizeCase
    print >> sys.stderr, "width=", width
    print >> sys.stderr, "height=", height

    scene = mySvgDrawer.Scene(name='homology_matrix', width=width, height=height)

    closeColorsGenerator = levelIdxGenerator()

    # draw lines of chromosomes
    # offset_genes : corresponds to the chromosome lines positions passing through the middle of the genes
    offset_genes_x = sizeCase + sizeCase/2
    offset_genes_y = sizeCase + sizeCase/2
    # scene.add(mySvgDrawer.Line(Point(offset_genes_x, height - offset_genes_y),
    #                       Point(width-sizeCase, height - offset_genes_y), width=0.1*sizeCase))
    # scene.add(mySvgDrawer.Line(Point(offset_genes_x, sizeCase),
    #                       Point(offset_genes_x, height - (offset_genes_y)), width=0.1*sizeCase))
    offset_matrix_x = 2 * sizeCase
    offset_matrix_y = 2 * sizeCase

    listOfMatrixItems = drawMatrix(nx, ny, (begC1, endC1), (begC2, endC2), hpSigns, diagsIndices, sizeCase, width, height,
                                   diagColorGenerator=None, scaleFactorRectangles=scaleFactorRectangles)
    listOfMatrixItems = mySvgDrawer.tanslateItems(listOfMatrixItems, 2 * sizeCase, - 2 * sizeCase)

    listOfItems = []
    if nx < 300 and ny < 300:
        (chromosome1, chromosome2) = drawChromosomes(genesStrandsC1, tbWithNoHomologyInWindowC1, genesRemovedDuringFilteringC1,
                                                     genesStrandsC2, tbWithNoHomologyInWindowC2, genesRemovedDuringFilteringC2,
                                                     homologyGroupsInWindow, closeColorsGenerator,
                                                     symbolsInGenes, sizeCase,
                                                     # offset_genes_x, offset_genes_y,
                                                     # offset_matrix_x, offset_matrix_y,
                                                     height)

        chromosome1 = mySvgDrawer.tanslateItems(chromosome1, 2 * sizeCase, height)
        for item in chromosome2:
            if isinstance(item, mySvgDrawer.Gene):
                item.start = Point(0, height - item.start.x - 2 * sizeCase)
                item.end = Point(0, height - item.end.x - 2 * sizeCase)
            elif isinstance(item, mySvgDrawer.Line):
                item.start = Point(0, height - item.start.x - 2 * sizeCase)
                item.end = Point(0, height - item.end.x - 2 * sizeCase)
            else:
                raise TypeError
        listOfItems = chromosome1 + chromosome2

    listOfItems += listOfMatrixItems
    listOfItems = mySvgDrawer.tanslateItems(listOfItems, sizeCase, -sizeCase)

    for item in listOfItems:
        scene.add(item)
    if outputFileName is not None:
        scene.write_svg(filename=str(outputFileName))
    #scene.display()
    return scene.strarray()

FilterType = myDiags.FilterType

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
                         identifyMonoGenicInversion=False,
                         identifyBreakpointsWithinGaps=True,
                         overlapMax=None,
                         consistentSwDType=True,
                         validateImpossToCalc_mThreshold=3,
                         nbHpsRecommendedGap=2,
                         targetProbaRecommendedGap=0.01,
                         chromosomesRewrittenInTbs=False,
                         scaleFactorRectangles=2.0,
                         considerAllPairComps=True,
                         switchOnDirectView=False,
                         optimisation=None,
                         inSbsInPairComp=None,
                         outSyntenyBlocksFileName="./syntenyBlocksDrawer.txt",
                         outImageFileName="./homologyMatrix.svg",
                         verbose=True):

    # if True, this opens the output image in firefox at the end of the computation
    assert isinstance(genome1, myLightGenomes.LightGenome)
    assert isinstance(genome2, myLightGenomes.LightGenome)
    assert isinstance(families, myLightGenomes.Families)

    kwargs = {'gapMax': gapMax,
              'distinguishMonoGenicDiags': distinguishMonoGenicDiags,
              'gapMaxMicroInv': gapMaxMicroInv,
              'distanceMetric': distanceMetric,
              'identifyMonoGenicInversion': identifyMonoGenicInversion,
              'identifyBreakpointsWithinGaps': identifyBreakpointsWithinGaps,
              'overlapMax': overlapMax,
              'consistentSwDType': consistentSwDType,
              'pThreshold': pThreshold,
              'nbHpsRecommendedGap': nbHpsRecommendedGap,
              'targetProbaRecommendedGap': targetProbaRecommendedGap,
              'validateImpossToCalc_mThreshold': validateImpossToCalc_mThreshold,
              'optimisation': optimisation,
              'verbose': True}

    assert distanceMetric == 'DPD' or distanceMetric == 'MD' or distanceMetric == 'CD' or distanceMetric == 'ED'
    assert (convertGenicToTbCoordinates and chromosomesRewrittenInTbs) or not convertGenicToTbCoordinates
    # Change genome format
    genome1Name = genome1.name
    genome2Name = genome2.name

    if not chromosomesRewrittenInTbs:

        # define the ROI (Region Of Interest)
        (chr1, range1) = parseChrRange(CDF1, genome1)
        (chr2, range2) = parseChrRange(CDF2, genome2)
        chrom1 = myLightGenomes.LightGenome()
        chrom2 = myLightGenomes.LightGenome()
        chrom1[chr1] = genome1[chr1][range1[0]:range1[1]]
        chrom2[chr2] = genome2[chr2][range2[0]:range2[1]]
        nbSpeciesSpecificGenes1 = len([gn for gn in chrom1.getGeneNames(asA=list, checkNoDuplicates=False)
                                       if families.getFamilyByName(gn, default=None) is None])
        print >> sys.stderr, "the ROI1 contains %s genes (%s species specific genes)" % (len(chrom1[chr1]), nbSpeciesSpecificGenes1)
        nbSpeciesSpecificGenes2 = len([gn for gn in chrom2.getGeneNames(asA=list, checkNoDuplicates=False) if families.getFamilyByName(gn, default=None) is None])
        print >> sys.stderr, "the ROI2 contains %s genes (%s species specific genes)" % (len(chrom2[chr2]), nbSpeciesSpecificGenes2)

        if inSbsInPairComp is None:
            if considerAllPairComps:
                comparedGenome1 = genome1
                comparedGenome2 = genome2
            else:
                # extract diagonals in the ROI without considering other pairwise comparisons
                comparedGenome1 = chrom1
                comparedGenome2 = chrom2
            sbsInPairComp = myDiags.extractSbsInPairCompGenomes(comparedGenome1,
                                                                comparedGenome2,
                                                                families,
                                                                tandemGapMax=tandemGapMax,
                                                                minChromLength=minChromLength,
                                                                filterType=filterType,
                                                                **kwargs)
        else:
            sbsInPairComp = inSbsInPairComp

        new_sbsInPairComp = myTools.Dict2d(list)
        for sb in sbsInPairComp[chr1][chr2]:
            newl1 = []
            newl2 = []
            newla = []
            for (idxHp, aG) in enumerate(sb.la):
                tb1 = []
                tb2 = []
                for i1g in sb.l1[idxHp]:
                    if (range1[0] <= i1g and i1g <= range1[1]):
                        tb1.append(i1g)
                for i2g in sb.l2[idxHp]:
                    if (range2[0] <= i2g and i2g <= range2[1]):
                        tb2.append(i2g)
                if len(tb1) > 0 and len(tb2) > 0:
                    newl1.append(tb1)
                    newl2.append(tb2)
                    newla.append(aG)
            if len(newla) > 0:
                new_sbsInPairComp[chr1][chr2].append(myDiags.SyntenyBlock(myDiags.Diagonal(sb.dt, newl1, newl2, newla),
                                                                          sb.pVal))
        sbsInPairComp = new_sbsInPairComp

        genesDiagIndices = []
        for sb in sbsInPairComp[chr1][chr2]:
            genesDiagIndices.append([])
            for idxHp, aG in enumerate(sb.la):
                for (gi1, gi2) in itertools.product(sb.l1[idxHp], sb.l2[idxHp]):
                    assert range1[0] <= gi1 <= range1[1]
                    assert range2[0] <= gi2 <= range2[1]
                    genesDiagIndices[-1].append((gi1 - range1[0], gi2 - range2[0]))

        ###
        # Build Genes Strands
        ###
        genesStrandsC1 = [s for (_, s) in chrom1[chr1]]
        genesStrandsC2 = [s for (_, s) in chrom2[chr2]]

        ((genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2), genesHomologiesHpSign,
         (genesNoHomologiesInWindowC1, genesNoHomologiesInWindowC2), genesHomologyGroupsInWindow) = \
        genesComputeHomologyInformations(chr1, chr2, chrom1, chrom2,
                                         families,
                                         filterType,
                                         minChromLength,
                                         tandemGapMax)

        strArray = drawHomologyMatrix((range1, range2),
                                      (genesStrandsC1, genesStrandsC2),
                                      (genesRemovedDuringFilteringC1, genesRemovedDuringFilteringC2),
                                      (genesNoHomologiesInWindowC1, genesNoHomologiesInWindowC2),
                                      genesHomologiesHpSign,
                                      genesHomologyGroupsInWindow,
                                      genesDiagIndices,
                                      outputFileName=outImageFileName,
                                      maxWidth=100,
                                      maxHeight=100,
                                      scaleFactorRectangles=scaleFactorRectangles)

    else:
        assert chromosomesRewrittenInTbs
        ((g1_tb, mtb2g1, (nCL1, nGL1)), (g2_tb, mtb2g2, (nCL2, nGL2))) = \
            myDiags.editGenomes(genome1, genome2, families,
                                filterType=filterType, labelWith='FamID', tandemGapMax=tandemGapMax,
                                minChromLength=minChromLength, keepOriginal=True)
        if not convertGenicToTbCoordinates:
            (chr1, range1) = parseChrRange(CDF1, g1_tb)
            (chr2, range2) = parseChrRange(CDF2, g2_tb)
        else:
            mg2tb1 = dict((c, m.old) for (c, m) in mtb2g1.iteritems())
            mg2tb2 = dict((c, m.old) for (c, m) in mtb2g2.iteritems())
            (chr1, range1) = parseChrRange(CDF1, genome1, g2gtb=mg2tb1)
            (chr2, range2) = parseChrRange(CDF2, genome2, g2gtb=mg2tb2)

        chrom1_tb = myLightGenomes.LightGenome()
        chrom2_tb = myLightGenomes.LightGenome()
        chrom1_tb[chr1] = g1_tb[chr1][range1[0]:range1[1]]
        chrom2_tb[chr2] = g2_tb[chr2][range2[0]:range2[1]]
        print >> sys.stderr, "the ROI1 contains %s genes (%s genes deleted during edition)" % \
                             (sum(len(mtb2g1[chr1][itb]) for (itb, _) in enumerate(chrom1_tb[chr1])), nGL1)
        print >> sys.stderr, "the ROI2 contains %s genes (%s genes deleted during edition)" % \
                             (sum(len(mtb2g2[chr2][itb]) for (itb, _) in enumerate(chrom2_tb[chr2])), nGL2)
        #Focus on the chromosome of the window, just give simple name to the chromosome of interest
        tb2g1 = mtb2g1[chr1]
        g2tb1 = mtb2g1[chr1].old
        tb2g2 = mtb2g2[chr2]
        g2tb2 = mtb2g2[chr2].old

        # load precomputed sbs if any
        if inSbsInPairComp is None:
            if considerAllPairComps:
                comparedGenome1 = g1_tb
                comparedGenome2 = g2_tb
            else:
                # extract diagonals in the ROI without considering other pairwise comparisons
                comparedGenome1 = chrom1_tb
                comparedGenome2 = chrom2_tb
            print >> sys.stderr, kwargs
            inSbsInPairComp = myDiags.extractSbsInPairCompGenomesInTbs(comparedGenome1,
                                                                     comparedGenome2,
                                                                     **kwargs)
        else:
            for sb in inSbsInPairComp[chr1][chr2]:
                # change the sb.lX structure from list of lists to list of ints
                sb.l1 = [g2tb1[tb[0]] for tb in sb.l1]
                sb.l2 = [g2tb2[tb[0]] for tb in sb.l2]

        # truncate sbs to fit the ROI
        new_sbsInPairComp = myTools.Dict2d(list)
        for sb in inSbsInPairComp[chr1][chr2]:
            if (range1[0] <= sb.minOnG(1) and sb.maxOnG(1) <= range1[1]) \
                    and (range2[0] <= sb.minOnG(2) and sb.maxOnG(2) <= range2[1]):
                # sb is perfectly included in the ROI
                new_sbsInPairComp[chr1][chr2].append(sb)
            elif (sb.maxOnG(1) < range1[0] or range1[1] < sb.minOnG(1)) \
                    or (sb.maxOnG(2) < range2[0] or range2[1] < sb.minOnG(2)):
                # sb is not in the ROI
                continue
            else:
                # sb is partially included in the ROI
                sb.truncate(range1, range2)
                if len(sb.la) > 0:
                    new_sbsInPairComp[chr1][chr2].append(sb)
        sbsInPairComp = new_sbsInPairComp

        ###
        # Build TbNumberOfGenesInEachTbC1 : [ 4,5,1,1,6,2, ...] number og genes in each TB of C1
        ###
        TbNumberOfGenesInEachTbC1 = [len(tb2g1[i1_tb]) for i1_tb in range(len(chrom1_tb[chr1]))]
        TbNumberOfGenesInEachTbC2 = [len(tb2g2[i2_tb]) for i2_tb in range(len(chrom2_tb[chr2]))]

        ###
        # Build TBStrands
        ###
        TbStrandsC1 = [s for (_, s) in chrom1_tb[chr1]]
        TbStrandsC2 = [s for (_, s) in chrom2_tb[chr2]]

        ###
        # Build rangeXTB
        ###
        (TbHpSign, (TbNoHomologiesInWindowC1, TbNoHomologiesInWindowC2), TbHomologyGroupsInWindow) = \
            TbComputeHomologyInformations(chrom1_tb[chr1], chrom2_tb[chr2])
        ###
        # Convert into the correct format for the function TbHomologyGroupsInWindow
        ###
        tmpTbHomologyGroupsInWindow = []
        for (tbs1, tbs2) in TbHomologyGroupsInWindow:
            tmpTbs1 = []
            tmpTbs2 = []
            for tb1 in tbs1:
                tmpTbs1.append([tb1])
            for tb2 in tbs2:
                tmpTbs2.append([tb2])
            tmpTbHomologyGroupsInWindow.append((tmpTbs1, tmpTbs2))
        TbHomologyGroupsInWindow = tmpTbHomologyGroupsInWindow
        TbNoHomologiesInWindowC1 = [[tb1] for tb1 in TbNoHomologiesInWindowC1]
        TbNoHomologiesInWindowC2 = [[tb2] for tb2 in TbNoHomologiesInWindowC2]

        TbDiagIndices = []
        for sb in sbsInPairComp[chr1][chr2]:
            TbDiagIndices.append([])
            for idxHp, aG in enumerate(sb.la):
                idxHp1 = sb.l1[idxHp]
                idxHp2 = sb.l2[idxHp]
                assert isinstance(idxHp1, int) and isinstance(idxHp2, int)
                assert range1[0] <= idxHp1 <= range1[1]
                assert range2[0] <= idxHp2 <= range2[1]
                TbDiagIndices[-1].append((idxHp1 - range1[0], idxHp2 - range2[0]))

        strArray = \
            drawHomologyMatrix((range1, range2),
                               (TbStrandsC1, TbStrandsC2),
                               ([], []),
                               (TbNoHomologiesInWindowC1, TbNoHomologiesInWindowC2),
                               TbHpSign,
                               TbHomologyGroupsInWindow,
                               TbDiagIndices,
                               outputFileName=outImageFileName,
                               maxWidth=100,
                               maxHeight=100,
                               symbolsInGenes=(TbNumberOfGenesInEachTbC1, TbNumberOfGenesInEachTbC2),
                               scaleFactorRectangles=scaleFactorRectangles
                               )

    #copy the css style sheet
    dirNameImage = os.path.dirname(outImageFileName)
    dirNameImage = dirNameImage if dirNameImage != "" else "."
    print >> sys.stderr, "cp %s/styleForHomologyMatrixWithSBs.css %s/" % (os.path.dirname(os.path.realpath(sys.argv[0])), dirNameImage)
    os.system("cp %s/styleForHomologyMatrixWithSBs.css %s/" % (os.path.dirname(os.path.realpath(sys.argv[0])), dirNameImage))

    # write a simple file with all diagonals into output file
    with open(outSyntenyBlocksFileName, 'w') as f:
        print >> f, "Mode : %s" % 'Genic scale' if chromosomesRewrittenInTbs is False else 'Tandem Blocks scale'
        print >> f, "chromosome %s de %s\t%s\t%s\tchromosome %s de %s\t%s\t%s\t%s" % (chr1, genome1Name, 'beginC1', 'endC1', chr2, genome2Name, 'beginC2', 'endC2', 'length in families')
        print >> f, "c1\tbeg1\tend1\tc2\tbeg2\tend2\thps\tpVal"

        for sb in sbsInPairComp[chr1][chr2]:
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

    # Add lengends and title to the ouput matrix
    height = 100
    width = 100
    var = ['<?xml version="1.0" encoding="utf-8" standalone="no"?>\n',
           #'<?xml-stylesheet type="text/css" href="styleForHomologyMatrixWithSBs.css" ?>\n', #Warning : request the css file
           '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n',
           "<svg height=\"100%%\" version=\"1.1\" viewBox=\"0 0 %s %s\" width=\"100%%\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n" % (width, height),
           '<defs>\n',
           '<style type="text/css">\n',
           '*{stroke-linecap:square;stroke-linejoin:round;}\n',
           '</style>\n',
           '</defs>\n'
           '<g style="fill-opacity:1.0; stroke:black;\n',
           'stroke-width:1;">\n']
    #Title
    # if not arguments['mode:chromosomesRewrittenInTbs']:

    title = \
        "%s, f=%s, tgm=%s tbs, gm=%s%s, gmmi=%s, ibwg=%s, om=%s, %s sbs" % \
        ('MHP' if chromosomesRewrittenInTbs else 'MH',
         filterType,
         tandemGapMax,
         gapMax,
         distanceMetric,
         gapMaxMicroInv,
         identifyBreakpointsWithinGaps,
         overlapMax,
         len(list(sbsInPairComp.iteritems2d())))
    # else:
    #     title =\
    #         "%s, tandemGapMax=%s tbs, gapMax=%s%s, overlapMax=%s tbs, %s sbs" %\
    #         ('MHP',
    #          arguments['tandemGapMax'],
    #          arguments['gapMax'],
    #          arguments['distanceMetric'],
    #          arguments['overlapMax'],
    #          len(list(sbsInPairComp.iteritems2d())))\

    var += [#'<svg x="5" y="0" viewBox="5 0 95 5" width="95" height="5" xmlns="http://www.w3.org/2000/svg" xmlns:xhtml="http://www.w3.org/1999/xhtml" xmlns:xlink="http://www.w3.org/1999/xlink">\n',
            '<svg x="5" y="0" viewBox="5 0 95 5" width="95" height="5">\n',
            ''.join(mySvgDrawer.Text(mySvgDrawer.Point(float(5 + 95)/2.0, 5.0/2.0), title,
                                     text_anchor='middle', fontWeight=300, size=2).strarray()),
            # '<foreignObject x="0" y="0" width="95" height="5">\n',
            # '<xhtml:div style="display:table; height:100%; width:100%; overflow:hidden; text-align:center; ">\n',
            #         '<xhtml:div style="display: table-cell; vertical-align: middle;">\n',
            #         '<xhtml:div style="color:black; word-wrap:break-word; font-size:1.5px; font-family:Arial" >' + title + '\n',
            #                         '</xhtml:div>\n',
            #                 '</xhtml:div>\n',
            #         '</xhtml:div>\n',
            # '</foreignObject>\n',
            '</svg>\n']

    #Add legends (genomes names and ranges)
    var += [#'<svg x="0" y="0" viewBox="0 0 5 95" width="5" height="95" xmlns="http://www.w3.org/2000/svg" xmlns:xhtml="http://www.w3.org/1999/xhtml" xmlns:xlink="http://www.w3.org/1999/xlink">\n',
            '<svg x="0" y="0" viewBox="0 0 5 95" width="5" height="95">\n',
            ''.join(mySvgDrawer.Text(mySvgDrawer.Point(5.0/2.0, float(5 + 95)/2.0), str(genome2.name) + " chr" + str(chr2) + ":" + str(range2[0]+1) + "-" + str(range2[1]),
                                     text_anchor='middle', size=2, fontWeight=300, transform='rotate(-90, %s, %s)' % (5.0/2.0, float(5 + 95)/2.0)).strarray()),
            # '<foreignObject x="0" y="0" width="95" height="5" transform="translate(5,0) rotate(90) translate(0,0)">\n',
            # '<xhtml:div style="display:table; height:100%; width:100%; overflow:hidden; text-align:center; ">\n',
            #         '<xhtml:div style="display: table-cell; vertical-align: middle;">\n',
            #         '<xhtml:div style="color:black; word-wrap:break-word; font-size:2px; font-family:Arial" >' + str(arguments["genome2"]) + " <br />  chr" + str(chr2) + ":" + str(range2[0]+1) + "-" + str(range2[1]) + '\n',
            #                         '</xhtml:div>\n',
            #                 '</xhtml:div>\n',
            #         '</xhtml:div>\n',
            # '</foreignObject>\n',
            '</svg>\n',
            '<svg x="5" y="95" viewBox="0 0 95 5" width="95" height="5">\n',
            ''.join(mySvgDrawer.Text(mySvgDrawer.Point(float(5 + 95)/2, 5.0/2.0),  str(genome1.name) + " chr" + str(chr1) + ":" + str(range1[0]+1) + "-" + str(range1[1]),
                                     text_anchor='middle', fontWeight=300, size=2).strarray()),
            # '<foreignObject x="0" y="0" width="95" height="5">\n',
            # '<xhtml:div style="display:table; height:100%; width:100%; overflow:hidden; text-align:center; ">\n',
            #         '<xhtml:div style="display: table-cell; vertical-align: middle;">\n',
            #         '<xhtml:div style="color:black; word-wrap:break-word; font-size:2px; font-family:Arial" >' + str(arguments["genome1"]) + " <br /> chr" + str(chr1) + ":" + str(range1[0]+1) + "-" + str(range1[1]) + '\n',
            #                         '</xhtml:div>\n',
            #                 '</xhtml:div>\n',
            #         '</xhtml:div>\n',
            # '</foreignObject>\n',
            '</svg>\n']

    var+=['<svg preserveAspectRatio="xMidYMid meet" x="5" y="5" viewBox="0 0 100 100" width="90" height="90" >\n'] # little transformation : viewBox = "the part of the forecoming images that we want to see", width and height = the width and height of the image that will be printed on the screen. This instructions takes a viewBox of the forecoming images and display it in a image of the specified width and height
    for line in strArray:
        if line.find("<?xml")>=0 or line.find("<!DOCTYPE")>=0: #or line.find("<svg")>=0 or line.find("</svg")>=0:
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
                print >> sys.stderr, "assymetric (on genome1 and genome2) SB of ", len([ anc for (anc,_,_,_) in la]), " HPs"
                print >> sys.stderr, "diag on G1 from %s to %s" % (min(l1[0][0]+range1[0], l1[-1][0]+range1[0]) , max(l1[0][0]+range1[0], l1[-1][0]+range1[0]))
                print >> sys.stderr, "diag on G2 from %s to %s" % (min(l2[0][0]+range2[0], l2[-1][0]+range2[0]) , max(l2[0][0]+range2[0], l2[-1][0]+range2[0]))
                if ecart < abs(len(l1) -  len(la)):
                    ecart = abs(len(l1) -  len(la))
                    diag_=diag
                if ecart < abs(len(l2) -  len(la)):
                    ecart = abs(len(l2) -  len(la))
                    diag_=diag
        diag = diag_ if diag_ != None else listOfDiags[0]
        l1 = diag[0][1]
        l2 = diag[1][1]
        la = diag[2]
        print >> sys.stderr, "The most assymetric diag"
        print >> sys.stderr, " SB of ", len([ anc for (anc,_,_,_) in la]), " HPs"
        print >> sys.stderr, "diag on G1 = ", [ gene[0]+range1[0] for gene in l1]
        print >> sys.stderr, "diag on G2 = ", [ gene[0]+range2[0] for gene in l2]

        minIndiceG1 = sys.maxint
        maxIndiceG1 = -sys.maxint-1
        minIndiceG2 = sys.maxint
        maxIndiceG2 = -sys.maxint-1
        for diag in listOfDiags :
            ((_,l1),(_,l2),_) = diag
            for (i1,_) in l1:
                minIndiceG1 = i1 if i1 < minIndiceG1 else minIndiceG1
                maxIndiceG1 = i1 if i1 > maxIndiceG1 else maxIndiceG1
            for (i2,_) in l2:
                minIndiceG2 = i2 if i2 < minIndiceG2 else minIndiceG2
                maxIndiceG2 = i2 if i2 > maxIndiceG2 else maxIndiceG2
        print >> sys.stderr, 'Indices entre lesquels se trouvents les diagonales sur le G1 = [%s,%s]' % (minIndiceG1+range1[0], maxIndiceG1+range1[0])
        print >> sys.stderr, 'Indices entre lesquels se trouvents les diagonales sur le G2 = [%s,%s]' % (minIndiceG2+range2[0], maxIndiceG2+range2[0])

        return

    # ask the user for the desired chromosome ranges
    @myTools.deprecated
    def chooseChrsAndRanges(genome1, genome2, families, distanceMetric = 'DPD'):
        while True:
            try:
                (chr1, range1) = parseChrRange(raw_input("chr1:deb1-fin1 = "), genome1)
                break
            except ValueError:
                print >> sys.stderr, "You need to write somtehing as chr1:deb1-fin1 with chr1 a chr of G1 and deb1 and fin1 indices of the first and last gene (indices start at 1)"

        while True:
            try:
                (chr2, range2) = parseChrRange(raw_input("chr2:deb2-fin2 = "), genome2)
                break
            except ValueError:
                print >> sys.stderr, "You need to write somtehing as chr2:deb2-fin2 with chr2 a chr of G2 and deb2 and fin2 indices of the first and last gene (indices start at 1)"

        chrom1 ={}
        chrom2 ={}
        chrom1[chr1] = genome1[chr1][range1[0]:range1[1]]
        chrom2[chr2] = genome2[chr2][range2[0]:range2[1]]
        listOfDiags = myDiags.extractSbsInPairCompGenomes(chrom1, chrom2, families,
                                                          gapMax=gapMax, consistentSwDType=consistentSwDType,
                                                          filterType=filterType, minChromLength=minChromLength,
                                                          distanceMetric=distanceMetric, verbose=verbose)
        listOfDiags = list(listOfDiags)
        print >> sys.stderr, "pairwise comparison of the two chromosomes yields" , len(listOfDiags), "diagonals."
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