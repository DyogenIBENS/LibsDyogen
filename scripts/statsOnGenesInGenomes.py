#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import utils.myTools
import utils.myGenomes

import sys, itertools, collections
from utils import myTools, myFile, myLightGenomes
from matplotlib.ticker import FuncFormatter
from scipy import stats

def intOrNone(string):
    return int(string) if string != 'None' else None

def formaterBasePairs(unit='bp'):
    'The two args are the value and tick position'
    def f(x, pos):
        if x == 0:
            res = '0'
        else:
            if unit == 'bp':
                factor = 1.0
            elif unit == 'kb':
                factor = 0.001
            elif unit == 'Mb':
                factor = 0.000001
            else:
                assert unit is None
                factor = 1
            res = '%1.0f%s' % (x * factor, unit if unit is not None else '')
        return res
    return f

__doc__ = """ plot different stats about genes """

# arguments
arguments = myTools.checkArgs(
    [
        ("genome", file),
    ],
    [
        ("removeUnofficialChrNames", bool, False),
        ("orderChromosomesBy", bool, 'names'),
        ('mode', bool, 'distribOnChr')
    ],
    __doc__)
assert arguments['orderChromosomesBy'] in {'decreasingNbOfGenes', 'names'}
assert arguments['mode'] in {'geneLengths', 'distribOnChr', 'distribOnChrs', 'overlap',
                             'correlationChromNbGenes', 'longestIntergene', 'minGeneLength'}
# longestIntergene computes the longer intergene, i.e. space between two genes.
# This gives the length of the longer rearrangement possibly unseen, except telomeres.

def readerDependingOnFileWithDebAndEnd(fileName):
        flb = myFile.firstLineBuffer(myFile.openFile(fileName, 'r'))
        c = flb.firstLine.split("\t")
        if len(c) == 6:
            print >> sys.stderr, "(c, beg, end, s, gName, transcriptName) -> (c, s, gName)"
            # c, beg, end, s,  gName, transcriptName
            reader = myFile.myTSV.readTabular(fileName, [str, str, str, str, str, str])
            reader = ((c, intOrNone(beg), intOrNone(end), intOrNone(strand), gName) for (c, beg, end, strand, gName, tName) in reader)
        elif len(c) == 5:
            print >> sys.stderr, "(c, beg, end, s, gName) -> (c, s, gName)"
            # c, beg, end, s,  gName
            tmpReader = myFile.myTSV.readTabular(fileName, [str, str, str, str, str])
            # check, with the first line, if there are several gene names (the format genome of Matthieu contains several gene names)
            (c, beg, end, strand, gNames) = tmpReader.next()
            severalNames = True if len(gNames.split(' ')) > 0 else False
            reader = itertools.chain([(c, beg, end, strand, gNames)], tmpReader)
            if severalNames:
                # if gNames contains more than one gene name, only take the first gene name
                reader = ((c, intOrNone(beg), intOrNone(end), intOrNone(strand), gNames.split(' ')[0]) for (c, beg, end, strand, gNames) in reader)
            else:
                reader = ((c, intOrNone(beg), intOrNone(end), intOrNone(strand), gName) for (c, beg, end, strand, gName) in reader)
        else:
            raise ValueError("%s file is badly formatted" % fileName)
        return reader

reader = readerDependingOnFileWithDebAndEnd(arguments['genome'])
unofficialChromRemoved = False
genomeListByChr = myTools.DefaultOrderedDict(list)
for (chr, x1, x2, s, gNames) in reader:
    if arguments['removeUnofficialChrNames'] and myLightGenomes.contigType(chr) != myLightGenomes.ContigType.Chromosome:
        unofficialChromRemoved = True
        continue
    if x1 < x2 and s == -1:
        print >> sys.stderr, 'Warning, genome not formated by the 5\' rule of x1=beg x2=end'
    elif x2 < x1 and s == +1:
        print >> sys.stderr, 'Warning, genome not formated by the 5\' rule of x1=beg x2=end'
    genomeListByChr[chr].append((x1, x2, s, gNames))
if unofficialChromRemoved:
    print >> sys.stderr, 'At least one unofficial chromosome (int(chr) > 100 or discarded name in myLightGenome.contigType()) has been removed'

rankOfChrHasChanged = False
if arguments['orderChromosomesBy'] == 'decreasingNbOfGenes':
    if not myTools.isSorted(genomeListByChr.items(), key=lambda x: len(x[1]), stricly=False, increasingOrder=False):
        genomeListByChr = collections.OrderedDict(sorted(genomeListByChr.items(), key=lambda x: len(x[1]), reverse=True))
        rankOfChrHasChanged = True
elif arguments['orderChromosomesBy'] == 'names':
    if not myTools.isSorted(genomeListByChr.items(), key=lambda x: myTools.keyNaturalSort(x[0])):
        genomeListByChr = collections.OrderedDict(sorted(genomeListByChr.items(), key=lambda x: myTools.keyNaturalSort(x[0]), reverse=False))
        rankOfChrHasChanged = True
if rankOfChrHasChanged:
    print >> sys.stderr, 'The rank of at least one chromosome has changed while sorting chrNames using the length of chromosomes'

print >> sys.stderr,  genomeListByChr.keys()

lengthsOfGenes = []
for chr, chrom in genomeListByChr.iteritems():
    for (idxGA, (begA, endA, sA, gNamesA)) in enumerate(chrom):
        lengthsOfGenes.append(abs(begA - endA))

print >> sys.stderr, "total original nb of genes = %s" % len(lengthsOfGenes)
print >> sys.stderr, "minimal length = %s" % min(lengthsOfGenes)
print >> sys.stderr, "maximal length = %s" % max(lengthsOfGenes)
print >> sys.stderr, "mean length = %.2f" % (float(sum(lengthsOfGenes)) / len(lengthsOfGenes))
print >> sys.stderr, "#genes with a length < 1kb  : %s" % len([gene for gene in lengthsOfGenes if gene < 1000])
print >> sys.stderr, "#genes with a length < 10kb : %s" % len([gene for gene in lengthsOfGenes if gene < 10000])
print >> sys.stderr, "#genes with a length < 100kb: %s" % len([gene for gene in lengthsOfGenes if gene < 100000])

if arguments['mode'] == 'distribOnChr':
    # plot the density along a given chromosome
    chr = 'X'
    geneCoordinates = []
    # import random
    # for _ in range(2000):
    #     geneCoordinates.append(random.random()*1000)
    for (idxGA, (begA, endA, sA, gNamesA)) in enumerate(genomeListByChr[chr]):
        geneCoordinates.append(begA)
    plt.figure()
    formatter = FuncFormatter(formaterBasePairs(unit='kb'))
    plt.gca().xaxis.set_major_formatter(formatter)
    plt.hist(geneCoordinates, bins=500, color='k')
    plt.title('Distribution of genes along the chromosome %s' % chr)
    plt.xlabel('coordinates of the genes')
    plt.ylabel('number of genes')
    plt.show(block=True)

if arguments['mode'] == 'longestIntergene':
    longestIntergenePerChromosome = collections.OrderedDict()
    for chr, chrom in genomeListByChr.iteritems():
        longestIntergene = 0
        for ((begA, endA, sA, gNamesA), (begB, endB, sB, gNamesB)) in myTools.myIterator.slidingTuple(chrom, 2):
            assert begA <= begB
            if sA == +1:
                assert begA < endA
            elif sA == -1:
                assert endA < begA
            if sB == +1:
                assert begB < endB
            elif sB == -1:
                assert endB < begB
            maxA = max(begA, endA)
            minB = min(begB, endB)
            tmp_longestIntergene = minB - maxA
            # if sA == +1:
            #     assert begA <= endA
            #     longestIntergene = minB - endA
            # elif sA == -1:
            #     assert endA <= begA
            #     longestIntergene = minB - begA
            # else:
            #     assert sA == None
            #     longestIntergene = minB - maxA
            if tmp_longestIntergene < 0:
                tmp_longestIntergene = 0
            longestIntergene = max(longestIntergene, tmp_longestIntergene)
        longestIntergenePerChromosome[chr] = longestIntergene
    print >> sys.stderr, 'the longest intergene length =%sbp' % max(longestIntergenePerChromosome.values())
    # print >> sys.stderr, ''
    # table = [['chromosome'],['longestIntergene (Mbp)']]
    # for chr, longestIntergene in longestIntergenePerChromosome.iteritems():
    #     table[0].append(chr)
    #     longestIntergene = longestIntergene / 1000000.0
    #     table[1].append('%.0f' % longestIntergene)
    # print >> sys.stderr, myTools.tableIntoLatex(table)
    plt.figure()
    formatter = FuncFormatter(formaterBasePairs(unit='Mb'))
    plt.gca().yaxis.set_major_formatter(formatter)
    X = [x for x in range(len(longestIntergenePerChromosome))]
    width = 0.8
    margin = (1.0-width)/2.0
    plt.bar([x + margin for x in X], longestIntergenePerChromosome.values(), width=width, color='k')
    plt.ticklabel_format()
    middle = (width + 2 * margin) / 2.0
    plt.xticks([x + middle for x in X], longestIntergenePerChromosome.keys())
    plt.title('Longest intergene length per chromosome of %s' % arguments['genome'])
    plt.xlabel('chromosome names')
    plt.ylabel('longest intergene length')
    plt.show()

if arguments['mode'] == 'minGeneLength':
    minLengthsOfGenesPerChr = {}
    for chr, chrom in genomeListByChr.iteritems():
        minLengthsOfGenesPerChr[chr] = sys.maxint
        for (beg, end, s, gNames) in chrom:
            minLengthsOfGenesPerChr[chr] = min(minLengthsOfGenesPerChr[chr], abs(beg - end))
    plt.figure()
    formatter = FuncFormatter(formaterBasePairs(unit='bp'))
    plt.gca().yaxis.set_major_formatter(formatter)
    X = [x for x in range(len(minLengthsOfGenesPerChr))]
    width = 0.8
    margin = (1.0-width)/2.0
    plt.bar([x + margin for x in X], minLengthsOfGenesPerChr.values(), width=width, color='k')
    plt.ticklabel_format()
    middle = (width + 2 * margin) / 2.0
    plt.xticks([x + middle for x in X], minLengthsOfGenesPerChr.keys())
    plt.title('Minimum gene length per chromosome of %s' % arguments['genome'])
    plt.xlabel('chromosome names')
    plt.ylabel('minimum gene length')
    plt.show()

if arguments['mode'] == 'correlationChromNbGenes':
    lengthsChrom = {}
    for chr, chrom in genomeListByChr.iteritems():
        currMaxLen = 0
        for (beg, end, s, gNames) in chrom:
            currMaxLen = max(max(beg, end), currMaxLen)
        lengthsChrom[chr] = currMaxLen
    X = []
    Y = []
    for chr in lengthsChrom:
        nbGenes = len(genomeListByChr[chr])
        X.append(lengthsChrom[chr])
        Y.append(nbGenes)
    slope, intercept, r_value, p_value, std_err = stats.linregress(X, Y)
    coeffOfPol1d = [slope, intercept]
    # print >> sys.stderr, slope, intercept, r_value, p_value, std_err
    #coeffOfPol1d = np.polyfit(X, Y, 1)
    #print >> sys.stderr, coeffOfPol1d
    ploy1d = np.poly1d(coeffOfPol1d)
    plt.figure()
    plt.plot(X, Y, linestyle='', marker='o')
    plt.plot([min(X), max(X)], [ploy1d(min(X)), ploy1d(max(X))], color='k')
    plt.gca().annotate('f(x bp)=%.2e x + %1.2f genes, R=%1.2f' % (coeffOfPol1d[0], coeffOfPol1d[1], r_value),(0.1,0.9), textcoords='axes fraction')
    #plt.text(0.1, 0.1, , textcoords='axes fraction')
    formatter = FuncFormatter(formaterBasePairs)
    plt.gca().xaxis.set_major_formatter(formatter)
    # http://stackoverflow.com/questions/5147112/matplotlib-how-to-put-individual-tags-for-a-scatter-plot
    # http://stackoverflow.com/questions/21140385/matplotlib-annotate-doesnt-work-on-log-scale
    for chr, x, y in itertools.izip(list(lengthsChrom.keys()), X ,Y):
        print chr
        xy = (x, y)
        print xy
        plt.gca().annotate(chr, xy=xy, xycoords='data', xytext=(-17, -7), textcoords='offset points')#, arrowprops=dict(arrowstyle="->"))
    plt.title('Correlation between the number of genes and the length of the chromosome')
    plt.xlabel('length of the chromosome')
    plt.ylabel('number of genes')
    plt.show(block=True)

if arguments['mode'] == 'geneLengths':

    # remove all genes that have a length > 100 000
    nbGenesBeforeRemoval = len(lengthsOfGenes)
    removedGenes = [gene for gene in lengthsOfGenes if gene >= 100000]
    lengthsOfGenes = [gene for gene in lengthsOfGenes if gene < 100000]
    nbGenesAfterRemoval = len(lengthsOfGenes)
    nbGenesRemoved = nbGenesBeforeRemoval - nbGenesAfterRemoval
    print >> sys.stderr,    \
        "For a better visualisation %s genes have been removed because they have a length > 100 kbp" % nbGenesRemoved

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    n, bins, patches = ax1.hist(lengthsOfGenes, 200, histtype='bar', color='grey')
    handle1, label1 = ax1.get_legend_handles_labels()
    ax1.yaxis.set_label_position("left")
    ax1.yaxis.tick_left()
    plt.ylabel('Absolute nb of genes')
    plt.ylim(0,max(n))
    ax2 = fig.add_subplot(111, sharex=ax1, frameon=False)
    # cumulative
    c_n, c_bins, c_patches = ax2.hist(lengthsOfGenes, 200, histtype='step', color='k', cumulative=True)
    handle2, label2 = ax2.get_legend_handles_labels()
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    plt.xlabel("Cumulated nb of genes")
    plt.ylim(0,max(c_n))
    #plt.legend(handle1 + handle2, label1 + label2, loc='upper left')
    #
    #plt.xticks(np.arange(min(0), max(lengthsOfGenes)+1, 200))
    def mjrFormatter(x, pos):
        if x > 1000:
            return "%.fkb" % (float(x) / 1000)
        elif x > 1000000:
            return "%.fMb" % (float(x) / 1000000)
        else:
            return x
    #def mjrFormatter_no_TeX(x, pos):
    #        return "2^{0}".format(x)
    #    )"}"))"}}}")
    ax = plt.gca()
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjrFormatter))
    #plt.locator_params(axis = 'x', nbins = 10)

    l = [float(e-s) for (s, e) in utils.myTools.myIterator.slidingTuple(bins)]
    bin_size = float(sum(l))/len(l)
    plt.title("Distribution of gene lengths (bp) in \n %s" % arguments["genome"])
    plt.xlabel("Gene lengths in bp, bins_size = %s" % bin_size)
    plt.show()
    #plt.savefig(sys.stdout, format='svg')



    # Figure with log scale
    print >> sys.stderr, "Be aware that the figure in log scale contains all the genes (even the %s genes previously removed )" % len(removedGenes)
    lengthsOfGenes = lengthsOfGenes + removedGenes

    fig2 = plt.figure()
    ax1 = fig2.add_subplot(111)
    n, bins, patches = ax1.hist(np.log10(lengthsOfGenes), 100, histtype='bar', color='grey')
    handle1, label1 = ax1.get_legend_handles_labels()
    ax1.yaxis.set_label_position("left")
    ax1.yaxis.tick_left()
    #ax1.set_xscale('log')
    plt.ylabel('Absolute nb of genes')
    plt.ylim(0,max(n))
    ax2 = fig2.add_subplot(111, sharex=ax1, frameon=False)
    # cumulative
    c_n, c_bins, c_patches = ax2.hist(np.log10(lengthsOfGenes), 200, histtype='step', color='black', cumulative=True)
    handle2, label2 = ax2.get_legend_handles_labels()
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    #ax2.set_xscale('log')
    plt.xlabel("Cumulated nb of genes")
    plt.ylim(0,max(c_n))
    #plt.legend(handle1 + handle2, label1 + label2, loc='upper left')

    def mjrFormatter(x, pos):
        return "$10^{%s}$" % x
    #plt.xticks(np.arange(min(0), max(lengthsOfGenes)+1, 200))
    ax = plt.gca()
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(mjrFormatter))
    #plt.locator_params(axis = 'x', nbins = 10)

    l = [float(e-s) for (s, e) in utils.myTools.myIterator.slidingTuple(bins)]
    bin_size = float(sum(l))/len(l)
    plt.title("Distribution of gene lengths (bp) in \n %s" % arguments["genome"])
    plt.xlabel("Gene lengths in bp, bins_size = $10^{%.3f}$" % bin_size)
    plt.xlim(xmin=1)
    plt.show()
    #plt.savefig(sys.stdout, format='svg')

if arguments['mode'] == 'distribOnChrs':
    lengthsChrom = {}
    for chr, chrom in genomeListByChr.iteritems():
        currMaxLen = 0
        for (beg, end, s, gNames) in chrom:
            currMaxLen = max(max(beg, end), currMaxLen)
        lengthsChrom[chr] = currMaxLen
    densityChrom = collections.OrderedDict()
    for chr in sorted(list(lengthsChrom.keys()), key=myTools.keyNaturalSort):
        densityChrom[chr] = float(len(genomeListByChr[chr])) / lengthsChrom[chr]
        print >> sys.stderr, 'density of genes in chrom %s = %s' %  (chr, densityChrom[chr])

    # plot the densities in all chromosomes
    def formaterDensity(unit='Mb'):
        'The two args are the value and tick position'
        def f (x, pos):
            if x == 0:
                res = 0
            else:
                if unit == 'kb':
                    res = '%0.0f' % (x*1000)
                elif unit == 'Mb':
                    res = '%0.0f' % (x*1000000)
                else:
                    res = '$%0.0f' % (x)
            return res
        return f
    plt.figure()
    unit = 'Mb'
    formatter = FuncFormatter(formaterDensity(unit=unit))
    plt.gca().yaxis.set_major_formatter(formatter)
    X = [x for x in range(len(densityChrom))]
    width = 0.8
    margin = (1.0-width)/2.0
    plt.bar([x + margin for x in X], densityChrom.values(), width=width, color='k')
    plt.ticklabel_format()
    middle = (width + 2 * margin) / 2.0
    plt.xticks([x + middle for x in X], densityChrom.keys())
    plt.title('Density of genes in %s' % arguments['genome'])
    plt.xlabel('chromosome names')
    plt.ylabel('density of genes (#genes/length of the chr in %sp)' % unit)
    plt.show(block=True)


if arguments['mode'] == 'overlap':
    nbOverlappingGenes = 0
    overlapPerGene = []
    for chr, chrom in genomeListByChr.iteritems():
        for (idxGA, (begA, endA, sA, gNamesA)) in enumerate(chrom):
            overlapA = 0
            for (begB, endB, sB, gNamesB) in chrom[0:idxGA] + chrom[idxGA+1:]:
                overlapAB = 0
                assert begA != endA
                assert begB != endB
                if endA < begA:
                    (begA, endA) = (endA, begA)
                if endB < begB:
                    (begB, endB) = (endB, begB)
                assert begA < endA and begB < endB
                assert (begA, endA) != (begB, endB)
                if begB < begA:
                    (begA, endA, sA, gNamesA), (begB, endB, sB, gNamesB) = (begB, endB, sB, gNamesB), (begA, endA, sA, gNamesA)
                assert begA <= begB
                if endA < begB:
                    # no overlap
                    pass
                else:
                    if endB <= endA:
                        # B is included in A
                        overlapAB = (endB - begB)
                        assert overlapAB > 0
                    else:
                        overlapAB = endA - begB + 1
                        assert overlapAB > 0
                overlapA += overlapAB
            overlapPerGene.append(overlapA)
    meanLengthOfGene = float(sum(lengthsOfGenes)) / len(lengthsOfGenes)
    assert len(lengthsOfGenes) == len(overlapPerGene)
    meanOverlapPerGene = float(sum(overlapPerGene)) / len(overlapPerGene)
    print >> sys.stderr, 'mean length of gene = %s' %  meanLengthOfGene
    ratioOverlap = float(meanOverlapPerGene) / meanLengthOfGene
    print >> sys.stderr, 'mean overlap per gene is %s (%s%% of mean gene length)' % (meanOverlapPerGene, 100 * ratioOverlap)
