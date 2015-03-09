# -*- coding: utf-8 -*-
# LibsDyogen
# python 2.7
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

import sys
import copy
import collections
import myFile
import myGenomes
import myTools

# TODO : unoriented genes
# TODO : unoriented adjacencies

# n = name (ENSG00001 for a human gene)
# s = strand (+1 ou -1)
class OGene(collections.namedtuple("Gene", ['n', 's'])):
    # need to use the __new__ magick since the namedtuple is an immutable class
    # (cf http://stackoverflow.com/questions/1565374/subclassing-python-tuple-with-multiple-init-arguments)
    def __new__(cls, *args):
        if len(args) == 2:
            n = args[0]
            s = args[1]
            return super(OGene, cls).__new__(cls, n, s)
        else:
            # args contain only one tuple
            assert len(args) == 1
            return super(OGene, cls).__new__(cls, args[0][0], args[0][1])
    def reversed(cls):
        try:
            return super(OGene, cls).__new__(type(cls), cls.n, - cls.s)
        except TypeError:
            assert cls.s == None
            return super(OGene, cls).__new__(type(cls), cls.n, None)
    def __repr__(self):
        #return self.__class__.__name__ + tuple.__repr__(self)
        return 'G' + tuple.__repr__(self)
# fn = family name (it may be the ancestral gene name ENSGT0001.a.b.a for
# instance)
# ns = names of homologs, a tuple
Family = collections.namedtuple("Family", ['fn', 'dns'])

# Gene Position
# c = chromosome
# idx = index of the gene in that chromosome
GeneP = collections.namedtuple("GeneP", ['c', 'idx'])

# It might be more interesting for certain applications to only use a classical
# collections.defaultdict instead of the less specialised class
# myTools.DefaultOrderedDict.
# However, myTools.DefaultOrderedDict has been a polyvalent class for all
# usages. It might be used to shuffle chromosomes (since it is ordered) and
# adding a chromosome is simple since the underlying list object is automatically
# built.
class LightGenome(myTools.DefaultOrderedDict):

    def __init__(self, *args, **kwargs):
        self.name = None
        myTools.DefaultOrderedDict.__init__(self, default_factory=list)
        self.withDict = kwargs.get("withDict", False)
        if self.withDict:
            self.g2p = {}
        if len(args) == 0:
            return
        else:
            assert len(args) == 1, args
            arg = args[0]

        if isinstance(arg, str):
            fileName = arg
            self.name = fileName
            print >> sys.stderr, "Loading LightGenome from", fileName,
            # FIXME use myFile.firstLineBuffer to choose which format is in
            # input.
            # choice of the loading function
            flb = myFile.firstLineBuffer(myFile.openFile(fileName, 'r'))
            c = flb.firstLine.split("\t")
            if len(c) == 6:
                print >> sys.stderr, "(c, beg, end, s,  gName, transcriptName) -> (c, s, gName)",
                # c, beg, end, s,  gName, transcriptName
                reader = myFile.myTSV.readTabular(fileName, [str, int, int, int, str, str])
                reader = ((c, strand, gName) for (c, beg, end, strand, gName, tName) in reader)
            elif len(c) == 3:
                print >> sys.stderr, "(c, s, gName)",
                # c, s, gName
                reader = myFile.myTSV.readTabular(fileName, [str, int, str])
            elif len(c) == 5:
                print >> sys.stderr, "(c, beg, end, s,  gName) -> (c, s, gName)",
                # c, beg, end, s,  gName
                reader = myFile.myTSV.readTabular(fileName, [str, int, int, int, str])
                reader = ((c, strand, gName) for (c, beg, end, strand, gName) in reader)
            else:
                raise ValueError("%s file is badly formatted" % fileName)
            print >> sys.stderr, "...",
            # FIXME do not need beg, end and tName
            # c = chromosome name
            # beg = coordinate in nucleotides of the beginning of
            # transcription of the shortest transcript
            # end = coordinate in nucleotides of the ending of
            # transcription of the shortest transcript
            # gName = gene name
            # tName = transcript name
            idx = -1
            c_old = None
            for (c, strand, gName) in reader:
                idx = (idx + 1) if c == c_old else 0
                self[c].append(OGene(gName, strand))
                if self.withDict:
                    # dict 'G'ene to (pronounced '2') 'P'osition
                    self.g2p[gName] = GeneP(c, idx)
                    c_old = c
            print >> sys.stderr, 'OK'
        elif isinstance(arg, myGenomes.Genome):
            genome = arg
            self.name = genome.name
            for c in genome.lstGenes.keys():
                for (idx, g) in enumerate(genome.lstGenes[c]):
                    self[str(c)].append(OGene(g.names[0], g.strand))
                    if self.withDict:
                        self.g2p[g.names[0]] = GeneP(str(c), idx)
        elif isinstance(arg, LightGenome):
            self.name = arg.name
            self.withDict = arg.withDict
            for c in arg:
                self[c] = [OGene(gene.n, gene.s) for gene in arg[c]]
            if self.withDict:
                self.g2p = dict((gn, GeneP(gp.c, gp.idx)) for (gn, gp) in arg.g2p.iteritems())
        elif isinstance(arg, dict):
            genome = arg
            for c in genome:
                for (idx, (gName, strand)) in enumerate(genome[c]):
                    self[c].append(OGene(gName, strand))
                    if self.withDict:
                        # dict 'G'ene to (pronounced '2') 'P'osition
                        self.g2p[gName] = GeneP(c, idx)
        else:
            raise ValueError('Constructor needs a file')

    def computeDictG2P(self):
        self.withDict = True
        self.g2p = {}
        for c in self:
            for (idx, g) in enumerate(self[c]):
                self.g2p[g.n] = GeneP(c, idx)

    # in certain rare utilisations there might be several positions for the same gene name
    def computeDictG2Ps(self):
        # dict gene name to position's'
        self.g2ps = collections.defaultdict(list)
        for c in self:
            for (idx, g) in enumerate(self[c]):
                self.g2ps[g.n].append(GeneP(c, idx))

    def getPosition(self, name, default=None):
        try:
            # if this raise an error, g2p may need to be initialised
            return self.g2p.get(name, default)
        except:
            raise ValueError("self.g2p needs to be initialised")

    def __repr__(self):
        res = []
        res.append("LightGenome(")
        # id(self) is the uniq python identifier of the current self object
        res.append("id=%s" % id(self))
        if self.name is not None:
            res.append("name=%s" % self.name)
        res.append("withDict=%s" % self.withDict)
        for c in self:
            for gene in self[c]:
                res.append("\t%s\t%s" % (c, gene))
        res.append(")")
        return "\n".join(res)

    def __copy__(self):
        newLightGenome = LightGenome(withDict=self.withDict)
        newLightGenome.name = self.name
        newLightGenome.withDict = self.withDict
        for c in self.keys():
            newLightGenome[c] = self[c]
        if self.withDict:
            newLightGenome.g2p = self.g2p
        return newLightGenome

    def __deepcopy__(self, memo):
        # inspired from http://pymotw.com/2/copy/
        objectAlreadyBuilt = memo.get(id(self), None)
        if objectAlreadyBuilt is not None:
            # already copied
            return objectAlreadyBuilt
        newLightGenome = LightGenome(withDict=self.withDict)
        newLightGenome.name = copy.deepcopy(self.name, {})
        for c in self.keys():
            for gene in self[c]:
                newLightGenome[c].append(OGene(gene.n, gene.s))
        if self.withDict:
            newLightGenome.g2p = copy.deepcopy(self.g2p, {})
        memo[id(self)] = newLightGenome
        return newLightGenome

    def printIn(self, stream, format='Ensembl'):
        if format == 'Ensembl':
            for c, chrom in self.iteritems():
                for (i, g) in enumerate(chrom):
                    print >> stream, myFile.myTSV.printLine([c, i, i+1, g.s, g.n])

    # Returns the gene names around the intergene
    def getIntervStr(self, c, x):
        # DEBUG assertion
        assert 0<= x and x <= len(self[c]), \
            "x=%s and len(self[c])=%s" % (x, len(self[c]))
        if x == 0:
            return "End-" + str(self[c][x][0])
        elif x == len(self[c]):
            return str(self[c][x-1][0]) + "-End"
        else:
            return "%s-%s" % (self[c][x-1][0], self[c][x][0])


class Families(list):

    def __init__(self, *args):
        list.__init__(self)
        self.fidMax = 0
        self.g2fid = {}
        if len(args) == 0:
            # null constructor
            return
        if len(args) == 1 and isinstance(args[0], str):
            fileName = args[0]
            self.name = fileName
            print >> sys.stderr, "Loading Families from", fileName, "...",
            # FIXME use myFile.firstLineBuffer to choose which format is in
            # input.
            # A more synthetic format would have only 3 columns:
            # c, s and gName
            file = myFile.openFile(fileName, 'r')
            for l in file:
                names = l.replace('\n', '').split(' ')
                # family name
                fn = names[0]
                # modern names
                # FIXME: dns might be only a set
                dns = names[1:]
                self.append(Family(fn, dns))
                for n in [fn] + dns:
                    fID = self.fidMax
                    # Each (fID + 1) corresponds to the line number in the output file
                    # of families obtained with self.printIn(file)
                    # Be carefull, a fID number is equal to the line number - 1
                    # (if line number begins from 1)
                    self.g2fid[n] = fID
                self.fidMax += 1
            file.close()
            assert self.fidMax == len(self)
        else:
            raise ValueError('Constructor needs a file')
        print >> sys.stderr, 'OK'

    def addFamily(self, family):
        assert isinstance(family, Family)
        self.append(family)
        for n in [family.fn] + family.dns:
            self.g2fid[n] = self.fidMax
        self.fidMax += 1

    # gene name is either a family name or homologs names
    def getFamID(self, name, default=None):
        return self.g2fid.get(name, default)

    def getFamilyByID(self, famID, default=None):
        try:
            return self[famID]
        except:
            return default

    def getFamilyByName(self, name, default=None):
        famID = self.getFamID(name)
        return self.getFamilyByID(famID, default=default)

    def getFamNameByName(self, name, default=None):
        famID = self.getFamID(name)
        if famID is not None:
            return self.getFamilyByID(famID, default=default).fn
        else:
            return default

    def getFamNameByID(self, famID, default=None):
        try:
            return self[famID].fn
        except:
            return default

    def printIn(self, stream):
        for (famID, family) in enumerate(self):
            line = [str(family.fn)]
            line.extend([str(name) for name in family.dns])
            print >> stream, myFile.myTSV.printLine([" ".join(line)])

    def __repr__(self):
        res = []
        res.append('Family:')
        for family in self:
            res.append(' '.join([family.fn] + family.dns))
        return '\n'.join(res)
