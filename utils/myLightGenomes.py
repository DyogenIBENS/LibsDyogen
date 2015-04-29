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


def newChromName(genome):
    assert isinstance(genome, LightGenome)
    chromNames = set([int(c) for c in genome.keys()])
    maxChromName = max(chromNames)
    chromNamesToFillGaps = set(range(maxChromName + 1)) - chromNames
    if len(chromNamesToFillGaps) > 0:
        chosenChromName = str(chromNamesToFillGaps.pop())
    else:
        chosenChromName = str(maxChromName + 1)
    return chosenChromName

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
        # kwargs.get('name', default=None)
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
                print >> sys.stderr, "(c, beg, end, s, gName, transcriptName) -> (c, s, gName)",
                # c, beg, end, s,  gName, transcriptName
                reader = myFile.myTSV.readTabular(fileName, [str, int, int, int, str, str])
                reader = ((c, strand, gName) for (c, beg, end, strand, gName, tName) in reader)
            elif len(c) == 3:
                print >> sys.stderr, "(c, s, gName)",
                # c, s, gName
                reader = myFile.myTSV.readTabular(fileName, [str, int, str])
            elif len(c) == 5:
                print >> sys.stderr, "(c, beg, end, s, gName) -> (c, s, gName)",
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

    def getPosition(self, name, default=None):
        try:
            # if this raise an error, g2p may need to be initialised
            return self.g2p.get(name, default)
        except:
            raise ValueError("self.g2p needs to be initialised")

    # in certain rare utilisations there might be several positions for the same gene name
    def computeDictG2Ps(self):
        # dict gene name to position's'
        self.g2ps = collections.defaultdict(set)
        for c in self:
            for (idx, g) in enumerate(self[c]):
                self.g2ps[g.n].add(GeneP(c, idx))

    def getPositions(self, name, default=None):
        """ return a set of positions of genes with the gene.n = name
        """
        try:
            # if this raise an error, g2p may need to be initialised
            return self.g2ps.get(name, default)
        except:
            raise ValueError("self.g2ps needs to be initialised")

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

    # see also myMapping.remapFilterGeneContent(genome, removedNames, mOld=None)
    def removeGenes(self, setRemovedGeneNames):
        nbRemovedGenes = 0
        nbRemovedChrs = 0
        for (c, chrom) in self.iteritems():
            newChrom = []
            lenNewChrom = 0
            thereWasARemoval = False
            for ogene in chrom:
                if ogene.n in setRemovedGeneNames:
                    if self.withDict:
                        del self.g2p[ogene.n]
                        try:
                            del self.g2ps[ogene.n]
                        except:
                            pass
                    thereWasARemoval = True
                    nbRemovedGenes += 1
                else:
                    newChrom.append(ogene)
                    if self.withDict and thereWasARemoval:
                        # if there was no removal, no need to change the g2p, it stays at it was for these first genes
                        # of the chrom until a gene is removed
                        self.g2p[ogene.n] = GeneP(c, lenNewChrom)
                        try:
                            self.g2ps[ogene.n].add(GeneP(c, lenNewChrom))
                        except:
                            pass
                    lenNewChrom += 1
            if lenNewChrom > 0:
                assert lenNewChrom == len(newChrom)
                self[c] = newChrom
            else:
                del self[c]
                nbRemovedChrs += 1
        return (nbRemovedGenes, nbRemovedChrs)

    # TODO change name to getGeneNames (no ambiguity any more with a possible setter)
    def getGeneNames(self, asA=set, checkNoDuplicates=True):
        res = asA()
        for chrom in self.values():
            for gene in chrom:
                if checkNoDuplicates and gene.n in res:
                    raise ValueError("%s contains two times the same gene name %s" % (self.name, gene.n))
                if asA == set:
                    res.add(gene.n)
                else:
                    assert asA == list
                    res.append(gene.n)
        return res

    def getOwnedFamilyNames(self, families, asA=set):
        assert asA in [set, list]
        assert isinstance(families, Families)
        res = asA()
        for chrom in self.values():
            for gene in chrom:
                family = families.getFamilyByName(gene.n, default=None)
                if family:
                    if asA == set:
                        res.add(family.fn)
                    else:
                        assert asA == list
                        res.append(family.fn)
        return res

    def removeChrsStrictlySmallerThan(self, minChrLen):
        sCs = self.keys()
        nbRemovedChrs = 0
        for c in sCs:
            if len(self[c]) < minChrLen:
                if self.withDict:
                    for gene in self[c]:
                        del self.g2p[gene.n]
                        try:
                            del self.g2ps[gene.n]
                        except:
                            pass
                del self[c]
                nbRemovedChrs += 1
        return nbRemovedChrs

    def sort(self):
        """ Sort sbs by decreasing sizes """
        l = self.items()
        for c in self.keys():
            del self[c]
        for (c, chrom) in sorted(l, key=lambda x: len(x[1]), reverse=True):
            self[c] = chrom


# FIXME, it could also be easier to use a dict here and no family IDs, but IDs are processed faster than strings when
# FIXME managing rewritten genomes
# FIXME families should be organised this way:
# 1st element of the line: name of the family (or ancestral gene)
# next elements of the line: the names of the last descendants. Most of them are in modern species, but some may be in
# intermediary ancestors
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
                # FIXME: dns cannot be a set, since the first element is the positional ortholog as much as possible
                dns = set(names[1:])
                self.append(Family(fn, dns))
                for n in {fn} | dns:
                    # Each (fID + 1) corresponds to the line number in the output file
                    # of families obtained with self.printIn(file)
                    # Be careful, a fID number is equal to the line number - 1
                    # (if line number begins from 1)
                    self.g2fid[n] = self.fidMax
                self.fidMax += 1
            file.close()
            assert self.fidMax == len(self)
        else:
            raise ValueError('Constructor needs a file')
        print >> sys.stderr, 'OK'

    def addFamily(self, family):
        assert isinstance(family, Family)
        selfFamId = self.getFamID(family.fn, default=None)
        if selfFamId is None:
            self.append(family)
            for n in {family.fn} | family.dns:
                self.g2fid[n] = self.fidMax
            self.fidMax += 1
        else:
            selfFam = self.getFamilyByID(selfFamId)
            assert family.fn == selfFam.fn
            self[selfFamId] = Family(selfFam.fn, selfFam.dns | family.dns)
            for n in family.dns:
                self.g2fid[n] = selfFamId

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
            line = [str(family.fn)] + [str(name) for name in sorted(family.dns)]
            print >> stream, myFile.myTSV.printLine([" ".join(line)])

    def __repr__(self):
        res = []
        res.append('Family:')
        for family in self:
            res.append(' '.join([family.fn] + sorted(family.dns)))
        return '\n'.join(res)

    def __eq__(self, other):
        return set((fn, dns) for (fn, dns) in self) == set((fn, dns) for (fn, dns) in other)

    def __ne__(self, other):
        return not self == other

    def __or__(self, other):
        """ self + other """
        res = Families()
        assert isinstance(other, Family)
        for family in self:
            res.addFamily(family)
        for family in other:
            self.addFamily(family)
        return res

    def __ior__(self, other):
        """ self += other """
        for family in other:
            self.addFamily(family)

def f_A0_A1_from_f_A0_Ds_and_f_A1_Ds(f_A0_D, f_A1_D):
    """
    Return a new Families with fn of A0 and dns of A1

    Warning: ancestor A0 must be older than A1 and A1 should be in the evolutive path from A0 to D

    :param f_A0_D: Families, fn of A0 and dns of D (Descendant)
    :param f_A1_D: Families, fn of A1 and dns of D (same Descendant)
    :return f_A0_A1: Families, fn of A0 and dns of A1
    """
    assert isinstance(f_A0_D, Families) and isinstance(f_A1_D, Families)
    f_A0_A1 = Families()
    tmp = collections.defaultdict(set)
    for i, (gn_A1, gns_D) in enumerate(f_A1_D):
        # FIXME, what if the gene is kept between A0 and A1 and lost after ? Use Xaft
        if len(gns_D) > 0:
            # If the ancestral gene has descendant genes, preferentially genes in extant species...
            gns_D_tmp = [gn for gn in gns_D if 'Xaft' not in gn]
            gns_D = gns_D_tmp if len(gns_D_tmp) > 0 else gns_D
            dn = gns_D[0]
            fn_A0 = f_A0_D.getFamilyByName(dn, default=None)
            # assert all(f_A0_D.getFamilyByName(dn, default=None) == fn_A0 for dn in gns_D)
            if fn_A0 is not None:
                tmp[fn_A0.fn].add(gn_A1)
    for (gn_A0, gns_A1) in tmp.iteritems():
        f_A0_A1.addFamily(Family(gn_A0, sorted(list(gns_A1))))
    return f_A0_A1