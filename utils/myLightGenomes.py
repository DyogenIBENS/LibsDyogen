# -*- coding: utf-8 -*-

# LibsDyogen version 1.0 (6/11/2015)
# python v2.7 at least is needed
# Copyright © 2015 IBENS/Dyogen : Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : jlucas@ens.fr
# Licences GLP v3 and CeCILL v2

import sys
import copy
import collections
import itertools
import re
import myFile
import myGenomes
import myTools
import enum


def _keyFuncNaturalSort((c, chrom)):
    '''
    alist.sort(key=natural_keys) sorts in human wished order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    res = [myTools.atoi(c) for c in re.split('([0-9]+)', str(c))]
    return res

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
    chromNames = set([])
    for c in genome.keys():
        try:
            c = int(c)
            chromNames.add(c)
        except:
            pass
    if len(chromNames) == 0:
        chosenChromName = 1
    else:
        maxChromName = max(chromNames)
        chromNamesToFillGaps = set(range(1, maxChromName + 1)) - chromNames
        if len(chromNamesToFillGaps) > 0:
            chosenChromName = str(chromNamesToFillGaps.pop())
        else:
            chosenChromName = str(maxChromName + 1)
    return chosenChromName

ContigType = enum.Enum('Chromosome', 'Mitochondrial', 'Scaffold', 'None', 'Random')

# return the type of contig depending on its name
# FIXME, add sexual chromosomes
def contigType(chrName):
    if chrName in [None, "Un_random", "UNKN", "Un"]:
        # chromosome named "Un" in Monodelphis.domestica corresponds to a scaffold that concatenates all unassembled
        # fragments
        return ContigType.None
    try:
        x = int(chrName)
        if x < 100:
            return ContigType.Chromosome
        else:
            return ContigType.Scaffold
    except:
        chrNameLow = chrName.lower()
        if "rand" in chrNameLow:
            # FIXME : what is a random contig ? It seems that a random contig is always associated to a chromosome number.
            # Thus we know that this random chromosome is on a specific chromosome number but we do not know where.
            return ContigType.Random
        # mitochondrion in Drosophila melanogaster
        for x in ['mt', 'mitochondrion']:
            if x in chrNameLow:
                return ContigType.Mitochondrial
        # ki in Homo.sapiens data81
        # jh in mus.musculus, Canis.lupus.familiaris
        # aaex in Canis.lupus.familiaris
        # aadn in Gallus.gallus
        keys = ["cont", "scaff", "ultra", "reftig", "_", "un", "gl", "ki", "jh", "aaex", "aadn"]
        for x in keys:
            if x in chrNameLow:
                return ContigType.Scaffold
        if (chrName in ["U", "E64", "2-micron"]) or chrName.endswith("Het"):
            return ContigType.Scaffold
        else:
            # sex chromosomes of the chicken as (Z, W) or sex chromosomes of the human (X, Y)
            return ContigType.Chromosome

# TODO: implement all the basic classes of myLigthGenomes.LightGenome into myLightGenome.Chromosome for more modularity
# then in myLigthGenomes.LightGenome, call the classes of myLigthGenomes.Chromosome whenever it is possible
class Chromosome(list):
    def __init__(self, seq=(), name=None):
        if name is not None:
            self.name = name
        list.__init__(self, seq)

    # def __getitem__(self, key): # Methode appellee par obj[3] <=> obj.__getitem__(3)
    #     return list.__getitem__(self, key)
    #

    # FIXME: is it optimal ?
    def __getslice__(self, i,j): # If user enters A[a:b] <=> A.__getslice__(a,b) # Impossible to restrict to __getitem__ ...
        return Chromosome(self.__getitem__(slice(i,j,1)))

    def __getitem__(self, item):
        result = list.__getitem__(self, item)
        return result
        # try:
        #
        # except TypeError:
        #     return result

    def __add__(self, rhs):
        return Chromosome(list.__add__(self, rhs))

    def getGeneNamesNotInFamilies(self, families, asA=set):
        assert asA in [set, list]
        assert isinstance(families, Families)
        res = asA()
        for gene in self:
            family = families.getFamilyByName(gene.n, default=None)
            if not family:
                if asA == set:
                    res.add(gene.n)
                else:
                    assert asA == list
                    res.append(gene.n)
        return res

    def getGeneNamesInFamilies(self, families, asA=set):
        assert asA in [set, list]
        assert isinstance(families, Families)
        res = asA()
        for gene in self:
            family = families.getFamilyByName(gene.n, default=None)
            if family is not None:
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
        for gene in self:
            family = families.getFamilyByName(gene.n, default=None)
            if family:
                if asA == set:
                    res.add(family.fn)
                else:
                    assert asA == list
                    res.append(family.fn)
        return res

# It might be more interesting for certain applications to only use a classical
# collections.defaultdict instead of the less specialised class
# myTools.DefaultOrderedDict.
# However, myTools.DefaultOrderedDict has been a polyvalent class for all
# usages. It might be used to shuffle chromosomes (since it is ordered) and
# adding a chromosome is simple since the underlying list object is automatically
# built.
class LightGenome(myTools.DefaultOrderedDict):

    @staticmethod
    def readerDependingOnFile(fileName):
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
            tmpReader = myFile.myTSV.readTabular(fileName, [str, int, int, int, str])
            # check, with the first line, if there are several gene names (the format genome of Matthieu contains several gene names)
            (c, beg, end, strand, gNames) = tmpReader.next()
            severalNames = True if len(gNames.split(' ')) > 0 else False
            reader = itertools.chain([(c, beg, end, strand, gNames)], tmpReader)
            if severalNames:
                # if gNames contains more than one gene name, only take the first gene name
                reader = ((c, strand, gNames.split(' ')[0]) for (c, beg, end, strand, gNames) in reader)
            else:
                reader = ((c, strand, gName) for (c, beg, end, strand, gName) in reader)
        else:
            raise ValueError("%s file is badly formatted" % fileName)
        return reader

    def __init__(self, *args, **kwargs):
        self.name = None
        # this dict contains the sets chromosomes per type of contig (cf class ContigType)
        self.chrSet = collections.defaultdict(set)
        # kwargs.get('name', default=None)
        myTools.DefaultOrderedDict.__init__(self, default_factory=Chromosome)
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
            # FIXME use myFile.firstLineBuffer to choose which format is in input.
            # choice of the loading function
            reader = self.readerDependingOnFile(fileName)
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
                self.chrSet[contigType(c)].add(c)
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
            self.chrSet = arg.chrSet
            for c in genome.lstGenes.keys():
                for (idx, g) in enumerate(genome.lstGenes[c]):
                    self[str(c)].append(OGene(g.names[0], g.strand))
                    if self.withDict:
                        self.g2p[g.names[0]] = GeneP(str(c), idx)
        elif isinstance(arg, LightGenome):
            self.name = arg.name
            self.chrSet = arg.chrSet
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
        except Exception as e:
            raise ValueError("maybe self.g2p needs to be initialised, otherwise %s" % e.message)

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

    def printIn(self, stream, format='Ensembl', begIsTranscriptionStart=True):
        if format == 'Ensembl':
            for c, chrom in self.iteritems():
                for (i, g) in enumerate(chrom):
                    if begIsTranscriptionStart:
                        if g.s in {+1, None}:
                            (beg, end) = (i, i+1)
                        else:
                            assert g.s == -1
                            (beg, end) =(i+1, i)
                    else:
                        (beg, end) = (i, i+1)
                    print >> stream, myFile.myTSV.printLine([c, beg, end, g.s, g.n])


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
    def getGeneNames(self, asA=set, checkNoDuplicates=False):
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

    def getOwnedFamilies(self, families):
        res = []
        for fn in self.getOwnedFamilyNames(families, asA=set):
            res.append(families.getFamilyByName(fn))
        return res

    def getGeneNamesInFamilies(self, families, asA=set):
        assert asA in [set, list]
        assert isinstance(families, Families)
        res = asA()
        for chrom in self.values():
            for gene in chrom:
                family = families.getFamilyByName(gene.n, default=None)
                if family:
                    if asA == set:
                        res.add(gene.n)
                    else:
                        assert asA == list
                        res.append(gene.n)
        return res

    def getGeneNamesNotInFamilies(self, families, asA=set):
        assert asA in [set, list]
        assert isinstance(families, Families)
        res = asA()
        for chrom in self.values():
            for gene in chrom:
                family = families.getFamilyByName(gene.n, default=None)
                if not family:
                    if asA == set:
                        res.add(gene.n)
                    else:
                        assert asA == list
                        res.append(gene.n)
        return res

    def removeChrsStrictlySmallerThan(self, minChrLen):
        sCs = self.keys()
        nbRemovedChrs = 0
        nbRemovedGenes = 0
        for c in sCs:
            if len(self[c]) < minChrLen:
                nbRemovedGenes += len(self[c])
                if self.withDict:
                    for gene in self[c]:
                        del self.g2p[gene.n]
                        try:
                            del self.g2ps[gene.n]
                        except:
                            pass
                del self[c]
                nbRemovedChrs += 1
        return (nbRemovedChrs, nbRemovedGenes)

    def removeUnofficialChromosomes(self):
        if len(self.chrSet.keys()) == 0:
            for c in self.keys():
                self.chrSet[contigType(c)].add(c)
        sCs = self.keys()
        nbRemovedChrs = 0
        nbRemovedGenes = 0
        for c in sCs:
            if c not in self.chrSet[ContigType.Chromosome]:
                nbRemovedGenes += len(self[c])
                if self.withDict:
                    for gene in self[c]:
                        del self.g2p[gene.n]
                        try:
                            del self.g2ps[gene.n]
                        except:
                            pass
                del self[c]
                nbRemovedChrs += 1
        return (nbRemovedChrs, nbRemovedGenes)

    def sort(self, byName=False):
        """ Sort chrs by decreasing sizes (default), or byName """
        l = self.items()
        for c in self.keys():
            del self[c]
        if byName is False:
            # sort by decreasing size
            for (c, chrom) in sorted(l, key=lambda x: len(x[1]), reverse=True):
                self[c] = chrom
        else:
            assert byName is True
            for (c, chrom) in sorted(l, key=_keyFuncNaturalSort):
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
        if isinstance(args, Families):
            families = args
            self.name = families.name
            for family in families:
                self.append(Family(copy.deepcopy(family.fn), copy.deepcopy(family.dns)))
            self.g2fid = copy.deepcopy(families.g2fid)
            self.fidMax = families.fidMax
            assert self.fidMax == len(self)
            # DEBUG assert
            assert set(self.g2fid.key()) == set().union(*[{fn} | dns for (fn, dns) in self])
        elif len(args) == 1 and isinstance(args[0], str):
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
        # "name is None" in the case of a chromosome rewritten with famNames that contain lineage specific genes without famNames
        assert isinstance(name, str) or isinstance(name, int) or name is None
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

    def genesAreInSameFamily(self, geneName1, geneName2):
        famName1 = self.getFamNameByName(geneName1, default=None)
        if famName1 is None:
            return False
        famName2 = self.getFamNameByName(geneName2)
        if famName2 is None:
            return False
        if famName1 == famName2:
            return True

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

# if __name__ == '__main__':
#     genome = LightGenome('/home/jlucas/Libs/PhylDiag/data/genesST.Homo.sapiens.list.bz2')
#     genome = LightGenome('/home/jlucas/Libs/PhylDiag/data/genesST.Mus.musculus.list.bz2')
#     genome = LightGenome('/home/jlucas/Libs/PhylDiag/data/genesST.Monodelphis.domestica.list.bz2')
#     genome = LightGenome('/home/jlucas/Libs/PhylDiag/data/genesST.Canis.lupus.familiaris.list.bz2')
#     genome = LightGenome('/home/jlucas/Libs/PhylDiag/data/genesST.Gallus.gallus.list.bz2')
#     genome = LightGenome('/home/jlucas/Libs/PhylDiag/data/genesST.Drosophila.melanogaster.list.bz2')
#     genome = LightGenome('/home/jlucas/Libs/PhylDiag/data/genesST.Drosophila.melanogaster.list.bz2')
#     print >> sys.stderr, genome.chrSet
#     genome.removeUnofficialChromosomes()
#     print >> sys.stderr, genome.keys()