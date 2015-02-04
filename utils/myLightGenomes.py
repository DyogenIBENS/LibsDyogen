# -*- coding: utf-8 -*-
# LibsDyogen
# python 2.7
# Copyright © 2013 IBENS/Dyogen Joseph LUCAS, Matthieu MUFFATO and Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or jlucas@ens.fr
# This is free software, you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 (GPL v3) or later and the CeCiLL v2 license in France

import sys
import collections
import myFile
import myGenomes

# n = name (ENSG00001 for a human gene)
# s = strand (+1 ou -1)
Gene = collections.namedtuple("Gene", ['n', 's'])

# fn = family name (it may be the ancestral gene name ENSGT0001.a.b.a for
# instance)
# ns = names of homologs, a tuple
Family = collections.namedtuple("Family", ['fn', 'ns'])

# c = chromosome
# idx = index of the gene in that chromosome
GeneP = collections.namedtuple("GeneP", ['c', 'idx'])


class LightGenome(collections.defaultdict):

    def __init__(self, *args, **kwargs):
        self.name = None
        collections.defaultdict.__init__(self, list)
        if len(args) == 0:
            return
        else:
            assert len(args) == 1
            arg = args[0]
        withDict = kwargs.get("withDict", True)
        if withDict:
            self.g2p = collections.defaultdict(list)

        if isinstance(arg, str):
            fileName = arg
            self.name = fileName
            print >> sys.stderr, "Loading LightGenome from", fileName, "...",
            # FIXME use myFile.firstLineBuffer to choose which format is in
            # input.
            # A more synthetic format would have only 3 columns:
            # c, s and gName
            reader = myFile.myTSV.readTabular(arg, [str, int, int, int, str, str])
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
            for (c, beg, end, strand, gName, tName) in reader:
                idx = (idx + 1) if c == c_old else 0
                self[c].append(Gene(gName, strand))
                if withDict:
                    # dict 'G'ene to (pronounced '2') 'P'osition
                    self.g2p[gName] = GeneP(c, idx)
                    c_old = c
            print >> sys.stderr, 'OK'
        elif isinstance(arg, myGenomes.Genome):
            genome = arg
            self.name = genome.name
            for c in genome.lstGenes.keys():
                for (idx, g) in enumerate(genome.lstGenes[c]):
                    self[str(c)].append(Gene(g.names[0], g.strand))
                    if withDict:
                        self.g2p[g.names[0]] = GeneP(str(c), idx)
        elif isinstance(arg, LightGenome):
            self.name = arg.name
            collections.defaultdict.__init__(arg, list)
            for c in arg:
                self[c] = [Gene(gene.n, gene.s) for gene in arg[c]]
            if withDict:
                try:
                    self.g2p = arg.g2p
                except:
                    pass
        elif isinstance(arg, dict):
            genome = arg
            self = genome
            idx = -1
            c_old = None
            for c in genome:
                for (idx, (gName, strand)) in enumerate(genome[c]):
                    self[c].append(Gene(gName, strand))
                    if withDict:
                        # dict 'G'ene to (pronounced '2') 'P'osition
                        self.g2p[gName] = GeneP(c, idx)
        else:
            raise ValueError('Constructor needs a file')

    def getPosition(self, name):
        # if this raise an error, g2p may need to be initialised
        return self.g2p[name] if name in self.g2p else None

    def __repr__(self):
        res = []
        if self.name is not None:
            res.append("LightGenome %s" % self.name)
        for c in self:
            for gene in self[c]:
                res.append("%s\t%s" % (c, gene))
        return "\n".join(res)


class Families(list):

    def __init__(self, arg):
        if isinstance(arg, str):
            fileName = arg
            self.name = fileName
            print >> sys.stderr, "Loading Families from", fileName, "...",
            # FIXME use myFile.firstLineBuffer to choose which format is in
            # input.
            # A more synthetic format would have only 3 columns:
            # c, s and gName
            f = myFile.openFile(fileName, 'r')
            fiD = 0
            self.g2fid = {}
            for l in f:
                names = l.replace('\n', '').split(' ')
                # names[0] is the family name
                # names[:1] are the modern names
                self.append(Family(names[0], names[1:]))
                for name in names:
                    self.g2fid[name] = fiD
                fiD += 1
            f.close()
            assert fiD == len(self)
        else:
            raise ValueError('Constructor needs a file')
        print >> sys.stderr, 'OK'

    # gene name is either a family name or homologs names
    def getFamID(self, name):
        return self.g2fid[name] if name in self.g2fid else None

    def getAncGene(self, famID):
        return self[famID] if famID in self else None
