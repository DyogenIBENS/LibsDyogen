# LibsDyogen version 1.0 (6/11/2015)
# python 2.7
# Copyright © 2015 IBENS/Dyogen : Nga THI THUY NGUYEN,  Joseph LUCAS and Hugues ROEST CROLLIUS
# mail : jlucas@ens.fr
# Licences GLP v3 and CeCILL v2

import itertools
import collections
import sys

import myTools
import myDiags
import myProbas
import myLightGenomes

def extractDiagsInPairCompChr(g1,g2,consistentSwDType,distanceMetric):
	cdef float c11, c21, c111, c211, c1N, c2N
	res = collections.defaultdict(lambda: collections.defaultdict(dict))
	diagsInPairComp = myTools.Dict2d(list)
	N12s = myTools.Dict2d(int)

	(stb1, stb1_g, _) = statsTbOrientation(g1, loc=False)
	(stb2, stb2_g, locG2) = statsTbOrientation(g2, loc=True)

	totalNbComps = len(g1.keys()) * len(g2.keys())
	progressBar = myTools.ProgressBar(totalNbComps)
	currCompNb = 0
	for (i, (c1,c2)) in enumerate(itertools.product(g1.keys(), g2.keys())):
		c11 = stb1[c1][+1]
		c21 = stb2[c2][+1]
		c111 = stb1[c1][-1]
		c211 = stb2[c2][-1]
		c1N = stb1[c1][None]
		c2N = stb2[c2][None]

		res[c1][c2][+1] = c11*c21 + c111*c211
		res[c1][c2][-1] = c11*c211 + c111*c21
		res[c1][c2][None] = c11*c2N + c1N*c21 + c111*c2N + c1N*c211 + c1N*c2N

		#####
		# Find list of diags
		#####
		listOfDiags = []
		nbHomo = 0

		if not locG2[c2]:
			N12s[c1][c2] = 0
			continue

		M = homologyMatrix(g1[c1], locG2[c2])
		nbHomo = sum([len(M[i]) for i in M])

		la = []
		l1 = []
		l2 = []
		diagType = None
		i1_old = None

		for (i1,(f,_)) in enumerate(g1[c1]):
			if f != None and f in locG2[c2]:
				i1_old = i1
				for (i2,_) in locG2[c2][f]:
					i1 = i1_old
					while i1 in M and i2 in M[i1]:
						f = g1[c1][i1][0]
						if len(la) == 0:
							diagType = myDiags.findDiagType(i1,i2,M,consistentSwDType)
						if g1[c1][i1][1] != None:
							ancestralStrand = g1[c1][i1][1]
						elif diagType == '/' and g2[c2][i2][1] != None:
							ancestralStrand = g2[c2][i2][1]
						elif diagType == '\\' and g2[c2][i2][1] != None:
							ancestralStrand = -g2[c2][i2][1]
						else:
							ancestralStrand = None
						la.append((f,ancestralStrand,len(la)+1))
						l1.append(i1)
						l2.append(i2)

						del M[i1][i2]
						if len(M[i1].keys()) == 0:
							del M[i1]
						if diagType == '/' and i1+1 in M and i2+1 in M[i1+1] and ((M[i1+1][i2+1] in [+1,None]) if consistentSwDType else True):
							i1 = i1+1
							i2 = i2+1
						elif diagType == '\\' and i1+1 in M and i2-1 in M[i1+1] and ((M[i1+1][i2-1] in [-1,None]) if consistentSwDType else True):
							i1 = i1+1
							i2 = i2-1
						else:
							listOfDiags.append(myDiags.Diagonal(diagType,l1,l2,la))
							l1=[]
							l2=[]
							la=[]
							diagType=None
							break
		if len(listOfDiags) > 0:
			(N12s[c1][c2], diagsInPairComp[c1][c2]) = (nbHomo, listOfDiags)
		currCompNb += 1
		progressBar.printProgressIn(sys.stderr, currCompNb)

	N12s_g = sum([nbH for nbH in N12s.values2d()])
	p_hpSign_g = {}
	#cdef float c11, c21, c111, c211, c1N, c2N
	c11 = stb1_g[+1]
	c21 = stb2_g[+1]
	c111 = stb1_g[-1]
	c211 = stb2_g[-1]
	c1N = stb1_g[None]
	c2N = stb2_g[None]
	p_hpSign_g[+1] = c11*c21 + c111*c211
	p_hpSign_g[-1] = c11*c211 + c111*c21
	p_hpSign_g[None] = c1N*(c21 + c211) + c2N*(c11+c111) + c1N*c2N

	return (res, p_hpSign_g, N12s, N12s_g, diagsInPairComp)


def statsTbOrientation(genome_tb, loc=False):
	p_tbO = {}
	locG = {}

	if loc:
		locG = collections.defaultdict(lambda: collections.defaultdict(list))
	cdef int nbTb_plus, nbTb_minus, nbTb_None, lenc, nbTb_plus_g, nbTb_minus_g, nbTb_None_g

	nbTb_plus_g, nbTb_minus_g, nbTb_None_g = 0, 0, 0
	for c in genome_tb:
		lenc = len(genome_tb[c])
		nbTb_plus, nbTb_minus, nbTb_None = 0, 0, 0
		for (i, (f,s)) in enumerate(genome_tb[c]):
			if s == +1:
				nbTb_plus += 1
			elif s == -1:
				nbTb_minus += 1
			elif s == None:
				nbTb_None += 1
			else:
				raise ValueError()
			if loc:
				if f != None:
					locG[c][f].append((i,s))

		p_tbO[c] = [float(v)/lenc for v in (nbTb_plus, nbTb_minus, nbTb_None)]
		p_tbO[c] = dict( ((+1,p_tbO[c][0]),(-1,p_tbO[c][1]),(None,p_tbO[c][2])))
		nbTb_plus_g += nbTb_plus
		nbTb_minus_g += nbTb_minus
		nbTb_None_g += nbTb_None

	nbTbs_g = sum([len(genome_tb[c]) for c in genome_tb])
	p_tbO_g = [float(v)/nbTbs_g for v in (nbTb_plus_g, nbTb_minus_g, nbTb_None_g)]
	p_tbO_g = dict( ((+1,p_tbO_g[0]),(-1,p_tbO_g[1]), (None,p_tbO_g[2]) ))

	return (p_tbO, p_tbO_g, locG)

def homologyMatrix(gc1, locG2):
	M = collections.defaultdict(lambda: collections.defaultdict(int))
	for (i1, (f,s1)) in enumerate(gc1):
		if f!=None and f in locG2:
			for (i2,s2) in locG2[f]:
				M[i1][i2] = strandProduct(s1,s2)
	return M

def strandProduct(sa, sb):
	if sa != None and sb != None:
		return sa*sb
	return None

