# coding=utf-8
from __future__ import division
import dendropy
import datetime
import time
import sys

#reload(sys)
#import xlrd
#import types
import sys
import os
#import shutil
#import string
#from random import choice
#import datetime
#import time
#from ete2 import Tree, NodeStyle, TreeStyle, phyloxml
#from Bio import Phylo
#from Bio.Phylo.PhyloXML import Phylogeny
#import pylab
#import sys
#import matplotlib
#from numpy.random import randn
#import matplotlib.pyplot as plt
#from matplotlib.ticker import FuncFormatter
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch
from reportlab.lib.colors import green, blue, yellow, red, maroon, purple, pink, gray, black, darkslategray
#import numpy as np
import math
#from datetime import datetime
#import random
#sys.setdefaultencoding('utf8')


def calMPD(pdc, internalNode):
	listLeaves = internalNode.leaf_nodes()
	numLeaves = len(listLeaves)
	sumD = 0
	count = 0
	for i in range(numLeaves):
		leafi = listLeaves[i]
		for j in range(i):
			leafj = listLeaves[j]
			dij = pdc(leafi.taxon, leafj.taxon)
			sumD = sumD + dij
			count = count + 1
	return sumD/count

def countMoreSubtypes(internalNode):
	listLeaves = internalNode.leaf_nodes()
	numLeaves = len(listLeaves)
	dictHtypeStrains = {}
	dictNtypeStrains = {}
	for i in range(numLeaves):
		leafi = listLeaves[i]
		gi = leafi.taxon.label
		subtype = gi.split(' ')[-1]
		# print(subtype)
		htype = subtype[:subtype.find('N')]
		ntype = subtype[subtype.find('N'):]
		genomeID = gi.split(' ')[0]
		
		if htype in dictHtypeStrains:
			dictHtypeStrains[htype].add(genomeID)
		else:
			dictHtypeStrains[htype] = set()
			dictHtypeStrains[htype].add(genomeID)
		
		if ntype in dictNtypeStrains:
			dictNtypeStrains[ntype].add(genomeID)
		else:
			dictNtypeStrains[ntype] = set()
			dictNtypeStrains[ntype].add(genomeID)
		
		

	return dictHtypeStrains,dictNtypeStrains

			
	

def get_median(data):
	data.sort()
	half = len(data) // 2
	return (data[half] + data[~half]) / 2
	

	
	
segs = ['HA','PB2','PB1','PA','NP','NA','MP','NS']
#segs = ['PA', 'NP']
#segs = ['NA', 'MP', 'NS']
segsHANA = ['HA','NA']
treeDir = 'iqtree'

prefixSeg = ''
suffixSeg = '.fasta.codon.fas.mafft.fas.correctSubtypes.fas.correctSubtypes2.fas.treefile'

dictHtypeMax = dict()
dictNtypeMax = dict()



fwrite = open('HANAcorrectSubtypes2_checkWrongSubtypes.txt', 'w')

for seg in segsHANA:
	print(seg)
	
	tree1 = dendropy.Tree.get(path=os.path.join(treeDir,prefixSeg+seg+suffixSeg), schema='newick')
	
	maxInternalEdgeLength = 0
	internalNodeToBeReRooted = None
	for internalNode in tree1.internal_nodes():
		if internalNode.edge_length != None and internalNode.edge_length > maxInternalEdgeLength :
			maxInternalEdgeLength = internalNode.edge_length
			internalNodeToBeReRooted = internalNode
	#sys.exit()
	treeRoot = internalNodeToBeReRooted
	tree1.reroot_at_node(internalNodeToBeReRooted)
	print(internalNodeToBeReRooted, 'is the new root.')
	
	tree1.ladderize(ascending=True)	
	internalNodes = tree1.internal_nodes()

	if seg == 'HA':
		for internalNode in internalNodes:
			numLeaves = len(internalNode.leaf_nodes())
			dictHtypeStrains, dictNtypeStrains = countMoreSubtypes(internalNode)
			# print(dictHtypeStrains)
			#print(numLeaves)
			for htype in dictHtypeStrains:
				# print(htype)
				htypeNum = len(dictHtypeStrains[htype])
				if htype in dictHtypeMax:
					if htypeNum > dictHtypeMax[htype][2]-3 and numLeaves < dictHtypeMax[htype][1]:
						dictHtypeMax[htype] = (internalNode,numLeaves,htypeNum)						
				else:
					dictHtypeMax[htype] = (internalNode,numLeaves,htypeNum)
	
		print(dictHtypeMax)
		
		numCountted = 0
		
		for i in range(17):
			htype = 'H'+str(i+1)
			(internalNode,numLeaves,htypeNum) = dictHtypeMax[htype]
			print(htype,':', numLeaves, htypeNum)
			numCountted = numCountted + htypeNum
			dictHtypeStrains, dictNtypeStrains = countMoreSubtypes(internalNode)
			if htypeNum/numLeaves > 0.8 and htypeNum != numLeaves:
				print('  exceptions:')
				for hhtype in dictHtypeStrains:
					if hhtype != htype:
						print(' ', hhtype, dictHtypeStrains[hhtype])
						numCountted = numCountted + len(dictHtypeStrains[hhtype])
						for genomeID in dictHtypeStrains[hhtype]:
							fwrite.write(genomeID+'\t'+seg+'\t'+hhtype+'\t'+htype+'\t'+str(numLeaves)+'\t'+str(htypeNum)+'\n')
		print('For all the ', len(tree1.leaf_nodes()), 'strains, we counted', numCountted, 'strains!')
		
		
		
	
	elif seg == 'NA':
		for internalNode in internalNodes:
			numLeaves = len(internalNode.leaf_nodes())
			dictHtypeStrains, dictNtypeStrains = countMoreSubtypes(internalNode)
			#print(numLeaves)
			for ntype in dictNtypeStrains:
				ntypeNum = len(dictNtypeStrains[ntype])
				if ntype in dictNtypeMax:
					if ntypeNum > dictNtypeMax[ntype][2]-3 and numLeaves < dictNtypeMax[ntype][1]:
						dictNtypeMax[ntype] = (internalNode,numLeaves,ntypeNum)						
				else:
					dictNtypeMax[ntype] = (internalNode,numLeaves,ntypeNum)
		print(dictNtypeMax)
		numCountted = 0
		for i in range(10):
			ntype = 'N'+str(i+1)
			(internalNode,numLeaves,ntypeNum) = dictNtypeMax[ntype]
			numCountted = numCountted + ntypeNum
			print(ntype,':', numLeaves, ntypeNum)
			dictHtypeStrains, dictNtypeStrains = countMoreSubtypes(internalNode)
			if ntypeNum/numLeaves > 0.8 and ntypeNum != numLeaves:
				print('  exceptions:')
				for nntype in dictNtypeStrains:
					if nntype != ntype:
						print(' ', nntype, dictNtypeStrains[nntype])
						for genomeID in dictNtypeStrains[nntype]:
							fwrite.write(genomeID+'\t'+seg+'\t'+nntype+'\t'+ntype+'\t'+str(numLeaves)+'\t'+str(ntypeNum)+'\n')
						numCountted = numCountted + len(dictNtypeStrains[nntype])
		print('For all the ', len(tree1.leaf_nodes()), 'strains, we counted', numCountted, 'strains!')
fwrite.close()
	
print('Finished!')


