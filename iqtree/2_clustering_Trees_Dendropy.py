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
#from reportlab.lib import colors
#from reportlab.lib.units import inch

#from reportlab.pdfgen import canvas
#from reportlab.lib.units import inch
#from reportlab.lib.colors import green, blue, yellow, red, maroon, purple, pink, gray, black, darkslategray
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

def get_median(data):
	data.sort()
	half = len(data) // 2
	return (data[half] + data[~half]) / 2

	
	
# segs = ['PB2','PB1','PA','HA','NP','NA','MP','NS']
segs = ['NS']
#segs = ['NA', 'MP', 'NS']
treeDir = ''

prefixSeg = ''
suffixSeg = '.fasta.codon.fas.correctSubtypes.fas.correctSubtypes2.fas.treefile'

minimumClusterSize = 100


for seg in segs:
	print(seg)
	
	mpdFile = seg+'_mpdLengths.txt'
	mpdList = []
	nodeNumList = []
	fread = open(mpdFile, 'r')
	print('Reading mpdfiles...', mpdFile, 'calculating median values...')
	for rline in fread:
		lrline = rline.strip().split('\t')
		mpd = float(lrline[0])
		nodeNum = int(lrline[1])
		mpdList.append(mpd)
		nodeNumList.append(nodeNum)
	mpdMedian = get_median(mpdList)
	nodeNumMedian = max(get_median(nodeNumList), minimumClusterSize)
	
	print('mpdMedian', mpdMedian)
	print('nodeNumMedian', nodeNumMedian)
	
	fread.close()
	
	tree1 = dendropy.Tree.get(path=os.path.join(treeDir,prefixSeg+seg+suffixSeg), schema='newick')

	# rootedFlag = 0
	# for leafNode in tree1.leaf_nodes():
		# taxonLeaf = leafNode.taxon
		# labelLeaf = taxonLeaf.label
		# if '105892' in labelLeaf and 'H17N10' in labelLeaf:
			# print(labelLeaf, 'is the new root.')
			# rootedFlag = 1
			# tree1.reroot_at_node(leafNode)
			# break
	# if rootedFlag == 0:
		# print('Not Rooted well.')
		# sys.exit()
	
	#tree1.ladderize(ascending=True)
	
	
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
	
	
	internalNodes = treeRoot.child_nodes()
	leafnodesOfIN = treeRoot.leaf_nodes()
	
	dictVisited = dict()
	
	
	
	pdc = tree1.phylogenetic_distance_matrix()
	
	
	queue = []
	queue.append(treeRoot)
	
	clusters = []
	
	#sys.exit()
	
	#fwrite = open(seg+'_mpdLengths.txt', 'w')
	
	countInternal = 0
	
	while True:
		nodeCurrent = queue.pop(0)
		if nodeCurrent.is_internal() :
			countInternal = countInternal + 1
			mpdIN = calMPD(pdc, nodeCurrent)
			#fwrite.write(str(mpdIN)+'\t'+str(len(nodeCurrent.leaf_nodes()))+'\t'+str(nodeCurrent)+'\n')
			#print(countInternal)
			numLeaves = len(nodeCurrent.leaf_nodes())
			if numLeaves < nodeNumMedian+1 or mpdIN < mpdMedian:
				clusters.append(nodeCurrent)
			else:				
				nodeChildren = nodeCurrent.child_nodes()
				for nodeChild in nodeChildren:
					queue.append(nodeChild)
		else:
			clusters.append(nodeCurrent)
		if len(queue) < 1:
			break
	print('clusters num is', len(clusters))
	
	
	fwrite = open(seg+'_clusters.txt'+str(minimumClusterSize), 'w')
	for node in clusters:
		if node.is_internal():
			mpdIN = calMPD(pdc, node)
			fwrite.write('>>>'+ str(node)+'\t'+str(len(node.leaf_nodes()))+'\t'+str(mpdIN)+'\n')
			for leaf in node.leaf_nodes():
				fwrite.write(str(leaf.taxon.label)+'\n')
		else:
			fwrite.write('>>>'+ str(node)+'\t'+'1'+'\t'+'0'+'\n')
			fwrite.write(str(node.taxon.label)+'\n')
	fwrite.close()
	#fwrite.close()
	# countClusters = 0
	# for node in clusters:
		# if node.is_internal():
			# countClusters = countClusters + 1
			# numLeaves = len(node.leaf_nodes())
			# print(numLeaves, node)
			# mpdIN = calMPD(pdc, node)
			# print mpdIN
			# # if numLeaves > 1000:
				# # sys.exit()
	# print(countClusters, 'internalNodes')

	# for i, t1 in enumerate(tree1.taxon_namespace[:-1]):
		# for t2 in tree1.taxon_namespace[i+1:]:
			# print("Distance between '%s' and '%s': %s" % (t1.label, t2.label, pdc(t1, t2)))
	
	#sys.exit()
	
print('Finished!')


