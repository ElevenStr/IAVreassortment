# coding=utf-8
from __future__ import division
import dendropy
import datetime
import time
import sys

#reload(sys)
#import xlrd
import types
import sys
import os
import shutil
import string
from random import choice
import datetime
import time
#from ete2 import Tree, NodeStyle, TreeStyle, phyloxml
# from Bio import Phylo
# from Bio.Phylo.PhyloXML import Phylogeny
import pylab
import sys
import matplotlib
from numpy.random import randn
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch
from reportlab.lib.colors import green, blue, yellow, red, maroon, purple, pink, gray, black, darkslategray
import numpy as np
import math
from datetime import datetime
import random

#reload(sys)
#sys.setdefaultencoding('utf8')

#colorList = [green, yellow, blue, red, maroon, purple, pink]
# colorList = [green, yellow, blue, red]
colorList = [green, blue, red]
#['Avian','Swine','Human']

#['Avian','Equine','Human','Swine']

colorJetList = ['#00008F','#00009F','#0000AF','#0000BF','#0000CF','#0000DF','#0000EF','#0000FF','#000FFF','#001FFF','#002FFF','#003FFF','#004FFF','#005FFF','#006FFF','#007FFF','#008FFF','#009FFF','#00AFFF','#00BFFF','#00CFFF','#00DFFF','#00EFFF','#00FFFF','#0FFFEF','#1FFFDF','#2FFFCF','#3FFFBF','#4FFFAF','#5FFF9F','#6FFF8F','#7FFF7F','#8FFF6F','#9FFF5F','#AFFF4F','#BFFF3F','#CFFF2F','#DFFF1F','#EFFF0F','#FFFF00','#FFEF00','#FFDF00','#FFCF00','#FFBF00','#FFAF00','#FF9F00','#FF8F00','#FF7F00','#FF6F00','#FF5F00','#FF4F00','#FF3F00','#FF2F00','#FF1F00','#FF0F00','#FF0000','#EF0000','#DF0000','#CF0000','#BF0000','#AF0000','#9F0000','#8F0000','#7F0000']

dictHTypeNumberColor = dict()

for i in range(18):
	colorHere = colorJetList[int(3.5*i)%64]
	dictHTypeNumberColor[i+1] = colorHere

dictNTypeNumberColor = dict()

for i in range(11):
	colorHere = colorJetList[int(6.3*i)%64]
	dictNTypeNumberColor[i+1] = colorHere



minimumClusterSize = 100

def get_median(data):
	data = sorted(data)
	size = len(data)
	if size % 2 == 0:   # 判断列表长度为偶数
		median = (data[size//2]+data[size//2-1])/2
		data[0] = median
	if size % 2 == 1:   # 判断列表长度为奇数
		median = data[(size-1)//2]
		data[0] = median
	return data[0]

def getYAngle(listChildPos):
	ymax = -9999999999
	ymin = 9999999999
	for (x,y) in listChildPos:
		if y > ymax:
			ymax = y
		if y < ymin:
			ymin = y
	return ((ymax+ymin)/2, ymax, ymin)


def getTopPercent(data, percent):
	data = sorted(data, reverse=True)
	size = len(data)
	posTop = int(size*percent)
	
	# if size % 2 == 0: # 判断列表长度为偶数
		# median = (data[size//2]+data[size//2-1])/2
		# data[0] = median
	# if size % 2 == 1: # 判断列表长度为奇数
		# median = data[(size-1)//2]
		# data[0] = median
	
	return data[posTop]

def getPosXY(xCircleCenter, yCircleCenter, radius, angle):
	xPosition = xCircleCenter + radius * math.cos(angle)
	yPosition = yCircleCenter + radius * math.sin(angle)
	return (xPosition, yPosition)

	
def arcs(canvas):
	# canvas.setLineWidth(4)
	# canvas.setStrokeColorRGB(0.8, 1, 0.6)
	# # draw rectangles enclosing the arcs
	# canvas.rect(inch, inch, 1.5*inch, inch)
	# canvas.rect(3*inch, inch, inch, 1.5*inch)
	# canvas.setStrokeColorRGB(0, 0.2, 0.4)
	# canvas.setFillColorRGB(1, 0.6, 0.8)
	p = canvas.beginPath()
	#p.moveTo(0.2*inch, 0.2*inch)
	#p.arcTo(inch, inch, 2.5*inch,2*inch, startAng=-30, extent=135)
	p.arc(3*inch, inch, 4*inch, 2.5*inch, startAng=90, extent=180)
	canvas.drawPath(p, fill=0, stroke=1)
	
def plotArc(c, xCircleCenter, yCircleCenter, radius, yAngleMax, yAngleMin):
	# x1 = xCircleCenter + radius * math.cos(yAngleMin)
	# y1 = yCircleCenter + radius * math.sin(yAngleMin)
	# x2 = xCircleCenter + radius * math.cos(yAngleMax)
	# y2 = yCircleCenter + radius * math.sin(yAngleMax)
	
	# x122 = xCircleCenter + radius * math.cos((yAngleMax+yAngleMin)*2/3)
	# y122 = xCircleCenter + radius * math.sin((yAngleMax+yAngleMin)*2/3)
	# x121 = xCircleCenter + radius * math.cos((yAngleMax+yAngleMin)/3)
	# y121 = xCircleCenter + radius * math.sin((yAngleMax+yAngleMin)/3)
	
	startAngArc = (360/(2*math.pi)) * yAngleMin
	extendArc = (360/(2*math.pi)) * (yAngleMax - yAngleMin)
	
	p = c.beginPath()
	x1rectN = xCircleCenter - radius
	y1rectN = yCircleCenter - radius
	x2rectN = xCircleCenter + radius
	y2rectN = yCircleCenter + radius
	p.arc(x1rectN, y1rectN, x2rectN, y2rectN, startAng=startAngArc, extent=extendArc)

	c.drawPath(p, stroke=1, fill=0)

segs = ['PB2','PB1','PA','HA','NP','NA','MP','NS']
# segs = ['PB2']
segsHANA = ['HA', 'NA']

treeDir = '../iqtree'

prefixSeg = ''
suffixSeg = '.fasta.codon.fas.correctSubtypes.fas.correctSubtypes2.fas.treefile'


# for seg in segs:
	# print seg

	# tree1 = dendropy.Tree.get(path=os.path.join(treeDir,prefixSeg+seg+suffixSeg), schema='newick')
	# rootedFlag = 0
	# for leafNode in tree1.leaf_nodes():
		# taxonLeaf = leafNode.taxon
		# labelLeaf = taxonLeaf.label
		# if '105892' in labelLeaf and 'H17N10' in labelLeaf:
			# print labelLeaf, 'is the new root.'
			# rootedFlag = 1
			# tree1.reroot_at_node(leafNode)
			# break
	# if rootedFlag == 0:
		# print 'Not Rooted well.'
		# sys.exit()

def plotLegend(c, posX, posY, lenLegend):
	lenLegend = lenLegend/10
	c.line(posX, posY, posX+lenLegend, posY)
	c.line(posX, posY, posX, posY+lenLegend/10)
	c.line(posX+lenLegend, posY, posX+lenLegend, posY+lenLegend/10)
	c.drawString(posX+lenLegend/3, posY+lenLegend/10,'0.1')
		
def plotTree(tree1, seg, pageWidth, pageLength, dictLeafCluster, dictIsolateInfo):
	
	dictLeafColor = dict()
	dictClusterLeaves = dict()
	for leaf in dictLeafCluster:
		cluster = dictLeafCluster[leaf]
		if cluster in dictClusterLeaves:
			dictClusterLeaves[cluster].add(leaf)
		else:
			s = set()
			s.add(leaf)
			dictClusterLeaves[cluster] = s
	
	# if seg != 'PB1':
	rootedFlag = 0
	maxInternalEdgeLength = 0
	internalNodeToBeReRooted = None
	for internalNode in tree1.internal_nodes():
		if internalNode.edge_length != None and internalNode.edge_length > maxInternalEdgeLength :
			maxInternalEdgeLength = internalNode.edge_length
			internalNodeToBeReRooted = internalNode
	#sys.exit()
	treeRoot = internalNodeToBeReRooted
	tree1.reroot_at_node(internalNodeToBeReRooted)
	# for leafNode in internalNodeToBeReRooted.leaf_nodes():
	# 	taxonLeaf = leafNode.taxon
	# 	labelLeaf = taxonLeaf.label
	# 	if 'H17N10' in labelLeaf:
	# 		print(labelLeaf, 'is the H17N10 new root.')
	print(internalNodeToBeReRooted, 'is the new root.')
	tree1.ladderize(ascending=True)
	
	



	countPloted = 0
	c = canvas.Canvas(seg+'_polar_clustered_withHostsNowaiquan'+str(minimumClusterSize)+'.pdf', pagesize=(pageWidth, pageLength))
	#c.translate(3*inch,3*inch)
	
	xCircleCenter = pageWidth / 2 - 80
	yCircleCenter = (pageLength / 2) * 0.8
	
	barLength = 16
	
	maxDistanceFromRoot = max([leaf.distance_from_root() for leaf in tree1.leaf_nodes()])
	minDistanceFromRoot = min([leaf.distance_from_root() for leaf in tree1.leaf_nodes()])
	#medianDistanceFromRoot = get_median([leaf.distance_from_root() for leaf in tree1.leaf_nodes()])

	#the Radius of the outermost layer
	dRadius = (pageWidth/2)*0.75/maxDistanceFromRoot
	print(maxDistanceFromRoot,dRadius)
	
	hostRadius = maxDistanceFromRoot*dRadius*1.05
	hTypeRadius = hostRadius + 2*barLength
	nTypeRadius = hTypeRadius + 1.5*barLength
	
	dAngle = (2*math.pi)*0.95/len(tree1.leaf_nodes())
	
	rHost = hostRadius*dAngle/3
	
	# lineWidth = r/3
	#initial = 0.1
	lineWidth = 0.01

	#initial = 0.5
	r = 0.05
	
	
	
	plotLegend(c, pageWidth/8, pageLength/1.2, dRadius)
	
	#plot host Color legend
	hostCount = 0
	for host in dictHostColor:
		hostCount = hostCount + 1
		colorHere = dictHostColor[host]
		c.setStrokeColor(colorHere)
		c.setFillColor(colorHere)
		posHost = [pageWidth/8,pageLength/1.2 - hostCount*10]
		#c.rect(posHost[0], posHost[1], r*3, r*3, fill=1)
		c.circle(posHost[0], posHost[1], r*5, fill=1)
		c.setFillColor(black)
		c.drawString(posHost[0]+5, posHost[1],host)
		
	if seg == 'HA':
		for i in range(18):
			typeNumber = i+1
			colorHere = dictHTypeNumberColor[typeNumber]			
			htype = 'H'+str(typeNumber)
			c.setStrokeColor(colorHere)
			c.setFillColor(colorHere)
			posType = [pageWidth/1.3,pageLength/1.2 - i*10]
			'''
			c.rect(posType[0], posType[1], r*3, r*3, fill=1)
			'''
			#c.circle(posType[0], posType[1], r*5, fill=1)
			c.setFillColor(black)
			c.drawString(posType[0]+5, posType[1],htype)		
	if seg == 'NA':
		for i in range(11):
			typeNumber = i+1
			colorHere = dictNTypeNumberColor[typeNumber]			
			ntype = 'N'+str(typeNumber)
			c.setStrokeColor(colorHere)
			c.setFillColor(colorHere)
			posType = [pageWidth/1.3,pageLength/1.2 - i*10]
			'''
			c.rect(posType[0], posType[1], r*5, r*5, fill=1)
			'''
			#c.circle(posType[0], posType[1], r*5, fill=1)
			c.setFillColor(black)
			c.drawString(posType[0]+5, posType[1],ntype)
	
	
	
	
	# dRadius = inch*2
	# r = inch/360
	# dAngle = inch/180

	# lineWidth = dAngle/3

	dictClusterColor = dict()
	clusterCount = 0
	
	leafNotInClusterCount = 0

	dictNodePos = dict()
	j = 0
	for leaf in tree1.leaf_nodes():
		
		leafDistFromRoot = leaf.distance_from_root()
		leafDist = leaf.edge_length
		# c.setLineWidth(r/3)
		# c.setStrokeColor(red)
		# c.line(xStart+(leafDistFromRoot-leafDist)*dRadius , yStart-j*dAngle, xStart+leafDistFromRoot*dRadius, yStart-j*dAngle)
		# c.setLineWidth(0.01)
		# c.setStrokeColor(green)
		# c.setFillColor(green)
		# c.circle(xStart+leafDistFromRoot*dRadius, yStart-j*dAngle, r, fill=1)
		
		xRadius = leafDistFromRoot * dRadius
		#yAngle = math.pi-j*dAngle
		yAngle = j*dAngle
		
		#xPosition = xCircleCenter + leafDistFromRoot * dRadius * math.cos(math.pi-i*dAngle)
		#yPosition = yCircleCenter + leafDistFromRoot * dRadius * math.sin(math.pi-i*dAngle)		
		dictNodePos[leaf] = (xRadius, yAngle)
		
		j = j+1
		if j%1000 == 0:
			print(seg, j, 'postioning leafs')
		
		
		
		leafName = str(leaf.taxon.label)
		# print(leafName,'leafName')
		if leafName in dictLeafCluster:
			leafCluster = dictLeafCluster[leafName]
			if len(dictClusterLeaves[leafCluster]) == 1:
				dictLeafColor[leafName] = gray
				# print('here')
			else:
				if leafCluster in dictClusterColor:
					dictLeafColor[leafName] = dictClusterColor[leafCluster]
				else:
					clusterCount += 1
					# dictClusterColor[leafCluster] = colorList[clusterCount%len(colorList)]
					dictClusterColor[leafCluster] = colorJetList[clusterCount%len(colorJetList)]
					dictLeafColor[leafName] = dictClusterColor[leafCluster]
		else:
			#print(leafName,'not in dictLeafCluster...')
			dictLeafColor[leafName] = gray
			leafNotInClusterCount += 1
	print(leafNotInClusterCount, 'leaves not in clusters...')
			
		
		
		
	print('positioning internal nodes')
	while True:
		needTodoFlag = 0
		for n in tree1.internal_nodes():
			if n not in dictNodePos and len(n.leaf_nodes()) > 1:
				flagToDoAgain = 0
				listChildPos = []
				for childNode in n.child_nodes():
					if childNode not in dictNodePos:
						flagToDoAgain = 1
					else:
						listChildPos.append(dictNodePos[childNode])
				if flagToDoAgain == 0:
					(posYAngleN, yAngleMax, yAngleMin) = getYAngle(listChildPos)
					posXRadiusN = dRadius*n.distance_from_root()
					dictNodePos[n] = (posXRadiusN, posYAngleN)
					for childNode in n.child_nodes():
						posChildNode = dictNodePos[childNode]
						c.setLineWidth(lineWidth)
						c.setStrokeColor(darkslategray)
						(posX0, posY0) = getPosXY(xCircleCenter, yCircleCenter, posXRadiusN, posChildNode[1])
						(posX1, posY1) = getPosXY(xCircleCenter, yCircleCenter, posChildNode[0], posChildNode[1])						
						c.line(posX0, posY0, posX1, posY1)
						countPloted = countPloted + 1
						if countPloted % 1000 == 0:
							print(seg, countPloted, 'internal_nodes positioned.')
					c.setLineWidth(lineWidth)
					c.setStrokeColor(darkslategray)
					# (posX0, posY0) = getPosXY(xCircleCenter, yCircleCenter, posXRadiusN, yAngleMax)
					# (posX1, posY1) = getPosXY(xCircleCenter, yCircleCenter, posXRadiusN, yAngleMin)
					# c.line(posX0, posY0, posX1, posY1)
					# if posXRadiusN != 0:
						# angleSaw = (lineWidth/2)*2*math.pi/posXRadiusN
					# else:
						# angleSaw = 0
					angleSaw = 0
					plotArc(c, xCircleCenter, yCircleCenter, posXRadiusN, yAngleMax+angleSaw, yAngleMin-angleSaw)
					#c.line(posXofN, ymax+(lineWidth/2), posXofN, ymin-(lineWidth/2))
				else:
					needTodoFlag = 1
		if needTodoFlag == 0:
			break
	j = 0
	colorHere = black
	for leaf in tree1.leaf_nodes():
		j = j+1
		# leafDistFromRoot = leaf.distance_from_root()
		# leafDist = leaf.edge_length
		c.setLineWidth(0.01)
		leafName = str(leaf.taxon.label)
		# print(leafName)
		isolateID = leafName.split(' ')[0]
		# print(isolateID)
		host = 'unkown'
		Htype = 'unkown'
		Ntype = 'unkown'

		if isolateID in dictIsolateInfo:
			host = dictIsolateInfo[isolateID][1]
			if 'avian' in host.lower():
				host = 'Avian'
			if host == 'human':
				host = 'Human'
			if host == 'swine':
				host = 'Swine'
			# print(host)
		hostColor = gray
		if host in dictHostColor:
			# print('here')
			hostColor = dictHostColor[host]
		if isolateID in dictIsolateInfo:
			subtype = dictIsolateInfo[isolateID][3]
			Htype = subtype[:subtype.find('N')]
			Ntype = subtype[subtype.find('N'):]
		

		htypeColor = gray
		ntypeColor = gray
		
		if Htype in dictHtypeColor:
			htypeColor = dictHtypeColor[Htype]
		if Ntype in dictNtypeColor:
			ntypeColor = dictNtypeColor[Ntype]
		

		#sys.exit()
		#fwrite.write(leafName+'\n')
		# colorHere = random.choice(colorList)

		clusterColor = dictLeafColor[leafName]
		#initial = hostColor
		c.setStrokeColor(clusterColor)
		c.setFillColor(clusterColor)
		#def getPosXY(xCircleCenter, yCircleCenter, radius, angle):
		(posX, posY) = getPosXY(xCircleCenter, yCircleCenter, dictNodePos[leaf][0], dictNodePos[leaf][1])
		#plot every leaf node circle
		c.circle(posX, posY, r, fill=1)
		#initial = 10*
		c.setLineWidth(100*lineWidth)
		
		#add hosts info
		
		
		clusterColor = dictLeafColor[leafName]
		c.setStrokeColor(clusterColor)
		c.setFillColor(clusterColor)
		
		#def getPosXY(xCircleCenter, yCircleCenter, radius, angle):
		# if seg not in ['HA']:
			# hostRadius = hostRadius / 2
		(posX1, posY1) = getPosXY(xCircleCenter, yCircleCenter, hostRadius-barLength/2, dictNodePos[leaf][1])
		(posX2, posY2) = getPosXY(xCircleCenter, yCircleCenter, hostRadius+barLength/2, dictNodePos[leaf][1])
		
		#c.circle(posX, posY, rHost, fill=1)
		#c.rect(posX-rHost/2, posY-rHost/2, rHost, rHost, fill=1)
		c.line(posX1, posY1, posX2, posY2)
		
		# plot subtypes? no.
		
		
		if seg == 'HA':		
			c.setStrokeColor(htypeColor)
			c.setFillColor(htypeColor)
			(posX1, posY1) = getPosXY(xCircleCenter, yCircleCenter, nTypeRadius - barLength/2, dictNodePos[leaf][1])
			(posX2, posY2) = getPosXY(xCircleCenter, yCircleCenter, nTypeRadius + barLength*1.5, dictNodePos[leaf][1])
			'''
			c.line(posX1, posY1, posX2, posY2)
			'''
		
		if seg == 'NA':
			c.setStrokeColor(ntypeColor)
			c.setFillColor(ntypeColor)
			(posX1, posY1) = getPosXY(xCircleCenter, yCircleCenter, nTypeRadius - barLength/2, dictNodePos[leaf][1])
			(posX2, posY2) = getPosXY(xCircleCenter, yCircleCenter, nTypeRadius + barLength*1.5, dictNodePos[leaf][1])
			'''
			c.line(posX1, posY1, posX2, posY2)
			'''
		
		
		
		
		
		
		
		
		c.setLineWidth(lineWidth)
	
	# put labels into the internal nodes...
	
	
	# c.setStrokeColor(black)
	# c.setFillColor(black)
	
	# edgeLengthGate = getTopPercent([internalNode.edge_length for internalNode in tree1.internal_nodes()], 0.01)
	
	# for internalNode in tree1.internal_nodes():
		# if internalNode.edge_length >= edgeLengthGate and internalNode.label:
			# (posX, posY) = getPosXY(xCircleCenter, yCircleCenter, dictNodePos[internalNode][0], dictNodePos[internalNode][1])
			# c.drawString(posX, posY, internalNode.label)
	
	colorHere = black
	c.setStrokeColor(colorHere)
	c.setFillColor(colorHere)
	c.circle(xCircleCenter, yCircleCenter, r*3, fill=1)
	
	c.showPage()
	c.save()

dictIsolateInfo = dict()
fread = open('isolateHostClassLocationRankedByTimeHostDetail-6.txt', 'r')
for rline in fread:
	if 'isolateID' not in rline:
		isolateID = rline[:rline.find('\t')]
		#isolateID	host_Classification	Continet_detail	subtype	year	Collecting Date
		dictIsolateInfo[isolateID] = rline.strip().split('\t')

fread.close()

dictHostColor = dict()
# hosts = ['Avian','Equine','Swine','Human']
hosts = ['Avian','Swine','Human']
i = 0
for host in hosts:
	dictHostColor[host] = colorList[i]
	i = i + 1

dictHtypeColor = dict()
dictNtypeColor = dict()
for i in range(18):
	typeNumber = i+1
	colorHere = dictHTypeNumberColor[typeNumber]
	dictHtypeColor['H'+str(typeNumber)] = colorHere
	
for i in range(11):
	typeNumber = i+1
	colorHere = dictNTypeNumberColor[typeNumber]
	dictNtypeColor['N'+str(typeNumber)] = colorHere	
	

	
for seg in segs:
	print(seg)	
	
	fread = open(seg+'_clusters.txt'+str(minimumClusterSize), 'r')
	dictLeafCluster = dict()
	for rline in fread:
		if rline.startswith('>>>'):
			cluster = rline.strip()
		else:
			leafName = rline.strip()
			if leafName in dictLeafCluster:
				print('repeated leafName...')
				sys.exit()
			dictLeafCluster[leafName] = cluster
	fread.close()
	
	tree1 = dendropy.Tree.get(path=os.path.join(treeDir,prefixSeg+seg+suffixSeg), schema='newick')
	rootedFlag = 0
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

	
	(pageWidth, pageLength) = (2100, 2970)
	
	plotTree(tree1, seg, pageWidth, pageLength, dictLeafCluster, dictIsolateInfo)
	
print('Finished!')


