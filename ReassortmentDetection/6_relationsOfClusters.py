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
from reportlab.lib.colors import Color, green, blue, yellow, red, maroon, purple, pink, gray, black, darkslategray
import numpy as np
import math
from datetime import datetime
import random

#reload(sys)
#sys.setdefaultencoding('utf8')

green50transparent = Color(0, 0, 100, alpha=0.5)



colorList = [green,blue,red]

dictHostColor = dict()
hosts = ['avian','swine','human']
i = 0
for host in hosts:
	dictHostColor[host] = colorList[i]
	i = i + 1

locList = ["East Africa","North Africa","South Africa","West Africa","Middle Africa","Central Asia","East Asia","South Asia","Southeast Asia","West Asia","Central Europe","Eastern Europe","Northern Europe","Southern Europe","Western Europe","Middle America","North America","The Caribbean","Oceania","Eastern South America","Midwest South American","Northern South America","Southern South America","Antarctica"]
colorJetList = ['#00008F','#00009F','#0000AF','#0000BF','#0000CF','#0000DF','#0000EF','#0000FF','#000FFF','#001FFF','#002FFF','#003FFF','#004FFF','#005FFF','#006FFF','#007FFF','#008FFF','#009FFF','#00AFFF','#00BFFF','#00CFFF','#00DFFF','#00EFFF','#00FFFF','#0FFFEF','#1FFFDF','#2FFFCF','#3FFFBF','#4FFFAF','#5FFF9F','#6FFF8F','#7FFF7F','#8FFF6F','#9FFF5F','#AFFF4F','#BFFF3F','#CFFF2F','#DFFF1F','#EFFF0F','#FFFF00','#FFEF00','#FFDF00','#FFCF00','#FFBF00','#FFAF00','#FF9F00','#FF8F00','#FF7F00','#FF6F00','#FF5F00','#FF4F00','#FF3F00','#FF2F00','#FF1F00','#FF0F00','#FF0000','#EF0000','#DF0000','#CF0000','#BF0000','#AF0000','#9F0000','#8F0000','#7F0000']

dictLocColor = {}

for i in range(len(locList)):
	loc = locList[i]
	colorIndex = int(64/len(locList)*i)%64
	# print(colorIndex)
	colorHere = colorJetList[colorIndex]
	dictLocColor[loc] = colorHere

print(os.getcwd())



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


def plotLegend(c, posX, posY, lenLegend):
	lenLegend = lenLegend/10
	c.line(posX, posY, posX+lenLegend, posY)
	c.line(posX, posY, posX, posY+lenLegend/10)
	c.line(posX+lenLegend, posY, posX+lenLegend, posY+lenLegend/10)
	c.drawString(posX+lenLegend/3, posY+lenLegend/10,'0.1')




		
def plotTree(tree1, seg, pageWidth, pageLength, dictLeafCluster, dictClusterLeaves, dictIsolateInfo, dictClusterInfo):
	
	dictLeafColor = dict()
	
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
	print(internalNodeToBeReRooted, 'is the new root.')
	tree1.ladderize(ascending=True)

	countPloted = 0
	c = canvas.Canvas(seg+'_relationsOfclustersHostDetail_202102241711_'+str(minimumClusterSize)+'.pdf', pagesize=(pageWidth, pageLength))
	#c.translate(3*inch,3*inch)
	
	xCircleCenter = pageWidth / 2
	yCircleCenter = (pageLength / 2) * 0.8
	
	barLength = 16

	leafDistFromRootList = [leaf.distance_from_root() for leaf in tree1.leaf_nodes()]
	leafDistFromRootList.sort()

	
	maxDistanceFromRoot = max(leafDistFromRootList)
	minDistanceFromRoot = min(leafDistFromRootList)
	#medianDistanceFromRoot = get_median([leaf.distance_from_root() for leaf in tree1.leaf_nodes()])
	
	dRadius = (pageWidth/2)*0.95/maxDistanceFromRoot
	# if seg == 'HA':
	# 	dRadius = 128

	dRadiusString = (leafDistFromRootList[-4])*1.01*dRadius
	print(leafDistFromRootList[-4],leafDistFromRootList[-1],leafDistFromRootList[0])
	
	hostRadius = maxDistanceFromRoot*dRadius*1.05
	hTypeRadius = hostRadius + 2*barLength
	nTypeRadius = hTypeRadius + 1.5*barLength

	print(seg, maxDistanceFromRoot, dRadius, dRadiusString, hostRadius, hTypeRadius, nTypeRadius)
	
	dAngle = (2*math.pi)*0.95/len(dictClusterLeaves)
	
	rHost = hostRadius*dAngle/3
	
	# lineWidth = r/3
	lineWidth = 0.1
	
	r = 0.5
	
	
	
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
			c.rect(posType[0], posType[1], r*3, r*3, fill=1)
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
			c.rect(posType[0], posType[1], r*5, r*5, fill=1)
			#c.circle(posType[0], posType[1], r*5, fill=1)
			c.setFillColor(black)
			c.drawString(posType[0]+5, posType[1],ntype)
	
	dictClusterNodeClusterNo = {}

	for clusterNo in dictClusterLeaves:
		taxon_labels = dictClusterLeaves[clusterNo]
		clusterNode = tree1.mrca(taxon_labels=taxon_labels)
		dictClusterNodeClusterNo[clusterNode] = clusterNo	  
	
	listClusterNode = []
	for nodeX in tree1.nodes():
		if nodeX in dictClusterNodeClusterNo:
			listClusterNode.append(nodeX)



	# dictClusterColor = dict()
	# clusterCount = 0
	
	# leafNotInClusterCount = 0

	dictNodePos = dict()
	j = 0
	# for leaf in tree1.leaf_nodes():
	for clusterNode in listClusterNode:
		
		clusterDistFromRoot = clusterNode.distance_from_root()
		clusterDist = clusterNode.edge_length
		# c.setLineWidth(r/3)
		# c.setStrokeColor(red)
		# c.line(xStart+(leafDistFromRoot-leafDist)*dRadius , yStart-j*dAngle, xStart+leafDistFromRoot*dRadius, yStart-j*dAngle)
		# c.setLineWidth(0.01)
		# c.setStrokeColor(green)
		# c.setFillColor(green)
		# c.circle(xStart+leafDistFromRoot*dRadius, yStart-j*dAngle, r, fill=1)
		
		xRadius = clusterDistFromRoot * dRadius
		#yAngle = math.pi-j*dAngle
		yAngle = j*dAngle
		
		#xPosition = xCircleCenter + leafDistFromRoot * dRadius * math.cos(math.pi-i*dAngle)
		#yPosition = yCircleCenter + leafDistFromRoot * dRadius * math.sin(math.pi-i*dAngle)		
		dictNodePos[clusterNode] = (xRadius, yAngle)
		
		j = j+1
		if j%1000 == 0:
			print(seg, j, 'postioning leafs')
	
	
	internalNodeFiltered = []

	for n in tree1.internal_nodes():
		offSpringFlag = 0
		for parentN in n.ancestor_iter(inclusive=False):
			if parentN in dictClusterNodeClusterNo:
				offSpringFlag = 1
		if offSpringFlag == 0:
			internalNodeFiltered.append(n)


		
		
	print('plsitioning internal nodes')
	while True:
		needTodoFlag = 0
		for n in internalNodeFiltered:			
			if n not in dictNodePos:
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
	for clusterNode in listClusterNode:
		j = j+1
		# leafDistFromRoot = leaf.distance_from_root()
		# leafDist = leaf.edge_length
		c.setLineWidth(0.1)
		clusterNo = dictClusterNodeClusterNo[clusterNode]
		# isolateID = leafName[:leafName.find('|')].replace(' ', '_')
		# host = 'unkown'
		# Htype = 'unkown'
		# Ntype = 'unkown'
		
		# if isolateID in dictIsolateInfo:
		# 	host = dictIsolateInfo[isolateID][1]
		# 	if 'avian' in host.lower():
		# 		host = 'Avian'
		# hostColor = gray
		# if host in dictHostColor:
		# 	hostColor = dictHostColor[host]
		# if isolateID in dictIsolateInfo:
		# 	subtype = dictIsolateInfo[isolateID][3]
		# 	Htype = subtype[:subtype.find('N')]
		# 	Ntype = subtype[subtype.find('N'):]
		

		# htypeColor = gray
		# ntypeColor = gray
		
		# if Htype in dictHtypeColor:
		# 	htypeColor = dictHtypeColor[Htype]
		# if Ntype in dictNtypeColor:
		# 	ntypeColor = dictNtypeColor[Ntype]
		

		# #sys.exit()
		# #fwrite.write(leafName+'\n')
		# # colorHere = random.choice(colorList)
		
		#c.setStrokeColor(hostColor)
		#c.setFillColor(hostColor)
		#def getPosXY(xCircleCenter, yCircleCenter, radius, angle):
		(posX, posY) = getPosXY(xCircleCenter, yCircleCenter, dictNodePos[clusterNode][0], dictNodePos[clusterNode][1])
		
		c.circle(posX, posY, r, fill=1)

		(posXString, posYString) = getPosXY(xCircleCenter, yCircleCenter, dRadiusString, dictNodePos[clusterNode][1])

		c.saveState()
		c.translate(posXString, posYString)
		c.rotate(dictNodePos[clusterNode][1]*360/(2*math.pi))
		c.setFont("Helvetica", rHost/2)
		c.setFillColor(black)
		
		
		clusterHostRange = []
		clusterLocRange = []
		clusterSubtypeRange = []
		
		if clusterNo in dictClusterInfo:
			clusterYearRange = list(dictClusterInfo[clusterNo])
			clusterYearRange.sort()
			strYear=str(clusterYearRange[0])+'-'+str(clusterYearRange[-1])
			for year in dictClusterInfo[clusterNo]:
				for (host, loc, subtype) in dictClusterInfo[clusterNo][year]:
					if 'avian' in host.lower():
						host = 'Avian'
					if host not in clusterHostRange:
						clusterHostRange.append(host)
					if loc not in clusterLocRange:
						clusterLocRange.append(loc)
					if subtype not in clusterSubtypeRange:
						clusterSubtypeRange.append(subtype)
			strHost = ','.join(clusterHostRange)
			strLoc = ','.join(clusterLocRange)
			strSubtype = ','.join(clusterSubtypeRange)

			#strCluster = clusterNo+': '+strYear+' | '+strHost+' | '+strLoc+' | '+strSubtype
			strCluster = clusterNo+': '+strYear
		else:
			strCluster = clusterNo
			print(clusterNo)
			strSubtype = ''
		
		hostIndex = 0
		rHostHere = rHost/4
		#rHost = 2*rHost
		for host in clusterHostRange:
			if 'avian' in host.lower():
				host = 'avian'
			if host in dictHostColor:
				c.setFillColor(dictHostColor[host])
			else:
				c.setFillColor(gray)
			c.circle(rHostHere*(hostIndex+22), 0, rHostHere, fill=1)
			hostIndex += 2.1
		
		locIndex = 0
		for loc in clusterLocRange:
			c.setFillColor(dictLocColor[loc])
			c.rect(rHostHere*(locIndex+31), -rHostHere, rHostHere*2, rHostHere*2, fill=1)
			locIndex = locIndex + 2.1



		c.setFillColor(black)
		c.drawString(rHostHere, -rHostHere, strCluster)

		c.drawString(rHostHere*(locIndex+1+31), -rHostHere, strSubtype)

		c.restoreState()


		c.setDash(2,1)
		c.setLineWidth(0.01)
		c.setFillColor(green50transparent)

		c.line(posX, posY, posXString, posYString)

		c.setDash(1,0)
		
		#add hosts info
		
		
		# clusterColor = dictLeafColor[leafName]
		# c.setStrokeColor(clusterColor)
		# c.setFillColor(clusterColor)
		
		# #def getPosXY(xCircleCenter, yCircleCenter, radius, angle):
		# # if seg not in ['HA']:
		# 	# hostRadius = hostRadius / 2
		# (posX1, posY1) = getPosXY(xCircleCenter, yCircleCenter, hostRadius-barLength/2, dictNodePos[leaf][1])
		# (posX2, posY2) = getPosXY(xCircleCenter, yCircleCenter, hostRadius+barLength/2, dictNodePos[leaf][1])
		
		# #c.circle(posX, posY, rHost, fill=1)
		# #c.rect(posX-rHost/2, posY-rHost/2, rHost, rHost, fill=1)
		# c.line(posX1, posY1, posX2, posY2)
		
		# plot subtypes? no.
		
		
		# if seg == 'HA':		
		# 	c.setStrokeColor(htypeColor)
		# 	c.setFillColor(htypeColor)
		# 	(posX1, posY1) = getPosXY(xCircleCenter, yCircleCenter, nTypeRadius - barLength/2, dictNodePos[leaf][1])
		# 	(posX2, posY2) = getPosXY(xCircleCenter, yCircleCenter, nTypeRadius + barLength*1.5, dictNodePos[leaf][1])
		# 	c.line(posX1, posY1, posX2, posY2)
		
		# if seg == 'NA':
		# 	c.setStrokeColor(ntypeColor)
		# 	c.setFillColor(ntypeColor)
		# 	(posX1, posY1) = getPosXY(xCircleCenter, yCircleCenter, nTypeRadius - barLength/2, dictNodePos[leaf][1])
		# 	(posX2, posY2) = getPosXY(xCircleCenter, yCircleCenter, nTypeRadius + barLength*1.5, dictNodePos[leaf][1])
		# 	c.line(posX1, posY1, posX2, posY2)
		
		
		
		
		
		
		
		
	# 	c.setLineWidth(lineWidth)
	
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

def getClusterInfo():
	dictClusterInfo = {}
	dictIsolateSegCluster = {}
	dictIsolateInfo = {}
	fread = open('isolatesCombinationsWithNameHostDetail-6.txt', 'r')
	for rline in fread:
		if 'isolateID' not in rline:
			lrline = rline.strip().split('\t')
			#isolateID	isolateName	host_Classification	Continet_detail	subtype	year	Collecting Date	PB2	PB1	PA	HA	NP	NA	MP	NS
			isolateID = lrline[0]
			isolateName = lrline[1]
			host_Classification = lrline[2]
			loc = lrline[3]
			subtype = lrline[4]
			year = int(lrline[5])
			dictIsolateInfo[isolateID] = rline.strip()
			dictIsolateSegCluster[isolateID] = {}
			hostLocSubtype = (host_Classification, loc, subtype)
			clusters = lrline[7:]
			segIndex = 0
			for seg in segs:
				cluster = lrline[7+segIndex]
				dictIsolateSegCluster[isolateID][seg] = cluster
				if cluster in dictClusterInfo:
					if year in dictClusterInfo[cluster]:
						if hostLocSubtype in dictClusterInfo[cluster][year]:
							dictClusterInfo[cluster][year][hostLocSubtype] = dictClusterInfo[cluster][year][hostLocSubtype] + 1
						else:
							dictClusterInfo[cluster][year][hostLocSubtype] = 1
					else:
						dictClusterInfo[cluster][year] = {}
						dictClusterInfo[cluster][year][hostLocSubtype] = 1
				else:
					dictClusterInfo[cluster] = {}
					dictClusterInfo[cluster][year] = {}
					dictClusterInfo[cluster][year][hostLocSubtype] = 1			
				segIndex += 1
	fread.close()
	return dictClusterInfo



if __name__ == "__main__":

	#colorList = [green, yellow, blue, red, maroon, purple, pink]
	colorList = [green, blue, red]

	segs = ['PB2','PB1','PA','HA','NP','NA','MP','NS']

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

	dictInternalNodeClusterNo = {}
	fread = open('clusterNo.txt', 'r')
	lineNo = 0
	for rline in fread:
		lineNo += 1
		if lineNo > 1:
			lrline = rline.strip().split('\t')
			internalNodeDSP = lrline[4]
			clusterNo = lrline[0]
			dictInternalNodeClusterNo[internalNodeDSP] = clusterNo
	fread.close()



	dictIsolateInfo = dict()
	fread = open('isolateHostClassLocationRankedByTimeHostDetail-6.txt', 'r')
	for rline in fread:
		# if rline.startswith('EPI'):
		isolateID = rline[:rline.find('\t')]
		#isolateID	host_Classification	Continet_detail	subtype	year	Collecting Date
		dictIsolateInfo[isolateID] = rline.strip().split('\t')

	fread.close()

	dictHostColor = dict()
	hosts = ['avian','swine','human']
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
	

	dictClusterInfo = getClusterInfo()
		
	for seg in segs:
		print(seg)

		fread = open(seg+'_clusters.txt'+str(minimumClusterSize), 'r')
		dictLeafCluster = dict()
		dictClusterLeaves = {}
		for rline in fread:
			if rline.startswith('>>>'):
				internalNodeDSP = rline.strip().split('\t')[0]
				clusterNo = dictInternalNodeClusterNo[internalNodeDSP]
				dictClusterLeaves[clusterNo] = []
			else:
				leafName = rline.strip()
				if leafName in dictLeafCluster:
					print('repeated leafName...')
					sys.exit()
				dictLeafCluster[leafName] = clusterNo
				dictClusterLeaves[clusterNo].append(leafName)
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

		# sys.exit()
		
		plotTree(tree1, seg, pageWidth, pageLength, dictLeafCluster, dictClusterLeaves, dictIsolateInfo, dictClusterInfo)
		
	print('Finished!')


# taxon_labels=['Python sebae',
#			   'Python regius',
#			   'Python curtus',
#			   'Python molurus']
# mrca = tree.mrca(taxon_labels=taxon_labels)
