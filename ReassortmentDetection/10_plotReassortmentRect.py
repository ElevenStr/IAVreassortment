# coding=utf-8
from __future__ import division
import datetime
import time
import sys

#reload(sys)
#import xlrd
import types
import sys
import os
from datetime import datetime
import matplotlib
from numpy.random import randn
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch
from reportlab.lib.colors import green, blue, yellow, red, maroon, purple, pink, gray, black, darkslategray
from reportlab.lib.colors import Color, seagreen, forestgreen, limegreen, green, lawngreen, yellow, blue, red
from reportlab.lib.colors import pink,wheat,green,purple,springgreen, yellow, blue, red

import math

from functools import cmp_to_key

def distTypesIJ(typeI, typeJ):
	simCount = 0
	for k in range(len(typeI)):
		if typeI[k] == typeJ[k]:
			simCount += 1
	return (8-simCount)

def strType(typeI):
	strT = '+'.join(list(typeI))
	return(strT)

dictTypeInfo = {}
lTypesRankedByTime = []
dictStrainTypeInfo = {}
dictStrainType = {}
fread = open('isolatesCombinationsWithNameHostDetail-6.txt', 'r')
lineNo = 0
for rline in fread:
	if lineNo > 0:
		lrline = rline.strip().split('\t')
		strainID = lrline[0]
		dictStrainTypeInfo[strainID] = lrline
		typeHere = tuple(lrline[7:])
		if typeHere not in lTypesRankedByTime:
			lTypesRankedByTime.append(typeHere)
		typeInfo = tuple(lrline[:7])
		dictStrainType[strainID] = typeHere
		if typeHere in dictTypeInfo:
			dictTypeInfo[typeHere].add(typeInfo)
		else:			
			s = set()
			s.add(typeInfo)
			dictTypeInfo[typeHere] = s
	lineNo += 1
fread.close()

lTypes = lTypesRankedByTime

lenTypes = len(lTypes)

print('There are ',str(lenTypes), ' types in total...')



dictTypeIndex = {}
dictIndexType = {}
dictIndexTypeSet = {}
dictTypeYearStart = {}
dictTypeIndexInfo = {}

dictTypeInfoRanges = {}


fwrite = open('typeIndexInfoAvianDetailFinalHostDetal.txt', 'w')
fwrite.write('typeIndex\tstrainCount\tCombination\tyears\thosts\tlocs\tlocStarts\thostStart\tsubtype\tstrains\n')
typeIndex = 0
for typeK in lTypes:
	strTypeIndex = 'type'+str(typeIndex)
	dictTypeIndex[typeK] = strTypeIndex
	dictIndexType[strTypeIndex] = typeK
	dictIndexTypeSet[strTypeIndex] = set(typeK)
	fwrite.write(dictTypeIndex[typeK]+'\t'+str(len(dictTypeInfo[typeK]))+'\t'+strType(typeK)+'\t')

	typeInfo = []

	typeHostRange = []
	typeLocRange = []
	typeSubtypeRange = []
	typeYearRange = []

	for typeKinfo in dictTypeInfo[typeK]:
		host = typeKinfo[2]
		loc = typeKinfo[3]
		subtype = typeKinfo[4]
		year = typeKinfo[5]

		# if 'avian' in host.lower():
		#	 host = 'Avian'
		if host not in typeHostRange:
			typeHostRange.append(host)
		if loc not in typeLocRange:
			typeLocRange.append(loc)
		if year not in typeYearRange:
			typeYearRange.append(year)
		if subtype not in typeSubtypeRange:
			typeSubtypeRange.append(subtype)

		typeYearRange.sort()

	fwrite.write(typeYearRange[0]+' '+typeYearRange[-1]+'\t')

	typeInfo.append((typeYearRange[0],typeYearRange[-1]))

	yearStart = typeYearRange[0]

	dictTypeYearStart[strTypeIndex] = yearStart

	locStart = set()
	hostStart = set()
	for typeKinfo in dictTypeInfo[typeK]:
		loc = typeKinfo[3]
		host = typeKinfo[2]
		year = typeKinfo[5]
		if year == yearStart:
			locStart.add(loc)
			hostStart.add(host)   



	
	for host in typeHostRange:
		fwrite.write(host+' ')
	fwrite.write('\t')
	typeInfo.append(typeHostRange)

	for loc in typeLocRange:
		fwrite.write(loc+' ')
	fwrite.write('\t')
	typeInfo.append(typeLocRange)

	for loc in locStart:
		fwrite.write(loc+' ')
	fwrite.write('\t')
	
	for host in hostStart:
		fwrite.write(host+' ')
	fwrite.write('\t')

	for subtype in typeSubtypeRange:
		fwrite.write(subtype+' ')
	fwrite.write('\t')
	typeInfo.append(typeSubtypeRange)



	for typeKinfo in dictTypeInfo[typeK]:
		fwrite.write(typeKinfo[1]+'_'+typeKinfo[4]+'\t')
	fwrite.write('\n')

	dictTypeIndexInfo[strTypeIndex] = typeInfo

	typeIndex += 1
	

fwrite.close()




dictTypeDist = {}

for i in range(lenTypes):
	for j in range(i):
		typeI = lTypes[i]
		typeJ = lTypes[j]
		distIJ = distTypesIJ(typeI,typeJ)
		#print(distIJ)
		if distIJ < 8:
			strTypeI = strType(typeI)
			strTypeJ = strType(typeJ)
			typeIndexI = dictTypeIndex[typeI]
			typeIndexJ = dictTypeIndex[typeJ]
			dictTypeDist[(typeIndexI,typeIndexJ)] = distIJ
			dictTypeDist[(typeIndexJ,typeIndexI)] = distIJ
			
#			 fwrite.write(typeIndexI+'\t'+typeIndexJ+'\t'+str(distIJ)+'\n')

# fwrite.close()

dictVertexEdgeDist = {}
dictVertexEdgeDiff = {}

for (typeIndexI,typeIndexJ) in dictTypeDist:
	if typeIndexI in dictVertexEdgeDist:
		dictVertexEdgeDist[typeIndexI][(typeIndexI,typeIndexJ)] = dictTypeDist[(typeIndexI,typeIndexJ)]
	else:
		dictVertexEdgeDist[typeIndexI] = {}
		dictVertexEdgeDist[typeIndexI][(typeIndexI,typeIndexJ)] = dictTypeDist[(typeIndexI,typeIndexJ)]


def calTypeDiff(typeIndexI, typeIndexJ):
	typeI = dictIndexType[typeIndexI]
	typeJ = dictIndexType[typeIndexJ]

	strainsI = dictTypeInfo[typeI]
	strainsJ = dictTypeInfo[typeJ]

	minDiff = 999999

	for straini in strainsI:
		hosti = straini[2]
		loci = straini[3]
		yeari = int(straini[5])
		for strainj in strainsJ:
			hostj = strainj[2]
			locj = strainj[3]
			yearj = int(strainj[5])
			timeGap = abs(yeari-yearj)

			sumSame = 0

			if hosti == hostj:
				sumSame += 1
			elif 'avian' in hosti.lower() and 'avian' in hostj.lower():
				sumSame += 0.5
			if loci == locj:
				sumSame += 1
			if timeGap < 4:
				if timeGap < 2:
					sumSame += 1
				else:
					sumSame += 0.5
			
			sumDiff = 3-sumSame
			if minDiff > sumDiff:
				minDiff = sumDiff
	return minDiff


print('cal dictVertexEdgeDiff...')
dictTypeDiff = {}

for (typeIndexI,typeIndexJ) in dictTypeDist:
	typeDiff = calTypeDiff(typeIndexI, typeIndexJ)
	dictTypeDiff[(typeIndexI, typeIndexJ)] = typeDiff
	dictTypeDiff[(typeIndexJ, typeIndexI)] = typeDiff
	if typeIndexI in dictVertexEdgeDiff:
		dictVertexEdgeDiff[typeIndexI][(typeIndexI,typeIndexJ)] = typeDiff
	else:
		dictVertexEdgeDiff[typeIndexI] = {}
		dictVertexEdgeDiff[typeIndexI][(typeIndexI,typeIndexJ)] = typeDiff



ltypeIndexRankedByTime = []
for typeK in lTypesRankedByTime:
	ltypeIndexRankedByTime.append(dictTypeIndex[typeK])

fread = open('reassortmentHistroy3HostDetail.txt', 'r')
lines = fread.readlines()
fwrite = open('typeReassortmentCostHostDetail.txt', 'w')
fwrite.write('typeID\tReassortNum\tReassortCost\n')
dictReassortment = {}
for i in range(len(lines)):
	rline = lines[i].rstrip()
	if rline[0] == '>':
		typeIDTarget = rline[1:rline.find('|')]
		if typeIDTarget in dictReassortment:
			print('repeated typeIDTarget!')
			sys.exit()
		else:
			dictReassortment[typeIDTarget] = set()
		lrline = rline.split('|')
		reassortNum = lrline[1]
		reassortCost = lrline[2]
		fwrite.write(typeIDTarget+'\t'+reassortNum+'\t'+reassortCost+'\n')
		for j in range(int(reassortNum)):
			typeIDSourceStrain = lines[i+j+2].rstrip().split('\t')[0]
			# print(typeIDSourceStrain,dictStrainType[typeIDSourceStrain])
			typeIDSource = dictTypeIndex[dictStrainType[typeIDSourceStrain]]
			# print(typeIDSource)
			dictReassortment[typeIDTarget].add(typeIDSource)
	# else:
	# 	typeIDSource = rline[:rline.find('\t')]
	# 	dictReassortment[typeIDTarget].add(typeIDSource)
fwrite.close()
fread.close()






fwrite = open('typeNetworkReassortmentHostDetail.txt', 'w')

fwrite.write('typeSource\ttypeTarget\tdist\tdiff\n')

dictTypeDist = {}
print('output the typeNetworkReassortment network')

for typeIndexTarget in ltypeIndexRankedByTime:
	if typeIndexTarget in dictReassortment:
		for typeIndexSource in dictReassortment[typeIndexTarget]:
			typeSource = dictIndexType[typeIndexSource]
			typeTarget = dictIndexType[typeIndexTarget]
			dist = distTypesIJ(typeSource,typeTarget)
			diff = dictVertexEdgeDiff[typeIndexSource][(typeIndexSource,typeIndexTarget)]
			fwrite.write(typeIndexSource+'\t'+typeIndexTarget+'\t'+str(dist)+'\t'+str(diff)+'\n')
fwrite.close()




# def getSpecificReassortment(typeIndexSList, typeIndexT, yearStartTarget):
	

	# typeT = dictIndexType[typeIndexT]
	# strainsT = dictTypeInfo[typeT]
	
	# sumMinDiff = 999999
	# strainTarget = ''
	
	
	# for straint in strainsT:
		# hostt = straint[2]
		# loct = straint[3]
		# yeart = int(straint[5])
		# if yeart == yearStartTarget:
			# strainSList = []
			# sumMinDiffLocal = 0
			# for typeS in typeIndexSList:
				# minDiffStrain = 999999
				# strainsS = dictIndexType[typeS]
				# strainS = ''
				# for strains in strainsS:
					# hosts = strains[2]
					# locs = strains[3]
					# years = int(strains[5])
					# timeGap = abs(years-yeart)
					
					# sumSame = 0

					# if hosts == hostt:
						# sumSame += 1
					# elif 'avian' in hosts.lower() and 'avian' in hostt.lower():
						# sumSame += 0.5
					# if locs == loct:
						# sumSame += 1
					# if timeGap < 4:
						# if timeGap < 2:
							# sumSame += 1
						# else:
							# sumSame += 0.5
					
					# sumDiff = 3-sumSame
					# if minDiffStrain > sumDiff:
						# minDiffStrain = sumDiff
				# sumMinDiffLocal += minDiffStrain
			# if sumMinDiffLocal < sumMinDiff:
				# sumMinDiff = sumMinDiffLocal
				# strainTarget = straint


	# strainsS = dictTypeInfo[typeS]

	# typeS = dictIndexType[typeIndexS]

	# minDiff = 999999

	# for straini in strainsI:
		# hosti = straini[2]
		# loci = straini[3]
		# yeari = int(straini[5])
		# for strainj in strainsJ:
			# hostj = strainj[2]
			# locj = strainj[3]
			# yearj = int(strainj[5])
			# timeGap = abs(yeari-yearj)

			# sumSame = 0

			# if hosti == hostj:
				# sumSame += 1
			# elif 'avian' in hosti.lower() and 'avian' in hostj.lower():
				# sumSame += 0.5
			# if loci == locj:
				# sumSame += 1
			# if timeGap < 4:
				# if timeGap < 2:
					# sumSame += 1
				# else:
					# sumSame += 0.5
			
			# sumDiff = 3-sumSame
			# if minDiff > sumDiff:
				# minDiff = sumDiff
	# return minDiff

def getSpecificReassortment(typeIndexSList, typeIndexT, yearStartTarget):
	typeT = dictIndexType[typeIndexT]
	strainsT = dictTypeInfo[typeT]
	minDiffAll = 999999
	(t,s) = ('', [])
	for straint in strainsT:
		hostt = straint[2]
		loct = straint[3]
		yeart = int(straint[5])
		if yeart == yearStartTarget:
			sumMinDiffST = 0
			strainsSelectedList = []
			for typeIndexS in typeIndexSList:
				typeS = dictIndexType[typeIndexS]
				strainsS = dictTypeInfo[typeS]
				minDiff = 999999
				strainsSelected = ''
				for strains in strainsS:
					hosts = strains[2]
					locs = strains[3]
					years = int(strains[5])
					#if years < yeart+1:
					if years < yeart:
						timeGap = abs(years-yeart)
						sumSame = 0
						if hosts == hostt:
							sumSame += 1
						elif 'avian' in hosts.lower() and 'avian' in hostt.lower():
							sumSame += 0.5
						if locs == loct:
							sumSame += 1
						if timeGap < 4:
							if timeGap < 2:
								sumSame += 1
							else:
								sumSame += 0.5
						sumDiff = 3 - sumSame
						if minDiff > sumDiff:
							minDiff = sumDiff
							strainsSelected = strains
				strainsSelectedList.append((dictTypeIndex[typeS],strainsSelected))
				sumMinDiffST += minDiff
			if minDiffAll > sumMinDiffST:
				minDiffAll = sumMinDiffST
				(t, s) = (straint, strainsSelectedList)
	return (t, s, minDiffAll)
						

			# sumSame = 0

			# if hosti == hostj:
				# sumSame += 1
			# elif 'avian' in hosti.lower() and 'avian' in hostj.lower():
				# sumSame += 0.5
			# if loci == locj:
				# sumSame += 1
			# if timeGap < 4:
				# if timeGap < 2:
					# sumSame += 1
				# else:
					# sumSame += 0.5
			
			# sumDiff = 3-sumSame
			# if minDiff > sumDiff:
				# minDiff = sumDiff
						

segs = ['PB2','PB1','PA','HA','NP','NA','MP','NS']
						
def calStrainTypeSame(idS, idT, dictStrainTypeInfo):
	typeS = tuple(dictStrainTypeInfo[idS][7:])
	typeT = tuple(dictStrainTypeInfo[idT][7:])
	sameList = []
	for i in range(len(typeS)):
		if typeS[i] == typeT[i]:
			sameList.append(segs[i])
	return sameList
		
						
						
dictYearLocStrains = {}
dictStrainST = dict()
dictStrainTS = {}
dictCostStrainTarget = dict()						
						
	
fwrite = open('reassortmentHistroySpecificHostDetail.txt', 'w')

for typeIndexTarget in ltypeIndexRankedByTime:
	if typeIndexTarget in dictReassortment:
		typeIndexSList = list(dictReassortment[typeIndexTarget])
		typeTargetInfo = dictTypeIndexInfo[typeIndexTarget]
		yearStartTarget = int(typeTargetInfo[0][0])
		#print(typeIndexTarget, yearStartTarget)
		(targetStrain, strainListSource, minDiffAll) = getSpecificReassortment(typeIndexSList, typeIndexTarget, yearStartTarget)		
		fwrite.write('>'+typeIndexTarget+'|'+str(minDiffAll)+'|')
		dictCostStrainTarget[targetStrain[0]] = minDiffAll
		for item in targetStrain:
			fwrite.write(item+'\t')
		fwrite.write('\n')
		for (types, strains) in strainListSource:
			fwrite.write(types+'|')
			for item in strains:
				fwrite.write(item+'\t')
			fwrite.write('\n')	
		#EPI_ISL_179403	A/Ohio/09/2015	Human	North America	H1N1	2015	2015-04-21	
		idT = targetStrain[0]
		hostT = targetStrain[2]
		locT = targetStrain[3]
		yearT = int(targetStrain[5])
		subtypeT = targetStrain[4]
		
		if yearT in dictYearLocStrains:
			if locT in dictYearLocStrains[yearT]:
				dictYearLocStrains[yearT][locT].add(idT)
			else:
				dictYearLocStrains[yearT][locT] = set()
				dictYearLocStrains[yearT][locT].add(idT)
		else:
			dictYearLocStrains[yearT] = {}
			dictYearLocStrains[yearT][locT] = set()
			dictYearLocStrains[yearT][locT].add(idT)
		
		for (types, strains) in strainListSource:
			idS = strains[0]
			hostS = strains[2]
			locS = strains[3]
			yearS = int(strains[5])
			subtypeS = strains[4]
			
			dictStrainST[(idS, idT)] = calStrainTypeSame(idS, idT, dictStrainTypeInfo)
			if idT in dictStrainTS:
				dictStrainTS[idT].add(idS)
			else:
				dictStrainTS[idT] = set()
				dictStrainTS[idT].add(idS)
			
			if yearS in dictYearLocStrains:
				if locS in dictYearLocStrains[yearS]:
					dictYearLocStrains[yearS][locS].add(idS)
				else:
					dictYearLocStrains[yearS][locS] = set()
					dictYearLocStrains[yearS][locS].add(idS)
			else:
				dictYearLocStrains[yearS] = {}
				dictYearLocStrains[yearS][locS] = set()
				dictYearLocStrains[yearS][locS].add(idS)
				
		
	
fwrite.close()	

# dictYearLocTypeIndex = dict()
# for typeIndex in ltypeIndexRankedByTime:
	

#locListOld = ["East Africa","North Africa","South Africa","West Africa","Middle Africa","Central Asia","East Asia","South Asia","Southeast Asia","West Asia","Central Europe","Eastern Europe","Northern Europe","Southern Europe","Western Europe","Middle America","North America","The Caribbean","Oceania","Eastern South America","Midwest South American","Northern South America","Southern South America"]
# locList = ["East Africa","North Africa","South Africa","West Africa","Middle Africa","Central Europe","Eastern Europe","Northern Europe","Southern Europe","Western Europe","Central Asia","East Asia","South Asia","Southeast Asia","West Asia","Middle America","North America","The Caribbean","Eastern South America","Midwest South American","Northern South America","Southern South America","Oceania"]
locList = ["East Africa","North Africa","South Africa","West Africa","Central Asia","East Asia","South Asia","Southeast Asia","West Asia","Central Europe","Eastern Europe","Northern Europe","Southern Europe","Western Europe","Middle America","North America","The Caribbean","Oceania","Eastern South America","Midwest South American","Northern South America","Southern South America"]
colorJetList = ['#00008F','#00009F','#0000AF','#0000BF','#0000CF','#0000DF','#0000EF','#0000FF','#000FFF','#001FFF','#002FFF','#003FFF','#004FFF','#005FFF','#006FFF','#007FFF','#008FFF','#009FFF','#00AFFF','#00BFFF','#00CFFF','#00DFFF','#00EFFF','#00FFFF','#0FFFEF','#1FFFDF','#2FFFCF','#3FFFBF','#4FFFAF','#5FFF9F','#6FFF8F','#7FFF7F','#8FFF6F','#9FFF5F','#AFFF4F','#BFFF3F','#CFFF2F','#DFFF1F','#EFFF0F','#FFFF00','#FFEF00','#FFDF00','#FFCF00','#FFBF00','#FFAF00','#FF9F00','#FF8F00','#FF7F00','#FF6F00','#FF5F00','#FF4F00','#FF3F00','#FF2F00','#FF1F00','#FF0F00','#FF0000','#EF0000','#DF0000','#CF0000','#BF0000','#AF0000','#9F0000','#8F0000','#7F0000']

yearList = list(dictYearLocStrains.keys())
yearList.sort()

print('begin to plot the network...')

colorList = [green, yellow, blue, red]

dictHostColor = dict()
hosts = ['Avian','Equine','Swine','Human']
hostsDetail = ['avian_shorebirds', 'avian_land birds', 'avian_waterfowl', 'avian_domestic birds', 'avian_unknown', 'equine', 'swine','human']
#colorListAvian = [seagreen, forestgreen, limegreen, green, lawngreen, yellow, blue, red]
colorListAvian = [pink,wheat,green,purple,springgreen, yellow, blue, red]
i = 0
for host in hostsDetail:
	dictHostColor[host] = colorListAvian[i]
	i = i + 1


dictLocColor = {}

for i in range(len(locList)):
	loc = locList[i]
	colorIndex = int(64/len(locList)*i)%64
	# print(colorIndex)
	colorHere = colorJetList[colorIndex]
	dictLocColor[loc] = colorHere

(pageWidth, pageLength) = (2970, 2970)


dictLocMaxNumStrains = {}
for year in dictYearLocStrains:
	for loc in dictYearLocStrains[year]:
		numStrains = len(dictYearLocStrains[year][loc])
		if loc in dictLocMaxNumStrains:
			if numStrains > dictLocMaxNumStrains[loc]:
				dictLocMaxNumStrains[loc] = numStrains
		else:
			dictLocMaxNumStrains[loc] = numStrains
	
yStep = (pageLength/(yearList[-1]-yearList[0]))*0.8

xStep = (6*pageWidth/sum(dictLocMaxNumStrains.values()))*0.8

dictLocStartX = {}
dictLocEndX = {}

xStart = 0.2*pageWidth
locGap = xStep*4

posIndex = xStart + locGap



for loc in locList:
	if loc in dictLocMaxNumStrains:
		dictLocStartX[loc] = posIndex
		posRange = dictLocMaxNumStrains[loc]
		posIndex = posIndex + posRange*xStep + locGap
		dictLocEndX[loc] = posIndex-locGap


# def mycmp(strainIDx,strainIDy):
	# subtypex = dictStrainTypeInfo[strainIDx][4]
	# subtypey = dictStrainTypeInfo[strainIDy][4]
	# return (subtypex > subtypey)
		
def orderStrains(strains):
	#strains.sort(cmp = mycmp)
	strains.sort(key=cmp_to_key(lambda a, b: dictStrainTypeInfo[a][4] < dictStrainTypeInfo[b][4]))

	return strains

		
		

c = canvas.Canvas('recombinationStrainsRectHostDetail.pdf', pagesize=(pageWidth*7, pageLength*6))


grayTransparent = Color(0, 0, 0, alpha=0.3)
greenTransparent = Color(0, 255, 0, alpha=0.3)
redTransparent = Color(255, 0, 0, alpha=0.3)

r = xStep * 0.3

yStart = 6 * pageLength * 0.9

dictStrainPos = {}



print('cal the pos of the strains...')
lineWidth = 0.5
c.setLineWidth(lineWidth)
dictStrainRadius = {}

for year in yearList:
	#yPos = yStart - (year - yearList[0])*yStep
	yPos = yStart - (year - yearList[0])*yStep*((1.03)**(year - yearList[0]))
	for loc in locList:
		if loc in dictLocMaxNumStrains:
			if year in dictYearLocStrains and loc in dictYearLocStrains[year]:
				strains = dictYearLocStrains[year][loc]
				xPosStart = dictLocStartX[loc]
				
				strainsOrdered = orderStrains(list(strains))
				
				xrange = dictLocMaxNumStrains[loc]
				numStrains = len(strainsOrdered)
				xStepHere = (xrange/(numStrains+1))*xStep
				
				r = min(xStepHere*0.4, yStep*0.8)
				#xPosStartHere = xPosStart + ((xrange - numStrains)/2)*xStep
				xPosStartHere = xPosStart + xStepHere
				xIndex = 0
				for strainX in strainsOrdered:
					dictStrainRadius[strainX] = r
					host = dictStrainTypeInfo[strainX][2]
					xPos = xPosStartHere + xIndex*xStepHere
					dictStrainPos[strainX] = (xPos, yPos)
					xIndex = xIndex + 1
					if host in dictHostColor:
						colorHere = dictHostColor[host]
					else:
						colorHere = gray
					c.setStrokeColor(colorHere)
					c.setFillColor(colorHere)
					c.circle(xPos, yPos, r, fill=1)

lineWidth25 = 0.25
lineWidth10 = 0.1
lineWidth = 0.5
c.setLineWidth(lineWidth)

radiusArrow = xStep/10
					
for (sourceStrain, targetStrain) in dictStrainST:
	locS = dictStrainTypeInfo[sourceStrain][3]
	locT = dictStrainTypeInfo[targetStrain][3]
	yearS = dictStrainTypeInfo[sourceStrain][5]
	yearT = dictStrainTypeInfo[targetStrain][5]
	
	
	if locS == locT:
		c.setStrokeColor(greenTransparent)
		c.setFillColor(green)
	else:
		c.setStrokeColor(redTransparent)
		c.setFillColor(red)
	
	x1 = dictStrainPos[sourceStrain][0]
	y1 = dictStrainPos[sourceStrain][1]
	x2 = dictStrainPos[targetStrain][0]
	y2 = dictStrainPos[targetStrain][1]
	
	if yearS == yearT:
		print('This is impossible...')
		xDist = (x2-x1)/3
		x121 = x1+xDist
		y121 = y1-xStep
		x122 = x1+2*xDist
		y122 = y1-xStep
	
		p = c.beginPath()
		p.moveTo(x1,y1)
		p.curveTo(x121,y121,x122,y122,x2,y2)
		c.drawPath(p, stroke=1, fill=0)
	else:
		
		radiusT = dictStrainRadius[targetStrain]
		distAT = radiusArrow+radiusT+2*lineWidth
		distST = math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
		xa = x2-(distAT/distST)*(x2-x1)
		ya = y2+(distAT/distST)*(y1-y2)
		c.setLineWidth(lineWidth25)
		c.line(x1,y1,xa,ya)
		if locS == locT:
			c.setStrokeColor(green)
		else:
			c.setStrokeColor(red)
			
		c.setLineWidth(lineWidth10)
		c.circle(xa,ya,radiusArrow,fill=1)
		c.setFillColor(black)
		c.circle(x1,y1,radiusArrow,fill=1)

print('redraw the circles...')
lineWidth = 0.5
fontSize = 100
c.setLineWidth(lineWidth)


for year in range(yearList[0], yearList[-1]+1):
	yPos = yStart - (year - yearList[0])*yStep*((1.03)**(year - yearList[0]))
	
	c.setFillColor(black)
	c.setFont("Helvetica", fontSize)
	c.setLineWidth(1)
	c.setStrokeColor(black)
	if year % 5 == 0:
		c.drawString(xStart + (-18)*xStep, yPos-yStep*0.3, str(year))
		c.line(xStart + (-3)*xStep, yPos, xStart + xStep, yPos)
	else:
		c.line(xStart + (-3)*xStep, yPos, xStart + (-1)*xStep, yPos)
	#c.setLineWidth(lineWidth)
	c.setLineWidth(lineWidth10)
	
	for loc in locList:
		if loc in dictLocMaxNumStrains:
			if year in dictYearLocStrains and loc in dictYearLocStrains[year]:
				strains = dictYearLocStrains[year][loc]
				xPosStart = dictLocStartX[loc]
				
				strainsOrdered = orderStrains(list(strains))
				
				xrange = dictLocMaxNumStrains[loc]
				numStrains = len(strainsOrdered)
				xStepHere = (xrange/(numStrains+1))*xStep
				
				r = min(xStepHere*0.4, yStep*0.8)
				#xPosStartHere = xPosStart + ((xrange - numStrains)/2)*xStep
				xPosStartHere = xPosStart + xStepHere
				xIndex = 0
				for strainX in strainsOrdered:
					dictStrainRadius[strainX] = r
					host = dictStrainTypeInfo[strainX][2]
					xPos = xPosStartHere + xIndex*xStepHere
					dictStrainPos[strainX] = (xPos, yPos)
					xIndex = xIndex + 1
					if host in dictHostColor:
						colorHere = dictHostColor[host]
					else:
						colorHere = gray
					c.setStrokeColor(colorHere)
					c.setFillColor(colorHere)
					c.circle(xPos, yPos, r, fill=0)

c.setLineWidth(1)
c.setStrokeColor(black)
c.line(xStart + (-3)*xStep, yPos-yStep, xStart + (-3)*xStep, yStart+yStep)

yPosLocLegend = yPos - yStep * 3
lineWidth = 1.0
c.setStrokeColor(black)
c.setLineWidth(lineWidth)
c.setFillColor(black)
for loc in locList:
	if loc in dictLocMaxNumStrains:
		xPosStart = dictLocStartX[loc]
		xPosEnd = dictLocEndX[loc]
		c.setStrokeColor(black)
		c.setLineWidth(lineWidth)
		c.line(xPosStart, yPosLocLegend, xPosEnd, yPosLocLegend)
		c.setStrokeColor(gray)
		c.setLineWidth(lineWidth/4)
		c.line(xPosEnd+locGap/2, yPosLocLegend, xPosEnd+locGap/2, yStart)
		
		c.saveState()
		# c.translate(posXString, posYString)
		# c.rotate(dictNodePos[clusterNode][1]*360/(2*math.pi))
		c.setFont("Helvetica", fontSize)
		# c.setFillColor(black)
		
		posXString = (xPosStart+xPosEnd)/2
		posYString = yPosLocLegend-yStep/2
		
		c.translate(posXString, posYString)
		c.rotate(270)
		
		
		c.drawString(yStep, 0, loc)
		
		c.restoreState()
		
	
xPosHostLegend = xPosEnd-yStep*5
yPosHostLegend = yStart
r = yStep*0.8
lineWidth = 0.5
c.setLineWidth(lineWidth)
c.setFont("Helvetica", fontSize)
for host in hostsDetail:
	colorHere = dictHostColor[host]
	c.setFillColor(colorHere)
	c.circle(xPosHostLegend, yPosHostLegend, r, fill=1)
	
	c.setFillColor(black)
	c.drawString(xPosHostLegend+r*2, yPosHostLegend,host)
	
	yPosHostLegend = yPosHostLegend - 2*yStep


c.showPage()
c.save()

dictLocYearReassortNum = {}
dictLocYearLongReassortNum = {}
dictLocYearHostDiffReassortNum = {}

dictYearHLsHLt = {}


for year in dictYearLocStrains:
	for loc in dictYearLocStrains[year]:
		strains = dictYearLocStrains[year][loc]
		for strain in strains:
			if strain in dictStrainTS:
				if loc in dictLocYearReassortNum:
					if year in dictLocYearReassortNum[loc]:
						dictLocYearReassortNum[loc][year] = dictLocYearReassortNum[loc][year] + 1
					else:
						dictLocYearReassortNum[loc][year] = 1
				else:
					dictLocYearReassortNum[loc] = {}
					dictLocYearReassortNum[loc][year] = 1
				strainT = strain
				strainS = dictStrainTS[strainT]
				locST = set()
				hostST = set()
				#strainS.add(strainT)
				ss = set()
				ss.add(strainT)
				for strainTemp in strainS:
					ss.add(strainTemp)
				HLs = []
				HLt = (dictStrainTypeInfo[strainT][2],dictStrainTypeInfo[strainT][3],dictStrainTypeInfo[strainT][4])
				
				for strain in strainS:
					host = dictStrainTypeInfo[strain][2]
					locLocal = dictStrainTypeInfo[strain][3]
					subtype = dictStrainTypeInfo[strain][4]
					HLs.append((host, locLocal,subtype))
				HLs.sort()
				HLs = tuple(HLs)
				
				if year in dictYearHLsHLt:
					if (HLs, HLt) in dictYearHLsHLt[year]:
						dictYearHLsHLt[year][(HLs, HLt)] = dictYearHLsHLt[year][(HLs, HLt)] + 1
					else:
						dictYearHLsHLt[year][(HLs, HLt)] = 1
				else:
					dictYearHLsHLt[year] = {}
					dictYearHLsHLt[year][(HLs, HLt)] = 1			
				
				
				for strain in ss:
					host = dictStrainTypeInfo[strain][2]
					locLocal = dictStrainTypeInfo[strain][3]
					hostST.add(host)
					locST.add(locLocal)
				if len(hostST) > 1:
					if loc in dictLocYearHostDiffReassortNum:
						if year in dictLocYearHostDiffReassortNum[loc]:
							dictLocYearHostDiffReassortNum[loc][year] = dictLocYearHostDiffReassortNum[loc][year] + 1
						else:
							dictLocYearHostDiffReassortNum[loc][year] = 1
					else:
						dictLocYearHostDiffReassortNum[loc] = {}
						dictLocYearHostDiffReassortNum[loc][year] = 1
					
				if len(locST) > 1:
					if loc in dictLocYearLongReassortNum:
						if year in dictLocYearLongReassortNum[loc]:
							dictLocYearLongReassortNum[loc][year] = dictLocYearLongReassortNum[loc][year] + 1
						else:
							dictLocYearLongReassortNum[loc][year] = 1
					else:
						dictLocYearLongReassortNum[loc] = {}
						dictLocYearLongReassortNum[loc][year] = 1
					
					
fwrite = open('locYearReassortNumHostDetail.txt', 'w')
fwrite.write('locs\t')
for year in yearList:
	fwrite.write(str(year)+'all\t'+str(year)+'longDist\t'+str(year)+'diffHost\t')
fwrite.write('\n')
for loc in locList:
	if loc in dictLocMaxNumStrains:
		fwrite.write(loc+'\t')
		for year in yearList:
			if loc in dictLocYearReassortNum:
				if year in dictLocYearReassortNum[loc]:
					fwrite.write(str(dictLocYearReassortNum[loc][year])+'\t')
				else:
					fwrite.write('0\t')
			else:
				#print(loc, 'not in dictLocYearReassortNum')
				fwrite.write('0\t')
			if loc in dictLocYearLongReassortNum:
				if year in dictLocYearLongReassortNum[loc]:
					fwrite.write(str(dictLocYearLongReassortNum[loc][year])+'\t')
				else:
					fwrite.write('0\t')
			else:
				fwrite.write('0\t')
			
			if loc in dictLocYearHostDiffReassortNum:
				if year in dictLocYearHostDiffReassortNum[loc]:
					fwrite.write(str(dictLocYearHostDiffReassortNum[loc][year])+'\t')
				else:
					fwrite.write('0\t')
			else:
				fwrite.write('0\t')
		fwrite.write('\n')
				



fwrite.close()


(pageWidth, pageLength) = (2970*2, 2970)

c = canvas.Canvas('reassortNumHostDetail.pdf', pagesize=(pageWidth, pageLength))
stepLoc = pageLength/400
stepYear = pageWidth/100

bottomY = pageLength/20
bottomX = pageWidth/10

barWidth = stepYear/5

yearStart = 1952

for year in yearList:
	xPosYear = bottomX + (year-yearStart)*stepYear
	heightAcc = bottomY
	heightAccLong = bottomY
	heightAccHostDiff = bottomY
	for loc in locList:
		if loc in dictLocYearReassortNum and year in dictLocYearReassortNum[loc]:
			reassortNum = dictLocYearReassortNum[loc][year]
			colorHere = dictLocColor[loc]
			c.setFillColor(colorHere)
			c.setStrokeColor(colorHere)
			c.rect(xPosYear, heightAcc, barWidth, reassortNum*stepLoc, fill=1)
			heightAcc = heightAcc + reassortNum*stepLoc
		if loc in dictLocYearLongReassortNum and year in dictLocYearLongReassortNum[loc]:
			reassortNumLong = dictLocYearLongReassortNum[loc][year]
			colorHere = dictLocColor[loc]
			c.setFillColor(colorHere)
			c.setStrokeColor(colorHere)
			c.rect(xPosYear+barWidth*1.1, heightAccLong, barWidth, reassortNumLong*stepLoc, fill=1)
			heightAccLong = heightAccLong + reassortNumLong*stepLoc
		if loc in dictLocYearHostDiffReassortNum and year in dictLocYearHostDiffReassortNum[loc]:
			reassortNumHostDiff = dictLocYearHostDiffReassortNum[loc][year]
			colorHere = dictLocColor[loc]
			c.setFillColor(colorHere)
			c.setStrokeColor(colorHere)
			c.rect(xPosYear+barWidth*2.2, heightAccHostDiff, barWidth, reassortNumHostDiff*stepLoc, fill=1)
			heightAccHostDiff = heightAccHostDiff + reassortNumHostDiff*stepLoc

			
c.setFillColor(black)
c.setStrokeColor(black)

c.line(bottomX, bottomY, xPosYear+stepYear, bottomY)
for year in yearList:
	xPosYear = bottomX + (year-yearStart)*stepYear
	xBarPos = xPosYear + 1.5*barWidth
	c.line(xBarPos, bottomY, xBarPos, bottomY - barWidth)
	c.drawString(xPosYear, bottomY - 2*barWidth, str(year))

c.line(bottomX, bottomY, bottomX, bottomY+300*stepLoc)
for i in range(320):
	if i % 50 == 0:
		yBarPos = bottomY+i*stepLoc
		c.line(bottomX, yBarPos, bottomX+barWidth, yBarPos)
		c.drawString(bottomX-4*barWidth, yBarPos, str(i))

posXLegend = bottomX + 6*barWidth
edgeLength = 3*barWidth
posYLegend = bottomY + 10*stepLoc
gapLegend = edgeLength+barWidth


for loc in locList:
	colorHere = dictLocColor[loc]
	c.setFillColor(colorHere)
	c.setStrokeColor(colorHere)
	c.rect(posXLegend, posYLegend, edgeLength, edgeLength, fill=1)
	c.setFillColor(black)
	c.drawString(posXLegend+edgeLength+barWidth, posYLegend, loc)
	posYLegend = posYLegend + edgeLength + barWidth
	

		
c.showPage()
c.save()

fwrite = open('hostLocYearSubtypeReassortHostDetail.txt', 'w')
for year in dictYearHLsHLt:
	for (HLs, HLt) in dictYearHLsHLt[year]:
		fwrite.write(str(year)+'\t')
		for (hosts, locs, subtype) in HLs:
			fwrite.write(hosts+'|'+locs+'|'+subtype+';')
		fwrite.write('\t')
		fwrite.write(HLt[0]+'|'+HLt[1]+'|'+HLt[2]+'\t')
		fwrite.write(str(dictYearHLsHLt[year][(HLs, HLt)]))
		fwrite.write('\n')


fwrite.close()


modeDir = os.path.join('.','specificReassortmentModeHostDetail')

dictYearReassortModeHost = dict()
dictYearReassortModeLoc = dict()
dictYearReassortModeSubtype = dict()
dictReassortModeHost = dict()
dictReassortModeLoc = dict()
dictReassortModeSubtype = dict()

dictReassortModeHostLoc = dict()
dictReassortModeHostSubtype = dict()
dictReassortModeLocSubtype = dict()

for year in dictYearHLsHLt:
	for (HLs, HLt) in dictYearHLsHLt[year]:
		# fwrite.write(str(year)+'\t')
		setHost = set()
		setLoc = set()
		setSubtype = set()
		
		for (hosts, locs, subtype) in HLs:
			setHost.add(hosts)
			setLoc.add(locs)
			setSubtype.add(subtype)
		hostt = HLt[0]
		loct = HLt[1]
		subtypet = HLt[2]
		sourceHostTuple = tuple(setHost)
		sourceLocTuple = tuple(setLoc)
		sourceSubtypeTuple = tuple(setSubtype)
		
		modeHostTuple = (sourceHostTuple, hostt)
		modeLocTuple = (sourceLocTuple, loct)
		modeSubtypeTuple = (sourceSubtypeTuple, subtypet)
		
		if modeHostTuple in dictReassortModeHost:
			dictReassortModeHost[modeHostTuple] = dictReassortModeHost[modeHostTuple] + 1
		else:
			dictReassortModeHost[modeHostTuple] = 1
		if modeLocTuple in dictReassortModeLoc:
			dictReassortModeLoc[modeLocTuple] = dictReassortModeLoc[modeLocTuple] + 1
		else:
			dictReassortModeLoc[modeLocTuple] = 1
		if modeSubtypeTuple in dictReassortModeSubtype:
			dictReassortModeSubtype[modeSubtypeTuple] = dictReassortModeSubtype[modeSubtypeTuple] + 1
		else:
			dictReassortModeSubtype[modeSubtypeTuple] = 1
		
		if modeHostTuple in dictReassortModeHostLoc:
			if modeLocTuple in dictReassortModeHostLoc[modeHostTuple]:
				dictReassortModeHostLoc[modeHostTuple][modeLocTuple] = dictReassortModeHostLoc[modeHostTuple][modeLocTuple] + 1
			else:
				dictReassortModeHostLoc[modeHostTuple][modeLocTuple] = 1
		else:
			dictReassortModeHostLoc[modeHostTuple] = {}
			dictReassortModeHostLoc[modeHostTuple][modeLocTuple] = 1
		
		if modeHostTuple in dictReassortModeHostSubtype:
			if modeSubtypeTuple in dictReassortModeHostSubtype[modeHostTuple]:
				dictReassortModeHostSubtype[modeHostTuple][modeSubtypeTuple] = dictReassortModeHostSubtype[modeHostTuple][modeSubtypeTuple] + 1
			else:
				dictReassortModeHostSubtype[modeHostTuple][modeSubtypeTuple] = 1
		else:
			dictReassortModeHostSubtype[modeHostTuple] = {}
			dictReassortModeHostSubtype[modeHostTuple][modeSubtypeTuple] = 1
			
		
		
		if modeLocTuple in dictReassortModeLocSubtype:
			if modeSubtypeTuple in dictReassortModeLocSubtype[modeLocTuple]:
				dictReassortModeLocSubtype[modeLocTuple][modeSubtypeTuple] = dictReassortModeLocSubtype[modeLocTuple][modeSubtypeTuple] + 1
			else:
				dictReassortModeLocSubtype[modeLocTuple][modeSubtypeTuple] = 1
		else:
			dictReassortModeLocSubtype[modeLocTuple] = {}
			dictReassortModeLocSubtype[modeLocTuple][modeSubtypeTuple] = 1
		
		if year in dictYearReassortModeHost:
			if modeHostTuple in dictYearReassortModeHost[year]:
				dictYearReassortModeHost[year][modeHostTuple] = dictYearReassortModeHost[year][modeHostTuple] + 1
			else:
				dictYearReassortModeHost[year][modeHostTuple] = 1
		else:
			dictYearReassortModeHost[year] = {}
			dictYearReassortModeHost[year][modeHostTuple] = 1
		
		if year in dictYearReassortModeLoc:
			if modeLocTuple in dictYearReassortModeLoc[year]:
				dictYearReassortModeLoc[year][modeLocTuple] = dictYearReassortModeLoc[year][modeLocTuple] + 1
			else:
				dictYearReassortModeLoc[year][modeLocTuple] = 1
		else:
			dictYearReassortModeLoc[year] = {}
			dictYearReassortModeLoc[year][modeLocTuple] = 1
		
		if year in dictYearReassortModeSubtype:
			if modeSubtypeTuple in dictYearReassortModeSubtype[year]:
				dictYearReassortModeSubtype[year][modeSubtypeTuple] = dictYearReassortModeSubtype[year][modeSubtypeTuple] + 1
			else:
				dictYearReassortModeSubtype[year][modeSubtypeTuple] = 1
		else:
			dictYearReassortModeSubtype[year] = {}
			dictYearReassortModeSubtype[year][modeSubtypeTuple] = 1
		
		
			# fwrite.write(hosts+'|'+locs+'|'+subtype+';')
		# fwrite.write('\t')
		# fwrite.write(HLt[0]+'|'+HLt[1]+'|'+HLt[2]+'\t')
		# fwrite.write(str(dictYearHLsHLt[year][(HLs, HLt)]))
		# fwrite.write('\n')

fwrite = open(os.path.join(modeDir,'yearReassortModeHostHostDetail.txt'), 'w')
for year in dictYearReassortModeHost:
	for (hosts,hostt) in dictYearReassortModeHost[year]:
		fwrite.write(str(year)+'\t')
		fwrite.write(','.join(list(hosts)))
		fwrite.write('--')
		fwrite.write(hostt+'\t')
		fwrite.write(str(dictYearReassortModeHost[year][(hosts,hostt)]))
		fwrite.write('\n')
fwrite.close()

fwrite = open(os.path.join(modeDir,'yearReassortModeLocHostDetail.txt'), 'w')
for year in dictYearReassortModeLoc:
	for (Locs,Loct) in dictYearReassortModeLoc[year]:
		fwrite.write(str(year)+'\t')
		fwrite.write(','.join(list(Locs)))
		fwrite.write('--')
		fwrite.write(Loct+'\t')
		fwrite.write(str(dictYearReassortModeLoc[year][(Locs,Loct)]))
		fwrite.write('\n')
fwrite.close()

fwrite = open(os.path.join(modeDir,'yearReassortModeSubtypeHostDetail.txt'), 'w')
for year in dictYearReassortModeSubtype:
	for (Subtypes,Subtypet) in dictYearReassortModeSubtype[year]:
		fwrite.write(str(year)+'\t')
		fwrite.write(','.join(list(Subtypes)))
		fwrite.write('--')
		fwrite.write(Subtypet+'\t')
		fwrite.write(str(dictYearReassortModeSubtype[year][(Subtypes,Subtypet)]))
		fwrite.write('\n')
fwrite.close()


fwrite = open(os.path.join(modeDir,'reassortModeHostHostDetail.txt'), 'w')
for (Hosts, Hostt) in dictReassortModeHost:
	fwrite.write(','.join(list(Hosts)))
	fwrite.write('--')
	fwrite.write(Hostt+'\t')
	fwrite.write(str(dictReassortModeHost[(Hosts,Hostt)]))
	fwrite.write('\n')
fwrite.close()

fwrite = open(os.path.join(modeDir,'reassortModeLocHostDetail.txt'), 'w')
for (Locs, Loct) in dictReassortModeLoc:
	fwrite.write(','.join(list(Locs)))
	fwrite.write('--')
	fwrite.write(Loct+'\t')
	fwrite.write(str(dictReassortModeLoc[(Locs,Loct)]))
	fwrite.write('\n')
fwrite.close()


fwrite = open(os.path.join(modeDir,'reassortModeSubtypeHostDetail.txt'), 'w')
for (Subtypes, Subtypet) in dictReassortModeSubtype:
	fwrite.write(','.join(list(Subtypes)))
	fwrite.write('--')
	fwrite.write(Subtypet+'\t')
	fwrite.write(str(dictReassortModeSubtype[(Subtypes,Subtypet)]))
	fwrite.write('\n')
fwrite.close()

# dictReassortModeHo-stLoc
# dictReassortModeHostSubtype

fwrite = open(os.path.join(modeDir,'reassortModeHostLocHostDetail.txt'), 'w')
for (Hosts, Hostt) in dictReassortModeHostLoc:
	for (Locs, Loct) in dictReassortModeHostLoc[(Hosts, Hostt)]:
		fwrite.write(','.join(list(Hosts)))
		fwrite.write('--')
		fwrite.write(Hostt+'\t')
		fwrite.write(','.join(list(Locs)))
		fwrite.write('--')
		fwrite.write(Loct+'\t')
		fwrite.write(str(dictReassortModeHostLoc[(Hosts, Hostt)][(Locs, Loct)]))
		fwrite.write('\n')
fwrite.close()

fwrite = open(os.path.join(modeDir,'reassortModeHostSubtypeHostDetail.txt'), 'w')
for (Hosts, Hostt) in dictReassortModeHostSubtype:
	for (Subtypes, Subtypet) in dictReassortModeHostSubtype[(Hosts, Hostt)]:
		fwrite.write(','.join(list(Hosts)))
		fwrite.write('--')
		fwrite.write(Hostt+'\t')
		fwrite.write(','.join(list(Subtypes)))
		fwrite.write('--')
		fwrite.write(Subtypet+'\t')
		fwrite.write(str(dictReassortModeHostSubtype[(Hosts, Hostt)][(Subtypes, Subtypet)]))
		fwrite.write('\n')
fwrite.close()


fwrite = open(os.path.join(modeDir,'reassortModeLocSubtypeHostDetail.txt'), 'w')
for (Locs, Loct) in dictReassortModeLocSubtype:
	for (Subtypes, Subtypet) in dictReassortModeLocSubtype[(Locs, Loct)]:
		fwrite.write(','.join(list(Locs)))
		fwrite.write('--')
		fwrite.write(Loct+'\t')
		fwrite.write(','.join(list(Subtypes)))
		fwrite.write('--')
		fwrite.write(Subtypet+'\t')
		fwrite.write(str(dictReassortModeLocSubtype[(Locs, Loct)][(Subtypes, Subtypet)]))
		fwrite.write('\n')
fwrite.close()





#接下来需要将关系进行筛选，因为有些重配应该是不太合理的，我觉得。。。




