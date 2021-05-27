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


fwrite = open('typeIndexInfoAvianDetailFinalHostDetal2.txt', 'w')
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
fwrite = open('typeReassortmentCostHostDetail2.txt', 'w')
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
			typeIDSourceStrain = lines[i + j + 2].rstrip().split('\t')[0]
			# print(typeIDSourceStrain,dictStrainType[typeIDSourceStrain])
			typeIDSource = dictTypeIndex[dictStrainType[typeIDSourceStrain]]
			# print(typeIDSource)
			dictReassortment[typeIDTarget].add(typeIDSource)
# else:
# 	typeIDSource = rline[:rline.find('\t')]
# 	dictReassortment[typeIDTarget].add(typeIDSource)
fwrite.close()
fread.close()






fwrite = open('typeNetworkReassortmentHostDetail2.txt', 'w')

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
						

segs = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']						
						
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
						
	
fwrite = open('reassortmentHistroySpecific4HostDetail.txt', 'w')

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

xCircleCenter = pageWidth * 3
yCircleCenter = pageLength * 3



maxDistanceFromRoot = (yearList[-1]-yearList[0])*1.3

dRadius = (pageWidth/2)*0.95/maxDistanceFromRoot/2
yearStart = yearList[0] - (yearList[-1]-yearList[0])*0.3


dAngle = (2*math.pi)*0.85/sum(dictLocMaxNumStrains.values())
			
yStep = (pageLength/(yearList[-1]-yearList[0]))*0.8

xStep = (pageWidth/sum(dictLocMaxNumStrains.values()))*0.8

dictLocStartAngle = {}
dictLocEndAngle = {}

angleStart = (2*math.pi)*0.03
locGap = dAngle*2

angleIndex = angleStart + locGap



for loc in locList:
	if loc in dictLocMaxNumStrains:
		dictLocStartAngle[loc] = angleIndex
		posRange = dictLocMaxNumStrains[loc]
		angleIndex = angleIndex + posRange*dAngle + locGap
		dictLocEndAngle[loc] = angleIndex-locGap


# def mycmp(strainIDx,strainIDy):
	# subtypex = dictStrainTypeInfo[strainIDx][4]
	# subtypey = dictStrainTypeInfo[strainIDy][4]
	# return (subtypex > subtypey)
		
def orderStrains(strains):
	#strains.sort(cmp = mycmp)
	strains.sort(key=cmp_to_key(lambda a, b: dictStrainTypeInfo[a][4] < dictStrainTypeInfo[b][4]))

	return strains

def getPosXY(xCircleCenter, yCircleCenter, radius, angle):
	xPosition = xCircleCenter + radius * math.cos(angle)
	yPosition = yCircleCenter + radius * math.sin(angle)
	return (xPosition, yPosition)		
		

c = canvas.Canvas('recombinationStrainsPolarHostDetail.pdf', pagesize=(6*pageWidth, 6*pageLength))

r = xStep * 0.3
c.setStrokeColor(black)
c.setFillColor(black)
c.circle(xCircleCenter, yCircleCenter, r, fill=1)

grayTransparent = Color(0, 0, 0, alpha=0.3)
greenTransparent = Color(0, 255, 0, alpha=0.15)
redTransparent = Color(255, 0, 0, alpha=0.15)

fontSize = 100


yStart = pageLength * 0.9

dictStrainPos = {}



print('cal the pos of the strains...')
lineWidth = 0.5
c.setLineWidth(lineWidth)
dictStrainRadius = {}
dRadiusFib = [1,1]

for loc in locList:
	if loc in dictLocMaxNumStrains:
		angleGap = dictLocStartAngle[loc] - locGap/2
		# rGap = dRadius * maxDistanceFromRoot * 2
		rGap = ((yearList[-1] - yearStart) * dRadius * 2)*((1.02)**(yearList[-1] - yearStart)) + dRadius
		posGapEnd = getPosXY(xCircleCenter, yCircleCenter, rGap, angleGap)
		c.setLineWidth(1)
		c.setStrokeColor(grayTransparent)	
		c.setFillColor(grayTransparent)
		c.line(xCircleCenter,yCircleCenter,posGapEnd[0], posGapEnd[1])
		
		c.setFillColor(black)
		c.setFont("Helvetica", fontSize)
		
		angleString = (dictLocStartAngle[loc]+dictLocEndAngle[loc])/2
		posString = getPosXY(xCircleCenter, yCircleCenter, rGap, angleString)
		
		c.saveState()
		c.translate(posString[0], posString[1])
		c.rotate(angleString*360/(2*math.pi))
		c.setFont("Helvetica", 200)
		c.setFillColor(black)	
		c.drawString(5*dRadius, (-3)*(dRadius), loc)	
		c.restoreState()
		
		#c.drawString(posGapEnd[0], posGapEnd[1], loc)

		angleEnd = dictLocEndAngle[loc] + locGap/2
		posGapEnd = getPosXY(xCircleCenter, yCircleCenter, rGap, angleEnd)
		c.line(xCircleCenter,yCircleCenter,posGapEnd[0], posGapEnd[1])

		c.setStrokeColor(black)
		c.setFillColor(black)
		posYearEnd = getPosXY(xCircleCenter, yCircleCenter, rGap, 0)
		c.setLineWidth(2)
		c.line(xCircleCenter,yCircleCenter,posYearEnd[0], posYearEnd[1])


		lineWidth = 0.5
		c.setLineWidth(lineWidth)




for year in range(yearList[0], yearList[-1]+1):
	#dRadiusFibNum = dRadiusFib[-1]+dRadiusFib[-2]
	radius = ((year - yearStart) * dRadius * 2)*((1.02)**(year - yearStart))
	c.setStrokeColor(black)
	c.setFillColor(black)	
	c.setLineWidth(2)
	if year % 5 == 0:
		c.line(xCircleCenter+radius, yCircleCenter, xCircleCenter+radius, yCircleCenter+10*dRadius)
		
		c.saveState()
		c.translate(xCircleCenter+radius, yCircleCenter-dRadius)
		c.rotate(270)
		c.setFont("Helvetica", 100)
		c.setFillColor(black)	
		c.drawString(0, (-1)*(dRadius), str(year))	
		c.restoreState()
	else:
		c.line(xCircleCenter+radius, yCircleCenter, xCircleCenter+radius, yCircleCenter+5*dRadius)
	#dRadiusFib.append(dRadiusFibNum)
	if year in yearList:
		for loc in locList:
			if loc in dictLocMaxNumStrains:
				# angleGap = dictLocStartAngle[loc] - locGap/2
				# rGap = dRadius * maxDistanceFromRoot * 2
				# posGapEnd = getPosXY(xCircleCenter, yCircleCenter, rGap, angleGap)
				# c.setLineWidth(lineWidth)
				# c.setStrokeColor(gray)	
				# c.setFillColor(gray)
				# c.line(xCircleCenter,yCircleCenter,posGapEnd[0], posGapEnd[1])
				
				# c.drawString(posGapEnd[0], posGapEnd[1], loc)
				if year in dictYearLocStrains and loc in dictYearLocStrains[year]:
					strains = dictYearLocStrains[year][loc]
					anglePosStart = dictLocStartAngle[loc]
					
					strainsOrdered = orderStrains(list(strains))
					
					xrange = dictLocMaxNumStrains[loc]
					numStrains = len(strainsOrdered)
					dAngleHere = (xrange/(numStrains+1))*dAngle
					
					r = min(dAngleHere*radius*0.4, dRadius*4)
					#xPosStartHere = xPosStart + ((xrange - numStrains)/2)*xStep
					anglePosStartHere = anglePosStart + dAngleHere
					angleIndex = 0
					for strainX in strainsOrdered:
						dictStrainRadius[strainX] = r
						host = dictStrainTypeInfo[strainX][2]
						anglePos = anglePosStartHere + angleIndex*dAngleHere
						(xPos, yPos) = getPosXY(xCircleCenter, yCircleCenter, radius, anglePos)
						dictStrainPos[strainX] = (xPos, yPos)
						angleIndex = angleIndex + 1
						if host in dictHostColor:
							colorHere = dictHostColor[host]
						else:
							colorHere = gray
						c.setStrokeColor(colorHere)
						c.setFillColor(colorHere)
						c.setLineWidth(lineWidth)
						c.circle(xPos, yPos, r, fill=1)

lineWidth = 0.5
c.setLineWidth(lineWidth)


radiusArrow = dRadius/2
					
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
		c.line(x1,y1,xa,ya)
		
		if locS == locT:
			c.setStrokeColor(green)
		else:
			c.setStrokeColor(red)
		c.circle(xa,ya,radiusArrow,fill=1)
		c.setFillColor(black)
		c.circle(x1,y1,radiusArrow,fill=1)

print('redraw the circles...')
lineWidth = 0.5
c.setLineWidth(lineWidth)
for year in range(yearList[0], yearList[-1]+1):
	#dRadiusFibNum = dRadiusFib[-1]+dRadiusFib[-2]
	radius = ((year - yearStart) * dRadius * 2)*((1.02)**(year - yearStart))
	#dRadiusFib.append(dRadiusFibNum)
	if year in yearList:
		for loc in locList:
			if loc in dictLocMaxNumStrains:
				# angleGap = dictLocStartAngle[loc] - locGap/2
				# rGap = dRadius * maxDistanceFromRoot * 2
				# posGapEnd = getPosXY(xCircleCenter, yCircleCenter, rGap, angleGap)
				# c.setLineWidth(lineWidth)
				# c.setStrokeColor(gray)	
				# c.setFillColor(gray)
				# c.line(xCircleCenter,yCircleCenter,posGapEnd[0], posGapEnd[1])
				
				# c.drawString(posGapEnd[0], posGapEnd[1], loc)
				if year in dictYearLocStrains and loc in dictYearLocStrains[year]:
					strains = dictYearLocStrains[year][loc]
					anglePosStart = dictLocStartAngle[loc]
					
					strainsOrdered = orderStrains(list(strains))
					
					xrange = dictLocMaxNumStrains[loc]
					numStrains = len(strainsOrdered)
					dAngleHere = (xrange/(numStrains+1))*dAngle
					
					r = min(dAngleHere*radius*0.4, dRadius*4)
					#xPosStartHere = xPosStart + ((xrange - numStrains)/2)*xStep
					anglePosStartHere = anglePosStart + dAngleHere
					angleIndex = 0
					for strainX in strainsOrdered:
						dictStrainRadius[strainX] = r
						host = dictStrainTypeInfo[strainX][2]
						anglePos = anglePosStartHere + angleIndex*dAngleHere
						(xPos, yPos) = getPosXY(xCircleCenter, yCircleCenter, radius, anglePos)
						dictStrainPos[strainX] = (xPos, yPos)
						angleIndex = angleIndex + 1
						if host in dictHostColor:
							colorHere = dictHostColor[host]
						else:
							colorHere = gray
						c.setStrokeColor(colorHere)
						c.setFillColor(colorHere)
						c.circle(xPos, yPos, r, fill=0)



# lineWidth = 0.5

# c.setLineWidth(lineWidth)


# for year in yearList:
	# yPos = yStart - (year - yearList[0])*yStep
	
	# c.setFillColor(black)
	# c.setFont("Helvetica", fontSize)
	# c.drawString(xStart + (-15)*xStep, yPos-yStep*0.3, str(year))
	# for loc in locList:
		# if loc in dictLocMaxNumStrains:
			# if year in dictYearLocStrains and loc in dictYearLocStrains[year]:
				# strains = dictYearLocStrains[year][loc]
				# xPosStart = dictLocStartX[loc]
								
				# strainsOrdered = orderStrains(list(strains))
				
				# xrange = dictLocMaxNumStrains[loc]
				# numStrains = len(strainsOrdered)
				# xStepHere = (xrange/(numStrains+1))*xStep
				# r = min(xStepHere*0.4, yStep*0.3)
				# #xPosStartHere = xPosStart + ((xrange - numStrains)/2)*xStep
				# xPosStartHere = xPosStart + xStepHere
				# xIndex = 0
				# for strainX in strainsOrdered:					
					# host = dictStrainTypeInfo[strainX][2]
					# xPos = xPosStartHere + xIndex*xStepHere
					# #dictStrainPos[strainX] = (xPos, yPos)
					# xIndex = xIndex + 1
					# if host in dictHostColor:
						# colorHere = dictHostColor[host]
					# else:
						# colorHere = gray
					# c.setStrokeColor(colorHere)
					# c.setFillColor(colorHere)
					# # if year > 2005 and year < 2014:
						# # c.circle(xPos, yPos, r, fill=1)


# yPosLocLegend = yPos - yStep
# lineWidth = 1.0
# c.setStrokeColor(black)
# c.setLineWidth(lineWidth)
# c.setFillColor(black)
# for loc in locList:
	# if loc in dictLocMaxNumStrains:
		# xPosStart = dictLocStartX[loc]
		# xPosEnd = dictLocEndX[loc]
		# c.line(xPosStart, yPosLocLegend, xPosEnd, yPosLocLegend)
		
		# c.saveState()
		# # c.translate(posXString, posYString)
		# # c.rotate(dictNodePos[clusterNode][1]*360/(2*math.pi))
		# c.setFont("Helvetica", fontSize)
		# # c.setFillColor(black)
		
		# posXString = (xPosStart+xPosEnd)/2
		# posYString = yPosLocLegend-yStep/2
		
		# c.translate(posXString, posYString)
		# c.rotate(270)
		
		
		# c.drawString(0,0, loc)
		
		# c.restoreState()
		
	
# xPosHostLegend = xPosEnd-yStep*5
# yPosHostLegend = yStart
# r = yStep*0.3
# lineWidth = 0.5
# c.setLineWidth(lineWidth)
# c.setFont("Helvetica", fontSize)
# for host in hostsDetail:
	# colorHere = dictHostColor[host]
	# c.setFillColor(colorHere)
	# c.circle(xPosHostLegend, yPosHostLegend, r, fill=1)
	
	# c.setFillColor(black)
	# c.drawString(xPosHostLegend+r*2, yPosHostLegend,host)
	
	# yPosHostLegend = yPosHostLegend - yStep


c.showPage()
c.save()

# dictLocYearReassortNum = {}
# dictLocYearLongReassortNum = {}
# dictLocYearHostDiffReassortNum = {}

# for year in dictYearLocStrains:
	# for loc in dictYearLocStrains[year]:
		# strains = dictYearLocStrains[year][loc]
		# for strain in strains:
			# if strain in dictStrainTS:
				# if loc in dictLocYearReassortNum:
					# if year in dictLocYearReassortNum[loc]:
						# dictLocYearReassortNum[loc][year] = dictLocYearReassortNum[loc][year] + 1
					# else:
						# dictLocYearReassortNum[loc][year] = 1
				# else:
					# dictLocYearReassortNum[loc] = {}
					# dictLocYearReassortNum[loc][year] = 1
				# strainT = strain
				# strainS = dictStrainTS[strainT]
				# locST = set()
				# hostST = set()
				# strainS.add(strainT)
				# ss = strainS
				# for strain in ss:
					# host = dictStrainTypeInfo[strain][2]
					# loc = dictStrainTypeInfo[strain][3]
					# hostST.add(host)
					# locST.add(loc)
				# if len(hostST) > 1:
					# if loc in dictLocYearHostDiffReassortNum:
						# if year in dictLocYearHostDiffReassortNum[loc]:
							# dictLocYearHostDiffReassortNum[loc][year] = dictLocYearHostDiffReassortNum[loc][year] + 1
						# else:
							# dictLocYearHostDiffReassortNum[loc][year] = 1
					# else:
						# dictLocYearHostDiffReassortNum[loc] = {}
						# dictLocYearHostDiffReassortNum[loc][year] = 1
					
				# if len(locST) > 1:
					# if loc in dictLocYearLongReassortNum:
						# if year in dictLocYearLongReassortNum[loc]:
							# dictLocYearLongReassortNum[loc][year] = dictLocYearLongReassortNum[loc][year] + 1
						# else:
							# dictLocYearLongReassortNum[loc][year] = 1
					# else:
						# dictLocYearLongReassortNum[loc] = {}
						# dictLocYearLongReassortNum[loc][year] = 1
					
					
# fwrite = open('locYearReassortNumPolar.txt', 'w')
# fwrite.write('locs\t')
# for year in yearList:
	# fwrite.write(str(year)+'all\t'+str(year)+'longDist\t'+str(year)+'diffHost\t')
# fwrite.write('\n')
# for loc in locList:
	# if loc in dictLocMaxNumStrains:
		# fwrite.write(loc+'\t')
		# for year in yearList:
			# if loc in dictLocYearReassortNum:
				# if year in dictLocYearReassortNum[loc]:
					# fwrite.write(str(dictLocYearReassortNum[loc][year])+'\t')
				# else:
					# fwrite.write('0\t')
			# else:
				# #print(loc, 'not in dictLocYearReassortNum')
				# fwrite.write('0\t')
			# if loc in dictLocYearLongReassortNum:
				# if year in dictLocYearLongReassortNum[loc]:
					# fwrite.write(str(dictLocYearLongReassortNum[loc][year])+'\t')
				# else:
					# fwrite.write('0\t')
			# else:
				# fwrite.write('0\t')
			
			# if loc in dictLocYearHostDiffReassortNum:
				# if year in dictLocYearHostDiffReassortNum[loc]:
					# fwrite.write(str(dictLocYearHostDiffReassortNum[loc][year])+'\t')
				# else:
					# fwrite.write('0\t')
			# else:
				# fwrite.write('0\t')
		# fwrite.write('\n')
				



# fwrite.close()	
		
		






#接下来需要将关系进行筛选，因为有些重配应该是不太合理的，我觉得。。。




