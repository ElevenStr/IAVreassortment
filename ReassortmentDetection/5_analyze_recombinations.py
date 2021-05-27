from __future__ import division
import os
import sys
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch
from reportlab.lib.colors import green, blue, yellow, red, maroon, purple, pink, gray, black, darkslategray


colorList = [green,blue,red]

dictHostColor = dict()
hosts = ['avian','swine','human']
i = 0
for host in hosts:
	dictHostColor[host] = colorList[i]
	i = i + 1


segs = ['PB2','PB1','PA','HA','NP','NA','MP','NS']



locList = ["East Africa","North Africa","South Africa","West Africa","Middle Africa","Central Asia","East Asia","South Asia","Southeast Asia","West Asia","Central Europe","Eastern Europe","Northern Europe","Southern Europe","Western Europe","Middle America","North America","The Caribbean","Oceania","Eastern South America","Midwest South American","Northern South America","Southern South America"]
colorJetList = ['#00008F','#00009F','#0000AF','#0000BF','#0000CF','#0000DF','#0000EF','#0000FF','#000FFF','#001FFF','#002FFF','#003FFF','#004FFF','#005FFF','#006FFF','#007FFF','#008FFF','#009FFF','#00AFFF','#00BFFF','#00CFFF','#00DFFF','#00EFFF','#00FFFF','#0FFFEF','#1FFFDF','#2FFFCF','#3FFFBF','#4FFFAF','#5FFF9F','#6FFF8F','#7FFF7F','#8FFF6F','#9FFF5F','#AFFF4F','#BFFF3F','#CFFF2F','#DFFF1F','#EFFF0F','#FFFF00','#FFEF00','#FFDF00','#FFCF00','#FFBF00','#FFAF00','#FF9F00','#FF8F00','#FF7F00','#FF6F00','#FF5F00','#FF4F00','#FF3F00','#FF2F00','#FF1F00','#FF0F00','#FF0000','#EF0000','#DF0000','#CF0000','#BF0000','#AF0000','#9F0000','#8F0000','#7F0000']

dictLocColor = {}

for i in range(len(locList)):
	loc = locList[i]
	colorIndex = int(64/len(locList)*i)%64
	# print(colorIndex)
	colorHere = colorJetList[colorIndex]
	dictLocColor[loc] = colorHere

print(os.getcwd())

# fread = open('clusterNo.txt', 'r')
# lineCount = 0

# dictSegClusterSize = dict()
# dictClusterSize = dict()

# for rline in fread:
# 	lineCount += 1
# 	if lineCount > 1:
# 		lrline = rline.strip().split('\t')
# 		clusterName = lrline[0]
# 		clusterSize = int(lrline[1])
# 		seg = clusterName[:clusterName.find('_')]
# 		if clusterSize > 9:
# 			dictClusterSize[clusterName] = clusterSize
# 			if seg in dictSegClusterSize:
# 				if clusterName in dictSegClusterSize[seg]:
# 					print('duplicated clusterNames', clusterName)
# 					sys.exit()
# 				else:
# 					dictSegClusterSize[seg][clusterName] = clusterSize
# 			else:
# 				d = dict()
# 				d[clusterName] = clusterSize
# 				dictSegClusterSize[seg] = d
		
# fread.close()

dictClusterInfo = dict()
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

#test combination
#isolateID = 'EPI_ISL_965'  # H5N1 1997
#isolateID = 'EPI_ISL_12597'  #H9N2 
#isolateID = 'EPI_ISL_3547'  #H7N7
#isolateID = 'EPI_ISL_138737' #H7N9 2013

#isolateID = 'EPI_ISL_12401'  # H5N1 2003
#isolateID = 'EPI_ISL_20818'  # H3N2 1968
#isolateID = 'EPI_ISL_170560'  # H2N2 1957
#isolateID = 'EPI_ISL_123154'  # H1N1 1999
#isolateID = 'EPI_ISL_7047'  # H1N1 2000 US
#isolateID = 'EPI_ISL_5873'  # H7N3
#isolateID = 'EPI_ISL_32639'  # H1N1 mexico
# isolateID = '1068849'  # H1N1 Taiwan
# isolateID = '1067554'  # H5N6 Viet Nam
# isolateID = '1024970'  # H7N9 Shanghai
# isolateID = '1066732'  # H1N1 Anhui
isolateID = '31087'  # H7N9 WuXi






print(dictIsolateInfo[isolateID])

(pageWidth, pageLength) = (2100, 2970)



lrline = dictIsolateInfo[isolateID].split('\t')
isolateName = lrline[1]
host_Classification = lrline[2]
loc = lrline[3]
subtype = lrline[4]
yearThis = int(lrline[5])

#isolateNameOut = isolateName.replace('/', '_')
c = canvas.Canvas(isolateID+'_'+subtype+'_'+loc+'_'+host_Classification+'_'+str(yearThis)+'_recombination'+'.pdf', pagesize=(pageWidth, pageLength))


dictSegCluster = dictIsolateSegCluster[isolateID]

setYears = set()

for seg in segs:
	segCluster = dictSegCluster[seg]
	#print(seg+':'+segCluster)

	years = list(dictClusterInfo[segCluster].keys())
	#years.sort()

	for year in years:
		#if year < yearThis + 1 :
		setYears.add(year)
years = list(setYears)
years.sort()

widthYear = pageWidth / (2 * len(years))
widthGap = widthYear / 5

lengthSeg = pageLength / (2 * len(segs))
lengthGap = lengthSeg / 5

startYearX = pageWidth/20
startYearY = pageLength*0.95

c.drawString(pageWidth/3, pageLength*0.98, isolateName)


yearX = startYearX

for year in years:
	c.setLineWidth(2)
	c.setFillColor(black)
	
	if year == yearThis:
		c.setFillColor(red)
		c.setStrokeColor(red)
	c.line(yearX ,startYearY, yearX + widthYear, startYearY)
	c.drawString(yearX + widthYear/4, startYearY+5,str(year))
	c.setFillColor(black)
	c.setStrokeColor(black)
	yearX = yearX + widthYear + widthGap


startSegX = yearX
startSegY = startYearY - lengthGap


segY = startSegY

for seg in segs:
	segCluster = dictSegCluster[seg]

	c.setLineWidth(2)
	c.setStrokeColor(black)
	c.line(startSegX, segY, startSegX, segY - lengthSeg)
	c.drawString(startSegX + 5, segY - lengthSeg/2, seg+':'+segCluster)
	
	yearX = startYearX
	for year in years:
		if year in dictClusterInfo[segCluster]:
			setTuple = set()
			for (host, loc, subtype) in dictClusterInfo[segCluster][year]:
				if 'avian' in host.lower():
					setTuple.add(('avian', loc, subtype))
				else:
					setTuple.add((host, loc, subtype))
					
					
			ytupleStep = lengthSeg / 5

			if len(setTuple) > 4:
				ytupleStep = lengthSeg / (len(setTuple)+1)


			r = ytupleStep / 8
			xHost = widthYear / 10
			xLoc = xHost + r*2
			xSubtype = xLoc + r*3
			
			ytuple = 0

			for (host, loc, subtype) in setTuple:
				if host in dictHostColor:
					colorHost = dictHostColor[host]
				else:
					colorHost = gray
				if loc in dictLocColor:
					colorLoc = dictLocColor[loc]
				else:
					colorLoc = gray
				
				c.setStrokeColor(colorHost)
				c.setFillColor(colorHost)
				c.circle(yearX + xHost, segY - ytuple - r, r, fill=1)
				c.setFillColor(black)
				c.setFillColor(black)


				c.setStrokeColor(colorLoc)
				c.setFillColor(colorLoc)
				c.rect(yearX + xLoc, segY - ytuple - 2*r, 2*r, 2*r, fill=1)
				c.setFillColor(black)
				c.setFillColor(black)
			

				#c.drawString(yearX + xHost, segY - ytuple, host)
				#c.drawString(yearX + xLoc, segY - ytuple, loc)
				c.drawString(yearX + xSubtype, segY - ytuple - 2*r, subtype)
				# print(segY - ytuple)
				ytuple = ytuple + ytupleStep


		yearX = yearX + widthYear + widthGap
	segY = segY - lengthSeg - lengthGap
			

#plot host Color legend
hostCount = 0
r = lengthSeg / 20
for host in dictHostColor:
	hostCount = hostCount + 1
	colorHere = dictHostColor[host]
	c.setStrokeColor(colorHere)
	c.setFillColor(colorHere)
	posHost = [pageWidth/1.2,pageLength/1.2 - hostCount*( r * 4)]
	#c.rect(posHost[0], posHost[1], r*3, r*3, fill=1)
	c.circle(posHost[0], posHost[1], r, fill=1)
	c.setFillColor(black)
	c.drawString(posHost[0]+3*r, posHost[1] - r,host)



locCount = 0
r = lengthSeg / 20
for loc in locList:
	locCount = locCount + 1
	colorHere = dictLocColor[loc]
	c.setStrokeColor(colorHere)
	c.setFillColor(colorHere)
	posloc = [pageWidth/1.1,pageLength/1.2 - locCount*( r * 4)]
	#c.rect(posloc[0], posloc[1], r*3, r*3, fill=1)
	c.rect(posloc[0], posloc[1], 2*r, 2*r, fill=1)
	c.setFillColor(black)
	c.drawString(posloc[0]+3*r, posloc[1],loc)


# for seg in segs:
# 	segCluster = dictSegCluster[seg]
# 	print(seg+':'+segCluster)
	
# 	for year in years:
# 		if year in dictClusterInfo[segCluster]:
			
# 		if year < yearThis + 1:
# 			print(year,dictClusterInfo[segCluster][year])
			

# 	#print(dictClusterInfo[segCluster])

c.showPage()
c.save()























# fread = open('isolatesCombinations.txt', 'r')
# fwrite = open('isolatesCombinationsFiltered10.txt', 'w')
# lineCount = 0
# for rline in fread:
# 	lineCount += 1
# 	if lineCount > 1:
# 		lrline = rline.strip().split('\t')
# 		for j in range(len(lrline)):
# 			"EPI_ISL_112670	Avian_waterfowl	East Asia	H6N6	2007"
# 			host_Classification = lrline[1]
# 			loc = lrline[2]
# 			subtype = lrline[3]
# 			year = int(lrline[4])
# 			if j > 5:
# 				clusterName = lrline[j]
# 				seg = clusterName[:clusterName.find('_')]
# 				if seg in dictSegClusterYearLocHost:
# 					if clusterName in dictSegClusterYearLocHost[seg]:
# 						if year in dictSegClusterYearLocHost[seg][clusterName]:
# 							if (loc, host_Classification) in dictSegClusterYearLocHost[seg][clusterName][year]:
# 								dictSegClusterYearLocHost[seg][clusterName][year][(loc, host_Classification)] = dictSegClusterYearLocHost[seg][clusterName][year][(loc, host_Classification)] + 1
# 							else:
# 								dictSegClusterYearLocHost[seg][clusterName][year][(loc, host_Classification)] = 1
# 						else:
							


# 				if clusterName in dictClusterSize:
# 					fwrite.write(clusterName)
# 				else:
# 					fwrite.write('-')
# 			else:
# 				fwrite.write(lrline[j])
# 			fwrite.write('\t')
# 		fwrite.write('\n')
# 	else:
# 		fwrite.write(rline)


# fwrite.close()
# fread.close()






# PB2 157
# PB1 153
# PA 157
# HA 164
# NP 165
# NA 152
# MP 147
# NS 169


# fread = open('isolateHostClassLocationRankedByTime.txt', 'r')
# lineNo = 0
# dictIsolateHost = dict()
# dictIsolateLoc = dict()
# dictIsolateSubtype = dict()
# dictIsolateYear = dict()
# dictIsolateInfo = dict()
# dictIsolateTime = dict()

# for rline in fread:
	# lineNo += 1
	# #isolateID	host_Classification	Continet_detail	subtype	year
	# if lineNo > 1:
		# lrline = rline.strip().split('\t')
		# isolateID = lrline[0]
		# dictIsolateInfo[isolateID] = rline.strip()
		# dictIsolateHost[isolateID] = lrline[1]
		# dictIsolateLoc[isolateID] = lrline[2]
		# dictIsolateSubtype[isolateID] = lrline[3]
		# dictIsolateYear[isolateID] = int(lrline[4])
		# dictIsolateTime[isolateID] = lrline[5]
# fread.close()


# dictSegIsolateCluster = dict()
# dictSegClusterNo = dict()
# dictSegClusterPeriod = dict()

	
# fwriteClusterNo = open('clusterNo.txt', 'w')
# fwriteClusterNo.write('segCluster\tclusterSize\tyearStart,yearEnd\tancestor\n')
# dictCluterSize = dict()
# dictClusterPeriod = dict()

# for seg in segs:
	# print(seg)
	
	# dictClusterIsolates = dict()
	
	
	# fread = open(seg+'_all_inClusters.txt', 'r')
	
	# dictIsolateCluster = dict()
	# dictClusterNo = dict()
	# clusterNo = 0
	
	
	# for rline in fread:
		# if rline.startswith('>>>'):
			# clusterAncestor = rline.strip()
			# clusterAncestorShort = clusterAncestor[:clusterAncestor.find('\t')]
			# if clusterAncestorShort in dictClusterNo:
				# print('repeated...', clusterAncestorShort)
			# else:
				# clusterNo += 1
				# dictClusterNo[clusterAncestorShort] = clusterNo
			# s = set()
			# dictClusterIsolates[clusterAncestorShort] = s
		# else:
			# if '|' not in rline:
				# isolateIDs = rline.strip().split('\t')
				# for isolateID in isolateIDs:
					# dictClusterIsolates[clusterAncestorShort].add(isolateID)
					# dictIsolateCluster[isolateID] = clusterAncestorShort
		
	# print(seg, len(dictClusterIsolates), 'clusters.')
	
	
	# fwrite = open(seg+'_clusters_isolates.txt', 'w')
	
	# for clusterAncestorShort in dictClusterIsolates:
		# #clusterAncestorShort = clusterAncestor[:clusterAncestor.find('\t')]
		# fwrite.write(clusterAncestorShort+'\t'+str(len(dictClusterIsolates[clusterAncestorShort]))+'\t')
		# for isolateID in dictClusterIsolates[clusterAncestorShort]:
			# fwrite.write(isolateID+' ')
		# fwrite.write('\n')
	# fwrite.close()
	
	# dictSegIsolateCluster[seg] = dictIsolateCluster
	# dictSegClusterNo[seg] = dictClusterNo
	
	# fread.close()
	
	# for clusterAncestorShort in dictClusterIsolates:
		# yearStart = 9999
		# yearEnd = 0
		# for isolateID in dictClusterIsolates[clusterAncestorShort]:
			# year = dictIsolateYear[isolateID]
			# if year < yearStart:
				# yearStart = year
			# if year > yearEnd:
				# yearEnd = year
		# clusterName = seg+'_'+str(dictClusterNo[clusterAncestorShort])
		# dictClusterPeriod[clusterName] = (yearStart, yearEnd)
		
	
	

	
	# for clusterAncestorShort in dictClusterNo:
		# clusterName = seg+'_'+str(dictClusterNo[clusterAncestorShort])
		# clusterPeriod = dictClusterPeriod[clusterName]
		# clusterPeriod2Write = str(clusterPeriod[0])+'\t'+str(clusterPeriod[1])
		# fwriteClusterNo.write(clusterName+'\t'+str(len(dictClusterIsolates[clusterAncestorShort]))+'\t'+clusterPeriod2Write+'\t'+clusterAncestorShort+'\n')
		# dictCluterSize[clusterName] = len(dictClusterIsolates[clusterAncestorShort])
# fwriteClusterNo.close()





# fwrite = open('isolatesCombinations.txt', 'w')
# fread = open('isolateHostClassLocationRankedByTime.txt', 'r')
# lineNo = 0
# dictIsolateCombination = dict()

# for rline in fread:
	# lineNo += 1
	# #isolateID	host_Classification	Continet_detail	subtype	year
	# if lineNo > 1:
		# lrline = rline.strip().split('\t')
		# isolateID = lrline[0]
		# fwrite.write(rline.strip()+'\t')
		# combination = []
		# for seg in segs:
			# clusterAncestorShort = dictSegIsolateCluster[seg][isolateID]
			# clusterNo = dictSegClusterNo[seg][clusterAncestorShort]
			# clusterName = seg+'_'+str(clusterNo)
			# fwrite.write(clusterName+'\t')
			# combination.append(clusterName)
		# fwrite.write('\n')
		# dictIsolateCombination[isolateID] = tuple(combination)
			
	# else:
		# fwrite.write(rline.strip()+'\t'+'PB2\tPB1\tPA\tHA\tNP\tNA\tMP\tNS\n')
		
		
# fread.close()

# fwrite.close()

# s = set(dictIsolateCombination.values())
# print(len(s), 'combinations for ', len(dictIsolateCombination), 'strians...' )


	
	
	
