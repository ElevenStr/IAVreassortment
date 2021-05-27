from __future__ import division
import os
import sys
import datetime
import time
import numpy as np
from random import choice

aaLigal = 'ACDEFGHIKLMNPQRSTVWY'
ntLigal = 'ATCG'

similarityGate = 0.99


dictEDNAFULL = dict()
fread = open('EDNAFULL_4table.txt', 'r')
allLines = fread.readlines()
fread.close()
nts = allLines[0].strip().split(' ')
for i in range(len(nts)):
	curLine = allLines[i+1].strip().split(' ')
	curNT = nts[i]
	for j in range(len(nts)):
		nowNT = nts[j]
		dictEDNAFULL[(curNT, nowNT)] = int(curLine[j+1])



dictCountryContinent = dict()
dictCountryContinentDetail = dict()
dictContinentDetailContinent = dict()
fread = open('countryContinetContinentDetail.txt', 'r')
for rline in fread:
	if rline[0] != '>':
		lrline = rline.strip().split('\t')
		country = lrline[0]
		continent = lrline[1]
		continentDetail = lrline[2]
		
		if country in dictCountryContinent:
			print('repeated country!!!', country)
		else:
			dictCountryContinent[country] = continent
			dictCountryContinentDetail[country] = continentDetail
		
		if continentDetail in dictContinentDetailContinent:
			continentIN = dictContinentDetailContinent[continentDetail]
			if continentIN != continent:
				print('not matched continent!!!!')
				sys.exit()
		dictContinentDetailContinent[continentDetail] = continent

fread.close()


fread = open('genomeAvianHostsClassification.txt', 'r')
dictGenomeAvianClassification = dict()

for rline in fread:
	if rline[0] != '>':
		lrline = rline.strip().split('\t')
		genomeID = lrline[0]
		strainName = lrline[1]
		bird = lrline[2]
		avianClassification = lrline[3]
		dictGenomeAvianClassification[genomeID] = avianClassification
fread.close()

dictHost = dict()
dictCountry = dict()
dictDateStrain = dict()
dictYear = dict()
dictSubtype = dict()
dictStrainName = dict()
dictGenomeInfo = dict()

dictContinentDetail = dict()
dictContinentDetailYear = dict()
dictContinentDetailHost = dict()


dictGenomeHost = dict()
dictGenomeCountry = dict()
dictGenomeDate = dict()
dictGenomeSubtype = dict()
dictGenomeStrainName = dict()

lSegs = ['PB2','PB1','PA','HA','NP','NA','MP','NS']

fread = open('genomeProteinCDScomplete.dat', 'r')
for rline in fread:
	if rline[0] == '>':
		genomeID = rline[1:rline.find('|')]
		genomeInfo = rline[rline.find('|')+1:].strip()
		dictGenomeInfo[genomeID] = genomeInfo
		genomeInfo = genomeInfo.split('\t')
		host = genomeInfo[0]
		if host == '':
			print(genomeID, 'host empty!')
		subtype = genomeInfo[1]
		if subtype == '':
			print(genomeID, 'subtype empty!')
		country = genomeInfo[2]
		if country == '':
			print(genomeID, 'country empty!')
		dateStrain = genomeInfo[3]
		if dateStrain == '':
			print(genomeID, 'dateStrain empty!')
		year = dateStrain[:4]
		if year in dictYear:
			dictYear[year] = dictYear[year] + 1
		else:
			dictYear[year] = 1
		
		dictGenomeHost[genomeID] = host
		dictGenomeCountry[genomeID] = country
		dictGenomeDate[genomeID] = dateStrain
		dictGenomeSubtype[genomeID] = subtype
		
		if dateStrain in dictDateStrain:
			dictDateStrain[dateStrain] = dictDateStrain[dateStrain] + 1
		else:
			dictDateStrain[dateStrain] = 1
		
		strainName = genomeInfo[4]
		dictGenomeStrainName[genomeID] = strainName
		dictStrainName[strainName] = genomeID
				
		if host in dictHost:
			dictHost[host] = dictHost[host] + 1
		else:
			dictHost[host] = 1
		if subtype in dictSubtype:
			dictSubtype[subtype] = dictSubtype[subtype] + 1
		else:
			dictSubtype[subtype] = 1
		if country in dictCountry:
			dictCountry[country] = dictCountry[country] + 1
		else:
			dictCountry[country] = 1
		
		if country not in dictCountryContinentDetail:
			print(country, 'not found in dictCountryContinentDetail!')
		continentDetail = dictCountryContinentDetail[country]
		if continentDetail in dictContinentDetail:
			dictContinentDetail[continentDetail] = dictContinentDetail[continentDetail] + 1
		else:
			dictContinentDetail[continentDetail] = 1
		
		year = int(year)
		if continentDetail in dictContinentDetailYear:			
			if year in dictContinentDetailYear[continentDetail]:
				dictContinentDetailYear[continentDetail][year].add(genomeID)
			else:
				s = set()
				s.add(genomeID)
				dictContinentDetailYear[continentDetail][year] = s
		else:
			s = set()
			s.add(genomeID)
			d = dict()
			d[year] = s
			dictContinentDetailYear[continentDetail] = d
fread.close()

#fwrite = open('continentDetailViruses.txt', 'w')

for continentDetail in dictContinentDetailYear:
	#continentDetailSum = sum(dictContinentDetailYear[continentDetail].values())
	continentDetailSum = 0
	for year in dictContinentDetailYear[continentDetail]:
		continentDetailSum = continentDetailSum + len(dictContinentDetailYear[continentDetail][year])
	print(continentDetail, continentDetailSum, dictContinentDetailContinent[continentDetail])
#	fwrite.write(continentDetail+'\t'+str(continentDetailSum)+'\t'+dictContinentDetailContinent[continentDetail]+'\n')
	
#fwrite.close()

dictGenomeHostClass = dict()

#dictOthers = dict()


for genomeID in dictGenomeHost:
	host = dictGenomeHost[genomeID]
	if host.lower() == 'avian':
		hostClass = 'avian_'+dictGenomeAvianClassification[genomeID]
		dictGenomeHostClass[genomeID] = hostClass.lower()
	elif host.lower() == 'human' or host.lower() == 'swine':
		dictGenomeHostClass[genomeID] = host.lower()
	else:
		# if host in dictOthers:
			# dictOthers[host] = dictOthers[host] + 1
		# else:
			# dictOthers[host] = 1
		dictGenomeHostClass[genomeID] = host.lower()	



dictHostYearLoationGenome = dict()		

for genomeID in dictGenomeHostClass:
	hostClass = dictGenomeHostClass[genomeID]
	#location = dictCountryContinentDetail[dictGenomeCountry[genomeID]]
	location = dictGenomeCountry[genomeID]
	year = int(dictGenomeDate[genomeID][:4])
	subtype = dictGenomeSubtype[genomeID]
	
	if year > 0 and location.lower() != 'unknown' and subtype.lower() != 'mixed':
		if (hostClass, location, year, subtype) in dictHostYearLoationGenome:
			dictHostYearLoationGenome[(hostClass, location, year, subtype)].add(genomeID)
		else:
			s = set()
			s.add(genomeID)
			dictHostYearLoationGenome[(hostClass, location, year, subtype)] = s


print('**************')
print(len(dictHostYearLoationGenome), 'clusters in total.')
	
lSegs8 = ['PB2','PB1','PA','HA','NP','NA','MP','NS']
lSegsInner = ['PB2','PB1','PA','NP','MP','NS']

dictGenomeSeq8 = dict()
dictGenomeSegSeq = dict()

for i in range(17):
	htype = 'H'+str(i+1)
	print(htype)
	fread = open(os.path.join('HANA', htype+'.fasta.codon.fas'), 'r')
	allLines = fread.readlines()
	fread.close()
	for j in range(len(allLines)):
		rline = allLines[j]
		if rline[0] == '>':
			genomeID = rline.split('|')[0][1:]
			seq = allLines[j+1].strip()
			if genomeID in dictGenomeSeq8:
				dictGenomeSeq8[genomeID] = dictGenomeSeq8[genomeID] + seq
				print('double genomeID!!!')
				sys.exit()
			else:
				dictGenomeSeq8[genomeID] = seq
				dictGenomeSegSeq[genomeID] = {}
				dictGenomeSegSeq[genomeID]['HA'] = seq
				
				

for i in range(10):
	ntype = 'N'+str(i+1)
	print(ntype)
	fread = open(os.path.join('HANA', ntype+'.fasta.codon.fas'), 'r')
	allLines = fread.readlines()
	fread.close()
	for j in range(len(allLines)):
		rline = allLines[j]
		if rline[0] == '>':
			genomeID = rline.split('|')[0][1:]
			seq = allLines[j+1].strip()
			if genomeID in dictGenomeSeq8:
				dictGenomeSeq8[genomeID] = dictGenomeSeq8[genomeID] + seq	
				dictGenomeSegSeq[genomeID]['NA'] = seq
			else:
				print(genomeID, 'not complete...')
				sys.exit()
				dictGenomeSeq8[genomeID] = seq



for seg in lSegsInner:
	print(seg)
	fread = open(seg+'.fasta.codon.fas', 'r')
	allLines = fread.readlines()
	fread.close()
	for i in range(len(allLines)):
		rline = allLines[i]
		if rline[0] == '>':
			genomeID = rline.split('|')[0][1:]
			seq = allLines[i+1].strip()
			if genomeID in dictGenomeSeq8:
				dictGenomeSeq8[genomeID] = dictGenomeSeq8[genomeID] + seq
				dictGenomeSegSeq[genomeID][seg] = seq
			else:
				print(genomeID, 'not complete...')
				sys.exit()
				dictGenomeSeq8[genomeID] = seq

print('finished reading genomes...')


def computeDNASimilarity(seqi, seqj):
	if len(seqi) != len(seqj):
		sys.exit('Not equal seq lengths!')
	else:
		lenAll = 0
		lenIdentity = 0
		lenSimilarity = 0
		for k in range(len(seqi)):
			nti = seqi[k]
			ntj = seqj[k]
			if nti != '-' and ntj != '-':
				lenAll = lenAll + 1
				scoreDNAFull = dictEDNAFULL[(nti, ntj)]
				if scoreDNAFull > 0:
					lenSimilarity = lenSimilarity + 1
				if nti == ntj:
					lenIdentity = lenIdentity + 1
			elif nti == '-' and ntj == '-':
				lenAll = lenAll
			else:
				lenAll = lenAll + 1
		#print(lenIdentity, lenSimilarity, lenAll)
		return float(lenSimilarity/lenAll)
				
def getFilteredGenomeIDs(dictGenomeSeq8, genomeIDs):
	s = set()
	lGenomeIDs = list(genomeIDs)
	s.add(lGenomeIDs[0])
	lGenomeIDs.pop(0)
	while True:
		if len(lGenomeIDs) < 1:
			break
		s2remove = set()
		for genomeIDinL in lGenomeIDs:
			for genomeIDinS in s:
				# print(genomeIDinL,genomeIDinS)
				similarity = computeDNASimilarity(dictGenomeSeq8[genomeIDinL], dictGenomeSeq8[genomeIDinS])
				if similarity > similarityGate:
					s2remove.add(genomeIDinL)
		for g2remove in s2remove:
			lGenomeIDs.remove(g2remove)
		if len(lGenomeIDs) < 1:
			break
		s.add(lGenomeIDs[0])
	return s		
		
		
#fwrite = open('genomeSeq.fas', 'w')
dictHostYearLoationGenomeFiltered = dict()			
for (hostClass, location, year, subtype) in dictHostYearLoationGenome:
	# genomeIDs = dictHostYearLoationGenome[(hostClass, location, year, subtype)]
	# for genomeID in genomeIDs:
		# genomeIDName = ('>'+genomeID+'|'+str(year)+'|'+location+'|'+subtype+'|'+hostClass).replace(' ', '~')
		# fwrite.write(genomeIDName+'\n')
		# fwrite.write(dictGenomeSeq8[genomeID]+'\n')
	genomeIDs = dictHostYearLoationGenome[(hostClass, location, year, subtype)]
	if len(dictHostYearLoationGenome[(hostClass, location, year, subtype)]) > 1:
		print(hostClass, location, year, subtype)
		print(len(genomeIDs), 'genomes in total!')
		# genomeIDSelected = choice(list(genomeIDs))
		# genomeIDs.remove(genomeIDSelected)
		# minSimilarity = 1
		# for genomeID in genomeIDs:
			# similarity = computeDNASimilarity(dictGenomeSeq8[genomeIDSelected], dictGenomeSeq8[genomeID])
			# #print(similarity)
			# if similarity < minSimilarity:
				# minSimilarity = similarity
		# if minSimilarity < 0.95:
			# print(minSimilarity)
		#print(minSimilarity)
		filteredGenomeIDs = getFilteredGenomeIDs(dictGenomeSeq8, genomeIDs)
		dictHostYearLoationGenomeFiltered[(hostClass, location, year, subtype)] = filteredGenomeIDs
		print(len(genomeIDs), len(filteredGenomeIDs), 'changed...')
	else:
		dictHostYearLoationGenomeFiltered[(hostClass, location, year, subtype)] = genomeIDs

# fwrite.close()

fwrite = open('dictHostYearCountryGenomeFiltered'+str(similarityGate)+'.txt', 'w')
for (hostClass, location, year, subtype) in dictHostYearLoationGenome:
	print((hostClass, location, year, subtype), len(dictHostYearLoationGenome[(hostClass, location, year, subtype)]))
	fwrite.write(hostClass+'\t'+location+'\t'+str(year)+'\t'+subtype+'\t'+str(len(dictHostYearLoationGenome[(hostClass, location, year, subtype)]))+'\t'+str(len(dictHostYearLoationGenomeFiltered[(hostClass, location, year, subtype)]))+'\t')
	for genomeID in dictHostYearLoationGenomeFiltered[(hostClass, location, year, subtype)]:
		fwrite.write(genomeID+' ')
	fwrite.write('\t')
	for genomeID in dictHostYearLoationGenome[(hostClass, location, year, subtype)]:
		fwrite.write(genomeID+' ')
	fwrite.write('\n')
fwrite.close()










# import matplotlib.pyplot as plt


# x1 = np.linspace(0.0, 5.0)
# x2 = np.linspace(0.0, 2.0)

# y1 = np.cos(2 * np.pi * x1) * np.exp(-x1)
# y2 = np.cos(2 * np.pi * x2)

# lYearSegs = []

# for i in range(100):
	# if i < 1:
		# segStart = -1
		# segEnd = 1975
	# else:
		# segStart = (i-1)*1 + 1975
		# segEnd = i*1 + 1975
	# if segEnd < 2016:
		# lYearSegs.append((segStart, segEnd))

# lenContinentDetail = len(dictContinentDetailYear)

# lContinentDetail = ["East Africa","North Africa","South Africa","West Africa","Central Asia","East Asia","South Asia","Southeast Asia","West Asia","Central Europe","Eastern Europe","Northern Europe","Southern Europe","Western Europe","Middle America","North America","The Caribbean","Oceania","Eastern South America","Midwest South American","Northern South America","Southern South America","Unknown"]

# rank = 0

# font = {'family' : 'serif',
        # 'color'  : 'black',
        # 'weight' : 'normal',
        # 'size'   : 12,
        # }

# listColors = ['darkcyan', 'darkred', 'gray', 'purple', 'red', 'yellowgreen', 'blue', 'yellow']

# colorList = ['blue', 'yellow','red','maroon','purple','pink','gray','black']

# index = 0
# dictContinentColor = dict()
# for continent in set(dictContinentDetailContinent.values()):
	# dictContinentColor[continent] = listColors[index]
	# index = index + 1



# def getYearRangeIndex(lYearSegs, year):
	# for i in range(len(lYearSegs)):
		# segStart = lYearSegs[i][0]
		# segEnd = lYearSegs[i][1]
		# if year > segStart and year < segEnd + 1:
			# return i
	# return 'unknown'

# lYearSegsShow = []

# lHostClass = ['human','swine','avian_domestic birds','avian_waterfowl','avian_shorebirds','avian_land birds','avian_unknown','others']

# for seg in lYearSegs:
	# lYearSegsShow.append(seg[1])

	
"""test!"""







