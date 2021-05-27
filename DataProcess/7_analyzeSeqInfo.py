import os
import sys
import datetime
import time

aaLigal = 'ACDEFGHIKLMNPQRSTVWY'
ntLigal = 'ATCG'

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


dictAvianClass = dict()
fread = open('avianClassification.txt', 'r')
for rline in fread:
	if rline[0] != '>':
		lrline = rline.strip().split('\t')
		bird = lrline[0].lower()
		classification = lrline[1]
		if bird in dictAvianClass:
			print('repeated bird:', bird)
			sys.exit()
		dictAvianClass[bird] = classification
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

dictAvianCount = dict()

lSegs = ['PB2','PB1','PA','HA','NP','NA','MP','NS']

fwrite = open('genomeAvianHostsClassification.txt', 'w')

fwrite.write('>genomeID	strainName	bird	classification\n')

dictBirdCount = dict()

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
		
		if host.lower() == 'avian':
			# print(strainName)
			bird = strainName.split('/')[1].lower()
			if bird in dictBirdCount:
				dictBirdCount[bird] = dictBirdCount[bird] + 1
			else:
				dictBirdCount[bird] = 1
			
			if bird not in dictAvianClass:
				print(bird, 'not in dictAvianClass!')
				sys.exit()
			fwrite.write(genomeID+'\t'+strainName+'\t'+bird + '\t' + dictAvianClass[bird] +'\n')
			if bird in dictAvianClass:
				classification = dictAvianClass[bird]
				if classification in dictAvianCount:
					dictAvianCount[classification] = dictAvianCount[classification] + 1
				else:
					dictAvianCount[classification] = 1
			else:
				print(bird, 'not in dictAvianClass!')
				sys.exit()
		
fread.close()
fwrite.close()


print('dictAvianCount:')
print(dictAvianCount)

fwrite = open('birdsHere.txt', 'w')
for bird in dictBirdCount:
	fwrite.write(bird+'\t'+str(dictBirdCount[bird])+'\n')

fwrite.close()




