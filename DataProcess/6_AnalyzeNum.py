import os
import sys
import datetime
import time


aaLigal = 'ACDEFGHIKLMNPQRSTVWY'
ntLigal = 'ATCG'


dictHost = dict()
dictCountry = dict()
dictDateStrain = dict()
dictYear = dict()
dictSubtype = dict()
dictStrainName = dict()
dictGenomeInfo = dict()

dictGenomeHost = dict()
dictGenomeCountry = dict()
dictGenomeDate = dict()
dictGenomeSubtype = dict()
dictGenomeStrainName = dict()

fread = open('genomeProteinCDScomplete.dat', 'r')
for rline in fread:
	if rline[0] == '>':
		genomeID = rline[1:rline.find('|')]
		genomeInfo = rline[rline.find('|')+1:].rstrip()
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
fread.close()

fwrite = open('country.txt', 'w')
for country in sorted(dictCountry.keys()):
	print(country+'\t'+str(dictCountry[country]))
	fwrite.write(country+'\t'+str(dictCountry[country]))
	fwrite.write('\n')
fwrite.close()


fwrite = open('subtype.txt', 'w')
for subtype in sorted(dictSubtype.keys()):
	print(subtype+'\t'+str(dictSubtype[subtype]))
	fwrite.write(subtype+'\t'+str(dictSubtype[subtype]))
	fwrite.write('\n')
fwrite.close()
	

fwrite = open('host.txt', 'w')
for host in dictHost:
	print(host+'\t'+str(dictHost[host]))
	fwrite.write(host+'\t'+str(dictHost[host]))
	fwrite.write('\n')
fwrite.close()

fwrite = open('year.txt', 'w')
for year in sorted(dictYear.keys()):
	print(year+'\t'+str(dictYear[year]))
	fwrite.write(year+'\t'+str(dictYear[year]))
	fwrite.write('\n')
fwrite.close()

