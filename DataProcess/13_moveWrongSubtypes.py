import os

fread = open('HANAcorrectSubtypes_checkWrongSubtypes.txt', 'r')
setGenomeWrong = set()
for rline in fread:
	genomeID = rline.split('\t')[0]
	setGenomeWrong.add(genomeID)
fread.close()

print('wrong genomes...')
print(setGenomeWrong)	
	
fread = open('dictHostYearCountryGenomeFiltered0.99.txt', 'r')
for rline in fread:
	genomes = rline.strip().split('\t')[-2].strip().split(' ')
	hitFlag = 0
	for genomeID in genomes:
		if genomeID in setGenomeWrong:
			hitFlag = 1
	if hitFlag == 1:
		# print(genomes)
		# break
		setGenomeWrong = setGenomeWrong | set(genomes)
fread.close()

print('new wrong genomes...')
print(setGenomeWrong)


fwrite = open('wrongSubtypesAll2.txt', 'w')
for genomeID in setGenomeWrong:
	fwrite.write(genomeID+'\t')
fwrite.close()

# fread = open('genomeProteinCDScompleteUniqueVirNameGoodInfo.dat', 'r')
# allLines = fread.readlines()
# fread.close()

# fwrite = open('genomeProteinCDScompleteUniqueVirNameGoodInfoSubtypeRight.dat','w')
# for i in range(len(allLines)):
	# rline = allLines[i]
	# if rline[0] == '>':
		# genomeID = rline[1:rline.find('|')]
		# if genomeID not in setGenomeWrong:
			# for j in range(11):
				# fwrite.write(allLines[i+j])
		# else:
			# print(genomeID,'not subtype right...')
# fwrite.close()


# # lSegs = ['PB2','PB1','PA','HA','NP','NA','M1','M2','NS1','NS2']

# # for seg in lSegs:
	# # fasFile = seg+'.fasta.codon.fas'
	# # fread = open('')
	


