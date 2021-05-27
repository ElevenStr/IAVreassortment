import os

aaLigal = 'ACDEFGHIKLMNPQRSTVWY'
ntLigal = 'ATCG'

dictCodon = dict()
fread = open('codon.txt', 'r')
for rline in fread:
	lrline = rline.strip().split('\t')
	aa = lrline[1]
	codons = lrline[0].split(',')
	for codon in codons:
		if len(codon) != 3:
			print('error codon!')
		else:
			dictCodon[codon] = aa 

fread.close()


def translateCDS(seq):
	length = len(seq)
	if length%3 != 0:
		print('wrong codon number, cannot divided by 3!')
		sys.exit()
	pSeq = ''
	for i in range(length):
		start = i*3
		codon = seq[start:start+3]
		if codon in dictCodon:
			aa = dictCodon[codon]
			if aa not in aaLigal or aa == 'STOP':
				break
			else:
				pSeq = pSeq + aa
		else:
			break
	return pSeq

lSegs = ['PB2','PB1','PA','HA','NP','NA','MP','NS']

fread = open('genomeProteinCDS.dat', 'r')
dictGenomeInfo = dict()
dictGenomePAccCDS = dict()
sPAcc = set()
sCDS = set()

for rline in fread:
	if rline[0] == '>':
		genomeID = rline[1:rline.find('|')]
		genomeInfo = rline[rline.find('|')+1:].rstrip()
		dictGenomeInfo[genomeID] = genomeInfo
		dictGenomePAccCDS[genomeID] = dict()
	else:
		lrline = rline.strip().split('\t')
		seg = lrline[0]
		pAcc = lrline[1]
		sPAcc.add(pAcc)
		cds = lrline[2]
		sCDS.add(cds)
		dictGenomePAccCDS[genomeID][seg] = (pAcc, cds)
fread.close()

print(len(sPAcc), 'proteins.')
print(len(sCDS), 'cdss.')

dictPAccSeq = dict()
fread = open('influenza1Line.faa', 'r')
allLines = fread.readlines()
fread.close()
for i in range(len(allLines)):
	rline = allLines[i]
	if rline[0] == '>':
		pAcc = rline.split('|')[3]
		if pAcc in sPAcc:
			dictPAccSeq[pAcc] = allLines[i+1]

dictCDSSeq = dict()
fread = open('influenza1Line.cds', 'r')
allLines = fread.readlines()
fread.close()
countTemp = 0
for i in range(len(allLines)):
	rline = allLines[i]
	if rline[0] == '>':
		cds = rline[1:rline.rfind('|')]
		if cds in sCDS and 'complete' in rline.lower():
			seq = allLines[i+1].strip()
			if len(seq) % 3 != 0:
				countTemp = countTemp + 1
				print(seq)
			else:
				seqAA = translateCDS(seq)
				if len(seqAA) > len(seq)/3-2:
					dictCDSSeq[cds] = seq

print(len(dictCDSSeq), 'completed cdss.')

sUnCompletedGenome = set()
for genomeID in dictGenomePAccCDS:
	genomePAccCDS = dictGenomePAccCDS[genomeID]
	for seg in genomePAccCDS:
		cds = genomePAccCDS[seg][1]
		if cds not in dictCDSSeq:
			sUnCompletedGenome.add(genomeID)
			break
print(len(sUnCompletedGenome), 'genomes are not completed.')


for genomeID in sUnCompletedGenome:
	dictGenomePAccCDS.pop(genomeID)


sLab = set()
for genomeID in dictGenomePAccCDS:
	genomeInfo = dictGenomeInfo[genomeID].strip().split('\t')
	if 'lab' in genomeInfo[2].lower():
		sLab.add(genomeID)

print(len(sLab), 'genomes are in LAB.')


for genomeID in sLab:
	dictGenomePAccCDS.pop(genomeID)
print(len(dictGenomePAccCDS), 'genome are left may be completed.')

NotEnaughLengthlist = ['1043795','1022015','1023170','1021506','24711']
for genomeID in NotEnaughLengthlist:
	dictGenomePAccCDS.pop(genomeID)
print(len(dictGenomePAccCDS), 'genome are left may be completed and enaugh length.')

fwrite = open('genomeProteinCDScomplete.dat', 'w')
for genomeID in dictGenomePAccCDS:
	fwrite.write('>'+genomeID+'|'+dictGenomeInfo[genomeID]+'\n')
	genomePAccCDS = dictGenomePAccCDS[genomeID]
	for seg in lSegs:
		fwrite.write(seg+'\t'+genomePAccCDS[seg][0]+'\t'+genomePAccCDS[seg][1]+'\n')
		
		
fwrite.close()


for genomeID in dictGenomePAccCDS:
	genomePAccCDS = dictGenomePAccCDS[genomeID]
	genomeInfo = dictGenomeInfo[genomeID]
	strainName = genomeInfo.split('\t')[4]
	strainName = strainName.replace(')', '').replace('(','_')
	for seg in genomePAccCDS:
		fwrite = open(seg+'.fasta', 'a')
		fwrite.write('>'+genomeID+'|'+strainName+'\n')
		cds = genomePAccCDS[seg][1]
		seq = dictCDSSeq[cds]
		fwrite.write(seq+'\n')
		fwrite.close()


