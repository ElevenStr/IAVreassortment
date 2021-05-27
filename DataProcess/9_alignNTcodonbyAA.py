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


def readFastaFile(fastaFile):
	giList = []
	fread = open(fastaFile, 'r')
	dictGiSeq = dict()
	lineNo = 0
	seq = ''
	for rline in fread:
		lineNo = lineNo + 1
		if rline[0] == '>':
			if lineNo > 1:
				dictGiSeq[gi] = seq
				seq = ''
			gi = rline[1:].rstrip()
			giList.append(gi)
		else:
			seq = seq + rline.strip()
	dictGiSeq[gi] = seq
	fread.close()
	return giList,dictGiSeq


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

for seg in lSegs:
	print(seg, ':')

	giList,dictGiSeq = readFastaFile(seg+'.fa.mafft.fas')
	
	dictGenomeSeq = dict()
	
	for gi in giList:
		genomeID = gi[:gi.find('|')]
		dictGenomeSeq[genomeID] = dictGiSeq[gi]
	
	fread = open(seg+'.fasta', 'r')
	allLines = fread.readlines()
	fread.close()
	fwrite = open(seg+'.fasta.codon.fas', 'w')
	
	for i in range(len(allLines)):
		rline = allLines[i]
		if rline[0] == '>':
			genomeID = rline[1:rline.find('|')]
			fwrite.write(rline.replace(' ',''))
			seqNT = allLines[i+1].strip()
			seqAA = dictGenomeSeq[genomeID]
			posNT = 0
			for pos in range(len(seqAA)):
				aa = seqAA[pos].upper()
				if aa in aaLigal:
					codonHere = seqNT[3*posNT:3*posNT+3]
					if aa != dictCodon[codonHere]:
						print(gi,seg,pos,codonHere,aa)
					else:
						fwrite.write(codonHere)
					posNT = posNT + 1
				else:
					fwrite.write('---')
			fwrite.write(seqNT[3*posNT:3*posNT+3])
			fwrite.write('\n')
	fwrite.close()

print('for HAs and NAs:')

hanaDir = 'HANA'

for seg in ['HA', 'NA']:
	print(seg)
	
	giList,dictGiSeq = readFastaFile(seg+'.fasta')
	
	dictGenomeNTSeq = dict()	
	for gi in giList:
		genomeID = gi[:gi.find('|')]
		dictGenomeNTSeq[genomeID] = dictGiSeq[gi]
	
	for aaFile in os.listdir(hanaDir):
		if aaFile.startswith(seg[0]) and aaFile.endswith('.mafft.fas'):
			print(aaFile)
			hntype = aaFile[:aaFile.find('.')]
			aaFilePath = os.path.join(hanaDir,aaFile)
			giList,dictGiSeq = readFastaFile(aaFilePath)
			
			fwrite = open(os.path.join(hanaDir, hntype+'.fasta.codon.fas'), 'w')
			
			dictGenomeAASeq = dict()
			for gi in giList:
				genomeID = gi[:gi.find('|')]
				dictGenomeAASeq[genomeID] = dictGiSeq[gi]
				fwrite.write('>'+gi.replace(' ','')+'\n')
				seqNT = dictGenomeNTSeq[genomeID]
				seqAA = dictGiSeq[gi]
				
				posNT = 0
				for pos in range(len(seqAA)):
					aa = seqAA[pos].upper()
					if aa in aaLigal:
						codonHere = seqNT[3*posNT:3*posNT+3]
						if aa != dictCodon[codonHere]:
							print(gi,seg,pos,codonHere,aa)
						else:
							fwrite.write(codonHere)
						posNT = posNT + 1
					else:
						fwrite.write('---')
				fwrite.write(seqNT[3*posNT:3*posNT+3])
				fwrite.write('\n')
			fwrite.close()
				
				
			
			
			
			
			
		
	


