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

for seg in lSegs:
	print(seg, ':')
	fread = open(seg+'.fasta', 'r')
	allLines = fread.readlines()
	fread.close()
	fwrite = open(seg+'.fa', 'w')
	for i in range(len(allLines)):
		rline = allLines[i]
		if rline[0] == '>':
			fwrite.write(rline)
			seq = allLines[i+1].strip()
			seqAA = translateCDS(seq)
			fwrite.write(seqAA)
			fwrite.write('\n')
	fwrite.close()
			


