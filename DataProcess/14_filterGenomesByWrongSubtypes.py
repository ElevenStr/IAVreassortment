import os

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

wrongIsolatelist = []
fread = open('../wrongSubtypesAll2.txt')
for line in fread:
    isolateID = line.rstrip().split('\t')
    for id in isolateID:
        wrongIsolatelist.append(id)

print(wrongIsolatelist)

for file in os.listdir():
    if file.endswith('.fasta.codon.fas.mafft.fas.correctSubtypes.fas'):
        giList,dictGiSeq = readFastaFile(file)
        fwrite = open(file + '.correctSubtypes2.fas', 'w')
        for gi in giList:
            isolateID = gi[0:gi.find('|')]
            if isolateID not in wrongIsolatelist:
                fwrite.write('>' + gi + '\n')
                fwrite.write(dictGiSeq[gi] + '\n')
        fwrite.close()