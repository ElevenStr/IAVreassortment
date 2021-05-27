#genomeProteinCDS.dat
import  sys

SeqidToInfo = dict()
fread1 = open('genomeset.txt','r')
lines = fread1.readlines()
fread1.close()
for i in range(len(lines)):
    line = lines[i].rstrip()
    if i > 0:
        info = line.split('\t')
        isolateID = info[0]
        SeqidToInfo[isolateID] = line

GeneidToProteinInfodict = dict()
BadGeneidList = list()
fread2 = open('influenza.dat','r')
lines = fread2.readlines()
fread2.close()
for i in range(len(lines)):
	line = lines[i].rstrip()
	info = line.split('\t')
	geneID = info[0]
	if len(info) > 2:
		GeneidToProteinInfodict[geneID] = line
	else:
		BadGeneidList.append(geneID)
print(len(BadGeneidList), ' geneid do not hava proteinCDS')

BadIsolateIDlist = list()
ProteinWronglist = list()
fwrite = open('genomeProteinCDS.dat','w')
for isolateID in SeqidToInfo:
	isolateInfo = SeqidToInfo[isolateID].split('\t')
	host = isolateInfo[2]
	subtype = isolateInfo[4]
	continent = isolateInfo[3]
	year = isolateInfo[5]
	name = isolateInfo[1]
	PB2id = isolateInfo[6]
	PB1id = isolateInfo[7]
	PAid = isolateInfo[8]
	HAid = isolateInfo[9]
	NPid = isolateInfo[10]
	NAid = isolateInfo[11]
	MPid = isolateInfo[12]
	NSid = isolateInfo[13]
	if PB2id in GeneidToProteinInfodict and PB1id in GeneidToProteinInfodict and PAid in GeneidToProteinInfodict and HAid in GeneidToProteinInfodict and NPid in GeneidToProteinInfodict and NAid in GeneidToProteinInfodict and MPid in GeneidToProteinInfodict and NSid in GeneidToProteinInfodict:
		PB2ProteinInfo = GeneidToProteinInfodict[PB2id].split('\t')
		if len(PB2ProteinInfo) > 3:
			print('here1')
			sys.exit()
		PB1ProteinInfo = GeneidToProteinInfodict[PB1id].split('\t')
		if len(PB1ProteinInfo) > 3:
			if ',' in PB1ProteinInfo[2] or ',' in PB1ProteinInfo[4]:
				if ',' in PB1ProteinInfo[2]:
					PB1Flag = 2;
				else:
					PB1Flag = 1;
			else:
				# print(isolateID)
				PB1_1 = PB1ProteinInfo[2]
				PB1_2 = PB1ProteinInfo[4]
				pos1_1 = PB1_1.find(':')
				if PB1_1[pos1_1+1] == '<':
					pos1_1 = pos1_1 + 1
				pos1_2 = PB1_1.find('-')
				pos2_1 = PB1_2.find(':')
				if PB1_2[pos2_1 + 1] == '<':
					pos2_1 = pos2_1 + 1
				pos2_2 = PB1_2.find('-')
				length1 = int(PB1_1[pos1_2+1:]) - int(PB1_1[pos1_1+1:pos1_2])
				length2 = int(PB1_2[pos2_2+1:]) - int(PB1_2[pos2_1+1:pos2_2])
				if length1 > length2:
					PB1Flag = 1
				else:
					PB1Flag = 2
		else:
			if ',' in PB1ProteinInfo[2]:
				PB1Flag = 0
			else:
				PB1Flag = 1
		PAProteinInfo = GeneidToProteinInfodict[PAid].split('\t')
		if len(PAProteinInfo) > 3:
			if ',' in PAProteinInfo[2]:
				PAFlag = 2
			else:
				PAFlag = 1
		else:
			if ',' in PAProteinInfo[2]:
				PAFlag = 0
			else:
				PAFlag = 1
		HAProteinInfo = GeneidToProteinInfodict[HAid].split('\t')
		if len(HAProteinInfo) > 3:
			print('here2')
			sys.exit()
		NPProteinInfo = GeneidToProteinInfodict[NPid].split('\t')
		if len(NPProteinInfo) > 3:
			print('here3')
			sys.exit()
		NAProteinInfo = GeneidToProteinInfodict[NAid].split('\t')
		if len(NAProteinInfo) > 3:
			print('here4')
			sys.exit()
		MPProteinInfo = GeneidToProteinInfodict[MPid].split('\t')
		if len(MPProteinInfo) > 3:
			if ',' in MPProteinInfo[2]:
				MPFlag = 2
			else:
				MPFlag = 1
		else:
			if ',' in MPProteinInfo[2]:
				MPFlag = 0
			else:
				MPFlag = 1
		NSProteinInfo = GeneidToProteinInfodict[NSid].split('\t')
		if len(NSProteinInfo) > 3:
			if ',' in NSProteinInfo[2]:
				NSFlag = 2
			else:
				NSFlag = 1
		else:
			if ',' in NSProteinInfo[2]:
				NSFlag = 0
			else:
				NSFlag = 1
		if PAFlag > 0 and MPFlag > 0 and NSFlag > 0 and PB1Flag > 0:
			fwrite.write('>' + isolateID + '|' + host + '\t' + subtype + '\t' + continent + '\t' + year + '\t' + name + '\n')
			fwrite.write('PB2' + '\t' + PB2ProteinInfo[1] + '\t' + PB2ProteinInfo[2] + '\n')
			if PB1Flag == 1:
				fwrite.write('PB1' + '\t' + PB1ProteinInfo[1] + '\t' + PB1ProteinInfo[2] + '\n')
			else:
				fwrite.write('PB1' + '\t' + PB1ProteinInfo[3] + '\t' + PB1ProteinInfo[4] + '\n')
			if PAFlag == 1:
				fwrite.write('PA' + '\t' + PAProteinInfo[1] + '\t' + PAProteinInfo[2] + '\n')
			else:
				fwrite.write('PA' + '\t' + PAProteinInfo[3] + '\t' + PAProteinInfo[4] + '\n')
			fwrite.write('HA' + '\t' + HAProteinInfo[1] + '\t' + HAProteinInfo[2] + '\n')
			fwrite.write('NP' + '\t' + NPProteinInfo[1] + '\t' + NPProteinInfo[2] + '\n')
			fwrite.write('NA' + '\t' + NAProteinInfo[1] + '\t' + NAProteinInfo[2] + '\n')
			if MPFlag == 1:
				fwrite.write('MP' + '\t' + MPProteinInfo[1] + '\t' + MPProteinInfo[2] + '\n')
			else:
				fwrite.write('MP' + '\t' + MPProteinInfo[3] + '\t' + MPProteinInfo[4] + '\n')
			if NSFlag == 1:
				fwrite.write('NS' + '\t' + NSProteinInfo[1] + '\t' + NSProteinInfo[2] + '\n')
			else:
				fwrite.write('NS' + '\t' + NSProteinInfo[3] + '\t' + NSProteinInfo[4] + '\n')
		else:
			BadIsolateIDlist.append(isolateID)
			ProteinWronglist.append(isolateID)
	else:
		BadIsolateIDlist.append(isolateID)

fwrite.close()
print(len(BadIsolateIDlist), ' IsolateID some segment do not have proteinCDS')
print(len(ProteinWronglist), ' IsolateID do not have complete protein')
print(ProteinWronglist)
