inforead = open('genomeProteinCDScomplete.dat','r')
lines = inforead.readlines()
inforead.close()
HSubtypeToidDict = dict()
NSubtypeToidDict = dict()
for i in range(len(lines)):
    line = lines[i].rstrip()
    if line[0] == '>':
        id = line[1:line.find('|')]
        subtype = line.split('\t')[1]
        htype = subtype[0:subtype.find('N')]
        ntype = subtype[subtype.find('N'):]
        if htype in HSubtypeToidDict:
            HSubtypeToidDict[htype].append(id)
        else:
            HSubtypeToidDict[htype] = list()
            HSubtypeToidDict[htype].append(id)
        if ntype in NSubtypeToidDict:
            NSubtypeToidDict[ntype].append(id)
        else:
            NSubtypeToidDict[ntype] = list()
            NSubtypeToidDict[ntype].append(id)

HAfread = open('HA.fa','r')
lines = HAfread.readlines()
HAfread.close()
HAidToseqdict = dict()
idToInfodict = dict()
for i in range(len(lines)):
    line = lines[i].rstrip()
    if '>' in line:
        id = line[1:line.find('|')]
        HAidToseqdict[id] = lines[i+1].rstrip()
        idToInfodict[id] = line

NAfread = open('NA.fa','r')
lines = NAfread.readlines()
NAfread.close()
NAidToseqdict = dict()
for i in range(len(lines)):
    line = lines[i].rstrip()
    if '>' in line:
        id = line[1:line.find('|')]
        NAidToseqdict[id] = lines[i+1].rstrip()


for htype in HSubtypeToidDict:
    filename = htype + '.fa'
    fwrite = open(filename,'w')
    for id in HSubtypeToidDict[htype]:
        # print(id)
        fwrite.write(idToInfodict[id] + '\n')
        fwrite.write(HAidToseqdict[id] + '\n')
    fwrite.close()

for ntype in NSubtypeToidDict:
    filename = ntype + '.fa'
    fwrite = open(filename,'w')
    for id in NSubtypeToidDict[ntype]:
        fwrite.write(idToInfodict[id] + '\n')
        fwrite.write(NAidToseqdict[id] + '\n')
    fwrite.close()

