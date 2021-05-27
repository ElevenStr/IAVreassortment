Segdict={'1':'PB2',
         '2':'PB1',
         '3':'PA',
         '4':'HA',
         '5':'NP',
         '6':'NA',
         '7':'MP',
         '8':'NS'}



def IfSubtypeHxNy(Subtype):
    if 'H' in Subtype and 'N' in Subtype:
        htype = Subtype[1:Subtype.find('N')]
        ntype = Subtype[Subtype.find('N')+1:]
        if htype.isdigit() and htype[0] != '0' and ntype.isdigit() and ntype[0] != '0':
            return True
        else:
            return False
    else:
        return False


def getDate(dateStrain):
    if len(dateStrain) > 1:
        if '/' in dateStrain:
            dateStrain = dateStrain.replace('/', '')
            if len(dateStrain) % 2 != 0:
                print('wrong dateStrain', dateStrain)
                sys.exit()
            else:
                if len(dateStrain) == 6:
                    dateStrain = dateStrain + '00'
        else:
            dateStrain = '00000000'
            # if dateStrain.isdigit():
            #     dateStrain = dateStrain + '0000'
            # else:
            #     dateStrain = '00000000'
            #     # print('00000000')
    else:
        dateStrain = '00000000'
        # print('00000000')
    if len(dateStrain) != 8:
        print(dateStrain, 'not ligal')
        sys.exit()
    return dateStrain


fread = open('genomeset.dat','r')
fwrite = open('genomeset.txt','w')
fwrite.write('isolateID\tname\thost\tContinet\tsubtype\tyear\tPB2id\tPB1id\tPAid\tHAid\tNPid\tNAid\tMPid\tNSid\n')
BadIsolateNameSet = set()
NiceIsolateNameSet = set()
InfluenzaBSet = set()
IsolateNameToInfodict = dict()
lines = fread.readlines()
fread.close()
SamenameSameID = set()
for i in range(len(lines)):
    line = lines[i]
    lineSplit = line.strip().split('\t')
    if 'Influenza A virus' in line:
        linelen = len(lineSplit)
        SeqNumber = lineSplit[0]
        host = lineSplit[1]
        Segment = Segdict[lineSplit[2]]
        Subtype = lineSplit[3]
        Continet = lineSplit[4]
        year = getDate(lineSplit[5])
        # namepos1 = line.find('(')
        # namepos2 = line.rfind(')')
        # name = line[namepos1 + 1:namepos2]
        name = lineSplit[7]
        # print(name)
        name = name.replace('Influenza A virus (', '').replace('Influenza A virus(', '').replace('Influenza A virus ','').replace('Influenza A virus','').replace(')', '').replace('(','_')
        # print(name)
        if '.' in name:
            dotpos1 = name.find('.')
            dotpos2 = name.rfind('.')
            namereplace = name[0:dotpos1]+'/'+name[dotpos2+1:len(name)]
            name = namereplace
        isolateID = lineSplit[linelen - 1]
        # print(SeqNumber, host, Segment, Subtype, Continet, year, name, isolateID)
        if IfSubtypeHxNy(Subtype) and 'NON' not in line and host != '' and Continet != '' and year != '00000000':
            if name in IsolateNameToInfodict:
                Lastline = lines[i-1]
                LastlineSplit = Lastline.strip().split('\t')
                if lineSplit[2] == LastlineSplit[2]:
                    SamenameSameID.add(name+'_'+isolateID)
                    continue
                else:
                    IsolateNameToInfodict[name].append(SeqNumber)
            else:
                IsolateNameToInfodict[name] = []
                IsolateNameToInfodict[name].append(isolateID)
                IsolateNameToInfodict[name].append(name)
                IsolateNameToInfodict[name].append(host)
                IsolateNameToInfodict[name].append(Continet)
                IsolateNameToInfodict[name].append(Subtype)
                IsolateNameToInfodict[name].append(year)
                IsolateNameToInfodict[name].append(SeqNumber)
            if i == len(lines)-1:
                if len(IsolateNameToInfodict) == 1:
                    if len(IsolateNameToInfodict[name]) == 14:
                        NiceIsolateNameSet.add(name+'_'+isolateID)
                        fwrite.write(isolateID+'\t'+name+'\t'+host+'\t'+Continet+'\t'+Subtype+'\t'+year+'\t'+IsolateNameToInfodict[name][6]+'\t'+IsolateNameToInfodict[name][7]+'\t'+IsolateNameToInfodict[name][8]+'\t'+IsolateNameToInfodict[name][9]+'\t'+IsolateNameToInfodict[name][10]+'\t'+IsolateNameToInfodict[name][11]+'\t'+IsolateNameToInfodict[name][12]+'\t'+IsolateNameToInfodict[name][13]+'\n')
                    else:
                        BadIsolateNameSet.add(name+'_'+isolateID)
                else:
                    minYearIsolateName = ''
                    minYear = 21001231
                    for IsolateName in IsolateNameToInfodict:
                        if int(IsolateNameToInfodict[IsolateName][5]) < minYear and len(IsolateNameToInfodict[IsolateName]) == 14:
                            minYear = int(IsolateNameToInfodict[IsolateName][5])
                            minYearIsolateName = IsolateName
                    for IsolateName in IsolateNameToInfodict:
                        if IsolateName == minYearIsolateName:
                            NiceIsolateNameSet.add(IsolateName+'_'+isolateID)
                            fwrite.write(isolateID + '\t' + name + '\t' + host + '\t' + Continet + '\t' + Subtype + '\t' + year + '\t' +IsolateNameToInfodict[IsolateName][6] + '\t' + IsolateNameToInfodict[IsolateName][7] + '\t' +IsolateNameToInfodict[IsolateName][8] + '\t' + IsolateNameToInfodict[IsolateName][9] + '\t' +IsolateNameToInfodict[IsolateName][10] + '\t' + IsolateNameToInfodict[IsolateName][11] + '\t' +IsolateNameToInfodict[IsolateName][12] + '\t' + IsolateNameToInfodict[IsolateName][13] + '\n')
                        else:
                            BadIsolateNameSet.add(IsolateName+'_'+isolateID)
            elif lines[i+1] in ['\n', '\r\n']:
                if len(IsolateNameToInfodict) == 1:
                    if len(IsolateNameToInfodict[name]) == 14:
                        NiceIsolateNameSet.add(name+'_'+isolateID)
                        fwrite.write(isolateID + '\t' + name + '\t' + host + '\t' + Continet + '\t' + Subtype + '\t' + year + '\t' + IsolateNameToInfodict[name][6] + '\t' + IsolateNameToInfodict[name][7] + '\t' + IsolateNameToInfodict[name][8] + '\t' + IsolateNameToInfodict[name][9] + '\t' + IsolateNameToInfodict[name][10] + '\t' + IsolateNameToInfodict[name][11] + '\t' + IsolateNameToInfodict[name][12] + '\t' + IsolateNameToInfodict[name][13] + '\n')
                    else:
                        BadIsolateNameSet.add(name+'_'+isolateID)
                else:
                    minYearIsolateName = ''
                    minYear = 21001231
                    for IsolateName in IsolateNameToInfodict:
                        if int(IsolateNameToInfodict[IsolateName][5]) < minYear and len(IsolateNameToInfodict[IsolateName]) == 14:
                            minYear = int(IsolateNameToInfodict[IsolateName][5])
                            minYearIsolateName = IsolateName
                    for IsolateName in IsolateNameToInfodict:
                        if IsolateName == minYearIsolateName:
                            NiceIsolateNameSet.add(IsolateName+'_'+isolateID)
                            fwrite.write(isolateID + '\t' + name + '\t' + host + '\t' + Continet + '\t' + Subtype + '\t' + year + '\t' +IsolateNameToInfodict[IsolateName][6] + '\t' + IsolateNameToInfodict[IsolateName][7] + '\t' +IsolateNameToInfodict[IsolateName][8] + '\t' + IsolateNameToInfodict[IsolateName][9] + '\t' +IsolateNameToInfodict[IsolateName][10] + '\t' + IsolateNameToInfodict[IsolateName][11] + '\t' +IsolateNameToInfodict[IsolateName][12] + '\t' + IsolateNameToInfodict[IsolateName][13] + '\n')
                        else:
                            BadIsolateNameSet.add(IsolateName+'_'+isolateID)
        else:
            BadIsolateNameSet.add(name+'_'+isolateID)
    elif 'Influenza B virus' in line:
        lineSplit = line.strip().split('\t')
        isolateID = lineSplit[10]
        namepos1 = line.find('(')
        namepos2 = line.rfind(')')
        name = line[namepos1 + 1:namepos2]
        if '.' in name:
            dotpos1 = name.find('.')
            dotpos2 = name.rfind('.')
            namereplace = name[0:dotpos1]+'/'+name[dotpos2+1:len(name)]
            name = namereplace
        InfluenzaBSet.add(isolateID)
    elif line in ['\n', '\r\n']:
        IsolateNameToInfodict.clear()


fwrite.close()
print('the number of NiceIsolateName is:', len(NiceIsolateNameSet))
print('the number of SamenameSameID is', len(SamenameSameID))
print('the number of BadIsolateName is', len(BadIsolateNameSet))
print('the number of Influenza B is', len(InfluenzaBSet))

print(len(NiceIsolateNameSet & SamenameSameID))
print(len(SamenameSameID))

