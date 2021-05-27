import collections

fread = open('../DataProcess/isolateHostClassLocation.txt','r')
lines = fread.readlines()
fread.close()
dictStrainHost = dict()
for i in range(len(lines)):
     if i > 0:
          line = lines[i].rstrip()
          info = line.split('\t')
          strain = info[0]
          host = info[1]
          dictStrainHost[strain] = host


idread = open('../DataProcess/seqsFiltered99/NS.fasta.codon.fas.correctSubtypes.fas.correctSubtypes2.fas','r')
idlist = list()
for line in idread:
     if line[0] == '>':
          id = line[1:line.find('|')]
          idlist.append(id)
idread.close()

fread = open('../DataProcess/countryContinetContinentDetail.txt','r')
countryContinentDetailDict = dict()
for line in fread:
     if line[0] != '>':
          info = line.rstrip().split('\t')
          country = info[0]
          continent = info[1]
          continentdetail = info[2]
          countryContinentDetailDict[country] = continentdetail
fread.close()

HostClass = ['human','swine','avian_domestic birds','avian_waterfowl','avian_shorebirds','avian_land birds','avian_unknown']

f1 = open('../DataProcess/genomeProteinCDScomplete.dat')
lines1 = f1.readlines()
f1.close()
idsToinfodict = dict()
idsTodatedict = dict()
for line in lines1:
     if line[0] == '>':
          infs = line.rstrip().replace('>','').replace('|','\t').split('\t')
          ids = infs[0]
          host = infs[1]
          hostdetail = dictStrainHost[ids]
          if hostdetail not in HostClass:
               hostdetail = 'others'
          loc = infs[3]
          continentdetail = countryContinentDetailDict[loc]
          subtype = infs[2]
          year = infs[4][0:4]
          date = infs[4]
          inf = ids + '\t' + hostdetail + '\t' + continentdetail + '\t' + subtype + '\t' + year + '\t' + date + '\n'
          if ids in idlist:
               idsTodatedict[ids] = int(date)
               idsToinfodict[ids] = inf

f2 = open('isolateHostClassLocationRankedByTimeHostDetail-6.txt','w')
f2.write('isolateID\thost_Classification\tContinet_detail\tsubtype\tyear\tCollecting Date\n')
sortedids = zip(idsTodatedict.values(),idsTodatedict.keys())
# print(sorted(sortedids))
for date_ids in sorted(sortedids):
     ids = date_ids[1]
     # print(date_ids)
     f2.write(idsToinfodict[ids])
f2.close()
          
          
