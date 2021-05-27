f1 = open('isolatesCombinationsHostDetail.txt')
lines1 = f1.readlines()
f1.close()
f2 = open('../genomeProteinCDScomplete.dat')
lines2 = f2.readlines()
f2.close()

seqname = dict()
for line in lines2:
     if line[0] == '>':
          info = line.rstrip().replace('|','\t').split('\t')
          seqid = info[0][1:]
          name = info[5]
          seqname[seqid] = name

f3 = open('isolatesCombinationsWithNameHostDetail-6.txt','w')
for  line in lines1:
     if 'isolateID' in line:
          line3 = line.replace('isolateID','isolateID\tisolateName')
          f3.write(line3)
     else:
          infs = line.split('\t')
          seqid = infs[0]
          if seqid not in seqname:
               print('error found ' + seqid)
          else:
               name = seqname[seqid]
               line3 = line.replace(seqid, seqid + '\t' + name)
               f3.write(line3)
f3.close()
          
          
