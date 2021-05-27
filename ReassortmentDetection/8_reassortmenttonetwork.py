f1=open('reassortmentHistroy3HostDetail.txt')
f2=open('typeIndexInfoAvianDetailHostDetail.txt')
f3=open('typeNetworkHostDetail.txt')


###loading network
lines3=f3.readlines()
f3.close()
network={}
for line in lines3:
     infs=line.replace(' ','').strip('\n').split('\t')
     if len(infs)<2:
          continue
     key=infs[0]+'\t'+infs[1]
     network[key]=infs[0]+'\t'+infs[1]+'\t'+infs[2]+'\t'+infs[3]
     key=infs[1]+'\t'+infs[0]
     network[key]=infs[1]+'\t'+infs[0]+'\t'+infs[2]+'\t'+infs[3]
     
     

###loading types
lines2=f2.readlines()
f2.close()
typedict={}
for line in lines2:
     infs=line.replace(' ','').strip('\n').split('\t')
     if len(infs)<3:
          continue
     key=infs[2]
     typex=infs[0]
     typedict[key]=typex

###loading reassortment
lines1=f1.readlines()
f1.close()
reassortmentnetwork=[]
for i in range(len(lines1)):
     if '>' in lines1[i]:
          infs=lines1[i].split('|')
          if len(infs)<2:
               continue
          num=infs[1]
          targettype=infs[0].replace('>','')
          for j in range(int(num)):
               source=lines1[i+j+2].strip('\n').split('\t')[-9:-1]
               sourcekey='+'.join(source)
               sourcetype='null'
               if sourcekey not in typedict:
                    print(sourcekey)
               else:
                    sourcetype=typedict[sourcekey]
               networkkey=sourcetype+'\t'+targettype
               if networkkey in network:
                    reassortmentnetwork.append(network[networkkey])
               else:
                    print(networkkey)

f4=open('reassortmentnetworkHostDetail.txt','w')
f4.write('typeI\ttypeJ\tdist\tdiff\n')
for reassort in reassortmentnetwork:
     f4.write(reassort+'\n')
f4.close()
                    
          
     
