import os

def getDominantHost(dictHostCount):
	allHostCount = sum(dictHostCount.values())
	for host in dictHostCount:
		hostNum = dictHostCount[host]
		if hostNum/allHostCount > 0.9:
			return host
	return 'noDominant'


segs = ['PB2','PB1','PA','HA','NP','NA','MP','NS']

fread = open('isolateHostClassLocationRankedByTimeHostDetail-6.txt', 'r')
lineNo = 0
dictIsolateHost = dict()
dictIsolateLoc = dict()
dictIsolateSubtype = dict()
dictIsolateYear = dict()
dictIsolateInfo = dict()
dictIsolateTime = dict()

for rline in fread:
	lineNo += 1
	#isolateID	host_Classification	Continet_detail	subtype	year
	if lineNo > 1:
		lrline = rline.strip().split('\t')
		isolateID = lrline[0]
		dictIsolateInfo[isolateID] = rline.strip()
		dictIsolateHost[isolateID] = lrline[1]
		dictIsolateLoc[isolateID] = lrline[2]
		dictIsolateSubtype[isolateID] = lrline[3]
		dictIsolateYear[isolateID] = int(lrline[4])
		dictIsolateTime[isolateID] = lrline[5]
fread.close()


dictSegIsolateCluster = dict()
dictSegClusterNo = dict()
dictSegClusterPeriod = dict()

for seg in segs:
	fread = open(seg+'_clusters_isolates.txt', 'r')
	fwrite = open(seg+'_clusters_isolates_withhostsHostDetal.txt', 'w')
	fwrite4HostMarkers = open(seg+'_clusters_isolates_withhosts4HostMarkersHostDetail.txt', 'w')
	for rline in fread:
		lrline = rline.strip().split('\t')
		clusterSize = int(lrline[1])
		if clusterSize > 1: #select only the clusters with more than one seq.
			listIsolates = lrline[-1].strip().split(' ')
			dictHostCount = dict()
			for isolate in listIsolates:
				host = dictIsolateHost[isolate]
				if 'avian' in host.lower():
					host = 'Avian'
				if host in dictHostCount:
					dictHostCount[host] = dictHostCount[host] + 1
				else:
					dictHostCount[host] = 1
			fwrite.write(rline.strip()+'\t')
			fwrite.write(str(len(dictHostCount))+'\t')
			for host in dictHostCount:
				fwrite.write(host+':'+str(dictHostCount[host])+',')
			fwrite.write('\n')
			dominantHost = getDominantHost(dictHostCount)
			print(dominantHost)
			
			if dominantHost.lower() in ['human', 'swine', 'equine', 'avian']:
				fwrite4HostMarkers.write('>'+dominantHost+'|'+lrline[0]+'\n')
				for isolate in listIsolates:
					fwrite4HostMarkers.write(isolate+' ')
				fwrite4HostMarkers.write('\n')
			
			
			
		else:
			fwrite.write(rline.strip()+'\t'+'x\n')
			
			
	fread.close()
	fwrite.close()
	fwrite4HostMarkers.close()



# fwriteClusterNo = open('clusterNo.txt', 'w')
# fwriteClusterNo.write('segCluster\tclusterSize\tyearStart,yearEnd\tancestor\n')
# dictCluterSize = dict()
# dictClusterPeriod = dict()

# for seg in segs:
	# print(seg)
	
	# dictClusterIsolates = dict()
	
	
	# fread = open(seg+'_all_inClusters.txt', 'r')
	
	# dictIsolateCluster = dict()
	# dictClusterNo = dict()
	# clusterNo = 0
	
	
	# for rline in fread:
		# if rline.startswith('>>>'):
			# clusterAncestor = rline.strip()
			# clusterAncestorShort = clusterAncestor[:clusterAncestor.find('\t')]
			# if clusterAncestorShort in dictClusterNo:
				# print('repeated...', clusterAncestorShort)
			# else:
				# clusterNo += 1
				# dictClusterNo[clusterAncestorShort] = clusterNo
			# s = set()
			# dictClusterIsolates[clusterAncestorShort] = s
		# else:
			# if '|' not in rline:
				# isolateIDs = rline.strip().split('\t')
				# for isolateID in isolateIDs:
					# dictClusterIsolates[clusterAncestorShort].add(isolateID)
					# dictIsolateCluster[isolateID] = clusterAncestorShort
		
	# print(seg, len(dictClusterIsolates), 'clusters.')
	
	
	# fwrite = open(seg+'_clusters_isolates.txt', 'w')
	
	# for clusterAncestorShort in dictClusterIsolates:
		# #clusterAncestorShort = clusterAncestor[:clusterAncestor.find('\t')]
		# fwrite.write(clusterAncestorShort+'\t'+str(len(dictClusterIsolates[clusterAncestorShort]))+'\t')
		# for isolateID in dictClusterIsolates[clusterAncestorShort]:
			# fwrite.write(isolateID+' ')
		# fwrite.write('\n')
	# fwrite.close()
	
	# dictSegIsolateCluster[seg] = dictIsolateCluster
	# dictSegClusterNo[seg] = dictClusterNo
	
	# fread.close()
	
	# for clusterAncestorShort in dictClusterIsolates:
		# yearStart = 9999
		# yearEnd = 0
		# for isolateID in dictClusterIsolates[clusterAncestorShort]:
			# year = dictIsolateYear[isolateID]
			# if year < yearStart:
				# yearStart = year
			# if year > yearEnd:
				# yearEnd = year
		# clusterName = seg+'_'+str(dictClusterNo[clusterAncestorShort])
		# dictClusterPeriod[clusterName] = (yearStart, yearEnd)
		
	
	

	
	# for clusterAncestorShort in dictClusterNo:
		# clusterName = seg+'_'+str(dictClusterNo[clusterAncestorShort])
		# clusterPeriod = dictClusterPeriod[clusterName]
		# clusterPeriod2Write = str(clusterPeriod[0])+'\t'+str(clusterPeriod[1])
		# fwriteClusterNo.write(clusterName+'\t'+str(len(dictClusterIsolates[clusterAncestorShort]))+'\t'+clusterPeriod2Write+'\t'+clusterAncestorShort+'\n')
		# dictCluterSize[clusterName] = len(dictClusterIsolates[clusterAncestorShort])
# fwriteClusterNo.close()





# fwrite = open('isolatesCombinations.txt', 'w')
# fread = open('isolateHostClassLocationRankedByTime.txt', 'r')
# lineNo = 0
# dictIsolateCombination = dict()

# for rline in fread:
	# lineNo += 1
	# #isolateID	host_Classification	Continet_detail	subtype	year
	# if lineNo > 1:
		# lrline = rline.strip().split('\t')
		# isolateID = lrline[0]
		# fwrite.write(rline.strip()+'\t')
		# combination = []
		# for seg in segs:
			# clusterAncestorShort = dictSegIsolateCluster[seg][isolateID]
			# clusterNo = dictSegClusterNo[seg][clusterAncestorShort]
			# clusterName = seg+'_'+str(clusterNo)
			# fwrite.write(clusterName+'\t')
			# combination.append(clusterName)
		# fwrite.write('\n')
		# dictIsolateCombination[isolateID] = tuple(combination)
			
	# else:
		# fwrite.write(rline.strip()+'\t'+'PB2\tPB1\tPA\tHA\tNP\tNA\tMP\tNS\n')
		
		
# fread.close()

# fwrite.close()

# s = set(dictIsolateCombination.values())
# print(len(s), 'combinations for ', len(dictIsolateCombination), 'strians...' )


	
	
	
