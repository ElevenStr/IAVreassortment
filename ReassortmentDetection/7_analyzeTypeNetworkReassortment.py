# coding=utf-8
from __future__ import division
import datetime
import time
import sys

#reload(sys)
#import xlrd
import types
import sys
import os
from datetime import datetime

def distTypesIJ(typeI, typeJ):
    simCount = 0
    for k in range(len(typeI)):
        if typeI[k] == typeJ[k]:
            simCount += 1
    return (8-simCount)

def strType(typeI):
    strT = '+'.join(list(typeI))
    return(strT)




dictTypeInfo = {}

fread = open('isolatesCombinationsWithNameHostDetail-6.txt', 'r')
lineNo = 0
for rline in fread:
    if lineNo > 0:
        lrline = rline.strip().split('\t')
        typeHere = tuple(lrline[7:])
        typeInfo = tuple(lrline[:7])
        if typeHere in dictTypeInfo:
            dictTypeInfo[typeHere].add(typeInfo)
        else:            
            s = set()
            s.add(typeInfo)
            dictTypeInfo[typeHere] = s
    lineNo += 1
fread.close()

lTypes = list(dictTypeInfo.keys())

lenTypes = len(lTypes)




dictTypeIndex = {}
dictIndexType = {}

fwrite = open('typeIndexInfoAvianDetailHostDetail.txt', 'w')
fwrite.write('typeIndex\tstrainCount\tCombination\tyears\thosts\tlocs\tlocStarts\tsubtype\tstrains\n')
typeIndex = 0
for typeK in lTypes:
    dictTypeIndex[typeK] = 'type'+str(typeIndex)
    dictIndexType['type'+str(typeIndex)] = typeK
    fwrite.write(dictTypeIndex[typeK]+'\t'+str(len(dictTypeInfo[typeK]))+'\t'+strType(typeK)+'\t')


    typeHostRange = []
    typeLocRange = []
    typeSubtypeRange = []
    typeYearRange = []

    for typeKinfo in dictTypeInfo[typeK]:
        host = typeKinfo[2]
        loc = typeKinfo[3]
        subtype = typeKinfo[4]
        year = typeKinfo[5]

        # if 'avian' in host.lower():
        #     host = 'Avian'
        if host not in typeHostRange:
            typeHostRange.append(host)
        if loc not in typeLocRange:
            typeLocRange.append(loc)
        if year not in typeYearRange:
            typeYearRange.append(year)
        if subtype not in typeSubtypeRange:
            typeSubtypeRange.append(subtype)

        typeYearRange.sort()

    fwrite.write(typeYearRange[0]+' '+typeYearRange[-1]+'\t')

    yearStart = typeYearRange[0]

    locStart = set()
    for typeKinfo in dictTypeInfo[typeK]:
        loc = typeKinfo[3]
        year = typeKinfo[5]
        if year == yearStart:
            locStart.add(loc)
    



    
    for host in typeHostRange:
        fwrite.write(host+' ')
    fwrite.write('\t')

    for loc in typeLocRange:
        fwrite.write(loc+' ')
    fwrite.write('\t')

    for loc in locStart:
        fwrite.write(loc+' ')
    fwrite.write('\t')

    for subtype in typeSubtypeRange:
        fwrite.write(subtype+' ')
    fwrite.write('\t')



    for typeKinfo in dictTypeInfo[typeK]:
        fwrite.write(typeKinfo[1]+'_'+typeKinfo[4]+'\t')
    fwrite.write('\n')

    typeIndex += 1

fwrite.close()


# sys.exit()


# fwrite = open('typeNetwork.txt', 'w')

# fwrite.write('typeI\ttypeJ\tdist\n')

dictTypeDist = {}

for i in range(lenTypes):
    for j in range(i):
        typeI = lTypes[i]
        typeJ = lTypes[j]
        distIJ = distTypesIJ(typeI,typeJ)
        #print(distIJ)
        if distIJ < 8:
            strTypeI = strType(typeI)
            strTypeJ = strType(typeJ)
            typeIndexI = dictTypeIndex[typeI]
            typeIndexJ = dictTypeIndex[typeJ]
            dictTypeDist[(typeIndexI,typeIndexJ)] = distIJ
            dictTypeDist[(typeIndexJ,typeIndexI)] = distIJ
            
#             fwrite.write(typeIndexI+'\t'+typeIndexJ+'\t'+str(distIJ)+'\n')

# fwrite.close()

dictVertexEdgeDist = {}
dictVertexEdgeDiff = {}

for (typeIndexI,typeIndexJ) in dictTypeDist:
    if typeIndexI in dictVertexEdgeDist:
        dictVertexEdgeDist[typeIndexI][(typeIndexI,typeIndexJ)] = dictTypeDist[(typeIndexI,typeIndexJ)]
    else:
        dictVertexEdgeDist[typeIndexI] = {}
        dictVertexEdgeDist[typeIndexI][(typeIndexI,typeIndexJ)] = dictTypeDist[(typeIndexI,typeIndexJ)]

def calTypeDiff(typeIndexI, typeIndexJ):
    typeI = dictIndexType[typeIndexI]
    typeJ = dictIndexType[typeIndexJ]

    strainsI = dictTypeInfo[typeI]
    strainsJ = dictTypeInfo[typeJ]

    minDiff = 999999

    for straini in strainsI:
        hosti = straini[2]
        loci = straini[3]
        yeari = int(straini[5])
        for strainj in strainsJ:
            hostj = strainj[2]
            locj = strainj[3]
            yearj = int(strainj[5])
            timeGap = abs(yeari-yearj)

            sumSame = 0

            if hosti == hostj:
                sumSame += 1
            elif 'avian' in hosti.lower() and 'avian' in hostj.lower():
                sumSame += 0.5
            if loci == locj:
                sumSame += 1
            if timeGap < 4:
                if timeGap < 2:
                    sumSame += 1
                else:
                    sumSame += 0.5
            
            sumDiff = 3-sumSame
            if minDiff > sumDiff:
                minDiff = sumDiff
    return minDiff


print('cal dictVertexEdgeDiff...')

for (typeIndexI,typeIndexJ) in dictTypeDist:
    typeDiff = calTypeDiff(typeIndexI, typeIndexJ)
    if typeIndexI in dictVertexEdgeDiff:
        dictVertexEdgeDiff[typeIndexI][(typeIndexI,typeIndexJ)] = typeDiff
    else:
        dictVertexEdgeDiff[typeIndexI] = {}
        dictVertexEdgeDiff[typeIndexI][(typeIndexI,typeIndexJ)] = typeDiff


fwrite = open('typeNetworkHostDetail.txt', 'w')

fwrite.write('typeI\ttypeJ\tdist\tdiff\n')

dictTypeDist = {}

for i in range(lenTypes):
    for j in range(i):
        typeI = lTypes[i]
        typeJ = lTypes[j]
        distIJ = distTypesIJ(typeI,typeJ)
        #print(distIJ)
        if distIJ < 8:
            strTypeI = strType(typeI)
            strTypeJ = strType(typeJ)
            typeIndexI = dictTypeIndex[typeI]
            typeIndexJ = dictTypeIndex[typeJ]
            dictTypeDist[(typeIndexI,typeIndexJ)] = distIJ
            dictTypeDist[(typeIndexJ,typeIndexI)] = distIJ
            
            fwrite.write(typeIndexI+'\t'+typeIndexJ+'\t'+str(distIJ)+'\t'+str(dictVertexEdgeDiff[typeIndexI][(typeIndexI,typeIndexJ)])+'\n')

fwrite.close()

def calCost(allLinesInfo, indexList):
    setLoc = set()
    setHost = set()
    listYear = []
    for indexI in indexList:
        loc = allLinesInfo[indexI][3]
        host = allLinesInfo[indexI][2]
        year = int(allLinesInfo[indexI][5])
        setLoc.add(loc)
        setHost.add(host)
        listYear.append(year)

    reassortCost = len(setLoc)*(len(indexList)-1) + len(setHost)*(len(indexList)-1) + (max(listYear)-min(listYear))
    return reassortCost

def analyzeReassortment(allLinesInfo, iHere):
    reassortCost = 9999999999
    hitFlag = 0
    reassortTuple = ()
    infoHere = allLinesInfo[iHere]
    yearHere = int(infoHere[5])
    setClustersHere = set(infoHere[7:])
    for i in range(iHere):
        infoI = allLinesInfo[i]
        yearI = int(infoI[5])
        if yearI < yearHere and yearI > yearHere - 6:
            setClusterI = set(infoI[7:])
            if (setClusterI & setClustersHere):
                for j in range(i):
                    infoJ = allLinesInfo[j]
                    yearJ = int(infoJ[5])
                    if yearJ < yearHere and yearJ > yearHere - 6:
                        setClusterJ = set(infoJ[7:])
                        if (setClusterJ & setClustersHere):
                            if setClustersHere.issubset(setClusterI|setClusterJ):
                                hitFlag = 1
                                reassortCostHere = calCost(allLinesInfo, [i, j, iHere])
                                if reassortCostHere < reassortCost :
                                    reassortCost = reassortCostHere
                                    reassortTuple = (i,j)
    if hitFlag == 1:
        print('two Reassortments')
        print(infoHere)
        print(allLinesInfo[reassortTuple[0]])
        print(allLinesInfo[reassortTuple[1]])        
        return(reassortTuple, reassortCost)
    
    for i in range(iHere):
        infoI = allLinesInfo[i]
        yearI = int(infoI[5])
        if yearI < yearHere and yearI > yearHere - 6:
            setClusterI = set(infoI[7:])
            if (setClusterI & setClustersHere):
                for j in range(i):
                    infoJ = allLinesInfo[j]
                    yearJ = int(infoJ[5])
                    if yearJ < yearHere and yearJ > yearHere - 6:
                        setClusterJ = set(infoJ[7:])
                        if (setClusterJ & setClustersHere):
                            for k in range(j):
                                infoK = allLinesInfo[k]
                                yearK = int(infoK[5])
                                if yearK < yearHere and yearK > yearHere - 6:
                                    setClusterK = set(infoK[7:])
                                    if (setClusterK & setClustersHere):
                                        if setClustersHere.issubset(setClusterI|setClusterJ|setClusterK):
                                            hitFlag = 1
                                            reassortCostHere = calCost(allLinesInfo, [i, j, k, iHere])
                                            if reassortCostHere < reassortCost :
                                                reassortCost = reassortCostHere
                                                reassortTuple = (i,j,k)
    if hitFlag == 1:
        print('three Reassortments')
        print(infoHere)
        print(allLinesInfo[reassortTuple[0]])
        print(allLinesInfo[reassortTuple[1]])
        print(allLinesInfo[reassortTuple[2]])
        return(reassortTuple, reassortCost)
    
    # for i in range(iHere):
    #     infoI = allLinesInfo[i]
    #     yearI = int(infoI[5])
    #     if yearI < yearHere and yearI > yearHere - 6:
    #         setClusterI = set(infoI[7:])
    #         for j in range(i):
    #             infoJ = allLinesInfo[j]
    #             yearJ = int(infoJ[5])
    #             if yearJ < yearHere and yearJ > yearHere - 6:
    #                 setClusterJ = set(infoJ[7:])
    #                 for k in range(j):
    #                     infoK = allLinesInfo[k]
    #                     yearK = int(infoK[5])
    #                     if yearK < yearHere and yearK > yearHere - 6:
    #                         setClusterK = set(infoK[7:])
    #                         for m in range(k):
    #                             infoM = allLinesInfo[m]
    #                             yearM = int(infoM[5])
    #                             if yearM < yearHere and yearM > yearHere - 6:
    #                                 setClusterM = set(infoM[7:])
    #                                 if setClustersHere.issubset(setClusterI|setClusterJ|setClusterK|setClusterM):
    #                                     hitFlag = 1
    #                                     reassortCostHere = calCost(allLinesInfo, [i, j, k, m, iHere])
    #                                     if reassortCostHere < reassortCost :
    #                                         reassortCost = reassortCostHere
    #                                         reassortTuple = (i,j,k,m)
    # if hitFlag == 1:
    #     print('four Reassortments')
    #     print(infoHere)
    #     print(allLinesInfo[reassortTuple[0]])
    #     print(allLinesInfo[reassortTuple[1]])
    #     print(allLinesInfo[reassortTuple[2]])
    #     print(allLinesInfo[reassortTuple[3]])
    #     return(reassortTuple, reassortCost)
    
    else:
        print('3-reassortment model still not working...considering 4???')
        print(allLinesInfo[iHere])
        return(reassortTuple,reassortCost)


dictIsolateInfo = dict()
dictTypeStartYear = dict()

segs = ['PB2','PB1','PA','HA','NP','NA','MP','NS']

fread = open('isolatesCombinationsWithNameHostDetail-6.txt', 'r')

print('Reading isolatesCombinationsWithName.txt file ...')
print('Assuming the isolatesCombinationsWithName.txt file has been ranked by year...')

listTypes = []
allLines = fread.readlines()
fread.close()

allLinesInfo = []

for i in range(1,len(allLines)):
    allLinesInfo.append(allLines[i].strip().split('\t'))


dictHard2FindReassorment = {}

dictYearClusters = dict()

dictTypeReassortmentCost = {}
dictTypeReassortmentTuple = {}

dictBeforeYearClusters = {}

fwrite = open('failedIn3ReassortmentingModeHostDetail.txt', 'w')

dictSameTypeYearLocHost = dict()

for i in range(len(allLinesInfo)):
    print(i, 'isolate being handled to see the reassortments...')

    lrline = allLinesInfo[i]
    typeHere = tuple(lrline[7:])
    typeInfo = tuple(lrline[:7])
    year = int(lrline[5])
    loc = lrline[3]
    host = lrline[2]
    if year in dictYearClusters:
        dictYearClusters[year] = dictYearClusters[year]|set(typeHere)
    else:
        dictYearClusters[year] = set(typeHere)

    if typeHere not in dictTypeStartYear:
        dictTypeStartYear[typeHere] = year
    else:
        if year < dictTypeStartYear[typeHere]:
            dictTypeStartYear[typeHere] = year
    
    isolateID = lrline[0]
    dictIsolateInfo[isolateID] = lrline

    if year not in dictBeforeYearClusters:
        dictBeforeYearClusters[year] = set()
        for yearOld in dictYearClusters:
            if yearOld < year and yearOld > year - 6:  #看看是不是在前面5年之内能够组合成现在的这个型？
                dictBeforeYearClusters[year] = dictBeforeYearClusters[year] | dictYearClusters[yearOld]
    if (typeHere,year,loc,host) in dictSameTypeYearLocHost:
        print('repeated analysis...pass...')
    elif year < dictTypeStartYear[typeHere] + 1 :
        if (set(typeHere)).issubset(dictBeforeYearClusters[year]):
            print('yes....')
            if typeHere in dictHard2FindReassorment and dictHard2FindReassorment[typeHere] == 1:
                print('Repeated not found in 3 reassortments...')
                fwrite.write('\t'.join(allLinesInfo[i]))
                fwrite.write('\t'+'r'+'\n')
            else:
                dictSameTypeYearLocHost[(typeHere,year,loc,host)] = 1
                nowtime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                print(nowtime)
                (reassortTuple,reassortCost) = analyzeReassortment(allLinesInfo, i)
                if reassortCost > 10000:
                    print('not found reassorment in 3 types...')
                    fwrite.write('\t'.join(allLinesInfo[i]))
                    fwrite.write('\n')
                    dictHard2FindReassorment[typeHere] = 1
                else:
                    dictHard2FindReassorment[typeHere] = 0
                    if typeHere in dictTypeReassortmentCost:
                        if dictTypeReassortmentCost[typeHere] > reassortCost:
                            dictTypeReassortmentCost[typeHere] = reassortCost
                            dictTypeReassortmentTuple[typeHere] = (i, reassortTuple)
                    else:
                        dictTypeReassortmentCost[typeHere] = reassortCost
                        dictTypeReassortmentTuple[typeHere] = (i, reassortTuple)
    
    
    if i % 100 == 0:
        print('wrinting dictTypeReassortmentTuple into temp Files...')
        fwriteTemp = open('reassortment3Histroy3Temp.txt', 'w')
        for typeHere in dictTypeReassortmentTuple:
            (ii, reassortTuple) = dictTypeReassortmentTuple[typeHere]
            fwriteTemp.write('>'+dictTypeIndex[typeHere]+'|'+str(len(reassortTuple))+'|'+str(dictTypeReassortmentCost[typeHere])+'\n')
            for item in allLinesInfo[ii]:
                fwriteTemp.write(item+'\t')
            fwriteTemp.write('\n')
            for x in reassortTuple:
                for item in allLinesInfo[x]:
                    fwriteTemp.write(item+'\t')
                fwriteTemp.write('\n')
        fwriteTemp.close()

fwrite.close()

fwrite = open('reassortmentHistroy3HostDetail.txt', 'w')
for typeHere in dictTypeReassortmentTuple:
    (i, reassortTuple) = dictTypeReassortmentTuple[typeHere]
    fwrite.write('>'+dictTypeIndex[typeHere]+'|'+str(len(reassortTuple))+'|'+str(dictTypeReassortmentCost[typeHere])+'\n')
    for item in allLinesInfo[i]:
        fwrite.write(item+'\t')
    fwrite.write('\n')
    for x in reassortTuple:
        for item in allLinesInfo[x]:
            fwrite.write(item+'\t')
        fwrite.write('\n')


fwrite.close()





    




# from scipy.sparse import csr_matrix
# from scipy.sparse.csgraph import minimum_spanning_tree

# matrixContent = []

# for i in range(lenTypes):
#     #print(i)
#     typeIndexI = 'type'+str(i)
#     arrContent = []
#     for j in range(lenTypes):
#         typeIndexJ = 'type'+str(j)
#         if j < i + 1 :
#             arrContent.append(0)
#         elif (typeIndexI, typeIndexJ) in dictTypeDist:
#             arrContent.append(dictTypeDist[(typeIndexI, typeIndexJ)])
#         else:
#             arrContent.append(0)
#     matrixContent.append(arrContent)

# X = csr_matrix(matrixContent)
# Tcsr = minimum_spanning_tree(X)

# MST_Array = Tcsr.toarray()

# fwrite = open('MST.txt', 'w')

# fwrite.write('typeIndexI\ttypeIndexJ\tdist\n')

# for i in range(lenTypes):
#     typeIndexI = 'type'+str(i)
#     for j in range(lenTypes):
#         if j > i:
#             typeIndexJ = 'type'+str(j)
#             if MST_Array[i][j] > 0:
#                 fwrite.write(typeIndexI+'\t'+typeIndexJ+'\t'+str(MST_Array[i][j])+'\n')


# fwrite.close()

# print('prim to MST...')

# listV = lTypes

# dictMarkedV = {}
# for typeK in lTypes:
#     typeIndexK = dictTypeIndex[typeK]
#     dictMarkedV[typeIndexK] = 0

# edges = dictTypeDist
# mst = []
# queueEdges = []


# def insertPriorityQueue(queueEdges, e):
#     if len(queueEdges) == 0:
#         queueEdges = [e]
#     else:
#         edgeDist = dictVertexEdgeDist[e[0]][e]
#         edgeDiff = dictVertexEdgeDiff[e[0]][e]
#         pos2Insert = 0
#         for eINQ in queueEdges:
#             if dictVertexEdgeDist[eINQ[0]][eINQ] < edgeDist:
#                 pos2Insert += 1
#         for eINQ in queueEdges[pos2Insert:]:
#             if dictVertexEdgeDist[eINQ[0]][eINQ] == edgeDist and not(edgeDiff < dictVertexEdgeDiff[eINQ[0]][eINQ]):
#                 pos2Insert += 1
#         queueEdges.insert(pos2Insert, e)
    
#     # for (a,b) in queueEdges:
#     #     print(a,b,dictVertexEdgeDist[a][(a,b)], dictVertexEdgeDiff[a][(a,b)])
#     # print('********************')

#     return queueEdges

# edges = dictTypeDist
# mst = []
# queueEdges = []

# type0 = lTypes[0]
# typeIndex0 = dictTypeIndex[type0]

# typeIndex0 = 'type461'

# dictMarkedV[typeIndex0] = 1
# for (typeIndexI,typeIndexJ) in dictVertexEdgeDist[typeIndex0]:
#     if typeIndexI != typeIndex0:
#         print('Wrong dictVertexEdge...')
#         sys.exit()
#     else:
#         pass
#     queueEdges = insertPriorityQueue(queueEdges, (typeIndexI,typeIndexJ))
    
# while True:
#     print('mstLength length...', len(mst))
#     print('queueEdges length...', len(queueEdges))
#     (typeIndexI,typeIndexJ) = queueEdges[0]
#     queueEdges.remove((typeIndexI,typeIndexJ))
#     # if dictMarkedV[typeIndexI] == 0:
#     #     mst.append((typeIndexI,typeIndexJ))
#     #     print('I:',typeIndexI,typeIndexJ,dictVertexEdgeDist[typeIndexI][(typeIndexI,typeIndexJ)],dictVertexEdgeDiff[typeIndexI][(typeIndexI,typeIndexJ)])
#     #     dictMarkedV[typeIndexI] = 1
#     #     for (typeIndexI,typeIndexII) in dictVertexEdgeDist[typeIndexI]:
#     #         if dictMarkedV[typeIndexII] == 0:
#     #             #queueEdges = insertPriorityQueue(queueEdges, (typeIndexI,typeIndexII))
                
#     #             edgeDist = dictVertexEdgeDist[e[0]][e]
#     #             edgeDiff = dictVertexEdgeDiff[e[0]][e]
#     #             pos2Insert = 0
#     #             for eINQ in queueEdges:
#     #                 if dictVertexEdgeDist[eINQ[0]][eINQ] < edgeDist:
#     #                     pos2Insert += 1
#     #             for eINQ in queueEdges[pos2Insert:]:
#     #                 if dictVertexEdgeDist[eINQ[0]][eINQ] == edgeDist and not(edgeDiff < dictVertexEdgeDiff[eINQ[0]][eINQ]):
#     #                     pos2Insert += 1
#     #             queueEdges.insert(pos2Insert, e) 

#     #         # if not((typeIndexI,typeIndexII) in queueEdges or (typeIndexII,typeIndexI) in queueEdges):
#     #         #     queueEdges = insertPriorityQueue(queueEdges, (typeIndexI,typeIndexII))

#     if dictMarkedV[typeIndexJ] == 0:
#         mst.append((typeIndexI,typeIndexJ))
#         print('J:',typeIndexI,typeIndexJ,dictVertexEdgeDist[typeIndexI][(typeIndexI,typeIndexJ)],dictVertexEdgeDiff[typeIndexI][(typeIndexI,typeIndexJ)])
#         dictMarkedV[typeIndexJ] = 1
#         for eINQ in queueEdges:
#             if dictMarkedV[eINQ[1]] == 1:
#                 queueEdges.remove(eINQ)

#         for (typeIndexJ,typeIndexJJ) in dictVertexEdgeDist[typeIndexJ]:
#             if dictMarkedV[typeIndexJJ] == 0:
#                 #queueEdges = insertPriorityQueue(queueEdges, (typeIndexJ,typeIndexJJ))
#                 e = (typeIndexJ,typeIndexJJ)
#                 edgeDist = dictVertexEdgeDist[e[0]][e]
#                 edgeDiff = dictVertexEdgeDiff[e[0]][e]
#                 pos2Insert = 0
#                 for eINQ in queueEdges:
#                     if dictVertexEdgeDist[eINQ[0]][eINQ] < edgeDist:
#                         pos2Insert += 1
#                 for eINQ in queueEdges[pos2Insert:]:
#                     if dictVertexEdgeDist[eINQ[0]][eINQ] == edgeDist and not(edgeDiff < dictVertexEdgeDiff[eINQ[0]][eINQ]):
#                         pos2Insert += 1
#                 queueEdges.insert(pos2Insert, e)

#             # if not((typeIndexJ,typeIndexJJ) in queueEdges or (typeIndexJJ,typeIndexJ) in queueEdges):
#             #     queueEdges = insertPriorityQueue(queueEdges, (typeIndexJ,typeIndexJJ))
        
#     if len(queueEdges) == 0:
#         break

# fwrite = open('mst_prime_withDiff_orphan461.txt', 'w')
# for (a,b) in mst:
#     fwrite.write(a+'\t'+b+'\t'+str(dictVertexEdgeDist[a][(a,b)])+'\t'+str(dictVertexEdgeDiff[a][(a,b)])+'\n')
    
    
# fwrite.close()

# fwrite = open('queueTemp.txt', 'w')

# for (a,b) in queueEdges:
#     fwrite.write(a+'\t'+b+'\t'+str(dictVertexEdgeDist[a][(a,b)])+'\t'+str(dictVertexEdgeDiff[a][(a,b)])+'\n')


# fwrite.close()

# fwrite = open('dictMarkedV.txt', 'w')

# for a in dictMarkedV:
#     if dictMarkedV[a] == 1:
#         fwrite.write(a+'\n')
# fwrite.close()
































#接下来需要将关系进行筛选，因为有些重配应该是不太合理的，我觉得。。。




