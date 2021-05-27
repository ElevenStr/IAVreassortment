from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import sys


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

typeIndex = 0
for typeK in lTypes:
    dictTypeIndex[typeK] = 'type'+str(typeIndex)
    dictIndexType['type'+str(typeIndex)] = typeK
    typeIndex += 1

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
            timeGap = abs(yeari - yearj)

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

            sumDiff = 3 - sumSame
            if minDiff > sumDiff:
                minDiff = sumDiff
    return minDiff

print('cal dictVertexEdgeDiff...')

for (typeIndexI, typeIndexJ) in dictTypeDist:
    typeDiff = calTypeDiff(typeIndexI, typeIndexJ)
    if typeIndexI in dictVertexEdgeDiff:
        dictVertexEdgeDiff[typeIndexI][(typeIndexI, typeIndexJ)] = typeDiff
    else:
        dictVertexEdgeDiff[typeIndexI] = {}
        dictVertexEdgeDiff[typeIndexI][(typeIndexI, typeIndexJ)] = typeDiff

matrixContent = []

for i in range(lenTypes):
    #print(i)
    typeIndexI = 'type'+str(i)
    arrContent = []
    for j in range(lenTypes):
        typeIndexJ = 'type'+str(j)
        if j < i + 1 :
            arrContent.append(0)
        elif (typeIndexI, typeIndexJ) in dictTypeDist:
            arrContent.append(dictTypeDist[(typeIndexI, typeIndexJ)])
        else:
            arrContent.append(0)
    matrixContent.append(arrContent)

X = csr_matrix(matrixContent)
Tcsr = minimum_spanning_tree(X)

MST_Array = Tcsr.toarray()

fwrite = open('MSTHostDetail.txt', 'w')

fwrite.write('typeIndexI\ttypeIndexJ\tdist\n')

for i in range(lenTypes):
    typeIndexI = 'type'+str(i)
    for j in range(lenTypes):
        if j > i:
            typeIndexJ = 'type'+str(j)
            if MST_Array[i][j] > 0:
                fwrite.write(typeIndexI+'\t'+typeIndexJ+'\t'+str(MST_Array[i][j])+'\n')
fwrite.close()

print('prim to MST...')

listV = lTypes

dictMarkedV = {}
for typeK in lTypes:
    typeIndexK = dictTypeIndex[typeK]
    dictMarkedV[typeIndexK] = 0

edges = dictTypeDist
mst = []
queueEdges = []


def insertPriorityQueue(queueEdges, e):
    if len(queueEdges) == 0:
        queueEdges = [e]
    else:
        edgeDist = dictVertexEdgeDist[e[0]][e]
        edgeDiff = dictVertexEdgeDiff[e[0]][e]
        pos2Insert = 0
        for eINQ in queueEdges:
            if dictVertexEdgeDist[eINQ[0]][eINQ] < edgeDist:
                pos2Insert += 1
        for eINQ in queueEdges[pos2Insert:]:
            if dictVertexEdgeDist[eINQ[0]][eINQ] == edgeDist and not(edgeDiff < dictVertexEdgeDiff[eINQ[0]][eINQ]):
                pos2Insert += 1
        queueEdges.insert(pos2Insert, e)

    # for (a,b) in queueEdges:
    #     print(a,b,dictVertexEdgeDist[a][(a,b)], dictVertexEdgeDiff[a][(a,b)])
    # print('********************')

    return queueEdges

edges = dictTypeDist
mst = []
queueEdges = []

type0 = lTypes[0]
typeIndex0 = dictTypeIndex[type0]

type0 = lTypes[0]

dictMarkedV[typeIndex0] = 1
for (typeIndexI,typeIndexJ) in dictVertexEdgeDist[typeIndex0]:
    if typeIndexI != typeIndex0:
        print('Wrong dictVertexEdge...')
        sys.exit()
    else:
        pass
    queueEdges = insertPriorityQueue(queueEdges, (typeIndexI,typeIndexJ))

while True:
    print('mstLength length...', len(mst))
    print('queueEdges length...', len(queueEdges))
    (typeIndexI,typeIndexJ) = queueEdges[0]
    queueEdges.remove((typeIndexI,typeIndexJ))
    # if dictMarkedV[typeIndexI] == 0:
    #     mst.append((typeIndexI,typeIndexJ))
    #     print('I:',typeIndexI,typeIndexJ,dictVertexEdgeDist[typeIndexI][(typeIndexI,typeIndexJ)],dictVertexEdgeDiff[typeIndexI][(typeIndexI,typeIndexJ)])
    #     dictMarkedV[typeIndexI] = 1
    #     for (typeIndexI,typeIndexII) in dictVertexEdgeDist[typeIndexI]:
    #         if dictMarkedV[typeIndexII] == 0:
    #             #queueEdges = insertPriorityQueue(queueEdges, (typeIndexI,typeIndexII))

    #             edgeDist = dictVertexEdgeDist[e[0]][e]
    #             edgeDiff = dictVertexEdgeDiff[e[0]][e]
    #             pos2Insert = 0
    #             for eINQ in queueEdges:
    #                 if dictVertexEdgeDist[eINQ[0]][eINQ] < edgeDist:
    #                     pos2Insert += 1
    #             for eINQ in queueEdges[pos2Insert:]:
    #                 if dictVertexEdgeDist[eINQ[0]][eINQ] == edgeDist and not(edgeDiff < dictVertexEdgeDiff[eINQ[0]][eINQ]):
    #                     pos2Insert += 1
    #             queueEdges.insert(pos2Insert, e)

    #         # if not((typeIndexI,typeIndexII) in queueEdges or (typeIndexII,typeIndexI) in queueEdges):
    #         #     queueEdges = insertPriorityQueue(queueEdges, (typeIndexI,typeIndexII))

    if dictMarkedV[typeIndexJ] == 0:
        mst.append((typeIndexI,typeIndexJ))
        print('J:',typeIndexI,typeIndexJ,dictVertexEdgeDist[typeIndexI][(typeIndexI,typeIndexJ)],dictVertexEdgeDiff[typeIndexI][(typeIndexI,typeIndexJ)])
        dictMarkedV[typeIndexJ] = 1
        for eINQ in queueEdges:
            if dictMarkedV[eINQ[1]] == 1:
                queueEdges.remove(eINQ)

        for (typeIndexJ,typeIndexJJ) in dictVertexEdgeDist[typeIndexJ]:
            if dictMarkedV[typeIndexJJ] == 0:
                #queueEdges = insertPriorityQueue(queueEdges, (typeIndexJ,typeIndexJJ))
                e = (typeIndexJ,typeIndexJJ)
                edgeDist = dictVertexEdgeDist[e[0]][e]
                edgeDiff = dictVertexEdgeDiff[e[0]][e]
                pos2Insert = 0
                for eINQ in queueEdges:
                    if dictVertexEdgeDist[eINQ[0]][eINQ] < edgeDist:
                        pos2Insert += 1
                for eINQ in queueEdges[pos2Insert:]:
                    if dictVertexEdgeDist[eINQ[0]][eINQ] == edgeDist and not(edgeDiff < dictVertexEdgeDiff[eINQ[0]][eINQ]):
                        pos2Insert += 1
                queueEdges.insert(pos2Insert, e)

            # if not((typeIndexJ,typeIndexJJ) in queueEdges or (typeIndexJJ,typeIndexJ) in queueEdges):
            #     queueEdges = insertPriorityQueue(queueEdges, (typeIndexJ,typeIndexJJ))

    if len(queueEdges) == 0:
        break

fwrite = open('mst_prime_withDiff_orphan461HostDetail.txt', 'w')
for (a,b) in mst:
    fwrite.write(a+'\t'+b+'\t'+str(dictVertexEdgeDist[a][(a,b)])+'\t'+str(dictVertexEdgeDiff[a][(a,b)])+'\n')


fwrite.close()

fwrite = open('queueTempHostDetail.txt', 'w')

for (a,b) in queueEdges:
    fwrite.write(a+'\t'+b+'\t'+str(dictVertexEdgeDist[a][(a,b)])+'\t'+str(dictVertexEdgeDiff[a][(a,b)])+'\n')


fwrite.close()

fwrite = open('dictMarkedVHostDetail.txt', 'w')

for a in dictMarkedV:
    if dictMarkedV[a] == 1:
        fwrite.write(a+'\n')
fwrite.close()