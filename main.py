#================================
# Mitchell Aschmeyer and Riley Collins
#
#
#================================
import os
import numpy
from compbio import *

def dictBuilder():
    alignments = {}
    currentDirect = os.getcwd()
    fileNames = ['\England (DEC052020)\HCOV19-ENGLAND-051220-A.fasta', '\England (DEC082020)\HCOV19-ENGLAND-081220-A.fasta', '\England (NOV032020)\HCOV19-ENGLAND-031120-A.fasta', '\England (NOV102020)\HCOV19-ENGLAND-101120-A.fasta', '\England (NOV272020)\HCOV19-ENGLAND-271120-A.fasta']
    for name in fileNames:
        file = open(currentDirect + name, "r")
        file = file.readlines()
        header = file[0]
        key = ''
        for line in file[1:]:
            key = key + line.strip('\n')
        alignments[header] = key
    return(alignments)

def findR(D):
    rs = []
    for row in D:
        rs.append(sum(row))
    return rs

def distanceMatrix(D, i, j):
    dij = D[i][j]
    newR = []
    for col in range(len(D[0]) - 1):
        newR.append((D[i][col] + D[j][col] - dij/2))

    for row in range(len(D) - 1):
        del(D[row][j])
        del(D[row][i])
        if row < len(newR):
            D[row].append(newR[row])
    if i > j:
        del(D[i])
        del(D[j])
    else:
        del(D[j])
        del(D[i])

    D.append(newR)
    return D

def updateD(D, i, j):
    dij = D[i][j]

    newr = []
    if i > j:
        del(D[i])
        del(D[j])
        for row in range(len(D)):
            value = ((D[row][i] + D[row][j] -dij)/2)
            D[row].append(value)
            del(D[row][i])
            del(D[row][j])
            newr.append(value)
    else:
        del(D[j])
        del(D[i])
        for row in range(len(D)):
            value = ((D[row][i] + D[row][j] -dij)/2)
            D[row].append(value)
            del(D[row][j])
            del(D[row][i])
            newr.append(value)
    newr.append(0)
    D.append(newr)

    return D


def neighborJoining(D, n, names):
    if n == 2:
        return(names)
    M = numpy.zeros((n, n))
    rs = findR(D)
    for row in range(n):
        ra = rs[row]
        for col in range(n):
            if row < col:
                rb = rs[col]
                M[row][col] = (n - 2) * D[row][col] - ra - rb
    min_flat_index = numpy.argmin(M)
    i, j = numpy.unravel_index(min_flat_index, (n, n))

    #establish lengths
    lengthB = D[i][j]/2 - ((rs[i] - rs[j])/(2*(n - 2)))
    lengthA = D[i][j]/2 - ((rs[j] - rs[i])/(2*(n - 2)))

    #update names array
    nameA, nameB = names[i], names[j]
    if type(nameA) != tuple:
        nA = nameA + ":" + str(lengthA)
    else:
        nA = nameA
    if type(nameB) != tuple:
        nB = nameB + ":" + str(lengthB)
    else:
        nB = nameB
    names.append((nA, nB))
    names.remove(nameA)
    names.remove(nameB)


    D = updateD(D, i, j)


    #ERROR IS OCCURING WITH ADDING
    return(neighborJoining(D, n - 1, names))



def main():
    #alignments = dictBuilder()
    #D = (distanceMatrix(alignments))
    names = ['dog', 'bear', 'raccoon', 'weasel', 'seal', 'sea lion', 'cat', 'monkey']
    example1 = [[0, 32, 48, 51, 50, 48, 98, 148],
                    [32, 0, 26, 34, 29, 33, 84, 136],
                    [48, 26, 0, 42, 44, 44, 92, 152],
                    [51, 34, 42, 0, 44, 38, 86, 142],
                    [50, 29, 44, 44, 0, 24, 89, 142],
                    [48, 33, 44, 38, 24, 0, 90, 142],
                    [98, 84, 92, 86, 89, 90, 0, 148],
                    [148, 136, 152, 142, 142, 142, 148, 0]]
    example2 = [[0, 32, 48, 51, 50, 48, 49.0], [32, 0, 26, 34, 29, 33, 36.0], [48, 26, 0, 42, 44, 44, 48.0], [51, 34, 42, 0, 44, 38, 40.0], [50, 29, 44, 44, 0, 24, 41.5], [48, 33, 44, 38, 24, 0, 42.0], [49.0, 36.0, 48.0, 40.0, 41.5, 42.0, 0.0]]
    T = (neighborJoining(example1, 8, names))
    print("After: ")
    print(T)
    return 0

main()
