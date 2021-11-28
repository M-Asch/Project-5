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

def neighborJoining(D, n , names):
    if n == 2:
        return(names[0], names[1])
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
    lengthA = D[i][j]/2 - ((rs[i] - rs[j])/(2*(n - 2)))
    lengthB = D[i][j]/2 - ((rs[j] - rs[i])/(2*(n - 2)))

    #update names array
    nameA, nameB = names[i], names[j]
    names.append((nameA, nameB))
    names.remove(nameA)
    names.remove(nameB)

    newR = []

    #find correct way to build the new row and update values in other rows
    '''for col in range(n):
        newR.append((D[i][col] + D[j][col] - D[i][j])/2)'''
    dij = D[i][j]
    newR = []
    for col in range(n):
        newR.append((D[i][col] + D[j][col] - D[i][j])/2)
    print(newR)

    if i > j:
        del(D[i])
        del(newR[i])
        for row in range(n - 1):
            D[row][n - 1] = newR[row]
        D[j] = newR
    if j > i:
        del(D[j])
        del(newR[j])
        for row in range(n - 1):
            D[row][n - 1] = newR[row]
        D[i] = newR


    #ERROR IS OCCURING WITH ADDING
    T = (neighborJoining(D, n - 1, names))
    #return T


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
    D = (neighborJoining(example1, 8, names))
    print("After: ")
    for line in D:
        print(line)
    return 0

main()
