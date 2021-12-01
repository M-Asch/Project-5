#================================
# Mitchell Aschmeyer and Riley Collins
#
#
#================================
import os
import numpy
from compbio import *
from Bio import Phylo
from io import StringIO
import matplotlib.pyplot as plt

def dictBuilder(fileNames):
    alignments = {}
    currentDirect = os.getcwd()
    for name in fileNames:
        file = open(currentDirect + name, "r")
        file = file.readlines()
        header = file[0]
        key = ''
        for line in file[1:]:
            if line[0] == '>':
                alignments[header] = key
                header = line
                key = ''
            else:
                key = key + line.strip('\n')
        alignments[header] = key
    return(alignments)

def findR(D):
    rs = []
    for row in D:
        rs.append(sum(row))
    return rs

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
        return("(" + str(names[0]) + ", " + str(names[1]) + ":" + str(D[0][1]) + ")")
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
    names.append("(" + str(nameA) + ":" + str(lengthA) + ", "  + str(nameB) + ":" + str(lengthB) + ")")
    names.remove(nameA)
    names.remove(nameB)


    D = updateD(D, i, j)


    #ERROR IS OCCURING WITH ADDING
    return(neighborJoining(D, n - 1, names))

def buildD(filename):
    currentDirect = os.getcwd()
    file = open(currentDirect + filename, "r")
    rows = []
    names = []
    D = []
    for line in file:
        rows.append(line.split())
    for row in rows:
        names.append(row[1][-25:-11])
        temp = []
        for item in row[2:]:
            temp.append(100 - float(item))
        D.append(temp)
    return D, names


def barChart(fileNames):
    alignments = []
    height = [0, 0, 0, 0 ,0]
    alignments.append(dictBuilder(fileNames[1:2]))
    alignments.append(dictBuilder(fileNames[2:3]))
    for name in range(len(fileNames)):
        fileNames[name] = fileNames[name][10:19]
        bars = tuple(fileNames)

    for align in range(len(alignments)):
        for a in alignments[align]:
            if alignments[align][a][21764:21770] == '------':
                if alignments[align][a][23062] == 'T':
                    if alignments[align][a][23603] == 'A':
                        height[align + 1] += 1

    y_pos = numpy.arange(len(bars))
    plt.bar(y_pos, height)
    plt.xticks(y_pos, bars)
    plt.show()


def main():
    fileNames = ['\England (DEC052020)\HCOV19-ENGLAND-051220-A.fasta', '\England (DEC082020)\HCOV19-ENGLAND-081220-A.fasta', '\England (NOV032020)\HCOV19-ENGLAND-031120-A.fasta', '\England (NOV102020)\HCOV19-ENGLAND-101120-A.fasta', '\England (NOV272020)\HCOV19-ENGLAND-271120-A.fasta']
    fileNames2 = ['\England (DEC052020)\HCOV19-ENGLAND-051220-D.pim', '\England (DEC082020)\HCOV19-ENGLAND-081220-D.pim', '\England (NOV032020)\HCOV19-ENGLAND-031120-D.pim', '\England (NOV102020)\HCOV19-ENGLAND-101120-D.pim', '\England (NOV272020)\HCOV19-ENGLAND-271120-D.pim']
    alignments = dictBuilder(fileNames)
    barChart(fileNames)
    #D = (distanceMatrix(alignments))

    '''
    #Example Distance Matrix
    names = ['dog', 'bear', 'raccoon', 'weasel', 'seal', 'sea_lion', 'cat', 'monkey']
    example1 = [[0, 32, 48, 51, 50, 48, 98, 148],
                    [32, 0, 26, 34, 29, 33, 84, 136],
                    [48, 26, 0, 42, 44, 44, 92, 152],
                    [51, 34, 42, 0, 44, 38, 86, 142],
                    [50, 29, 44, 44, 0, 24, 89, 142],
                    [48, 33, 44, 38, 24, 0, 90, 142],
                    [98, 84, 92, 86, 89, 90, 0, 148],
                    [148, 136, 152, 142, 142, 142, 148, 0]]
    T = (neighborJoining(example1, 8, names))
    T = str(T)
    print(T)
    T = Phylo.read(StringIO(T), 'newick')
    Phylo.draw(T)
    '''

    '''
    #.pim files
    for file in fileNames2:
        D, names = buildD(file)
        T = (neighborJoining(D, len(D), names))
        T = Phylo.read(StringIO(T), 'newick')
        Phylo.draw(T)
    '''


    return 0

main()
