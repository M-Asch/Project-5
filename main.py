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
    '''
    This function takes in a list of file names and will return a dictionary where each key
    is a sequence header and the value is the sequences

    Parameters:
        fileNames: a list of fileNames for .fasta files containing multiple sequences
    '''
    alignments = {}
    currentDirect = os.getcwd()
    for name in fileNames:
        file = open(currentDirect + name, "r")  #combine the current working directory with file names so .fasta files
        file = file.readlines()                 #could be stored in seperate files within the same directory
        header = file[0]
        key = ''
        for line in file[1:]:
            if line[0] == '>':                  #each new sequence starts with a > so we search for these to know when we
                alignments[header] = key        #should begin creating a new sequence and append the current one to our dictionary
                header = line
                key = ''
            else:
                key = key + line.strip('\n')    #clean up the file lines
        alignments[header] = key                #will miss one alignment otherwise
    return(alignments)

def findR(D):
    '''
    This function creates a list containing the sum of each row in a distance Matrix

    Parameters:
        D: a distance matrix
    '''
    rs = []
    for row in D:
        rs.append(sum(row))         #sum the current row within the D matrix
    return rs

def updateD(D, i, j):
    '''
    This function updates the distance matrix in each iteration of the neighbor joining algorithm
    by removing the nodes with the lowest M value and adding in a new row

    Parameters:
        D: a distance Matrix
        i: row index for the lowest M value
        j: column index for the lowest M value
    '''
    dij = D[i][j]               #store our total length of i,j

    newr = []
    if i > j:                   #check to see if i or j is bigger as we need to delete the
        del(D[i])               #larger one first to avoid creating an index out of range error
        del(D[j])
        for row in range(len(D)):
            value = ((D[row][i] + D[row][j] -dij)/2)    #calcualte the new value
            D[row].append(value)                        #add new value to row
            del(D[row][i])
            del(D[row][j])                              #delete the i and j rows
            newr.append(value)                          #build new row
    else:
        del(D[j])
        del(D[i])
        for row in range(len(D)):
            value = ((D[row][i] + D[row][j] -dij)/2)
            D[row].append(value)
            del(D[row][j])
            del(D[row][i])
            newr.append(value)
    newr.append(0)                  #add 0 to row for the (n, n) location
    D.append(newr)                  #append the new row to the new D matrix

    return D


def neighborJoining(D, n, names):
    '''
    This function creates a phylogenetic tree by using a distance matrix to calucalte
    the M values and the new distance matricies

    Parameters:
        D: a distance Matrix
        n: lengthof the distance matrix
        names: a list of names where the index coresponds to the same index in the distance matrix
    '''
    if n == 2:
        return("(" + str(names[0]) + ", " + str(names[1]) + ":" + str(D[0][1]) + ")")       #base case for when there are only 2 nodes left
    M = numpy.zeros((n, n))                 #create our M matrix
    rs = findR(D)
    for row in range(n):
        ra = rs[row]                        #get the R value for the current row
        for col in range(n):
            if row < col:
                rb = rs[col]                #get the R value for our current column
                M[row][col] = (n - 2) * D[row][col] - ra - rb       #calcualate the M for value (row, col)
    min_flat_index = numpy.argmin(M)
    i, j = numpy.unravel_index(min_flat_index, (n, n))              #get smallest value in the M matrix

    #establish lengths
    lengthB = D[i][j]/2 - ((rs[i] - rs[j])/(2*(n - 2)))             #find the lengths for the new nodes
    lengthA = D[i][j]/2 - ((rs[j] - rs[i])/(2*(n - 2)))

    #update names array
    nameA, nameB = names[i], names[j]
    names.append("(" + str(nameA) + ":" + str(lengthA) + ", "  + str(nameB) + ":" + str(lengthB) + ")") #add a tuple to names that has the nodes and their length
    names.remove(nameA)     #remove the updated names from the list
    names.remove(nameB)

    D = updateD(D, i, j)

    #ERROR IS OCCURING WITH ADDING
    return(neighborJoining(D, n - 1, names))        #recursively call the function

def buildD(filename):
    '''
    A funciton used to build a distance matrix from .pim files

    Parameters:
        filename: a list of .pim files to be read in
    '''
    currentDirect = os.getcwd()
    file = open(currentDirect + filename, "r")      #read in files
    rows = []
    names = []
    D = []
    for line in file:
        rows.append(line.split())        #split each row at its spaces in order to grab the name and values easier
    for row in rows:
        names.append(row[1][-25:-11])     #find the name for the current row and trim it so it is easier to read
        temp = []
        for item in row[2:]:
            temp.append(100 - float(item))  #pull and adjust the distance value for each column in the current row
        D.append(temp)
    return D, names


def barChart(fileNames):
    '''
    This funcion will read in alignments and calculate how many if them contain the B.1.1.7 mutations

    Parameters:
        fileNames: a list of files containing alignments
    '''
    alignments = []
    height = [0, 0, 0, 0 ,0]
    bars = []
    for bar in range(len(fileNames)):
        bars.append(fileNames[bar][10:19])  #get names for the bars
    for file in range(len(fileNames)):
        dict = dictBuilder(fileNames[file:file+1])
        count = 0
        for align in dict:
            if dict[align][21764:21770] == 'TTAATG':        #check the 3 mutations
                if dict[align][23064] == 'T':
                    if (dict[align][23601]) == 'A':
                        count += 1
                        print(align)
        height[file] = count

    y_pos = numpy.arange(len(bars))     #build the bar chart
    plt.bar(y_pos, height)
    plt.xticks(y_pos, bars)
    plt.show()


def main():
    fileNames = ['\England (DEC052020)\HCOV19-ENGLAND-051220-A.fasta', '\England (DEC082020)\HCOV19-ENGLAND-081220-A.fasta', '\England (NOV032020)\HCOV19-ENGLAND-031120-A.fasta', '\England (NOV102020)\HCOV19-ENGLAND-101120-A.fasta', '\England (NOV272020)\HCOV19-ENGLAND-271120-A.fasta']
    fileNames2 = ['\England (DEC052020)\HCOV19-ENGLAND-051220-D.pim', '\England (DEC082020)\HCOV19-ENGLAND-081220-D.pim', '\England (NOV032020)\HCOV19-ENGLAND-031120-D.pim', '\England (NOV102020)\HCOV19-ENGLAND-101120-D.pim', '\England (NOV272020)\HCOV19-ENGLAND-271120-D.pim']
    fileNames3 = ['\England (DEC052020)\HCOV19-ENGLAND-051220.fasta', '\England (DEC082020)\HCOV19-ENGLAND-081220.fasta', '\England (NOV032020)\HCOV19-ENGLAND-031120.fasta', '\England (NOV102020)\HCOV19-ENGLAND-101120.fasta', '\England (NOV272020)\HCOV19-ENGLAND-271120.fasta']
    barChart(fileNames3)
    #D = (distanceMatrix(alignments))

    #Example Distance Matrix
    '''names = ['dog', 'bear', 'raccoon', 'weasel', 'seal', 'sea_lion', 'cat', 'monkey']
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
