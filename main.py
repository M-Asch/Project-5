#================================
# Mitchell Aschmeyer and Riley Collins
#
#
#================================
import os


def dictBuilder():
    alignments = {}
    currentDirect = os.getcwd()
    fileNames = ['\England (DEC052020)\HCOV19-ENGLAND-051220-A.fasta', '\England (DEC082020)\HCOV19-ENGLAND-081220-A.fasta', '\England (NOV032020)\HCOV19-ENGLAND-031120-A.fasta', '\England (NOV102020)\HCOV19-ENGLAND-101120-A.fasta', '\England (NOV272020)\HCOV19-ENGLAND-271120-A.fasta']
    for name in fileNames:
        file = open(currentDirect + name, "r")
        file = file.readlines()
        alignments[file[0]] = file[1:]
    return(alignments)


def main():
    dictBuilder()
    return 0

main()
