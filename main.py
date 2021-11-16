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
        header = file[0]
        key = ''
        for line in file[1:]:
            key = key + line.strip('\n')
        alignments[header] = key
    return(alignments)


def main():
    alignments = dictBuilder()
    print((alignments))
    return 0

main()
