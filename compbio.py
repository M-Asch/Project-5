"""
Riley Collins
October 2021
compbio.py
Contains recurring functions from previous projects to use in future porjects.
"""

import math
import random
from Bio import Entrez, SeqIO

"""
writeDNAFile
Purpose:
    will pull a given fasta file from NCBI and write it to a text file
Parameters:
    db: the database to pull the file from
    ids: the id of the genome to pull
    rettype: the type of file to pull
    retmode: determines the form of the returned output
Return Values:
    dna: a string containing the genome of the wanted organism
    record: information about the pulled genome
    file: A text file containing the new genome
"""
def writeDNAFile(db, ids, rettype, retmode):
    Entrez.email = "collin_r2@denison.edu"
    handle = Entrez.efetch(db = db, id =ids, rettype = rettype, retmode = retmode)
    record = SeqIO.read(handle, 'fasta')
    dna = str(record.seq)
    file = open(ids, "w")
    file.write(dna)
    file.close()
    return dna, record

"""
readDNA
Purpose:
    will pull a given fasta file from NCBI and return a string containing it
Parameters:
    db: the database to pull the file from
    ids: the id of the genome to pull
    rettype: the type of file to pull
    retmode: determines the form of the returned output
Return Values:
    dna: a string containing the genome of the wanted organism
"""
def readFromOnline(db, ids, rettype, retmode):
    Entrez.email = "collin_r2@denison.edu"
    handle = Entrez.efetch(db = db, id =ids, rettype = rettype, retmode = retmode)
    record = SeqIO.read(handle, 'fasta')
    dna = str(record.seq)
    return dna

"""
readDNAFile
Purpose:
    will contain a string of the given file
Parameters:
    filename: name of the file to be read
Return Values:
    dna: a string containing the genome of the wanted organism
"""
def readFromFile(filename):
  file = open(filename, "r")
  seq_object = SeqIO.read(filename, 'fasta')
  dna = str(seq_object.seq)
  return dna

def readFastaFileSequence(filename):
  sequenceDict = SeqIO.index(filename, 'fasta')
  sequences = []
  for id in sequenceDict:
    sequences.append(str(sequenceDict[id].seq))
  return sequences


#######################################################################################

def getCounts(motifs):
  """Return nt counts given a list of motif strings."""
  counts = []
  for i in range(len(motifs[0])):
    counts.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})

  for motif in motifs:
    motif = motif.upper()
    for i in range(len(motif)):
        counts[i][motif[i]] += 1

  return counts

def getProfile(motifs):
  """Get a profile from a set of motifs."""
  counts = getCounts(list(motifs))
  total = len(motifs)
  for i in range(len(counts)):
    for nt in 'ACGT':
      counts[i][nt] /= total
  return counts  # really now a profile

def getProfileLaplace(motifs):
  """Get a profile from a set of motifs, adding pseudocounts
    (Laplace's rule of succession) to prevent zero-probability events.
  """
  counts = getCounts(list(motifs))
  total = len(motifs) + 4
  for i in range(len(counts)):
    for nt in 'ACGT':
      counts[i][nt] += 1
      counts[i][nt] /= total
  return counts  # really now a profile

def getMostProbable(sequence, profile, k):
  """Return the most probable k-mer in sequence, given a profile."""
  maxProb = -1
  best_kmer = ''
  for i in range(len(sequence) - k + 1):
    kmer = sequence[i:i + k]
    prob = 1
    for j in range(k):
      prob *= profile[j][kmer[j]]
    if prob > maxProb:
      maxProb = prob
      best_kmer = kmer
  return best_kmer

def getMostProbableMotifs(sequences, profile, k):
  """Return the most probable set of k-mers from a set of sequences,
    given a profile."""

  motifs = []
  for sequence in sequences:
    motifs.append(getMostProbable(sequence, profile, k))
  return motifs

def getConsensus(motifs):
  """Return the consensus sequence for a set of motifs."""

  counts = getCounts(list(motifs))

  consensus = ''
  for i in range(len(counts)):
    majority = 'A'
    for nt in counts[i]:
      if counts[i][nt] > counts[i][majority]:
        majority = nt
      consensus += majority

  return consensus

def getScore(motifs):
  """Return the score for a set of motifs."""
  counts = getCounts(list(motifs))
  consensus = getConsensus(motifs)

  t = len(motifs)
  score = 0
  for i in range(len(consensus)):
    nt = consensus[i]
    diff = t - counts[i][nt]
    score += diff

  return score

def motifSearchBF(sequences, k, indices):
  if len(indices) == len(sequences):
    motifs = []
    for i in range(len(indices)):
      dna = sequences[i]
      index = indices[i]
      motifs.append(dna[index:index+k])
    return motifs

  n = len(sequences[0])
  bestScore = math.inf
  bestMotifs = []
  for i in range(n-k+1):
    motifs = motifSearchBF(sequences, k, indices + [i])
    score = getScore(motifs)
    if score < bestScore:
      bestScore = score
      bestMotifs = motifs
  return bestMotifs

def motifSearchGreedy(sequences, k, iterations):
  motifs = []
  for i in range(len(sequences)):
    motifs.append(sequences[i][:k])

  bestMotifs = motifs
  for i in range(iterations):
    profile = getProfileLaplace(motifs)
    motifs = getMostProbableMotifs(sequences, profile, k)
    if getScore(motifs) < getScore(bestMotifs):
      bestMotifs = motifs
  return bestMotifs