########
# A04
########
import sys
import random

# (1) the size of genome
def readGenome(filename):
    genome=''
    with open (filename,'r') as f:
        for line in f:
	          if not line[0]=='>':
		            genome += line.rstrip()
    return genome
# run
genome = readGenome(sys.argv[1])

# (2) naive exact matching algorithm
def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        match = True
	      for j in range(len(p)):
	          if not t[i+j] == p[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    return occurrences

# (3) generate Reads from genome
def generateReads(genome, numReads, readLen):
    ''' Gernerate reads from random position in the given genome. '''
    reads = []
    for _ in range(numReads):
        start = random.randint(0, len(genome) - readLen - 1)
        reads.append(genome[start : start+readLen])
    return reads
# run
reads = generateReads(genome, 100, 100)
for r in reads:
    matches = naive(r, genome)
    if len(matches) > 0:
        numMatched += 1
print "%d / %d reads matched exactly" % (numMatched, len(reads))
