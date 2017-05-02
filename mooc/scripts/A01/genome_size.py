########
# A01
########
# (1) the size of genome
def readGenome(filename):
    genome=''
    with open (filename,'r') as f:
        for line in f:
            if not line[0]=='>':
                genome += line.rstrip()
    return genome
# run
genome = readGenome('genome.fasta')
#print genome[:100]
print len(genome)

# (2) frequency of each base
# method 1
counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
for base in genome:
    counts[base] += 1
print counts

# method 2
import collections
counts = collections.Counter(genome)
print counts
