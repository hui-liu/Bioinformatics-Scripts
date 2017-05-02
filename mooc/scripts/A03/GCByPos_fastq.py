########
# A03
########
import collections
import sys
import matplotlib.pyplot as plt
import gzip

def readFastq(filename):
    sequences = []
    qualities = []
    with gzip.open(filename, 'r') as fh:
        while True:
            fh.readline()
            seq = fh.readline().rstrip()
            fh.readline()
            qual = fh.readline().rstrip()
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

# (1)
def findGCByPos(reads):
    gc = [0] * 101 # the length of reads
    totals = [0] * 101
	
    for read in reads:
	      for i in range(len(read)):
		        if read[i] == 'C' or read[i] == 'G':
			          gc[i] += 1
		        totals[i] += 1
		
    for i in range(len(gc)):
		    if totals[i] > 0:
		        gc[i] /= float(totals[i])
    return gc
# run
seqs, quals = readFastq(sys.argv[1])
gc = findGCByPos(seqs)
plt.plot(range(len(gc), gc))
plt.show()

# (2)
count = collections.Counter()
for seq in seqs:
    count.update(seq)
print count
