########
# A03
########

"""
GC content is different from species to species, different species will have
different characteristics GC contents.
We just using GC content as a way of trying to figure out whether the mix of
different bases is changing as we move along read. We expect that it won't
change very much.
"""

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
plt.plot(range(1, len(gc) + 1), gc)
plt.savefig(sys.argv[1].split('.gz')[0] + "_GC.png")
