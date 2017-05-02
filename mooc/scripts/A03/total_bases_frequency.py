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

seqs, quals = readFastq(sys.argv[1])
count = collections.Counter()
for seq in seqs:
    count.update(seq)
print count
