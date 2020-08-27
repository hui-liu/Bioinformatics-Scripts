import sys
from functions import merge_intervals

# read bed file, and convert into 1-based
region = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        chr, start, end = line.split()
        region.setdefault(chr, []).append([int(start)+1, int(end)])

# merge overlapping regions
merged_region = {}
for i in region:
    merged_region[i] = merge_intervals(region[i])

size = 0
for i in merged_region:
    for j in merged_region[i]:
        start, end = j
        if (int(end) - int(start)) < 0: print("data need to check since the start coordinate larger than the end coordinate")
        length = int(end) - int(start) + 1
        size += length

print size
