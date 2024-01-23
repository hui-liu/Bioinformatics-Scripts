import sys
from collections import OrderedDict

def parseFasta(filename):
    fas = OrderedDict()
    id = None
    with open(filename, 'r') as fh:
        for line in fh:
            if line[0] == '>':
                header = line[1:].rstrip()
                id = header.split()[0]
                fas[id] = []
            else:
                fas[id].append(line.rstrip().upper())
        for id, seq in fas.iteritems():
            fas[id] = ''.join(seq)
    return fas

def merge_intervals(intervals):
    # https://stackoverflow.com/questions/43600878/merging-overlapping-intervals/43600953
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]
    for current in intervals:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
            #merged[-1][1] = max(previous[1], current[1])
        else:
            merged.append(current)
    return merged


seq_dict = parseFasta(sys.argv[1])
bed = sys.argv[2]
OUT = open(sys.argv[3], 'w')
OUT.write("\t".join(["chr", "start", "end", "percent"]) + "\n")

res = 100000
seq_len = OrderedDict()
for chr in seq_dict:
    seq = seq_dict[chr]
    seq_len[chr] = len(seq_dict[chr])

#
slide_win = OrderedDict()
for chr in seq_len:
    size = seq_len[chr]
    for n in range(0, size/res*res+res, res):
        slide_win[(chr, str(n+1), str(n+res))] = 0

bed_records = OrderedDict()
chroms = set()
with open(bed) as f:
    h = f.readline()
    for line in f:
        chr, start, end = line.split()
        chroms.add(chr)
        start, end = int(start)+1, int(end)
        pos = start/res * res
        bed_records.setdefault((chr, str(pos+1), str(pos+res)), []).append([start, end])

for index in bed_records:
    merged_reg = merge_intervals(bed_records[index])
    merged_reg_size = sum([x[1]-x[0]+1 for x in merged_reg])
    slide_win[index] += merged_reg_size

for index in slide_win:
    t = list(index)
    if t[0] in chroms:
        OUT.write("\t".join([t[0], str(int(t[1])-1), t[2]] + [str(slide_win[index]/(res*1.0))]) + "\n")

OUT.close()
