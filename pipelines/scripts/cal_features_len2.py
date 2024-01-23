import sys
from collections import OrderedDict

def parseFasta(filename):
    fas = {}
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

regions = OrderedDict()
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line[0] == "#": continue
        lsp = line.rstrip().split("\t")
        chr,start,end = lsp[0], lsp[3], lsp[4]
        if int(start) > int(end):
            start, end = end, start
        regions.setdefault(chr, []).append([int(start), int(end)])

len_sum = []
for chr in regions:
    merged_reg = merge_intervals(regions[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    len_sum.append(length)


seq_dict = parseFasta(sys.argv[2])
size = sum([len(seq_dict[i]) for i in seq_dict])

print sys.argv[1].split(".")[0], str(int(round(sum(len_sum)/1000000.0))) +  "(" + str(round(sum(len_sum)/ float(size) * 100,2)) + "%)"
