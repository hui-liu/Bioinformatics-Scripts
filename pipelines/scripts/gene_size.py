import sys
from collections import OrderedDict

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
        if lsp[2] != "gene": continue
        chr,start,end = lsp[0], lsp[3], lsp[4]
        if int(start) > int(end):
            start, end = end, start
        regions.setdefault(chr, []).append([int(start), int(end)])

res = []
for chr in regions:
    merged_reg = merge_intervals(regions[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    res.append(length)


print round(sum(res)/1000000.0, 2), "Mb"
