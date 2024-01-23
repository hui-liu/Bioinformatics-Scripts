import sys

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

regions = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = line.rstrip().split("\t")
        #if lsp[2] == "match":
        chr,start,end = lsp[0], lsp[3], lsp[4]
        if int(start) > int(end):
            start, end = end, start
        regions.setdefault(chr, []).append([int(start)-1, int(end)])

for chr in regions:
    merged_reg = merge_intervals(regions[chr])
    for i in merged_reg:
        print "\t".join([chr] + map(str, i))
