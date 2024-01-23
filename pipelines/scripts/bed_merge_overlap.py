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

repeat = {}
for line in sys.stdin:
    if line[0] == "#": continue
    lsp = line.rstrip().split("\t")
    chr,start,end = line.rstrip().split("\t")
    if int(start) > int(end):
        start, end = end, start
    repeat.setdefault(chr, []).append([int(start), int(end)])

for chr in repeat:
    merged_reg = merge_intervals(repeat[chr])
    for x in merged_reg:
        print chr + "\t" + "\t".join(map(str,x))
