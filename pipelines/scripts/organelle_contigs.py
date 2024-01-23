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


#
eva_cutoff = 1e-10
ide_cutoff = 99.0
cov_cutoff = 60.0

#
seq_dict = parseFasta(sys.argv[1])

regions_mt = OrderedDict()
regions_pt = OrderedDict()
with open(sys.argv[2], 'r') as f:
    for line in f:
        if line[0] == "#": continue
        lsp = line.rstrip().split("\t")
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = lsp
        pident = float(pident)
        evalue = float(evalue)
        qstart, qend = int(qstart), int(qend)
        if qstart > qend:
            qstart, qend = qend, qstart
        if evalue < eva_cutoff and pident > ide_cutoff:
                if sseqid == "Mt":
                    regions_mt.setdefault(qseqid, []).append([qstart, qend])
                elif sseqid == "Pt":
                    regions_pt.setdefault(qseqid, []).append([qstart, qend])

mt = OrderedDict()
n = 0
for chr in regions_mt:
    merged_reg = merge_intervals(regions_mt[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    cov = length * 1.0 / len(seq_dict[chr]) * 100
    if cov > cov_cutoff:
        n += 1
        mt[chr] = "Mt" # + str(n)

pt =  OrderedDict()
n = 0
for chr in regions_pt:
    merged_reg = merge_intervals(regions_pt[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    cov = length * 1.0 / len(seq_dict[chr]) * 100
    if cov > cov_cutoff:
        n += 1
        mt[chr] = "Pt" # + str(n)

for i in mt:
    print i, mt[i]

for i in pt:
    print i, pt[i]



"""
n = 0
for chr in seq_dict:
    if chr in mt or chr in pt:
        continue
    seq = seq_dict[chr]
    if "Chr" in chr:
        print(">" + chr + "\n" + seq)
    else:
        n += 1
        print(">Contig" + str(n) + "\n" + seq)

for chr in mt:
    seq = seq_dict[chr]
    print(">" + mt[chr] + "\n" + seq)

for chr in pt:
    seq = seq_dict[chr]
    print(">" + pt[chr] + "\n" + seq)
"""
