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
cov_cutoff = 80.0

#
seq_dict = parseFasta(sys.argv[1])

regions = OrderedDict()
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
            regions.setdefault(qseqid, []).append([qstart, qend])

redun = set()
for chr in regions:
    merged_reg = merge_intervals(regions[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    cov = length * 1.0 / len(seq_dict[chr]) * 100
    if cov > cov_cutoff:
        redun.add(chr)
for i in redun:
    print i
regions2 = OrderedDict()
with open(sys.argv[3], 'r') as f:
    for line in f:
        if line[0] == "#": continue
        lsp = line.rstrip().split("\t")
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = lsp
        if qseqid == sseqid: continue
        if qseqid in redun or sseqid in redun: continue
        pident = float(pident)
        evalue = float(evalue)
        qstart, qend = int(qstart), int(qend)
        if qstart > qend:
            qstart, qend = qend, qstart
        if evalue < eva_cutoff and pident > ide_cutoff:
            if len(seq_dict[qseqid]) <= len(seq_dict[sseqid]):
                regions2.setdefault(qseqid, []).append([qstart, qend])

print "---"
redun2 = set()
for chr in regions2:
    merged_reg = merge_intervals(regions2[chr])
    length = sum([i[1] - i[0] + 1 for i in merged_reg])
    cov = length * 1.0 / len(seq_dict[chr]) * 100
    if cov > cov_cutoff:
        print  chr, cov
        redun2.add(chr)


"""
for chr in seq_dict:
    seq = seq_dict[chr]
    if chr not in redun:
        seq = seq_dict[chr]
        print(">" + chr + "\n" + seq)
"""
