import re
import sys
from collections import OrderedDict

"""
AAACCCT
 AACCCTA
  ACCCTAA
   CCCTAAA
    CCTAAAC
     CTAAACC
      TAAACCC
"""

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
OUT = open(sys.argv[2], 'w')
OUT.write("\t".join(["chr", "start", "end", "win_start", "win_end", "motif"]) + "\n")
OUT2 = open(sys.argv[3], 'w')
OUT2.write("\t".join(["chr", "start", "end", "size", "percent"]) + "\n")
OUT3 = open(sys.argv[4], 'w')
OUT3.write("\t".join(["chr", "chr_len", "start", "end", "telo_len"]) + "\n")

res = 100000
times = 5
telomere_win = OrderedDict()
seq_len = OrderedDict()
telomere_loc = OrderedDict()
connect = 1000

for chr in seq_dict:
    seq = seq_dict[chr]
    seq_len[chr] = len(seq_dict[chr])

    telomere_F = re.finditer("CCCTAAA"*times, seq)
    indices_F = [m.start(0) for m in telomere_F]
    telomere_R = re.finditer("TTTAGGG"*times, seq)
    indices_R = [m.start(0) for m in telomere_R]

    for x in indices_F:
        pos = (x+1)/res * res
        # to connect the adjacent regions, add 1 or 200 to the ends
        telomere_win.setdefault((chr, str(pos+1), str(pos+res)), []).append([x+1, x+7*times+1])
        telomere_loc.setdefault(chr, []).append([x+1, x+7*times+connect])
        OUT.write("\t".join([chr, str(x+1), str(x+7*times), str(pos+1), str(pos+res), "CCCTAAA"]) + "\n")

    for x in indices_R:
        pos = (x+1)/res * res
        telomere_win.setdefault((chr, str(pos+1), str(pos+res)), []).append([x+1, x+7*times+1])
        telomere_loc.setdefault(chr, []).append([x+1, x+7*times+connect])
        OUT.write("\t".join([chr, str(x+1), str(x+7*times), str(pos+1), str(pos+res), "TTTAGGG"]) + "\n")

OUT.close()

#
slide_win = OrderedDict()
for chr in seq_len:
    size = seq_len[chr]
    for n in range(0, size/res*res+res, res):
        slide_win[(chr, str(n+1), str(n+res))] = 0

for index in telomere_win:
    merged_reg = merge_intervals(telomere_win[index])
    merged_reg_size = sum([x[1]-x[0]+1 for x in merged_reg])
    slide_win[index] += merged_reg_size

for index in slide_win:
    OUT2.write("\t".join(list(index) + [str(slide_win[index]), str(slide_win[index]/(res*1.0))]) + "\n")

OUT2.close()


for chr in telomere_loc:
    m = merge_intervals(telomere_loc[chr])
    for n in m:
        OUT3.write("\t".join(map(str, [chr, seq_len[chr], n[0], n[1]-connect, n[1]-connect-n[0]+1])) + "\n")

OUT3.close()
