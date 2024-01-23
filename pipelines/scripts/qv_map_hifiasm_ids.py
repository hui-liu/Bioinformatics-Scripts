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


n = 0
ids = {}
id_map = OrderedDict()
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = line.split()
        if lsp[0][0] == ">":
            ids[lsp[1]] = lsp[0].lstrip(">")
        else:
            n += 1
            t = [x.lstrip("-") for x in lsp]
            chr = "HiC_scaffold_" + str(n)
            for y in t:
                id_map.setdefault(chr, []).append(ids[y])

seq_dict = parseFasta(sys.argv[2])
n = 0
new_seq_dict = OrderedDict()
for chr in id_map:
    seq = seq_dict[chr]
    ids = id_map[chr]
    n += 1
    if n <= 24:
        new_seq_dict[chr] = seq
    else:
        if "N" in seq:
            seq_list = [x for x in seq.split("N") if x]
            for x,y in enumerate(ids):
                new_seq_dict[y] = seq_list[x]
        else:
            new_seq_dict[ids[0]] = seq

for chr in new_seq_dict:
    seq = new_seq_dict[chr]
    print(">" + chr + "\n" + seq)

