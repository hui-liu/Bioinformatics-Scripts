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
        for id, seq in fas.items():
            fas[id] = ''.join(seq)
    return fas

chroms = ["Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10","Chr11","Chr12"]

res = OrderedDict()
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line[0] == "#": continue
        temp = line.rstrip().split('\t')
        chr, source, feat_type, start, end, score, strand, phase, info = temp
        d_info = dict([v.split('=') for v in info.split(';')])
        #if chr not in chroms: continue
        if feat_type == "mRNA":
            res.setdefault(chr, []).append([chr.replace("Chr0", "").replace("Chr", ""), d_info["ID"],
                                            start, end, strand, d_info["Parent"]])
seqs = parseFasta(sys.argv[2])
gffout = open(sys.argv[3], 'w')
lensout = open(sys.argv[4], 'w')
for i in chroms:
    lst = sorted(res[i], key = lambda x: int(x[2]))
    lensout.write("\t".join([i.replace("Chr0", "").replace("Chr", ""), str(len(seqs[i])), str(len(lst))]) + "\n")
    n = 0
    for i in lst:
        n += 1
        gffout.write("\t".join(i[:-1] + [str(n), i[-1]]) + "\n")

lensout.close()
gffout.close()
