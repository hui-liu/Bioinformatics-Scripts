import sys
from collections import OrderedDict

chroms = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
          "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19"]

res = OrderedDict()
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line[0] == "#": continue
        temp = line.rstrip().split('\t')
        chr, source, feat_type, start, end, score, strand, phase, info = temp
        d_info = dict([v.split('=') for v in info.split(';')])
        #if chr not in chroms: continue
        if feat_type == "mRNA":
            res.setdefault(chr, []).append([chr.replace("chr", ""), d_info["ID"].replace(".Genoscope12X", ""),
                                            start, end, strand, d_info["Parent"].replace(".Genoscope12X", "")])
OUT = open(sys.argv[2], 'w')
for i in chroms:
    lst = sorted(res[i], key = lambda x: int(x[2]))
    n = 0
    for i in lst:
        n += 1
        OUT.write("\t".join(i[:-1] + [str(n), i[-1]]) + "\n")

OUT.close()
