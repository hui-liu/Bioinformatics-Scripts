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

chroms = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
          "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19"]

pep = parseFasta(sys.argv[1])
res = OrderedDict()
with open(sys.argv[2], 'r') as f:
    for line in f:
        if line[0] == "#": continue
        temp = line.rstrip().split('\t')
        chr, source, feat_type, start, end, score, strand, phase, info = temp
        d_info = dict([v.split('=') for v in info.split(';')])
        #if chr not in chroms: continue
        if feat_type == "mRNA":
            id = d_info["ID"].replace(".v2.1", "")
            if id in pep:
                res.setdefault(chr, []).append([chr.replace("chr", ""), id, start, end, 
                               strand, d_info["Parent"].replace(".v2.1", "")])

seqs = parseFasta(sys.argv[3])
gffout = open(sys.argv[4], 'w')
lensout = open(sys.argv[5], 'w')
chr_pep = set()
for i in chroms:
    lst = sorted(res[i], key = lambda x: int(x[2]))
    lensout.write("\t".join([i.replace("chr", ""), str(len(seqs[i])), str(len(lst))]) + "\n")
    n = 0
    for s in lst:
        n += 1
        chr_pep.add(s[1])
        gffout.write("\t".join(s[:-1] + [str(n), s[-1]]) + "\n")

lensout.close()
gffout.close()

pepout = open(sys.argv[6], 'w')
for i in chr_pep:
    pepout.write(">" + i + "\n" + pep[i] + "\n")
pepout.close()
