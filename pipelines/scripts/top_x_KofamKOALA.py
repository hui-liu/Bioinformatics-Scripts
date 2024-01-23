import sys
from collections import OrderedDict
import gzip

ev_cutoff = 0.01

infile = sys.argv[1]
topx = sys.argv[2]

df = OrderedDict()
with gzip.open(infile, 'rt') as f:
    for line in f:
        if line[0] == "#": continue
        lsp = line.rstrip().split('\t')
        e_value = float(lsp[5])
        if e_value > ev_cutoff: continue
        df.setdefault(lsp[1], []).append(lsp)


out = open(sys.argv[3], 'w')
for i in df:
    t = df[i]
    top_x = sorted(t, key = lambda x: -float(x[4]))
    if len(top_x) > 5:
        top_x = top_x[:int(topx)]
    for j in top_x:
        out.write("\t".join(j) + "\n")

out.close()
