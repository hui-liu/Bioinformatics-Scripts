import sys
from collections import OrderedDict
import gzip

ev_cutoff = 1e-5
id_cutoff = 0.0
bs_cutoff = 0.0

infile = sys.argv[1]
topx = sys.argv[2]

df = OrderedDict()
with gzip.open(infile, 'rt') as f:
    for line in f:
        lsp = line.rstrip().split('\t')
        identity = float(lsp[2])
        e_value = float(lsp[10])
        bit_score = float(lsp[11])
        if identity < id_cutoff or e_value > ev_cutoff or bit_score < bs_cutoff: continue
        df.setdefault(lsp[0], []).append(lsp)


out = open(sys.argv[3], 'w')
for i in df:
    t = df[i]
    top_x = sorted(t, key = lambda x: -float(x[11]))
    if len(top_x) > 5:
        top_x = top_x[:int(topx)]
    for j in top_x:
        out.write("\t".join(j) + "\n")

out.close()
