import sys
import itertools

species = ['QAL', 'QDE', 'QFA', 'QGR', 'QLI', 'QMO', 'QSE']

matrix = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = [x for x in line.split() if x]
        ogid = lsp[0].rstrip(":")
        og = lsp[1:]
        #if len(og) == 1: continue
        sp = set([x[:3] for x in og])
        #if len(sp) == 1: continue
        m = []
        for s in species:
            if s in sp:
                m.append(1)
            else:
                m.append(0)
        matrix.append(m)

matrix.sort(key = lambda x: (-sum(x), -x[6], -x[5], -x[4], -x[3], -x[2], -x[1], -x[0]))

out = open(sys.argv[2], 'w')
out.write("\t".join(species) + "\n")
for i in matrix:
    out.write("\t".join(map(str, i)) + "\n")

out.close()
