import sys
import itertools

#core = {}
species = ['QAL', 'QDE', 'QFA', 'QGR', 'QLI', 'QMO', 'QSE']
matrix = {x: set() for x in species}
n = len(species)
genes = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = [x for x in line.split() if x]
        ogid = lsp[0].rstrip(":")
        og = lsp[1:]
        genes[ogid] = og
        #if len(og) == 1: continue
        sp = list(set([x[:3] for x in og]))
        for t in matrix:
            if t in sp:
                matrix[t].add(ogid)
        #[matrix[t].add(ogid) for t in matrix if t in sp]
        #m = 0
        #while m <= len(sp):
        #    m += 1
        #    comb = [x for x in itertools.combinations(sp, m)]
        #    for s in comb:
        #        s = tuple(sorted(list(s)))
        #        core.setdefault(s, []).append(ogid)

pan = {}
core = {}
m = 0
while m <= len(species):
    m += 1
    comb = [x for x in itertools.combinations(species, m)]
    for s in comb:
        s = tuple(sorted(list(s)))
        if len(s) == 1:
            pan[s] = matrix[s[0]]
            core[s] = matrix[s[0]]
        else:
            for t in s:
                pan.setdefault(s, set()).update(matrix[t])
                if s in core:
                    core[s] = core[s].intersection(matrix[t])
                else:
                    core[s] = matrix[t]
core_genes = {}
for i in core:
    for j in core[i]:
        t = [g for g in genes[j] if g[:3] in i]
        print i, t
