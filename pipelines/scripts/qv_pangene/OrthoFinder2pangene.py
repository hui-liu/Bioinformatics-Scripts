import sys
import random
import math

sp_gene = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = line.split()
        sp_gene[lsp[1]] = lsp[0]

species = []
with open(sys.argv[2], 'r') as f:
    for line in f:
        species.append(line.rstrip())

print ",".join(["groupID", "pangene", "class"] + species)

sp_num = len(species)
soft_num = int(math.ceil(0.9 * sp_num))

with open(sys.argv[3], 'r') as f:
    for line in f:
        lsp = [x for x in line.split() if x]
        ogid = lsp[0].rstrip(":")
        og = lsp[1:]
        sp = set()
        d = {}
        for x in og:
            sp.add(sp_gene[x])
            d.setdefault(sp_gene[x], []).append(x)
        t = []
        genes = []
        clu = None
        for s in species:
            if s in d:
                g = d[s]
                t.append(random.choice(g))
                genes.append(";".join(g))
            else:
                genes.append("NA")
        if len(sp) == sp_num:
            clu = "core"
        elif soft_num <= len(sp) < sp_num:
            clu = "softcore"
        elif 2 <= len(sp) < sp_num:
            clu = "shell"
        else:
            clu = "private"
        print ",".join([ogid, t[0], clu] + genes)
