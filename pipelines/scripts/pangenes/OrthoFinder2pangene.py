import sys
import random

species = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        species.append(line.rstrip())

print ",".join(["groupID", "pangene", "class"] + species)

sp_num = len(species)
with open(sys.argv[2], 'r') as f:
    for line in f:
        lsp = [x for x in line.split() if x]
        ogid = lsp[0].rstrip(":")
        og = lsp[1:]
        sp = set()
        d = {}
        for x in og:
            sp.add(x[:3])
            d.setdefault(x[:3], []).append(x)
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
        elif 2 <= len(sp) <= sp_num:
            clu = "shell"
        else:
            clu = "private"
        print ",".join([ogid, t[0], clu] + genes)
