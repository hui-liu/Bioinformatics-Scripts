import sys
import gzip


res = {}
n = int(sys.argv[2])
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = line.rstrip().split(",")
        pan_id = lsp[0]
        if pan_id == "Pan_gene_ID": continue
        genes = [x for x in lsp[2:] if x!= "NA"]
        sp_num = len(set([x[:3] for x in genes]))
        if sp_num == n:
            for g in genes:
                t = [k for k in g.split(";")]
                for v in t:
                    if "gmap_ID" in v: continue
                    res.setdefault("core", []).append(v)
        elif 2 <= sp_num < n:
            for g in genes:
                t = [k for k in g.split(";")]
                for v in t:
                    if "gmap_ID" in v: continue
                    res.setdefault("shell", []).append(v)
        else:
            for g in genes:
                t = [k for k in g.split(";")]
                for v in t:
                    if "gmap_ID" in v: continue
                    res.setdefault("private", []).append(v)
domains = set()
with gzip.open(sys.argv[3], 'r') as f:
    for line in f:
        lsp = line.split()
        domains.add(lsp[0])

for i in res:
    for j in res[i]:
        if j in domains:
            print "\t".join([j[:3], j, i, "with domain"])
        else:
            print "\t".join([j[:3], j, i, "without domain"])
