import sys
import gzip

res = {}
n = int(sys.argv[2])
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = line.split()
        og = lsp[0].rstrip(":")
        genes = lsp[1:]
        sp_num = len(set([x[:3]for x in genes]))
        if sp_num == n:
            for g in genes:
                res.setdefault("core", []).append(g)
        elif 2 <= sp_num < n:
            for g in genes:
                res.setdefault("shell", []).append(g)
        else:
            for g in genes:
                res.setdefault("private", []).append(g)

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

