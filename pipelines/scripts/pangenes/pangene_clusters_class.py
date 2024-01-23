import sys

res = []
n = int(sys.argv[2])
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = line.split()
        og = lsp[0].rstrip(":")
        genes = lsp[1:]
        sp_num = len(set([x[:3]for x in genes]))
        if sp_num == n:
            res.append([og, "core"])
        elif 2 <= sp_num < n:
            res.append([og, "shell"])
        else:
            res.append([og, "specific"])

for i in res:
    print "\t".join(i)
