import sys

res = []
n = int(sys.argv[2])
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = line.rstrip().split(",")
        pan_id = lsp[0]
        if pan_id == "Pan_gene_ID": continue
        genes = [x for x in lsp[2:] if x!= "NA"]
        sp_num = len(set([x[:3] for x in genes]))
        if sp_num == n:
            res.append([pan_id, "core"])
        elif 2 <= sp_num < n:
            res.append([pan_id, "shell"])
        else:
            res.append([pan_id, "specific"])

for i in res:
    print "\t".join(i)
