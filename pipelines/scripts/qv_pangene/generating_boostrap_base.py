import sys

sp_gene = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = line.split()
        sp_gene[lsp[1]] = lsp[0]

species = []
with open(sys.argv[2], 'r') as f:
    for line in f:
        species.append(line.rstrip())

df = []
with open(sys.argv[3], 'r') as f:
    for line in f:
        lsp = [x for x in line.split() if x]
        ogid = lsp[0].rstrip(":")
        og = lsp[1:]
        sp = set([sp_gene[x] for x in og])
        t = []
        for s in species:
            if s in sp:
                t.append(ogid)
            else:
                t.append("NA")
        df.append(t)

for x, y in enumerate(species):
    with open(y + "_pan_id.csv", 'w') as f:
        for i in df:
            if i[x] != "NA":
                f.write(i[x] + "\n")
