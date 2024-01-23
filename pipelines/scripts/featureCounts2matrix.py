import sys

OUT = open(sys.argv[2], 'w')
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line[0] == "#": continue
        lsp = line.rstrip().split("\t")
        if lsp[0] == "Geneid":
            Hline = [lsp[0]] + lsp[6:]
            Hline = ["_".join(i.split("_")[:2]) for i in Hline]
            OUT.write("\t".join(Hline) + "\n")
        else:
            geneid = lsp[0]
            exp = map(float, lsp[6:])
            if sum(exp) > 0:
                Eline = [lsp[0]] + map(str, exp)
                OUT.write("\t".join(Eline) + "\n")

OUT.close()
