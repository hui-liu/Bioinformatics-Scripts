import sys


OUT = open(sys.argv[2], 'w')

n=0
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line[0] == "#":
            continue
        lsp = line.rstrip().split("\t")
        strand = lsp[6]
        if strand in ["?", "."]:
            strand = "+"
        OUT.write("\t".join(lsp[:6] + [strand] + lsp[7:]) + "\n")

OUT.close()
