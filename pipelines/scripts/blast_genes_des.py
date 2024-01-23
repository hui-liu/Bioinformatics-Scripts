import sys
import gzip

des = {}
with gzip.open(sys.argv[1], 'rt') as f:
    for line in f:
        lsp = line.rstrip().split("\t")
        des[lsp[0]] = lsp[1]

OUT = open(sys.argv[3], 'w')
with gzip.open(sys.argv[2], 'rt') as f:
    for line in f:
        lsp = line.split()
        OUT.write("\t".join(lsp + [des[lsp[1]]]) + "\n")

OUT.close()
