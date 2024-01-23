import sys
import gzip

def get_represent(filename):
    ids = set()
    with open(filename, 'r') as fh:
        for line in fh:
            if line[0] == '>':
                header = line[1:].rstrip()
                ids.add(header.split()[0])
    return ids


gene_ids = get_represent(sys.argv[1])

OUT = open(sys.argv[3], 'w')

with gzip.open(sys.argv[2], 'rt') as f:
    for line in f:
        lsp = line.rstrip().split("\t")
        if lsp[0] == "#": continue
        if lsp[1] in gene_ids and float(lsp[5]) < 0.001:
            OUT.write("\t".join([lsp[1], lsp[2], lsp[5], lsp[6]]) + "\n")

OUT.close()
