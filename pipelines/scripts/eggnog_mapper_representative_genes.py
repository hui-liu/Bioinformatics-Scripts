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
        lsp = line.split()
        if lsp[0] == "##": continue
        if lsp[0] == "#query":
            OUT.write(line)
        if lsp[0] in gene_ids:
            OUT.write(line)

OUT.close()
