import sys
import gzip

sub_taxon = set()
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line.split():
            sub_taxon.add(line.split()[0])

sub_taxon_ids = set()
with gzip.open(sys.argv[2], 'rt') as f:
    for line in f:
        lsp = line.split()
        if lsp[1] in sub_taxon:
            sub_taxon_ids.add(lsp[0])

out = open(sys.argv[4], 'w')
keep = 0
with gzip.open(sys.argv[3], 'rt') as f:
    for line in f:
        if line[0] == '>':
            if line[1:].split()[0] in sub_taxon_ids:
                keep = 1
            else:
                keep = 0
        if keep == 1:
            out.write(line)

out.close()

