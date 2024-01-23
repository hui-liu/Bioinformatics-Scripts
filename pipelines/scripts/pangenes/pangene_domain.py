import sys
import gzip

domains = set()
with gzip.open(sys.argv[1], 'rt') as f:
    for line in f:
        lsp = line.split()
        domains.add(lsp[0])

with open(sys.argv[2], 'r') as f:
    for line in f:
        lsp = line.rstrip().split(",")
        if lsp[0] == "groupID": continue
        grID, geneid, cls = lsp[:3]
        if geneid in domains:
            print "\t".join([grID, geneid, cls, "with domain"])
        else:
            print "\t".join([grID, geneid, cls, "without domain"])
