import sys

def getSeqID(filename):
    id = set()
    with open(filename, 'r') as fh:
        for line in fh:
            if line[0] == '>':
                header = line[1:].rstrip()
                id.add(header.split()[0])
    return id

ids = getSeqID(sys.argv[1])

OUT = open(sys.argv[3], 'w')

with open(sys.argv[2], 'r') as f:
    for line in f:
        lsp = line.rstrip().split("\t")
        if lsp[2] == "gene":
            OUT.write(line)
        elif lsp[2] == "mRNA":
            x = lsp[8].split(";")[0].split("=")[1]
            if x in ids:
                OUT.write(line)
        else:
            x = lsp[8].split(";")[1].split("=")[1]
            if x in ids:
                OUT.write(line)

OUT.close()
