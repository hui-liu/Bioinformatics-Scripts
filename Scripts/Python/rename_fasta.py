import sys
from collections import OrderedDict

def parseFasta(filename):
    fas = OrderedDict()
    id = None
    with open(filename, 'r') as fh:
        for line in fh:
            if line[0] == '>':
                header = line[1:].rstrip()
                id = header.split()[0]
                fas[id] = []
            else:
                fas[id].append(line.rstrip().upper())
        for id, seq in fas.iteritems():
            fas[id] = ''.join(seq)
    return fas

idmap = OrderedDict()
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = line.split()
        idmap[lsp[0]] = lsp[1]

seq_dict = parseFasta(sys.argv[2])

OUT = open(sys.argv[3], 'w')
for i in idmap:
    id = idmap[i]
    OUT.write(">"+id+"\n"+seq_dict[i]+"\n")

OUT.close()

