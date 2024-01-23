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

def countgaps(seq):
    return len([x for x in seq.split("N") if x]) - 1

seq_dict = parseFasta(sys.argv[1])

num = 0
for chr in seq_dict:
    n = countgaps(seq_dict[chr])
    num += n
    print chr, n

print "Total gaps:", num
