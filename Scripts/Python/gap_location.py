import re
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

seq_dict = parseFasta(sys.argv[1])

for s in seq_dict:
    seq = seq_dict[s]
    matches = list(re.finditer('N+', seq))
    for region_number, match in enumerate(matches, 1):
        print "\t".join(map(str, [s, match.start(), match.end(), match.end() - match.start()]))
