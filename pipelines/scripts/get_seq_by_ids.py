import sys

def parseFasta(filename):
    fas = {}
    id = None
    with open(filename, 'r') as fh:
        for line in fh:
            if line[0] == '>':
                header = line[1:].rstrip()
                id = header.split()[0]
                fas[id] = []
            else:
                fas[id].append(line.rstrip().upper())
        for id, seq in fas.items():
            fas[id] = ''.join(seq)
    return fas

seq_dict = parseFasta(sys.argv[1])

ids = []
with open(sys.argv[2], 'r') as f:
    for line in f:
        ids.append(line.rstrip())

for i in ids:
    print ">" + i + "\n" + seq_dict[i]

