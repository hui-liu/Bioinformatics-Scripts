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
ids = set()
with open(sys.argv[2], 'r') as f:
    for line in f:
        ids.add(line.rstrip().split()[0])

OUT = open(sys.argv[3], 'w')
for i in seq_dict:
    if i not in ids:
        OUT.write(">" + i + "\n" + seq_dict[i] + "\n")
OUT.close()
