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
id = sys.argv[2]
seq = seq_dict[id]

with open(sys.argv[3], 'w') as fw:
    fw.write(">" + id + "\n" + seq + "\n")

