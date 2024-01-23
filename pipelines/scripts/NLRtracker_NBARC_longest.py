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
        for id, seq in fas.iteritems():
            fas[id] = ''.join(seq)
    return fas

seq_dict = parseFasta(sys.argv[1])
species = sys.argv[2]

seq_info = {}
for i in seq_dict:
    id = i.split("_(")[0]
    seq_info.setdefault(id, []).append([i, len(seq_dict[i])])

OUT = open(sys.argv[3], 'w')
for i in seq_info:
    if len(seq_info[i]) == 1:
        OUT.write(">" + species + "|" + i + "\n" + seq_dict[i] + "\n")
    else:
        tmp = sorted(seq_info[i], key=lambda x: -x[1])
        OUT.write(">" + species + "|" + i  + "\n" + seq_dict[tmp[0][0]] + "\n")
OUT.close()
