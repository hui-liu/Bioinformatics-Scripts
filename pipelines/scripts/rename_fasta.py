import sys
import textwrap

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
prefix = sys.argv[2]
#out = open(sys.argv[3], 'w')

seq_len = [[s, len(seq_dict[s])] for s in seq_dict]
seq_len.sort(key=lambda x: -x[1])

n = 0
for x, _ in seq_len:
    n += 1
    seq = seq_dict[x]
    #seq_chunks = textwrap.wrap(seq, 60)
    newid = prefix + str(n)
    #out.write(">" + newid + "\n" + "\n".join(seq_chunks) + "\n")
    print(">" + newid + "\n" + seq)
#out.close()
