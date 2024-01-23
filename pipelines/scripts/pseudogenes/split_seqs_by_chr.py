import sys

USAGE = "\nusage: python extract_seq.py genome.fas prefix\n"

if len(sys.argv) != 3:
    print USAGE
    sys.exit()

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

# split
for i in seq_dict:
    outfa = open(prefix + "." + i + ".fa", 'w')
    outfa.write(">" + i + "\n" + seq_dict[i] + "\n")
    outfa.close()
