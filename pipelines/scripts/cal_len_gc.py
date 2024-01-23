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

for chr in seq_dict:
    seq = seq_dict[chr]
    gc = round((seq.count("C") + seq.count("G")) / float(len(seq)), 2)
    print chr, len(seq), gc

