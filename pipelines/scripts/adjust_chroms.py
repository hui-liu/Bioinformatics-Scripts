import sys
from collections import OrderedDict


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

def reverse_comp(sequence):
    comp_dict = {
        'A': 'T',
        'B': 'V',
        'C': 'G',
        'D': 'H',
        'G': 'C',
        'H': 'D',
        'M': 'K',
        'N': 'N',
        'R': 'Y',
        'S': 'S',
        'T': 'A',
        'U': 'A',
        'V': 'B',
        'W': 'W',
        'X': 'X',
        'Y': 'R'}
    sequence = sequence.upper()
    sequence_rev = ''
    for i in range(1, len(sequence)+1):
        sequence_rev += comp_dict[sequence[-i]]
    return sequence_rev


seq_dict = parseFasta(sys.argv[1])
chr_map = OrderedDict()
with open(sys.argv[2], 'r') as f:
    for line in f:
        id, chr, strand = line.rstrip().split("\t")
        chr_map[id] = [chr, strand]


for s in chr_map:
    seq = seq_dict[s]
    chr, strand = chr_map[s]
    if strand == "-":
        seq = reverse_comp(seq)
    print(">" + chr + "\n" + seq)
