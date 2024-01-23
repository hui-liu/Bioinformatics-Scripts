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
    #comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', '-': '-', 'N': 'N'}
    sequence = sequence.upper()
    sequence_rev = ''
    for i in range(1, len(sequence)+1):
        sequence_rev += comp_dict[sequence[-i]]
    return sequence_rev

def getSeq(coord, reference):
    seq = ''
    chr, start, end, strand = coord
    start, end = int(start), int(end)
    seq += reference[(start-1):end]
    if strand == "-":
        seq = reverse_comp(seq)
    return seq

LTR = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line[0] == "#": continue
        chrom, source, feat_type, start, end, score, strand, frame, info = line.rstrip().split("\t")
        d_info = dict([v.split('=') for v in info.split(';')])
        if feat_type in ["LTR_retrotransposon", "Copia_LTR_retrotransposon", "Gypsy_LTR_retrotransposon"]:
            LTR.append([chrom, start, end, "+", d_info["ID"]])

genome = parseFasta(sys.argv[2])

OUT = open(sys.argv[3], 'w')
for i in LTR:
    chr = i[0]
    ref_seq = genome[chr]
    seq = getSeq(i[:-1], ref_seq)
    OUT.write(">"+i[-1]+"\n"+seq+"\n")
OUT.close()
