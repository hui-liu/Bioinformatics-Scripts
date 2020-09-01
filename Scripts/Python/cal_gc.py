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

def revcomp(seq):
    bases = 'ABCDGHKMNRSTUVWXYabcdghkmnrstuvwxyTVGHCDMKNYSAABWXRtvghcdmknysaabwxr'
    complement_dict = {bases[i]:bases[i+34] for i in range(34)}
    seq = reversed(seq)
    result = [complement_dict[base] for base in seq]
    return ''.join(result)

def parseBed(filename):
    bed = {}
    with open(filename, 'r') as fh:
        for line in fh:
            chr, start, end, strand = line.split()
            bed.setdefault(chr, []).append([int(start)+1, int(end), strand])
    return bed

seq_dict = parseFasta(sys.argv[1])
bed_dict = parseBed(sys.argv[2])

for chr in bed_dict:
    for item in bed_dict[chr]:
        start, end, strand = item
        if strand == "+":
           seq = seq_dict[chr][start:end]
           gc = (seq.count("C") + seq.count("G")) / float(len(seq))
           gc = round(gc, 2)
           print gc
        else:
           seq = revcomp(seq_dict[chr][start: end])
           gc = (seq.count("C") + seq.count("G")) / float(len(seq))
           gc = round(gc, 2)
           print gc
