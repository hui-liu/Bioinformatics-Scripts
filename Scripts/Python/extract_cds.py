import sys

def gff2dict(filename):
    """
    The phase of a CDS feature depends on the associated upstream CDS feature.
    If there is the length/3 of the previous CDS feature leaves a remainder of 1,
    your CDS feature requires a phase of 2 (the last base of the previous and the
    first two bases of this form a codon triplett). Standalone features and the
    first features always have a 0 phase value.
    """
    from collections import OrderedDict
    ids = {}
    cds_info = OrderedDict()
    with open(filename, 'r') as f:
        for line in f:
            lsp = line.split()
            chr, source, type, start, end, score, strand, phase, attibutes = lsp
            d_attr = dict([v.split('=') for v in attibutes.strip(';').split(';')])
            if type == 'mRNA':
                id = d_attr['ID']
                ids.setdefault(chr, []).append([id, strand])
            elif type == 'CDS':
                id = d_attr['Parent']
                cds_info.setdefault(id, []).append([chr, int(start), int(end), strand])
    return ids, cds_info

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


def getCDS(cds_records, reference):
    seq = ''
    for record in cds_records:
        chr, start, end, strand = record
        seq += reference[(start-1):end]
    if strand == "-":
        seq = reverse_comp(seq)
    return seq


ids_dict, cds_info_dict = gff2dict(sys.argv[1])
genome = parseFasta(sys.argv[2])


for chr in genome:
    ref_seq = genome[chr]
    if chr not in ids_dict: continue
    for item in ids_dict[chr]:
        id, strand = item
        cds_info = cds_info_dict[id]
        cds_seq = getCDS(cds_info, ref_seq)
        print ">" + id
        print cds_seq
