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

id_map = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line[0] == "#": continue
        if not line.strip(): continue
        chr, source, type, start, end, score, strand, phase, attibutes = line.rstrip().split("\t")
        d_attr = dict([v.split('=') for v in attibutes.strip(';').split(';')])
        if type == 'mRNA':
            id_map.setdefault(d_attr["Parent"], []).append(d_attr["ID"])

cds_dict = parseFasta(sys.argv[2])
pep_dict = parseFasta(sys.argv[3])

cds_out = open(sys.argv[4], 'w')
pep_out = open(sys.argv[5], 'w')

for i in id_map:
    ids_len = [[x, len(cds_dict[x])] for x in id_map[i]]
    if len(ids_len) == 1:
        id = ids_len[0][0]
        cds_out.write(">" + id + "\n" + cds_dict[id] + "\n")
        pep_out.write(">" + id + "\n" + pep_dict[id] + "\n")
    else:
        ids_sorted = sorted(ids_len, key=lambda x: -int(x[1]))
        id = ids_sorted[0][0]
        cds_out.write(">" + id + "\n" + cds_dict[id] + "\n")
        pep_out.write(">" + id + "\n" + pep_dict[id] + "\n")

cds_out.close()
pep_out.close()

