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

# The standard code currently allows initiation from UUG and CUG in addition to AUG (https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG1)
stop_codon = set(["TAA", "TAG", "TGA"])
start_codon = set(["ATG", "TTG", "CTG"])
seq_dict = parseFasta(sys.argv[1])
raw_gff = sys.argv[2]
out = open(sys.argv[3], 'w')

incomplete_ids = set()
for id in seq_dict:
    seq = seq_dict[id]
    split_seq = [seq[i:i+3] for i in range(0, len(seq), 3)]
    a = set(split_seq[1:-1])
    if stop_codon.intersection(a):
        incomplete_ids.add(id)
    if split_seq[0] not in start_codon:
        incomplete_ids.add(id)
    if split_seq[-1] not in stop_codon:
        incomplete_ids.add(id)

id_map = {}
with open(raw_gff, 'r') as f:
    for line in f:
        if line[0] == "#": continue
        if not line.strip(): continue
        chr, source, type, start, end, score, strand, phase, attibutes = line.rstrip().split("\t")
        d_attr = dict([v.split('=') for v in attibutes.strip(';').split(';')])
        if type == 'mRNA':
            if d_attr["ID"] not in incomplete_ids:
                id_map.setdefault(d_attr["Parent"], []).append(d_attr["ID"])

with open(raw_gff, 'r') as f:
    for line in f:
        if line[0] == "#": continue
        if not line.strip(): continue
        chr, source, type, start, end, score, strand, phase, attibutes = line.rstrip().split("\t")
        if source == ".":
            source = "PASA"
        d_attr = dict([v.split('=') for v in attibutes.strip(';').split(';')])
        if type == "gene":
           if d_attr["ID"] in id_map:
               out.write("\t".join([chr, source, type, start, end, score, strand, phase, attibutes]) + "\n")
        elif type == "mRNA":
           if d_attr["ID"] not in incomplete_ids:
               out.write("\t".join([chr, source, type, start, end, score, strand, phase, attibutes]) + "\n")
        else:
           if d_attr["Parent"] not in incomplete_ids:
               out.write("\t".join([chr, source, type, start, end, score, strand, phase, attibutes]) + "\n")

out.close()
