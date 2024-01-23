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

def parseBUSCOtable(filename):
    res = {}
    with open(filename, 'r') as fh:
        for line in fh:
            if line[0] == "#": continue
            lsp = line.rstrip().split("\t")
            if lsp[1] == "Missing": continue
            id = lsp[0]
            chr, start, end = lsp[2:5]
            res[chr + ":" + start + "-" + end] = id
    return res

seq_dict = parseFasta(sys.argv[1])
id_dict = parseBUSCOtable(sys.argv[2])
gff_dir = sys.argv[3]
OUT = open(sys.argv[4], 'w')

CDS = {}
for i in id_dict:
    if i in seq_dict:
        transcript_id = i.replace(":", "_").replace("-", "_")
        with open(gff_dir + "/" + id_dict[i] + ".gff", 'r') as f:
            for line in f:
                lsp = line.rstrip().split("\t")
                if lsp[2] == "CDS":
                    CDS.setdefault(transcript_id, []).append(lsp[:8] + ["transcript_id " + '"' + transcript_id + '"'])
for s in CDS:
    t = CDS[s]
    if t[0][6] == "-":
        t.sort(key = lambda x: -int(x[3]))
    for u in t:
        OUT.write("\t".join(u) + "\n")

OUT.close()
