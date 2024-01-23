import re
import sys

#annot = "interproscan.tsv"
annot = sys.argv[1]
go_annot = {}
with open(annot, 'r') as f:
    for line in f:
        if "GO:" in line:
            GO = re.findall('GO:\d+', line)
            gene_id = line.split("\t")[0]
            go_annot.setdefault(gene_id, set()).update(GO)

with open(sys.argv[2], 'w') as f:
#    f.write("go_accession\tgene" + "\n")
    for i in go_annot:
        for j in go_annot[i]:
            f.write(i + '\t' + j + "\n")
