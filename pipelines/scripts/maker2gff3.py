import sys

input = sys.argv[1]
source = sys.argv[2]
out = open(sys.argv[3], 'w')

if source != "genemark":
    source = source + "_masked"

source2 = source.split("_")[0]

with open(input, 'r') as f:
    for line in f:
        if line[0] == "#": continue
        lsp = line.rstrip().split("\t")
        if lsp[2] == "contig": continue
        if lsp[1] == source:
            if lsp[2] == "match":
                d = dict([x.split("=")for x in lsp[8].split(";")])
                out.write("\t".join([lsp[0], source2, "gene"] + lsp[3:8] + ["ID=" + d["ID"] + ".g;Name=" + d["Name"]]) + "\n")
                out.write("\t".join([lsp[0], source2, "mRNA"] + lsp[3:8] + ["ID=" + d["ID"] + ";Parent=" + d["ID"] + ".g"]) + "\n")
            else:
                d = dict([x.split("=")for x in lsp[8].split(";")])
                out.write("\t".join([lsp[0], source2, "exon"] + lsp[3:8] + ["ID=" + d["ID"] + ".exon;Parent=" + d["Parent"]]) + "\n")
                out.write("\t".join([lsp[0], source2, "CDS"] + lsp[3:8] + ["ID=" + d["ID"] + ".cds;Parent=" + d["Parent"]]) + "\n")
out.close()
