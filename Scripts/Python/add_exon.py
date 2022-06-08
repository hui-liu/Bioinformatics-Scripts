import sys
from collections import OrderedDict

mRNA = OrderedDict()
five_UTR = {}
three_UTR = {}
CDS_tmp = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        lsp = line.rstrip().split('\t')
        chr, source, type, start, end, score, strand, phase, attibutes = lsp
        d_attr = dict([v.split('=') for v in attibutes.strip(';').split(';')])
        if type == "mRNA":
            mRNA[d_attr["ID"]] = lsp
        elif type == "five_prime_UTR":
            five_UTR[d_attr["Parent"]] = lsp
        elif type == "three_prime_UTR":
            three_UTR[d_attr["Parent"]] = lsp
        elif type == "CDS":
            CDS_tmp.setdefault(d_attr["Parent"], []).append(lsp)

CDS = {}
for i in CDS_tmp:
    CDS[i] =  sorted(CDS_tmp[i], key = lambda x: int(x[3]))

OUT = open(sys.argv[2], 'w')
for i in mRNA:
    a = mRNA[i]
    b = a[:2] + ["gene"] + a[3:]
    c = a[-1].replace(";", ";Parent=" + i + ";")
    d = a[:-1] + [c]
    OUT.write("\t".join(b) + "\n" + "\t".join(d) + "\n")
    start, end = a[3:5]
    strand = a[6]
    if len(CDS[i]) == 1:
        e = a[:2] + ["exon"] + a[3:8] + [CDS[i][0][-1].replace("_cds_", "_exon_")]
        if strand == "+":
            if i in five_UTR:
                OUT.write("\t".join(five_UTR[i]) + "\n" + "\t".join(e) + "\n" + "\t".join(CDS[i][0]) + "\n")
            else:
                OUT.write("\t".join(e) + "\n" + "\t".join(CDS[i][0]) + "\n")
            if i in three_UTR:
                OUT.write("\t".join(three_UTR[i]) + "\n")
        else:
            if i in three_UTR:
                OUT.write("\t".join(three_UTR[i]) + "\n" + "\t".join(e) + "\n" + "\t".join(CDS[i][0]) + "\n")
            else:
                OUT.write("\t".join(e) + "\n" + "\t".join(CDS[i][0]) + "\n")
            if i in five_UTR:
                OUT.write("\t".join(five_UTR[i]) + "\n")
    elif len(CDS[i]) == 2:
        if strand == "+":
            first_CDS = CDS[i][0]
            last_CDS = CDS[i][-1]
            if i in five_UTR:
                first_exon = first_CDS[:2] + ["exon"] + [five_UTR[i][3]] + first_CDS[4:8] + [first_CDS[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(five_UTR[i]) + "\n" + "\t".join(first_exon) + "\n" + "\t".join(first_CDS) + "\n")
            else:
                first_exon = first_CDS[:2] + ["exon"] + first_CDS[3:8] + [first_CDS[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(first_exon) + "\n" + "\t".join(first_CDS) + "\n")
            if i in three_UTR:
                last_exon = last_CDS[:2] + ["exon"] + [last_CDS[3]] + [three_UTR[i][4]] + last_CDS[5:8] + [last_CDS[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(last_exon) + "\n" + "\t".join(last_CDS) + "\n" + "\t".join(three_UTR[i]) + "\n")
            else:
                last_exon = last_CDS[:2] + ["exon"] + last_CDS[3:8] + [last_CDS[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(last_exon) + "\n" + "\t".join(last_CDS) + "\n")
        else:
            first_CDS = CDS[i][0]
            last_CDS = CDS[i][-1]
            if i in three_UTR:
                first_exon = first_CDS[:2] + ["exon"] + [three_UTR[i][3]] + first_CDS[4:8] + [first_CDS[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(three_UTR[i]) + "\n" + "\t".join(first_exon) + "\n" + "\t".join(first_CDS) + "\n")
            else:
                first_exon = first_CDS[:2] + ["exon"] + first_CDS[3:8] + [first_CDS[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(first_exon) + "\n" + "\t".join(first_CDS) + "\n")
            if i in five_UTR:
                last_exon = last_CDS[:2] + ["exon"] + [last_CDS[3]] + [five_UTR[i][4]] + last_CDS[5:8] + [last_CDS[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(last_exon) + "\n" + "\t".join(last_CDS) + "\n" + "\t".join(five_UTR[i]) + "\n")
            else:
                last_exon = last_CDS[:2] + ["exon"] + last_CDS[3:8] + [last_CDS[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(last_exon) + "\n" + "\t".join(last_CDS) + "\n")

    else:
        if strand == "+":
            first_CDS = CDS[i][0]
            last_CDS = CDS[i][-1]
            if i in five_UTR:
                first_exon = first_CDS[:2] + ["exon"] + [five_UTR[i][3]] + first_CDS[4:8] + [first_CDS[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(five_UTR[i]) + "\n" + "\t".join(first_exon) + "\n" + "\t".join(first_CDS) + "\n")
            else:
                first_exon = first_CDS[:2] + ["exon"] + first_CDS[3:8] + [first_CDS[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(first_exon) + "\n" + "\t".join(first_CDS) + "\n")
            #
            for x in CDS[i][1:-1]:
                y = x[:2] + ["exon"] + x[3:8] + [x[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(y) + "\n" + "\t".join(x) + "\n")
            if i in three_UTR:
                last_exon = last_CDS[:2] + ["exon"] + [last_CDS[3]] + [three_UTR[i][4]] + last_CDS[5:8] + [last_CDS[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(last_exon) + "\n" + "\t".join(last_CDS) + "\n" + "\t".join(three_UTR[i]) + "\n")
            else:
                last_exon = last_CDS[:2] + ["exon"] + last_CDS[3:8] + [last_CDS[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(last_exon) + "\n" + "\t".join(last_CDS) + "\n")
        else:
            first_CDS = CDS[i][0]
            last_CDS = CDS[i][-1]
            if i in three_UTR:
                first_exon = first_CDS[:2] + ["exon"] + [three_UTR[i][3]] + first_CDS[4:8] + [first_CDS[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(three_UTR[i]) + "\n" + "\t".join(first_exon) + "\n" + "\t".join(first_CDS) + "\n")
            else:
                first_exon = first_CDS[:2] + ["exon"] + first_CDS[3:8] + [first_CDS[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(first_exon) + "\n" + "\t".join(first_CDS) + "\n")
            #
            for x in CDS[i][1:-1]:
                y = x[:2] + ["exon"] + x[3:8] + [x[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(y) + "\n" + "\t".join(x) + "\n")
            if i in five_UTR:
                last_exon = last_CDS[:2] + ["exon"] + [last_CDS[3]] + [five_UTR[i][4]] + last_CDS[5:8] + [last_CDS[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(last_exon) + "\n" + "\t".join(last_CDS) + "\n" + "\t".join(five_UTR[i]) + "\n")
            else:
                last_exon = last_CDS[:2] + ["exon"] + last_CDS[3:8] + [last_CDS[8].replace("_cds_", "_exon_")]
                OUT.write("\t".join(last_exon) + "\n" + "\t".join(last_CDS) + "\n")
OUT.close()
