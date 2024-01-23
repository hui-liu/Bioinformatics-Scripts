import sys

gene_len = []
exons = {}
CDS = {}
intron = {}
utr = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        chr, source, type, start, end, score, strand, phase, attibutes = line.rstrip().split("\t")
        if type == "gene":
            gene_len.append(int(end) - int(start) + 1)
        elif type == "exon":
            id = attibutes.split(";")[1].split("=")[1]
            exons.setdefault(id, []).append(int(end) - int(start) + 1)
        elif type == "CDS":
            id = attibutes.split(";")[1].split("=")[1]
            CDS.setdefault(id, []).append(int(end) - int(start) + 1)
        elif type == "intron":
            id = attibutes.split("=")[1]
            intron.setdefault(id, []).append(int(end) - int(start) + 1)
        elif type == "three_prime_UTR" or type == "five_prime_UTR":
            id = attibutes.split(";")[1].split("=")[1]
            utr.setdefault(id, []).append(int(end) - int(start) + 1)



exon_len = []
exon_num = []
for i in exons:
    exon_len.append(sum(exons[i]))
    exon_num.append(len(exons[i]))

CDS_len = []
CDS_num = []
for i in CDS:
    CDS_len.append(sum(CDS[i]))
    CDS_num.append(len(CDS[i]))

intron_len = []
intron_num = []
for i in intron:
    intron_len.append(sum(intron[i]))
    intron_num.append(len(intron[i]))

utr_len = []
utr_num = []
for i in utr:
    utr_len.append(sum(utr[i]))
    utr_num.append(len(utr[i]))

average_gene_len = round(sum(gene_len)/float(len(gene_len)), 1)
average_exon_len = round(sum(exon_len)/float(sum(exon_num)), 1)
exon_num_per_gene = round(sum(exon_num)/float(len(gene_len)), 1)
average_CDS_len = round(sum(CDS_len)/float(len(CDS_len)), 1)
average_PEP_len = round(sum(CDS_len)/float(len(CDS_len) * 3), 1)
#average_intron_len = round(sum(intron_len)/float(sum(intron_num)), 1)
#intron_num_per_gene = round(sum(intron_num)/float(len(gene_len)), 1)
#average_utr_len = round(sum(utr_len)/float(sum(utr_num)), 1)
#utr_num_per_gene = round(sum(utr_num)/float(len(gene_len)), 1)


print("Average gene length:", average_gene_len)
print("Average exon length:", average_exon_len)
print("Average exon number per gene:", exon_num_per_gene)
print("Average CDS length:", average_CDS_len)
print("Average protein length:", average_PEP_len)
#print("Average intron length:", average_intron_len)
#print("Average intron number per gene:", intron_num_per_gene)
#print("Average utr length:", average_utr_len)
#print("Average utr number per gene:", utr_num_per_gene)

