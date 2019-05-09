import sys
import gzip

# (1) gene
# (1.1) CDS
# (1.1.1) 0-fold nonsynonymous
# (1.1.2) fourfold synonymous
# (1.2) 3'UTR
# (1.3) 5'UTR
# (1.4) intron

# (2) intergenic

def gff2dict(filename):
    gene = {}
    cds = {}
    three_utr = {}
    five_utr = {}
    with open(filename, 'r') as f:
        for line in f:
            lsp = line.split()
            chr, source, type, start, end, score, strand, phase, attibutes = lsp
            d_attr = dict([v.split('=') for v in attibutes.strip(';').split(';')])
            if type == 'gene':
                gene.setdefault(chr, []).append([int(start), int(end)])
            elif type == 'CDS':
                cds.setdefault(chr, []).append([int(start), int(end), strand, int(phase)])
            elif type == 'three_prime_UTR':
                three_utr.setdefault(chr, []).append([int(start), int(end)])
            elif type == 'five_prime_UTR':
                five_utr.setdefault(chr, []).append([int(start), int(end)])

    return gene, cds, three_utr, five_utr

def getvcf(line):
    lsp = line.split()
    return [lsp[0], int(lsp[1])]

def locateBin(a, b):
    result = True
    if a < b[0] or b[1] < a:
        result = False
    return result

def is_feature(feature_dict, chr, pos):
    result = False
    for i in feature_dict[chr]:
        if locateBin(pos, i):
            result = True
    return result

def is_CDS(cds_dict, chr, pos):
    result = False
    for i in cds_dict[chr]:
        start, end, strand, phase = i
        if strand == "+":
            if locateBin(pos, [start+phase, end]):
                return True
        else:
            if locateBin(pos, [start, end-phase]):
                return True

#
gene_dict, cds_dict, three_utr_dict, five_utr_dict = gff2dict(sys.argv[1])

#
cds = set()
three_utr = []
five_utr = []
intron = []
intergenic = []

csd_out = gzip.open(sys.argv[3], 'w')
three_utr_out = gzip.open(sys.argv[4], 'w')
five_utr_out = gzip.open(sys.argv[5], 'w')
intron_out = gzip.open(sys.argv[6], 'w')
intergenic_out = gzip.open(sys.argv[7], 'w')

with gzip.open(sys.argv[2], 'r') as f:
    for line in f:
        if line[0] == "#":
            csd_out.write(line)
            three_utr_out.write(line)
            five_utr_out.write(line)
            intron_out.write(line)
            intergenic_out.write(line)
        else:
            chr, pos = getvcf(line)
            if chr in gene_dict:
                if is_feature(gene_dict, chr, pos):
                    # "gene"
                    if is_CDS(cds_dict, chr, pos):
                        cds.add(line)
                        csd_out.write(line)
                    elif chr in three_utr_dict and is_feature(three_utr_dict, chr, pos):
                        three_utr.append(line)
                    elif chr in five_utr_dict and is_feature(five_utr_dict, chr, pos):
                        five_utr.append(line)
                    else:
                        intron.append(line)
                        intron_out.write(line)
                else:
                    #"intergenic"
                    intergenic_out.write(line)
            else:
                # "intergenic"
                intergenic_out.write(line)

for i in three_utr:
    if i not in cds:
        three_utr_out.write(i)
for i in five_utr:
    if i not in cds and i not in three_utr:
        five_utr_out.write(i)


csd_out.close()
three_utr_out.close()
five_utr_out.close()
intron_out.close()
intergenic_out.close()
