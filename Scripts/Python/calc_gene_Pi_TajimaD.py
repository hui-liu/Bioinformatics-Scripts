
"""
http://codextechnicanum.blogspot.com/2016/11/nucleotide-diversity-pi-or-pi.html
https://github.com/vcftools/vcftools/blob/f7aee6d26885064d834199052573e31338240b76/src/cpp/variant_file_output.cpp
###################
# Measures nucleotide divergency on a per-site basis.
###################
AC: allele count
AN: allele number
pairs = (AN * (AN-1)
mismatches = AC * (AN - AC) + (AN - AC) * AC
site_Pi = mismatches / pairs

###################
# Measures the nucleotide diversity in windows.
###################
N_variant_sites: Number of sites in a window that have VCF entries
N_variant_site_pairs: Number of possible pairwise mismatches at polymorphic sites within a window
N_mismatches: Number of actual pairwise mismatches at polymorphic sites within a window
N_polymorphic_sites: number of sites within a window where there is at least 1 sample that is polymorphic with respect to the reference alleleimport allele
win_length: the window size
N_kept_chr: the number of samples

N_variant_site_pairs = sum(AN * (AN - 1))
N_monomorphic_sites = win_length - N_polymorphic_sites
N_comparisons = (N_kept_chr * (N_kept_chr - 1))

N_pairs = N_variant_site_pairs + (N_monomorphic_sites * N_comparisons)
N_mismatches = sum(AC * (AN - AC) + (AN-AC) * AC)
window_Pi = N_mismatches / N_pairs

###################
# Measures the nucleotide diversity on a single gene.
###################

gene_length: the gene length
N_kept_chr: the number of samples

N_variant_site_pairs = sum(AN * (AN - 1))
N_monomorphic_sites = gene_length - N_polymorphic_sites
N_comparisons = (N_kept_chr * (N_kept_chr - 1))

N_pairs = N_variant_site_pairs + (N_monomorphic_sites * N_comparisons)
N_mismatches = sum(AC * (AN - AC) + (AN-AC) * AC)
gene_Pi = N_mismatches / N_pairs
"""
# hui.liu@umu.se/liuhui.bfu@gmail.com
# calculate Pi and Tajima's D for a single gene

import sys
from math import sqrt
import gzip

def parseVCF(file):
    AC = [] # allele count for each SNP
    AN = [] # allele number for each SNP
    N_indvs = 0 # the number of samples
    pos = []
    with gzip.open(file, 'r') as f:
        for line in f:
            if line[:2] == "##": continue
            else:
                lsp = line.rstrip().split("\t")
                if lsp[0] == "#CHROM":
                    N_indvs = (len(lsp) - 9)
                else:
                    pos.append([lsp[0], lsp[1]])
                    alleles = [x.split(":")[0] for x in lsp[9:]]
                    ref = alleles.count("0/0") + alleles.count("0|0")
                    alt = alleles.count("1/1") + alleles.count("1|1")
                    het = alleles.count("0/1") + alleles.count("1/0") + alleles.count("0|1") + alleles.count("1|0")
                    ac = alt * 2 + het
                    an = 2 * (ref + alt + het)
                    AC.append(ac)
                    AN.append(an)
    return pos, AC, AN, N_indvs

def parseGFF(file):
    res = {}
    with open(file, 'r') as f:
        for line in f:
            if line[0] == "#": continue
            lsp = line.split()
            if lsp[2] == "gene":
                geneid = lsp[8].split(";")[0].split("=")[1]
                chr, start, end = lsp[0], int(lsp[3]), int(lsp[4])
                res.setdefault(chr, []).append([geneid, start, end])
    return res


def geneLen(file):
    res = {}
    with open(file) as f:
        for line in f:
            if line[0] == "#": continue
            lsp = line.split()
            if lsp[2] == "gene":
                gene = lsp[8].split(";")[0].split("=")[1]
                start, end = int(lsp[3]), int(lsp[4])
                res[gene] = end - start + 1
    return res

def TajimaD(pi, S, n):
    #https://github.com/vcftools/vcftools/blob/f7aee6d26885064d834199052573e31338240b76/src/cpp/variant_file_output.cpp
    #https://github.com/AndersenLab/VCF-kit/blob/d44b597381d40b84f752b9dd7b345b1b76fa2e00/vcfkit/tajima.py
    # Carlson et al. Genome Res (2005)
    a1 = sum([1.0 / i for i in xrange(1, n)])
    a2 = sum([1.0 / i**2 for i in xrange(1, n)])
    n  = float(n)
    b1 = (n+1) /(3 * (n-1))
    b2 = (2 * ((n*n)+n+3)) / (9*n*(n-1))
    c1 = b1 - (1/a1)
    c2 = b2 - ((n+2) / (a1*n)) + (a2 / (a1*a1))
    e1 = c1 / a1
    e2 = c2 / ((a1*a1) + a2)

    tw = (S/a1)
    var = float(e1*S) + (e2*S*(S-1))
    D  = (pi - tw) / sqrt(var)
    return D

_pos, AC, AN, N_indvs = parseVCF(sys.argv[1])
geneid = sys.argv[2]
gene_length_dict = geneLen(sys.argv[3])
gene_length = gene_length_dict[geneid]

"""
site_Pi = []
for i, x in enumerate(AC):
    y = AN[i]
    pairs = y * (y-1)
    mismatches = x * (y - x) + (y - x) * x
    pi = mismatches / float(pairs)
    site_Pi.append(_pos[i]+[pi])
"""

# calculate Pi for a sigle gene
N_variant_site_pairs = sum([x*(x-1) for x in AN])
N_polymorphic_sites = len(AN)
N_monomorphic_sites = gene_length - N_polymorphic_sites
N_kept_chr = N_indvs * 2
N_comparisons = (N_kept_chr * (N_kept_chr - 1))
N_pairs = N_variant_site_pairs + (N_monomorphic_sites * N_comparisons)

N_mismatches = 0
for i, x in enumerate(AC):
    y = AN[i]
    mismatches = x * (y - x) + (y - x) * x
    N_mismatches += mismatches

gene_Pi = N_mismatches / float(N_pairs)

# calculate Tajima's D
# https://genome.cshlp.org/content/15/11/1553.full#sec-4
piqi = 0.0 # pi is the derived (nonancestral) allele frequency of the ith SNP
for i, x in enumerate(AC):
    y = AN[i]
    pi = x/float(y)
    qi = 1.0-pi
    piqi += pi * qi

pi = 2 * piqi*N_kept_chr/(N_kept_chr-1.0)

gene_tajamdD = TajimaD(pi, N_polymorphic_sites, N_kept_chr)


out = open(geneid + ".tsv", "w")
out.write("\t".join(["geneid", "Pi", "TajimaD"]) + "\n")
out.write("\t".join([geneid, str(gene_Pi), str(gene_tajamdD)]) + "\n")

out.close()
