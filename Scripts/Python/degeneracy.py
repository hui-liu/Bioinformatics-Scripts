#!/usr/bin/python
#Author: liuhui
#EMail: liuhui@bjfu.edu.cn / hui.liu@umu.se
#Description: generate 0 and 4fold sites

import sys
import math


USAGE = "\nusage: python %s <gff3> <fasta> <0-fold output> <4-fold output>\n" % sys.argv[0]

if len(sys.argv) != 5:
    print USAGE
    sys.exit()


standardCodonTable = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

#https://github.com/russcd/vcf2MK/blob/master/src/fold_counts.h
def get_degeneracy(codon_table):
    nucs = ['A', 'T', 'G', 'C']
    codon_degen = {}
    for codon in codon_table:
        degen_list = []
        for i in range(len(codon)):
            degen = 0
            for j in nucs:
                if i == 0:
                    mod_codon = j + codon[1:]
                elif i == 1:
                    mod_codon = codon[0] + j + codon[2]
                elif i == 2:
                    mod_codon = codon[0:2] + j

                if codon_table[codon] == codon_table[mod_codon]:
                    degen += 1

            degen_list.append(degen)

            for site in range(len(degen_list)):
                if degen_list[site] == 1:
                    degen_list[site] = 0
        codon_degen[codon] = degen_list

    return codon_degen

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

def getcodon(cds_seq, n_bases, strand):
    number = int(math.ceil(n_bases/float(3)))
    if strand == "+":
        codon = cds_seq[(number-1)*3: number*3]
    else:
        cds_seq = cds_seq[::-1]
        codon = cds_seq[(number-1)*3: number*3]
        codon = codon[::-1]

    m = n_bases % 3
    if m == 0:
        codon_pos = 2
    elif m == 1:
        codon_pos = 0
    else:
        codon_pos = 1
    if strand == "-":
        codon_pos = 2 - codon_pos
    return codon, codon_pos

# input
codonDegenDict = get_degeneracy(standardCodonTable)
ids_dict, cds_info_dict = gff2dict(sys.argv[1])
genome = parseFasta(sys.argv[2])
chroms = [chr for chr in genome]

# output
zero_fold_out = open(sys.argv[3], 'w')
four_fold_out = open(sys.argv[4], 'w')

# loop for each chomosome
for chr in chroms:
    ref_seq = genome[chr]
    # scan each mRNA in a chorosome
    if chr not in ids_dict: continue
    for item in ids_dict[chr]:
        id, strand = item
        full_cds = cds_info_dict[id]
        # extract the whole CDS region
        cds_seq = getCDS(full_cds, ref_seq)
        # scan each part of the CDS
        n = 0
        for r in full_cds:
            start, end = r[1:3]
            for pos in range(start, end+1):
                n += 1
                codon, codon_pos = getcodon(cds_seq, n, strand)
                if codon not in codonDegenDict: continue
                i_th_fold = codonDegenDict[codon][codon_pos]
                #print id, pos, codon, codon_pos, i_th_fold, r, cds_seq
                if i_th_fold == 0:
                    zero_fold_out.write("\t".join([chr, str(pos-1), str(pos), id]) + "\n")
                elif i_th_fold == 4:
                    four_fold_out.write("\t".join([chr, str(pos-1), str(pos), id]) + "\n")
zero_fold_out.close()
four_fold_out.close()
