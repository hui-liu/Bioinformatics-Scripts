import sys
import gzip

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

def locateBin(a, b):
    result = True
    if a < b[0] or b[1] < a:
        result = False
    return result

def gff2dict(filename):
    cds = {}
    cds_info = {}
    with open(filename, 'r') as f:
        for line in f:
            lsp = line.split()
            chr, source, type, start, end, score, strand, phase, attibutes = lsp
            d_attr = dict([v.split('=') for v in attibutes.strip(';').split(';')])
            if type == 'CDS':
                id = d_attr['Parent']
                cds.setdefault(chr, []).append([int(start), int(end), strand, int(phase), id])
                cds_info.setdefault(id, []).append([int(start), int(end), int(phase)])
    return cds, cds_info

def snp_cds(cds_dict, chr, SNPpos):
    for i in cds_dict[chr]:
        if locateBin(SNPpos, i):
            result = [chr] + i
    return result

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
    comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', '-': '-', 'N': 'N'}
    sequence = sequence.upper()
    sequence_rev = ''
    for i in range(1, len(sequence)+1):
        sequence_rev += comp_dict[sequence[-i]]
    return sequence_rev

def lstIndex(targetlist, dblist):
    for i,j in enumerate(dblist):
        if j == targetlist:
            return i

def getsnp(line):
    lsp = line.split()
    return [lsp[0], int(lsp[1]), lsp[3], lsp[4]]

def getcodon(cds_seq, snp_block_pos):
    for i in range(0, len(cds_seq), 3):
        if i <= snp_block_pos < i+3:
            codon = cds_seq[i:i+3]
            snp_codon_pos = snp_block_pos - i
            return codon, snp_codon_pos

def snp_codon(snp_position, cds_record, reference):
    if cds_record[3] == '+':
        cds_seq = reference[cds_record[1]-1+cds_record[4]:cds_record[2]]
        snp_block_pos = snp_position - cds_record[1] - cds_record[4]
        codon, snp_codon_pos = getcodon(cds_seq, snp_block_pos)
    else:
        cds_seq = reverse_comp(reference[cds_record[1]-1:cds_record[2]-cds_record[4]])
        snp_block_pos = cds_record[2] - snp_position - cds_record[4]
        codon, snp_codon_pos = getcodon(cds_seq, snp_block_pos)
    return codon, snp_codon_pos

# input
codonDegenDict = get_degeneracy(standardCodonTable)
cds_dict, cds_info_dict = gff2dict(sys.argv[1])
genome = parseFasta(sys.argv[2])

# output
zero_fold_out = gzip.open(sys.argv[4], 'w')
four_fold_out = gzip.open(sys.argv[5], 'w')

with gzip.open(sys.argv[3], 'r') as f:
    for line in f:
        if line[0] == "#":
            zero_fold_out.write(line)
            four_fold_out.write(line)
        else:
            chr, pos, ref, alt = getsnp(line)
            if chr in cds_dict:
                ref_seq = genome[chr]
                cds_record = snp_cds(cds_dict, chr, pos)
                result = snp_codon(pos, cds_record, ref_seq)
                if result:
                    codon, snp_codon_pos = result
                    if len(codon) == 3:
                        i_th_fold = codonDegenDict[codon][snp_codon_pos]
                        if i_th_fold == 0:
                            zero_fold_out.write(line)
                        elif i_th_fold == 4:
                            four_fold_out.write(line)

zero_fold_out.close()
four_fold_out.close()
