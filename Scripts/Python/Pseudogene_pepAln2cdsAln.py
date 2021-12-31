import sys
import gzip
import re

#Author: liuhui
#EMail: liuhui@bjfu.edu.cn
#Description: generate cds alignment based on the corresponding protein alignment
#Date: 2021-12-31
#Python 2.7.18

# pgenes.align.gz: generated from PseudoPipe pipeline
# genome.fasta: genome file with FASTA format
# cds.fasta: cds sequences
# pgenes_cds_aln.axt: the output file with axt format

USAGE = "\nusage: python %s <pgenes.align.gz> <genome.fasta> <cds.fasta> <pgenes_cds_aln.axt>\n" % sys.argv[0]

if len(sys.argv) != 5:
    print USAGE
    sys.exit()

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
        for id, seq in fas.iteritems():
            fas[id] = ''.join(seq)
    return fas

complement_table = {
'A': 'T',
'B': 'V',
'C': 'G',
'D': 'H',
'G': 'C',
'H': 'D',
'K': 'M',
'M': 'K',
'N': 'N',
'R': 'Y',
'S': 'S',
'T': 'A',
'U': 'A',
'V': 'B',
'W': 'W',
'X': 'X',
'Y': 'R',
'a': 't',
'b': 'v',
'c': 'g',
'd': 'h',
'g': 'c',
'h': 'd',
'k': 'm',
'm': 'k',
'n': 'n',
'r': 'y',
's': 's',
't': 'a',
'u': 'a',
'v': 'b',
'w': 'w',
'x': 'x',
'y': 'r'
}


def rev_comp(seq):
    new_seq = []
    line = seq.rstrip()
    for letter in line:
        complement_letter = complement_table[letter]
        new_seq.append(complement_letter)
    new_seq.reverse()
    return "".join(new_seq)

def get_seq(seq_dict, chr, start, end, strand):
    # samtools faidx genome.fasta chr:start-end
    # samtools -i faidx genome.fasta chr:start-end
    seq = seq_dict[chr][start-1:end]
    if strand == "-":
        seq = rev_comp(seq)
    return seq

def pal2nal(pal, nuc_seq):
    nal = ''
    n = 0
    for i in range(len(pal)):
        if pal[i] == "-":
            nal += "---"
        else:
            nal += nuc_seq[n: n+3]
            n += 3
    return nal

res = {}
with gzip.open(sys.argv[1], 'r') as f:
    for line in f:
        if line[0] == ">":
            coors = re.findall(r'[^ ]*\(.*?\)', line)
            pseudogene_coors = [[int(j) for j in i.split("..")] for i in coors[0].split("(")[1].rstrip(")").split()]
            query_coors = coors[1]
            q_strand = "+"
            if query_coors[0:3] == "com": q_strand = "-"
            if len(coors) > 2:
                query_coors = coors[-2]
            query_coors = [[int(j) for j in i.split("..")] for i in query_coors.split("(")[1].rstrip(")").split()]
            lsp = line.split()
            # the query_start is the coors of protein sequence
            pseudogene, _, start, end, strand, query, query_start, query_end, query_len, frac = lsp[0:10]
            num_of_ins, num_of_dels, num_of_shifts, num_of_stops, expect, ident, polya = lsp[10:17]
            chr = pseudogene.split("===")[0].lstrip(">")
            id = chr + "_" + start + "_" + end
            res[id] = [chr, int(start), int(end), strand, query, int(query_start), int(query_end), int(query_len), float(frac),
                       int(num_of_ins), int(num_of_dels), int(num_of_shifts), int(num_of_stops), float(expect), float(ident), int(polya), pseudogene_coors, query_coors, q_strand]
        else:
            res[id].append(line.rstrip("\n"))

genome_dict = parseFasta(sys.argv[2])
cds_dict = parseFasta(sys.argv[3])

pep_aln = []
for i in res:
    chr, start, end, strand, query, query_start, query_end, query_len, frac, num_of_ins, num_of_dels, num_of_shifts, num_of_stops, expect, ident, polya, pseudogene_coors, query_coors, q_strand, seq1,  seq2= res[i]
    #pseudogene_seq = "".join([get_seq(genome_dict, chr, s, e, strand) for s,e in pseudogene_coors])
    #cds_seq = "".join([cds_dict[query][s*3-3: e*3] for s,e in query_coors])
    if num_of_shifts > 0: continue
    if num_of_stops > 0: continue
    if q_strand == "+":
        cds_seq = ""
        for index, x in enumerate(query_coors):
            codon_aln = pal2nal(seq1.split()[index], get_seq(cds_dict, query, x[0]*3-3+1, x[1]*3, q_strand))
            cds_seq += codon_aln
    else:
        cds_seq = ""
        for index, x in enumerate(query_coors[::-1]):
            codon_aln = pal2nal(seq1.split()[index], get_seq(cds_dict, query, x[0]*3-3+1, x[1]*3, q_strand))
            cds_seq += codon_aln
    if strand == "+":
        pseudogene_seq = ""
        for index, x in enumerate(pseudogene_coors):
            codon_aln = pal2nal(seq2.split()[index], get_seq(genome_dict, chr, x[0], x[1], strand))
            pseudogene_seq += codon_aln
    else:
        pseudogene_seq = ""
        for index, x in enumerate(pseudogene_coors[::-1]):
            codon_aln = pal2nal(seq2.split()[index], get_seq(genome_dict, chr, x[0], x[1], strand))
            pseudogene_seq += codon_aln
    pep_aln.append([i, cds_seq, pseudogene_seq])
    #print ">" + i
    #print seq1
    #print seq2
    #print cds_seq
    #print pseudogene_seq

# axt format
OUT = open(sys.argv[4], 'w')
for i in pep_aln:
    x, y, z = i
    OUT.write(x + "\n" + y + "\n" + z + "\n\n")
