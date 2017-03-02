#!/usr/bin/python
#Filename: MSA_triplets_gaps_removed.py
#Author: liuhui
#EMail: liuhui@bjfu.edu.cn	
#Description: remove triplets gaps "---" in the msa file

import sys
from Bio import SeqIO

# usage
USAGE = "\nusage: python convert_fasta2phylip.py [MSA of cds sequences fasta file] [output file]\n"

if len(sys.argv) !=3:
    print USAGE
    sys.exit()

def condon_compare(codon1, codon2):
     mismatches = []
     for pos in range(3):
         if codon1[pos] != codon2[pos]:
             mismatches.append(pos)
     return len(mismatches) / float(3)

# {"id1": ["ATG", "---"], "id2", ["AT-", "--A"]}
id_seq = {}
# id list -- to keep the order
id_lis = []

for rec in SeqIO.parse(open(sys.argv[1]), 'fasta'):
    id_lis.append(rec.id)
    id_seq[rec.id] = []
    count = 0
    for c in range(len(rec.seq) / 3):
        id_seq[rec.id].append(str(rec.seq)[count : count + 3])
        count += 3
        
# the number of codons = (cds sequences length / 3)
length = len(id_seq.values()[0])

# the sites of sequence
seq_sites = [s for s in range(length * 3)]

# {0: ["ATG", "ATG", "---". "ATC", ...], ...}
pos_codon = {}
for pos in range(length):
    pos_codon[pos] = []
    for id in id_seq:
        pos_codon[pos].append(id_seq[id][pos])

# if the frequency of gaps codon "---" larger than 0.6, we hold them 
# in the list of "gaps" from small to large in order.
gaps = []
for columns in pos_codon: # search all codons from start to end
    numbers = 0
    for p in range(len(pos_codon[columns])): # check for certain column codon
        if 0.0 == condon_compare("---",pos_codon[columns][p]): # gaps codon "---"
            numbers += 1
    if numbers / float(len(pos_codon[columns])) >= 0.6:
        gaps.append(columns)

# sites of gaps
gaps_sites = []
for i in gaps:
    for j in range(i * 3, i * 3 + 3):
        gaps_sites.append(j)

# remove gaps codon and then export the results
res = open(sys.argv[2], 'w')

for id in id_lis:
    num = 0
    for pos in range(len(id_seq[id])):
        if pos in gaps:
            id_seq[id].pop(pos - num)
            num += 1
    res.write(">" + id + "\n" + "".join(id_seq[id]) + "\n")
res.close()

# report the ColumnsMap
ColumnsMap = open(sys.argv[2] + ".cols", 'w')
ColumnsMap.write("#ColumnsMap" + "\t" + ", ".join(map(str, set(seq_sites) - set(gaps_sites))) + "\n")
ColumnsMap.close()
