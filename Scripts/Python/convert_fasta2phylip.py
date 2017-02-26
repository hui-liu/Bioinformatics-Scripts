#!/usr/bin/python

import sys
from Bio import SeqIO

# usage
USAGE = "\nusage: python convert_fasta2phylip.py [input fasta file] [output phy file]\n"

if len(sys.argv) !=3:
    print USAGE
    sys.exit()

infile = sys.argv[1]
outfile = sys.argv[2]

sequence_list = [] # To keep order of sequence

sequence_dict = {}
for record in SeqIO.parse(open(infile, "r"), "fasta"):
    sequence_list.append(record.id)
    sequence_dict[record.id] = str(record.seq)

# Test length of the alignment:
alignment_length = 0
for gene in sequence_dict:
    if (alignment_length != 0) and (len(sequence_dict[gene]) != alignment_length):
        print "Error in alignment length, exit on error !!!"
        sys.exit()
    else:
        alignment_length = len(sequence_dict[gene])

number_of_seq = len(sequence_dict)

longest_id = sorted(sequence_dict.keys())[-1]

# Write alignment in Phylip format
phyfile = open(outfile, "w")
phyfile.write(str(number_of_seq)+" "+str(alignment_length)+"\n")

for gene in sequence_list:
    phyfile.write(gene.ljust(len(longest_id), ' ') + "   " + sequence_dict[gene] + "\n")
phyfile.close()
