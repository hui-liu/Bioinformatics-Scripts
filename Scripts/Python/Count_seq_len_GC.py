
#!/usr/bin/python
#Filename: Count_seq_len_GC.py
#Author: liuhui
#EMail: liuhui@bjfu.edu.cn	
#Description: calculate the length an GC content of sequences from a fasta files

import sys
from Bio import SeqIO

out_file = open(sys.argv[2], 'w')

for rec in SeqIO.parse(open(sys.argv[1]), 'fasta'):
    print >> out_file, rec.id, len(rec.seq), (rec.seq.count("C") + rec.seq.count("G")) / float(len(rec.seq))

out_file.close()
