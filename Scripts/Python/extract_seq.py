#!/usr/bin/python
#Filename: extract_seq.py
#Author: liuhui
#EMail: liuhui@bjfu.edu.cn	
#Description: extract a pair of sequences from two input fasta files based on a pair of name list

import os
import sys
from Bio import SeqIO

USAGE = "\nusage: python extract_seq.py pair_list protein.fas genome.fas\n"

if len(sys.argv) < 3:
    print USAGE
    sys.exit()

# pairs dict
accids = open(sys.argv[1], 'r') # 读取 pair of list 文件

# 以 pair of list 文件的第一列为 keys，第二列为 values，创建字典
ACC_DICT = {}
for line in accids:
    line = line.rstrip()
    pdID, ptID = line.split(' ')
    ACC_DICT[pdID] = ptID

# id vs seq dict
pd_infa = SeqIO.parse(open(sys.argv[2]), 'fasta')
pt_infa = SeqIO.parse(open(sys.argv[3]), 'fasta')

# 创建以 FASTA 文件的 ID 为 keys， sequence 为 values 的字典
SEQ_DICT = {}
for rec in pd_infa:
    SEQ_DICT[rec.id] = str(rec.seq)

for rec in pt_infa:
    SEQ_DICT[rec.id] = str(rec.seq)

# extracting, pair sequences per directory
for pd, pt in ACC_DICT.items():
    os.mkdir(pd)
    pd_outfa = open(pd + "/" + pd + ".fa", 'w')
    pt_outfa = open(pd + "/" + pt + ".fa", 'w')
    pd_outfa.write(">" + pd + "\n" + SEQ_DICT[pd] + "\n")
    pt_outfa.write(">" + pt + "\n" + SEQ_DICT[pt] + "\n")
    pd_outfa.close()
    pt_outfa.close()

accids.close()
pd_infa.close()
pt_infa.close()
