#!/usr/bin/python
#Filename: rev_comp.py
#Author: liuhui
#EMail: liuhui@bjfu.edu.cn
#Description: generate the reverse complementary strand

import sys
"""python rev_comp.py input.fa out.fa"""

complement_table = {
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
'Y': 'R',
'a': 't',
'b': 'v',
'c': 'g',
'd': 'h',
'g': 'c',
'h': 'd',
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

def pqrse_fasta(seqs):
    new_seqs = {}
    for line in seqs:
        if line.startswith(">"):
            name = line.rstrip()
            new_seqs[name] = ""
        else:
            new_seqs[name] = new_seqs[name] + line.rstrip()
    return new_seqs

def rev_comp(seq):
    new_seq = []
    line = seq.rstrip()
    for letter in line:
        complement_letter = complement_table[letter]
        new_seq.append(complement_letter)
    new_seq.reverse()
    return "".join(new_seq)


in_file = open(sys.argv[1])
out_file = open(sys.argv[2], 'w')

seqs = pqrse_fasta(in_file)

for name in seqs.keys():
    if name.endswith("-"):
        print >> out_file, name + '\n' + rev_comp(seqs[name])
    elif name.endswith("+"):
        print >> out_file, name + '\n' + seqs[name]
    else:
        print >> out_file, name + '\n' + rev_comp(seqs[name]) # 如果文件没有 '+' 或 '-' 号标记正负链，则默认为负链。
