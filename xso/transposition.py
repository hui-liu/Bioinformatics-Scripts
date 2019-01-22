#!/usr/bin/python
#Filename: transposition_v2.py
#Author: liuhui
#EMail: liuhui@bjfu.edu.cn
#Description: row to line

import sys

USAGE = "usage example: python transposition.py xso_exp_VST_mean_test.txt xso_expr_mat.txt > xso_genes.txt"

if len(sys.argv) !=3:
    print USAGE
    sys.exit()

exp = []
genes = []
with open(sys.argv[1], 'r') as infile:
    header = infile.readline().rstrip()
    for line in infile:
        lsp = line.split()
        genes.append(lsp[0])
        exp.append(lsp[1:])

for i in genes:
    print i

with open(sys.argv[2], 'w') as outfile:
    for x in zip(*exp):
        outfile.write("\t".join(x) + "\n")
