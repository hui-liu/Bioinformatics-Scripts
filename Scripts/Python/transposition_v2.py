#!/usr/bin/python
#Filename: transposition_v2.py
#Author: liuhui
#EMail: liuhui@bjfu.edu.cn
#Description: row to line

import sys

USAGE = "usage: python transposition_v2.py inputfile outfile"

if len(sys.argv) !=3:
    print USAGE
    sys.exit()

with open(sys.argv[1], 'r') as infile:
    lis = [x.split() for x in infile]

with open(sys.argv[2], 'w') as outfile:
    for x in zip(*lis):
        outfile.write("\t".join(x) + "\n")
