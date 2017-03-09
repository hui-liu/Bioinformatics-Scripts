#!/usr/bin/python
#Filename: concatenate_multiple_rows_into_single_line.py
#Author: liuhui
#EMail: liuhui@bjfu.edu.cn	
#Description: concatenate multiple rows into single line
#Usage: python concatenate_multiple_rows_into_single_line.py inputfile outputfile
#example: python concatenate_multiple_rows_into_single_line.py pfam.txt pfam_out.txt 
import sys

dict = {}
with open(sys.argv[1], 'r') as IN:
    for line in IN:
        s = line.split()
        if s[0] not in dict:
            dict[s[0]] = s[1:]
        else:
            dict[s[0]] = dict[s[0]] + ['|'] + s[1:]

with open(sys.argv[2], 'w') as OUT:
    for k in dict:
        OUT.write("%s %s\n" % (k, " ".join(dict[k])))
