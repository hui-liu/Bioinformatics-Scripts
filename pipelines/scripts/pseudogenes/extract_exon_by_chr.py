#!/usr/bin/env python

import os, sys

#USAGE = "\npython %s [pta.gff]\n" % sys.argv[0]

#if len(sys.argv) != 2:
#    print USAGE
#    sys.exit()

chrFiles = {}
exon_dict = {}
for l in open(sys.argv[1]):
    if l[0] == "#": continue
    if not l.rstrip(): continue
    f = l[:-1].split('\t')
    if f[2] != "exon": continue
    exon_dict.setdefault(f[0], []).append(f[3:5])

#
#with open(sys.argv[2],'r') as f:
#    alist = [line.rstrip() for line in f]

for k in exon_dict:
    #if len(exon_dict[k]) == 1: continue
    #if k in alist:
    lis_sorted = sorted(exon_dict[k], key=lambda x: int(x[1]))
    with open('./%s_exLocs' % k, 'w') as OUT:
        for exon in lis_sorted:
            OUT.write(k + "\t" + "-" + "\t" + "\t".join(exon) + "\n")
