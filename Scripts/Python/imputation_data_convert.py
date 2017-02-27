#!/usr/bin/python
#Filename: imputation_data_convert.py
#Author: liuhui
#EMail: liuhui@bjfu.edu.cn
#Description:

import sys
import os

USAGE = "\nusage: python imputation_data_convert.py [vcf file] [ord file] [linkage file] [linkage group] [output file]\n"

if len(sys.argv) !=6:
    print USAGE
    sys.exit()

vcfFile = open(sys.argv[1], 'r')
ordFile = open(sys.argv[2], 'r')
linkage = open(sys.argv[3], 'r')
lg_num = sys.argv[4]
res = open(sys.argv[5], 'w')

# vcf file
# to keep the order of snp locus
locus_list = []
for line in vcfFile:
    if '#' == line[0]: continue
    f = line[:-1].split('\t')[2]
    locus_list.append(f)

vcfFile.close()

# ord file
# line number vs (position and phase) -- {line num: [position, phase]}
num_pos = {}

for line in ordFile:
    if '#' == line[0]: continue # ingore the header
    l = line[:-1].split('\t')
    if len(l) == 4: continue # ingore the "duplicate"
    num_pos[l[0]] = [l[1]] + [l[4].rstrip("-")]

ordFile.close()

# linkage file
# holding the population markers
lis = []
# hoding the population name
pop_name = []
for line in linkage:
    if '#' == line[0]: continue
    f = line[:-1].split('\t')
    if 'F' == f[1][-1]: continue # ingore the 'F' family
    if 'M' == f[1][-1]: continue # ingore the 'M' family
    pop_name.append(f[1])
    lis.append(f[6:])

linkage.close()

header = ["No", "Nr", "Locus", "Segregation", "Phase", "Classification", "LG", "Position"] + pop_name
res.write("\t".join(header) + "\n")

# transposition of linkage file
linkage_trans = []
for x in zip(*lis):
    linkage_trans.append(list(x))

# a list holding the locus and markers of the population
locus_markers = []
for num1 in xrange(len(locus_list)):
     locus_markers.append([locus_list[int(num1)]] + linkage_trans[int(num1)])


# To keep the order according the positions
num_pos_sorted = sorted(num_pos.items(), key = lambda x: float(x[1][0]))
numbers = [i[0] for i in num_pos_sorted]


count = 1
for num2 in numbers:
    res.write(str(count) + "\t" + str(count) + "\t" + locus_markers[int(num2)-1][0] + "\t" + "<lmxll>" + "\t" + num_pos[num2][1] + "\t" +
              "(ll,lm)" + "\t" + lg_num + "\t" + num_pos[num2][0] + "\t" + "\t".join(locus_markers[int(num2)-1][1:]) + "\n")
    count += 1

res.close()

os.system("sed -i 's/0 0/--/g;s/1 2/lm/g;s/2 2/ll/g' %s" % sys.argv[5])
