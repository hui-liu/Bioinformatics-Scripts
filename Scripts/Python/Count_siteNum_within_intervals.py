#!/usr/bin/env python
import sys

def locateBin(Start, End, site):
    return site >= Start and site <= End

# usage
USAGE = "\nusage: python %s [snp sites] [genome regions] [output]\n" % sys.argv[0]

if len(sys.argv) !=4:
    print USAGE
    sys.exit()
#
IN1 = open(sys.argv[1],'r')
IN2 = open(sys.argv[2],'r')
OUT = open(sys.argv[3],'w')

# {'Chr01': [...]}
snp = {}
for eachline in IN1:
    if "CHROM" == eachline.split()[0]: continue
    split=eachline.rstrip().split()
    snp.setdefault(split[0],[]).append(split[1])
#
temp_list = []
for line in IN2:
    l = line.split()
    temp_list.append(l)

# sort the temp_list with the start of region
bin_list = sorted(temp_list, key = lambda k: int(k[2]))

# warn: the position values should be sorted
count_list = [] # partial gene list 
#gene_list = [] # all gene list
index = 0

while index < len(bin_list):
    num = 0
    el = bin_list[index] # ['Potri.001G000100part1', 'Chr01', '2081', '2502']
    for site in snp[el[1]][:]:
        if int(site) < int(el[2]):
            snp[el[1]].remove(site)
        elif locateBin(int(el[2]), int(el[3]), int(site)):
            num += 1
        else:
#            count_list.append([el[0], num]) 
            count_list.append(el + [num]) 
            break
    if not len(snp[el[1]]): break
    index += 1

# the lis hold the gene list representing in the count_list
lis = []
for l in count_list:
    lis.append(l[0])

# if genes did not represent in the count_list, then append them into the count_list
for gene in bin_list:
    if gene[0] not in lis:
        count_list.append(gene + [0])
#
for line in count_list:
    OUT.write("%s\t%d\n" % ("\t".join(line[:4]), line[4]))

IN1.close()
IN2.close()
OUT.close()
