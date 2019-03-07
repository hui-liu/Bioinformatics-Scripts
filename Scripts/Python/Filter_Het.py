import sys
from cyvcf2 import VCF
import numpy as np
import gzip
# Author: Hui Liu, liuhui@bjfu.edu.cn

"""
citation:
cyvcf2: fast, flexible variant analysis with Python
"""

inFile = sys.argv[1]
cutoff = float(sys.argv[2])

HetThreshod = 0.7
if HetThreshod != cutoff:
    HetThreshod = cutoff

# header
for line in gzip.open(inFile, 'r'):
    if line[0] == "#":
        print line.rstrip("\n")
    else: continue

# filter high het rate locis
vcf = VCF(inFile)
for variant in vcf:
    HetGenos = 0
    TotGenos = 0
    if variant.is_indel:
        print str(variant).rstrip("\n")
    # gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
    HetGenos = np.sum(variant.gt_types == 1)
    TotGenos = np.sum(variant.gt_types != 2)
    if TotGenos == 0:
        print str(variant).rstrip("\n")
    else:
        HetRate = HetGenos / float(TotGenos)
        if HetRate <= HetThreshod:
            print str(variant).rstrip("\n")
