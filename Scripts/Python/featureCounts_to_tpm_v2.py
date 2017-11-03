import sys
import os

# example input data
"""
# Program:featureCounts v1.5.3; Command:"featureCounts" "-p" "-Q" "10" "-t" "exon" "-g" "gene_id" "-s" "0" "-T" "10" "-a" "../mapping/illumina/test/pta.gtf" "-o" "read.counts" "../mapping/illumina/SRR1200298_sort.bam" "../mapping/illumina/SRR1200343_sort.bam"
Geneid  Chr     Start   End     Strand  Length  ../mapping/illumina/SRR1200298_sort.bam ../mapping/illumina/SRR1200343_sort.bam
PTA00000014     C27702832       281     363     -       83      0       0
PTA00000015     C32554308       14858   16516   +       1659    0       0
PTA00000016     scaffold789550  4351    4444    -       94      0       1
PTA00000017     C32300050       378     481     -       104     5       2
"""

# usage
USAGE ="\nTHis python script Convert raw feature counts generated from *featureCounts* to transcripts per million (TPM). \
        \n\nExample usage: python %s [featureCount result file] [output file]\n\n" % os.path.basename(sys.argv[0])

if len(sys.argv) != 3:
    print USAGE
    sys.exit()

"""
# RPKM
rpkm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(counts) * 1e6
}

# TPM
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
"""

# counts to tpm
def CountstoTPM(counts, lengths):
    tpm = {}
    rates = {i: float(counts[i]) / float(lengths[i]) for i in counts}
    Sum = sum([rates[i] for i in rates])
    tpm = {i: [rates[i] / float(Sum) * 1000000] for i in rates}
    return tpm

Counts = {}
Lengths = {}
sample_names = None
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line[0] == "#":
            continue
        elif line.split()[0] == "Geneid":
            sample_names = [os.path.basename(i) for i in line.split()[6:]]
        else:
            lsplit = line.split()
            Lengths[lsplit[0]] = float(lsplit[5])
            Counts[lsplit[0]] = lsplit[6:]

Gene_lis = [gene_id for gene_id in Counts]
TPM = None
temp_dict = {}
# calculate TPM value for each sample
for i in range(len(sample_names)):
    for j in Counts:
        temp_dict[j] = Counts[j][i]
    if not TPM:
        TPM = CountstoTPM(temp_dict, Lengths)
    elif TPM:
        temp_TPM = CountstoTPM(temp_dict, Lengths)
        for geneid in Gene_lis:
            TPM[geneid] = TPM[geneid] + temp_TPM[geneid]


OUT = open(sys.argv[2], 'w')
# print header
OUT.write("gene_id" + "\t" + "\t".join(sample_names) + "\n")

for i in TPM:
    OUT.write(i + "\t" + "\t".join([str(float('%.2f' % float(j))) for j in TPM[i]]) + "\n")
OUT.close()
