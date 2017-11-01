import sys
import os

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
    rates = {i: counts[i][0] / float(lengths[i]) for i in counts}
    Sum = sum([rates[i] for i in rates])
    tpm = {i: rates[i] / float(Sum) * 1000000 for i in rates}
    return tpm

Counts = {}
Lengths = {}
with open(sys.argv[1], 'r') as f:
    for line in f:
        if line[0] == "#" or line.split()[0] == "Geneid":
            continue
        else:
            lsplit = line.split()
            Lengths[lsplit[0]] = float(lsplit[5])
            Counts.setdefault(lsplit[0], []).append(float(lsplit[6]))

TPM = CountstoTPM(Counts, Lengths)
OUT = open(sys.argv[2], 'w')

for i in TPM:
    OUT.write(i + "\t" + str(float('%.2f' % TPM[i])) + "\n")
OUT.close()
