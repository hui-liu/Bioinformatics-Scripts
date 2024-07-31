import sys
import os
from Bio import SeqIO
from Bio import AlignIO
import argparse

def perc_identity(aln):
    i = 0
    s1, s2 = aln[0].seq, aln[1].seq
    s1, s2 = s1.replace("-", ""), s2.replace("-", "")
    for a in range(0,len(aln[0])):
        s = aln[:,a]
        if s[0] == s[1]:
            i += 1
    return 100*i/float(min(len(s1), len(s2)))

parser = argparse.ArgumentParser()
parser.add_argument('--input_file', type=str, help='The input fasta file', required=True)
args = parser.parse_args()

infa = args.input_file
#infa = sys.argv[1]
aln_fa = os.path.basename(infa) + "_aln"
os.system('mafft --auto --quiet --preservecase %s > %s' % (infa, aln_fa))
align = AlignIO.read(aln_fa, "fasta")
print(perc_identity(align))
